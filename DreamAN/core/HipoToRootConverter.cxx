#include "HipoToRootConverter.h"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <vector>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "TInterpreter.h"
#include "TSystem.h"
#include "hipo4/RHipoDS.hxx"
#include "Compression.h"
#include "RVersion.h"

namespace fs = std::filesystem;

// ---- Static helper functions (preload_hipo, makeRDFFromSlice, etc.) remain unchanged ----

static inline void preload_hipo() {
  if (const char* h = std::getenv("HIPO_HOME")) {
    std::string inc = std::string(h) + "/include";
    gInterpreter->AddIncludePath(inc.c_str());
  } else {
    gInterpreter->AddIncludePath(
      "/cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/"
      "almalinux9-gcc11/local/hipo/4.2.0/include");
  }
  gSystem->Load("libhipo4");
}

static ROOT::RDataFrame makeRDFFromSlice(std::vector<std::string> slice) {
  auto ds = std::make_unique<RHipoDS>(slice);
  return ROOT::RDataFrame(std::move(ds));
}

static std::unordered_map<std::string, std::string>
inspectSliceSchema(const std::vector<std::string>& slice) {
  ROOT::RDataFrame rdf =
      makeRDFFromSlice(std::vector<std::string>(slice.begin(), slice.end()));
  std::unordered_map<std::string, std::string> types;
  for (const auto& c : rdf.GetColumnNames()) types.emplace(c, rdf.GetColumnType(c));
  return types;
}

static std::vector<std::string>
computeStableColumns(const std::vector<std::vector<std::string>>& slices) {
  std::unordered_map<std::string, std::string> common;
  bool first = true;
  for (const auto& s : slices) {
    auto types = inspectSliceSchema(s);
    if (first) { common = std::move(types); first = false; continue; }
    for (auto it = common.begin(); it != common.end();) {
      auto jt = types.find(it->first);
      if (jt == types.end() || jt->second != it->second) it = common.erase(it);
      else ++it;
    }
    if (common.empty()) break;
  }
  std::vector<std::string> cols; cols.reserve(common.size());
  for (auto& kv : common) cols.push_back(kv.first);
  std::sort(cols.begin(), cols.end());
  return cols;
}

static void snapshotSliceToRoot(std::vector<std::string> sliceFiles,
                                const std::string& outPath,
                                const std::string& treeName,
                                const std::vector<std::string>& cols) {
  ROOT::RDataFrame rdf = makeRDFFromSlice(std::move(sliceFiles));
  ROOT::RDF::RSnapshotOptions opts;
  opts.fMode = "RECREATE";
  #if defined(ROOT_VERSION_CODE) && ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
  opts.fCompressionAlgorithm = ROOT::ECompressionAlgorithm::kLZ4;
  #else
  opts.fCompressionAlgorithm = ROOT::ECompressionAlgorithm::kZLIB;
  #endif
  opts.fCompressionLevel = 4;
  opts.fAutoFlush = 50'000;
  rdf.Snapshot(treeName, outPath, cols, opts);
}

// ---- Constructor and getHipoFilesInPath remain unchanged ----

HipoToRootConverter::HipoToRootConverter(const std::string& inputDir,
                                         const std::string& outputDir,
                                         int nFiles, int nThreads)
  : fInputDir_(inputDir), fOutputDir_(outputDir),
    fnFiles_(nFiles), fnThreads_(nThreads) {
  if (!fs::exists(fOutputDir_)) fs::create_directories(fOutputDir_);
}

std::vector<std::string>
HipoToRootConverter::getHipoFilesInPath(const std::string& directory, int nfiles) const {
    // ... implementation is unchanged ...
    std::vector<std::string> files;
    std::size_t count = 0;
    try {
        for (const auto& entry : fs::recursive_directory_iterator(directory)) {
            if (entry.is_regular_file() && entry.path().extension() == ".hipo") {
                files.push_back(entry.path().string());
                ++count;
                if (nfiles > 0 && static_cast<int>(count) >= nfiles) break;
            }
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "[Converter] Filesystem error: " << e.what() << std::endl;
    }
    std::cout << "[Converter] Found " << files.size() << " .hipo files\n";
    return files;
}


// --- CHANGED: `convertAndMerge` is now `convert` and returns the file list ---
std::vector<std::string> HipoToRootConverter::convert(const std::string& tempFilePrefix) {
  preload_hipo();

  // Keep conversion single-threaded at the RDF level per thread
  ROOT::DisableImplicitMT();
  ROOT::EnableThreadSafety();

  std::vector<std::string> hipoFiles = getHipoFilesInPath(fInputDir_, fnFiles_);
  lastInputCount_ = hipoFiles.size();
  if (hipoFiles.empty()) {
    std::cerr << "[Converter] No .hipo files found in: " << fInputDir_ << std::endl;
    return {};
  }
  std::sort(hipoFiles.begin(), hipoFiles.end());

  const std::size_t nWorkers = std::max<std::size_t>(
      1, std::min<std::size_t>(fnThreads_ > 0 ? (std::size_t)fnThreads_ : 1, hipoFiles.size()));
  const std::size_t filesPerThread = (hipoFiles.size() + nWorkers - 1) / nWorkers;

  std::vector<std::vector<std::string>> slices;
  slices.reserve(nWorkers);
  for (std::size_t i = 0; i < nWorkers; ++i) {
    const auto start = i * filesPerThread;
    const auto end   = std::min(start + filesPerThread, hipoFiles.size());
    if (start < end) slices.emplace_back(hipoFiles.begin()+start, hipoFiles.begin()+end);
  }

  auto stableCols = computeStableColumns(slices);
  if (stableCols.empty())
    throw std::runtime_error("[Converter] No common, type-consistent columns across slices.");

  std::vector<std::string> tempRoots;
  tempRoots.reserve(slices.size());
  std::vector<std::thread> threads;
  threads.reserve(slices.size());

  for (std::size_t i = 0; i < slices.size(); ++i) {
    std::string tmp = (fs::path(fOutputDir_) / (tempFilePrefix + std::to_string(i) + ".root")).string();
    tempRoots.push_back(tmp);
    threads.emplace_back(snapshotSliceToRoot, std::move(slices[i]), tmp,
                         std::string(kSnapshotTreeName), stableCols);
  }
  for (auto& t : threads) if (t.joinable()) t.join();

  // --- REMOVED: All TFileMerger logic is gone ---

  std::cout << "[Converter] Parallel conversion complete. Generated "
            << tempRoots.size() << " temporary ROOT files." << std::endl;
            
  return tempRoots;
}