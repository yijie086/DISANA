#include "HipoToRootConverter.h"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <unordered_set>

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TFileMerger.h"
#include "hipo4/RHipoDS.hxx"

namespace fs = std::filesystem;

// Build an RDF over a slice (by value so vector is mutable for RHipoDS)
static ROOT::RDataFrame makeRDFFromSlice(std::vector<std::string> slice) {
  auto ds = std::make_unique<RHipoDS>(slice);  // RHipoDS takes vector<string>&
  return ROOT::RDataFrame(std::move(ds));
}

// Inspect one slice: returns map<columnName, cxxTypeName>
static std::unordered_map<std::string, std::string> inspectSliceSchema(const std::vector<std::string>& slice) {
  ROOT::RDataFrame rdf = makeRDFFromSlice(std::vector<std::string>(slice.begin(), slice.end()));
  std::unordered_map<std::string, std::string> types;
  for (const auto& c : rdf.GetColumnNames()) {
    types.emplace(c, rdf.GetColumnType(c));  // e.g. "float", "ROOT::VecOps::RVec<float>"
  }
  return types;
}

// Compute intersection of columns whose type matches across all slices.
static std::vector<std::string> computeStableColumns(const std::vector<std::vector<std::string>>& slices) {
  std::unordered_map<std::string, std::string> common;  // name -> type
  bool first = true;

  for (const auto& s : slices) {
    auto types = inspectSliceSchema(s);
    if (first) {
      common = std::move(types);
      first = false;
      continue;
    }
    for (auto it = common.begin(); it != common.end();) {
      auto jt = types.find(it->first);
      if (jt == types.end() || jt->second != it->second)
        it = common.erase(it);
      else
        ++it;
    }
    if (common.empty()) break;
  }

  std::vector<std::string> cols;
  cols.reserve(common.size());
  for (auto& kv : common) cols.push_back(kv.first);
  std::sort(cols.begin(), cols.end());
  return cols;
}

// Snapshot a slice to ROOT with a fixed column list
static void snapshotSliceToRoot(std::vector<std::string> sliceFiles, const std::string& outPath, const std::string& treeName, const std::vector<std::string>& cols) {
  ROOT::RDataFrame rdf = makeRDFFromSlice(std::move(sliceFiles));

  ROOT::RDF::RSnapshotOptions opts;
  opts.fMode = "RECREATE";  // overwrite temp if exists
  rdf.Snapshot(treeName, outPath, cols, opts);
}

HipoToRootConverter::HipoToRootConverter(const std::string& inputDir, const std::string& outputDir, int nFiles, int nThreads)
    : fInputDir_(inputDir), fOutputDir_(outputDir), fnFiles_(nFiles), fnThreads_(nThreads) {
  if (!fs::exists(fOutputDir_)) fs::create_directories(fOutputDir_);
}

std::string HipoToRootConverter::convertAndMerge(const std::string& outputFileName) {
  std::vector<std::string> hipoFiles = getHipoFilesInPath(fInputDir_, fnFiles_);
  lastInputCount_ = hipoFiles.size();
  if (hipoFiles.empty()) {
    std::cerr << "[Converter] No .hipo files found in: " << fInputDir_ << std::endl;
    return {};
  }

  std::sort(hipoFiles.begin(), hipoFiles.end());  // determinism

  const std::size_t nWorkers = std::max<std::size_t>(1, std::min<std::size_t>(fnThreads_ > 0 ? static_cast<std::size_t>(fnThreads_) : 1, hipoFiles.size()));
  const std::size_t filesPerThread = (hipoFiles.size() + nWorkers - 1) / nWorkers;

  // Build slices
  std::vector<std::vector<std::string>> slices;
  slices.reserve(nWorkers);
  for (std::size_t i = 0; i < nWorkers; ++i) {
    const std::size_t start = i * filesPerThread;
    const std::size_t end = std::min(start + filesPerThread, hipoFiles.size());
    if (start >= end) continue;
    slices.emplace_back(hipoFiles.begin() + start, hipoFiles.begin() + end);
  }

  // NEW: determine a stable, type-consistent column set across all slices
  auto stableCols = computeStableColumns(slices);
  if (stableCols.empty()) {
    throw std::runtime_error(
        "[Converter] No common, type-consistent columns across slices. "
        "Consider reducing to a known-safe subset.");
  }
  // Optional: if you *know* some problematic columns, drop them here by name.

  // Fast path: if only one slice/temp, skip TFileMerger entirely
  if (slices.size() == 1) {
    const std::string merged = (fs::path(fOutputDir_) / outputFileName).string();
    const std::string tmp = (fs::path(fOutputDir_) / "temp_single.root").string();

    snapshotSliceToRoot(std::move(slices[0]), tmp, std::string(kSnapshotTreeName), stableCols);

    std::error_code ec;
    fs::remove(merged, ec);  // best-effort remove old
    fs::rename(tmp, merged, ec);
    if (ec) {
      // fallback: copy then remove
      fs::copy_file(tmp, merged, fs::copy_options::overwrite_existing, ec);
      fs::remove(tmp, ec);
    }
    std::cout << "[Converter] Conversion complete. Single slice -> " << merged << std::endl;
    return merged;
  }

  // Multi-slice: snapshot per slice in parallel, then merge trees
  std::vector<std::string> tempRoots;
  tempRoots.reserve(slices.size());
  std::vector<std::thread> threads;
  threads.reserve(slices.size());

  for (std::size_t i = 0; i < slices.size(); ++i) {
    const std::string tmpPath = (fs::path(fOutputDir_) / ("temp_" + std::to_string(i) + ".root")).string();
    tempRoots.push_back(tmpPath);

    threads.emplace_back(snapshotSliceToRoot, std::move(slices[i]), tmpPath, std::string(kSnapshotTreeName), stableCols);
  }
  for (auto& t : threads)
    if (t.joinable()) t.join();

  const std::string merged = (fs::path(fOutputDir_) / outputFileName).string();

  TFileMerger merger;
  merger.OutputFile(merged.c_str(), "RECREATE");
  for (const auto& f : tempRoots) merger.AddFile(f.c_str());
  if (!merger.Merge()) {
    for (const auto& f : tempRoots) {
      std::error_code ec;
      fs::remove(f, ec);
    }
    throw std::runtime_error("TFileMerger failed");
  }
  for (const auto& f : tempRoots) {
    std::error_code ec;
    fs::remove(f, ec);
  }

  std::cout << "[Converter] Conversion complete. Merged file: " << merged << std::endl;
  return merged;
}

std::vector<std::string> HipoToRootConverter::getHipoFilesInPath(const std::string& directory, int nfiles) const {
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
