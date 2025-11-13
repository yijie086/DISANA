#include "HipoToRootConverter.h"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <vector>
#include <fnmatch.h>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "TInterpreter.h"
#include "TSystem.h"
#include "hipo4/RHipoDS.hxx"
#include "Compression.h"
#include "RVersion.h"

namespace fs = std::filesystem;

// ---- Static helper functions (preload_hipo, makeRDFFromSlice, etc.) ----

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

// Inspect a slice's schema; if the slice has zero entries, return empty (caller will ignore)
static std::unordered_map<std::string, std::string>
inspectSliceSchema(const std::vector<std::string>& slice) {
  ROOT::RDataFrame rdf =
      makeRDFFromSlice(std::vector<std::string>(slice.begin(), slice.end()));
  // Count is lazy; force evaluation
  auto n = *rdf.Count();
  if (n == 0) {
    std::cerr << "[Converter] Schema discovery: slice has 0 rows, ignoring in intersection\n";
    return {};
  }
  std::unordered_map<std::string, std::string> types;
  for (const auto& c : rdf.GetColumnNames()) types.emplace(c, rdf.GetColumnType(c));
  return types;
}

// Compute columns common to all *non-empty* slices. Empty slices are ignored.
static std::vector<std::string>
computeStableColumns(const std::vector<std::vector<std::string>>& slices) {
  std::unordered_map<std::string, std::string> common;
  bool haveSeed = false;

  for (const auto& s : slices) {
    auto types = inspectSliceSchema(s);
    if (types.empty()) {
      // ignore empty/no-tree slices for the intersection
      continue;
    }
    if (!haveSeed) {
      common = std::move(types);
      haveSeed = true;
      continue;
    }
    for (auto it = common.begin(); it != common.end();) {
      auto jt = types.find(it->first);
      if (jt == types.end() || jt->second != it->second) it = common.erase(it);
      else ++it;
    }
    if (common.empty()) break;
  }

  if (!haveSeed) {
    throw std::runtime_error("[Converter] All slices appear empty; cannot determine stable schema.");
  }

  std::vector<std::string> cols; cols.reserve(common.size());
  for (auto& kv : common) cols.push_back(kv.first);
  std::sort(cols.begin(), cols.end());
  return cols;
}

// Snapshot a slice to ROOT only if it has rows; otherwise, skip writing the file.
static void snapshotSliceToRoot(std::vector<std::string> sliceFiles,
                                const std::string& outPath,
                                const std::string& treeName,
                                const std::vector<std::string>& cols) {
  ROOT::RDataFrame rdf = makeRDFFromSlice(std::move(sliceFiles));

  // Skip empty slices to avoid producing empty/no-tree files that break the chain
  auto n = *rdf.Count();
  if (n == 0) {
    std::cerr << "[Converter] Slice " << outPath
              << " is empty -> not writing a ROOT file\n";
    return;
  }

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

// ---- Constructor and getHipoFilesInPath ----

HipoToRootConverter::HipoToRootConverter(const std::string& inputDir,
                                         const std::string& outputDir,
                                         int nFiles, int nThreads)
  : fInputDir_(inputDir), fOutputDir_(outputDir),
    fnFiles_(nFiles), fnThreads_(nThreads) {
  if (!fs::exists(fOutputDir_)) fs::create_directories(fOutputDir_);
}

// Minimal '*' and '?' glob matcher (case-sensitive)
static bool globMatch(const std::string& pat, const std::string& s) {
  size_t pi = 0, si = 0, star = std::string::npos, match = 0;
  while (si < s.size()) {
    if (pi < pat.size() && (pat[pi] == '?' || pat[pi] == s[si])) {
      ++pi; ++si;                   // single-char match
    } else if (pi < pat.size() && pat[pi] == '*') {
      star = pi++; match = si;      // remember star position
    } else if (star != std::string::npos) {
      pi = star + 1;                // backtrack: extend '*' to cover one more char
      si = ++match;
    } else {
      return false;
    }
  }
  while (pi < pat.size() && pat[pi] == '*') ++pi;
  return pi == pat.size();
}

std::vector<std::string>
HipoToRootConverter::getHipoFilesInPath(const std::string& pathOrPattern, int nfiles) const {
  namespace fs = std::filesystem;

  std::vector<std::string> files;
  std::size_t count = 0;
  auto push = [&](const fs::path& p) {
    files.push_back(p.string());
    ++count;
    return !(nfiles > 0 && static_cast<int>(count) >= nfiles);
  };

  const fs::path p(pathOrPattern);
  std::error_code ec;

  // CASE 0: wildcard pattern
  const std::string base = p.filename().string();
  if (base.find_first_of("*?") != std::string::npos) {
    const fs::path parent = p.has_parent_path() ? p.parent_path() : fs::path(".");
    if (!fs::exists(parent, ec) || !fs::is_directory(parent, ec)) {
      std::cerr << "[Converter] Parent directory not found for pattern: " << parent << "\n";
      std::cout << "[Converter] Found " << files.size() << " .hipo files\n";
      return files;
    }

    for (fs::directory_iterator it(parent, ec), end; !ec && it != end; it.increment(ec)) {
      const fs::directory_entry& de = *it;
      if (!de.is_regular_file(ec)) continue;
      const fs::path& ep = de.path();
      if (ep.extension() != ".hipo") continue;
      if (fnmatch(base.c_str(), ep.filename().string().c_str(), 0) == 0) {
        if (!push(ep)) break;
      }
    }
    std::cout << "[Converter] Found " << files.size() << " .hipo files (pattern)\n";
    return files;
  }

  // CASE 1: existing file
  if (fs::is_regular_file(p, ec)) {
    if (p.extension() == ".hipo") push(p);
    std::cout << "[Converter] Found " << files.size() << " .hipo files (file)\n";
    return files;
  }

  // CASE 2: existing directory (recursive)
  if (fs::is_directory(p, ec)) {
    fs::recursive_directory_iterator it(p, fs::directory_options::skip_permission_denied, ec), end;
    for (; !ec && it != end; it.increment(ec)) {
      const fs::directory_entry& de = *it;
      if (de.is_regular_file(ec) && de.path().extension() == ".hipo") {
        if (!push(de.path())) break;
      }
    }
    if (ec) std::cerr << "[Converter] Iteration warning: " << ec.message() << "\n";
    std::cout << "[Converter] Found " << files.size() << " .hipo files (dir)\n";
    return files;
  }

  // CASE 3: not found and not a glob
  std::cerr << "[Converter] Path not found: " << p << "\n";
  std::cout << "[Converter] Found " << files.size() << " .hipo files\n";
  return files;
}

// --- convert(): parallel per-slice snapshots with stable schema & pruning ---
std::vector<std::string> HipoToRootConverter::convert(const std::string& tempFilePrefix) {
  preload_hipo();

  // Keep conversion single-threaded inside each worker (RDF-level)
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
    throw std::runtime_error("[Converter] No common, type-consistent columns across NON-empty slices.");

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

  // Prune any missing or empty files (slices that were empty won't have written a file)
  std::vector<std::string> pruned;
  pruned.reserve(tempRoots.size());
  for (auto& f : tempRoots) {
    std::error_code ec;
    if (fs::exists(f, ec) && !ec && fs::file_size(f, ec) > 0 && !ec) {
      pruned.push_back(f);
    } else {
      std::cerr << "[Converter] Dropping missing/empty temp file: " << f << "\n";
    }
  }
  tempRoots.swap(pruned);

  std::cout << "[Converter] Parallel conversion complete. Generated "
            << tempRoots.size() << " temporary ROOT files." << std::endl;

  return tempRoots;
}
