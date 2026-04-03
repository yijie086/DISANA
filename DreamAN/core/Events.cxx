#include "TSystem.h"
#include "TInterpreter.h"
#include "TROOT.h"

#include "Events.h"

#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace fs = std::filesystem;

Events::Events(const std::string& directory, const std::string& outputDirectory,
               bool fIsReprocessRootFile,
               const std::string& fInputROOTtreeName,
               const std::string& fOutputROOTfileName,
               int nfiles, int nthreads)
  : fOutputDir_(outputDirectory),
    fIsReprocessRootFile_(fIsReprocessRootFile),
    fnfiles_(nfiles),
    fnthreads_(nthreads),
    fInputROOTtreeName_(fInputROOTtreeName),
    fOutputROOTfileName_(fOutputROOTfileName)
{
  try {
    // ------------------------------------------------------------------------
    // REPROCESS MODE: read an existing ROOT file / tree
    // ------------------------------------------------------------------------
    if (fIsReprocessRootFile_) {
      finalInputPath_ = (fs::path(directory) / fOutputROOTfileName_).string();
      if (!fs::exists(finalInputPath_)) {
        throw std::runtime_error("Reprocess mode: ROOT file not found: " + finalInputPath_);
      }

      fileCount_ = 0;
      std::cout << "[Events] Reprocessing existing ROOT file: " << finalInputPath_ << "\n";

      if (fnthreads_ == 0) {
        ROOT::EnableImplicitMT();
      } else if (fnthreads_ > 1) {
        ROOT::EnableImplicitMT(fnthreads_);
      }

      const std::string treeName =
        fInputROOTtreeName_.empty()
          ? "hipo"
          : fInputROOTtreeName_;

      auto rdf = ROOT::RDataFrame(treeName, finalInputPath_);
      dfNode_.emplace(rdf);

      std::cout << "[Events] DataFrame initialized successfully.\n";
      return;
    }

    // ------------------------------------------------------------------------
    // HIPO MODE: use RHipoDS directly
    // ------------------------------------------------------------------------

    // Collect .hipo input files
    inputFiles = this->getHipoFilesInPath(directory, nfiles);

    if (inputFiles.empty()) {
      throw std::runtime_error("No .hipo files found in directory: " + directory);
    }
    for (const auto& file : inputFiles) {
        std::cout << "Input file: " << file << std::endl;
    }

    fileCount_ = inputFiles.size();
    finalInputPath_.clear();

    // Enable ROOT MT
    if (fnthreads_ == 1) {
      ROOT::EnableImplicitMT(1);
    } else if (fnthreads_ > 1) {
      ROOT::EnableImplicitMT(fnthreads_);
    }

    std::cout << "[Events] Creating RHipoDS from " << inputFiles.size() << " input file(s)...\n";
    dataSource = std::make_unique<RHipoDS>(inputFiles, 1000000);

    auto rdf = ROOT::RDataFrame(std::move(dataSource));
    dfNode_.emplace(rdf);

    std::cout << "[Events] DataFrame initialized successfully.\n";
  }
  catch (const std::exception& e) {
    std::cerr << "[Events] ERROR: " << e.what() << std::endl;
    throw;
  }
}

Events::~Events() {
  // No temporary files are created in RHipoDS mode.
}

std::optional<ROOT::RDF::RNode> Events::getNode() const { return dfNode_; }
std::size_t Events::getFileCount() const { return fileCount_; }
std::string Events::getFinalInputPath() const { return finalInputPath_; }
std::vector<std::string>Events::getHipoFilesInPath(const std::string& pathOrPattern,
                                        int nfiles) const {
  std::vector<std::string> files;
  std::size_t count = 0;

  auto push = [&](const fs::path& p) {
    files.push_back(p.string());
    ++count;
    return !(nfiles > 0 && static_cast<int>(count) >= nfiles);
  };

  const fs::path p(pathOrPattern);
  std::error_code ec;

  // CASE 0: wildcard pattern in filename
  const std::string base = p.filename().string();
  if (base.find_first_of("*?") != std::string::npos) {
    const fs::path parent =
        p.has_parent_path() ? p.parent_path() : fs::path(".");

    if (!fs::exists(parent, ec) || !fs::is_directory(parent, ec)) {
      std::cerr << "[Converter] Parent directory not found for pattern: "
                << parent << "\n";
      std::cout << "[Converter] Found " << files.size() << " .hipo files\n";
      return files;
    }

    for (fs::directory_iterator it(parent, ec), end;
         !ec && it != end; it.increment(ec)) {
      const fs::directory_entry& de = *it;
      if (!de.is_regular_file(ec)) continue;

      const fs::path& ep = de.path();
      if (ep.extension() != ".hipo") continue;

      if (fnmatch(base.c_str(), ep.filename().string().c_str(), 0) == 0) {
        if (!push(ep)) break;
      }
    }

    std::cout << "[Converter] Found " << files.size()
              << " .hipo files (pattern)\n";
    return files;
  }

  // CASE 1: existing file
  if (fs::is_regular_file(p, ec)) {
    if (p.extension() == ".hipo") push(p);
    std::cout << "[Converter] Found " << files.size()
              << " .hipo files (file)\n";
    return files;
  }

  // CASE 2: existing directory (recursive)
  if (fs::is_directory(p, ec)) {
    fs::recursive_directory_iterator it(
        p, fs::directory_options::skip_permission_denied, ec),
        end;

    for (; !ec && it != end; it.increment(ec)) {
      const fs::directory_entry& de = *it;
      if (de.is_regular_file(ec) && de.path().extension() == ".hipo") {
        if (!push(de.path())) break;
      }
    }

    if (ec) std::cerr << "[Converter] Iteration warning: " << ec.message()
                      << "\n";

    std::cout << "[Converter] Found " << files.size()
              << " .hipo files (dir)\n";
    return files;
  }

  // CASE 3: not found and not a glob
  std::cerr << "[Converter] Path not found: " << p << "\n";
  std::cout << "[Converter] Found " << files.size() << " .hipo files\n";
  return files;
}
