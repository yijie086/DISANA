#include "Events.h"

#include <filesystem>
#include <iostream>
#include <stdexcept>

#include "HipoToRootConverter.h"
#include "TROOT.h"

namespace fs = std::filesystem;

Events::Events(const std::string& directory, const std::string& outputDirectory, bool fIsReprocessRootFile,
               const std::string& fInputROOTtreeName, const std::string& fOutputROOTfileName,
               int nfiles, int nthreads)
    : fOutputDir_(outputDirectory),
      fIsReprocessRootFile_(fIsReprocessRootFile),
      fnfiles_(nfiles),
      fnthreads_(nthreads),
      fInputROOTtreeName_(fInputROOTtreeName),
      fOutputROOTfileName_(fOutputROOTfileName) {
  try {
    if (!fIsReprocessRootFile_) {
      HipoToRootConverter converter(directory, fOutputDir_, fnfiles_, fnthreads_);

      // --- CHANGED: Use the new `convert` method and get a vector of file paths ---
      fTempRootFiles_ = converter.convert("temp_hipo_conversion_");
      fileCount_ = converter.lastInputCount();
      finalInputPath_ = ""; // No single input path anymore

      if (fTempRootFiles_.empty()) {
        throw std::runtime_error("No .hipo files were converted from: " + directory);
      }
      std::cout << "[Events] Using " << fTempRootFiles_.size()
                << " temporary ROOT files for analysis.\n";

      if (fnthreads_ > 0) ROOT::EnableImplicitMT(fnthreads_);

      // --- CHANGED: Construct RDataFrame from the vector of file paths ---
      auto rdf = ROOT::RDataFrame(std::string(HipoToRootConverter::kSnapshotTreeName), fTempRootFiles_);
      dfNode_.emplace(rdf);

    } else {
      finalInputPath_ = (fs::path(directory) / fOutputROOTfileName_).string();
      if (!fs::exists(finalInputPath_)) {
        throw std::runtime_error("Reprocess mode: ROOT file not found: " + finalInputPath_);
      }
      fileCount_ = 0; // Or determine from the file if needed
      std::cout << "[Events] Reprocessing existing ROOT file: " << finalInputPath_ << "\n";

      if (fnthreads_ > 0) ROOT::EnableImplicitMT(fnthreads_);

      const std::string treeName = fInputROOTtreeName_.empty()
                                       ? std::string(HipoToRootConverter::kSnapshotTreeName)
                                       : fInputROOTtreeName_;
      auto rdf = ROOT::RDataFrame(treeName, finalInputPath_);
      dfNode_.emplace(rdf);
    }

    std::cout << "[Events] DataFrame initialized successfully.\n";
  } catch (const std::exception& e) {
    std::cerr << "[Events] ERROR: " << e.what() << std::endl;
    throw;
  }
}

// --- NEW: Destructor implementation for cleaning up temp files ---
Events::~Events() {
  if (!fTempRootFiles_.empty()) {
    std::cout << "[Events] Cleaning up " << fTempRootFiles_.size() << " temporary files...\n";
    for (const auto& path : fTempRootFiles_) {
      std::error_code ec;
      fs::remove(path, ec);
      if (ec) {
        std::cerr << "[Events] Warning: could not remove temporary file: " << path
                  << " (" << ec.message() << ")\n";
      }
    }
  }
}

std::optional<ROOT::RDF::RNode> Events::getNode() const { return dfNode_; }
std::size_t Events::getFileCount() const { return fileCount_; }
std::string Events::getFinalInputPath() const { return finalInputPath_; }