#include "Events.h"

#include <filesystem>
#include <iostream>
#include <stdexcept>

#include "HipoToRootConverter.h"
#include "TROOT.h"

namespace fs = std::filesystem;

Events::Events(const std::string& directory, const std::string& outputDirectory, bool fIsReprocessRootFile, const std::string& fInputROOTtreeName,
               const std::string& fOutputROOTfileName, int nfiles, int nthreads)
    : fOutputDir_(outputDirectory),
      fIsReprocessRootFile_(fIsReprocessRootFile),
      fnfiles_(nfiles),
      fnthreads_(nthreads),
      fInputROOTtreeName_(fInputROOTtreeName),
      fOutputROOTfileName_(fOutputROOTfileName) {
  try {
    if (!fIsReprocessRootFile_) {
      ROOT::DisableImplicitMT();
      ROOT::EnableThreadSafety();
      HipoToRootConverter converter(directory, fOutputDir_, fnfiles_, fnthreads_);
      finalInputPath_ = converter.convertAndMerge("converted_hipo_ROOT.root");
      fileCount_ = converter.lastInputCount();

      if (finalInputPath_.empty()) throw std::runtime_error("No .hipo files found in: " + directory);

      std::cout << "[Events] Converted/merged into: " << finalInputPath_ << "\n";

      ROOT::EnableImplicitMT();
      auto rdf = ROOT::RDataFrame(std::string(HipoToRootConverter::kSnapshotTreeName), finalInputPath_);
      dfNode_.emplace(rdf);
    } else {
      finalInputPath_ = (fs::path(directory) / fOutputROOTfileName_).string();
      if (!fs::exists(finalInputPath_)) throw std::runtime_error("Reprocess mode: ROOT file not found: " + finalInputPath_);
      fileCount_ = 0;
      std::cout << "[Events] Reprocessing existing ROOT file: " << finalInputPath_ << "\n";

      ROOT::EnableImplicitMT();
      const std::string treeName = fInputROOTtreeName_.empty() ? std::string(HipoToRootConverter::kSnapshotTreeName) : fInputROOTtreeName_;
      auto rdf = ROOT::RDataFrame(treeName, finalInputPath_);
      dfNode_.emplace(rdf);
    }

    std::cout << "[Events] DataFrame initialized on: " << finalInputPath_ << "\n";
  } catch (const std::exception& e) {
    std::cerr << "[Events] ERROR: " << e.what() << std::endl;
    throw;
  }
}

std::optional<ROOT::RDF::RNode> Events::getNode() const { return dfNode_; }
std::size_t Events::getFileCount() const { return fileCount_; }
std::string Events::getFinalInputPath() const { return finalInputPath_; }
