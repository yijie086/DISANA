#include "Events.h"

#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "HipoToRootConverter.h" // used only for file discovery helper
#include "TInterpreter.h"
#include "TROOT.h"
#include "TSystem.h"
#include "hipo4/RHipoDS.hxx"

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
    // ------------------------------------------------------------------------
    // REPROCESS MODE: read an existing ROOT file / tree
    // ------------------------------------------------------------------------
    if (fIsReprocessRootFile_) {
      finalInputPath_ = (fs::path(directory) / fOutputROOTfileName_).string();
      if (!fs::exists(finalInputPath_)) {
        throw std::runtime_error("Reprocess mode: ROOT file not found: " + finalInputPath_);
      }
      fileCount_ = 0; // can be filled by inspecting the file if desired
      std::cout << "[Events] Reprocessing existing ROOT file: " << finalInputPath_ << "\n";

      if (fnthreads_ == 0) {
        ROOT::EnableImplicitMT();          // ROOT decides
      } else if (fnthreads_ > 1) {
        ROOT::EnableImplicitMT(fnthreads_);// fixed threads
      }                                   // fnthreads_==1 -> single thread

      const std::string treeName = fInputROOTtreeName_.empty()
                                       ? std::string(HipoToRootConverter::kSnapshotTreeName)
                                       : fInputROOTtreeName_;
      auto rdf = ROOT::RDataFrame(treeName, finalInputPath_);
      dfNode_.emplace(rdf);

      std::cout << "[Events] DataFrame initialized successfully.\n";
      return;
    }

    // ------------------------------------------------------------------------
    // HIPO MODE: use RHipoDS directly (NO slicing / NO temp ROOT files)
    // ------------------------------------------------------------------------

    // 1) Collect .hipo input files
    HipoToRootConverter converter(".", fOutputDir_, nfiles, nthreads); // used only to locate files
    inputFiles = converter.getHipoFilesInPath(directory, nfiles);

    if (inputFiles.empty()) {
      throw std::runtime_error("No .hipo files found in directory: " + directory);
    }

    fileCount_ = inputFiles.size();
    finalInputPath_.clear();

    // 2) Load headers/libs BEFORE enabling IMT / creating the DataFrame
    /*if (const char* h = std::getenv("HIPO_HOME")) {
      std::string inc = std::string(h) + "/include";
      gInterpreter->AddIncludePath(inc.c_str());
    }
    gSystem->Load("libhipo4");
    gSystem->Load("libHipoDataFrame");*/

    // 3) Enable ROOT MT if requested
    //    - nthreads == 0 : ROOT decides
    //    - nthreads == 1 : single-thread
    //    - nthreads >  1 : fixed number of threads
    if (fnthreads_ == 0) {
      ROOT::EnableImplicitMT(4);
    } else if (fnthreads_ > 1) {
      ROOT::EnableImplicitMT(fnthreads_);
    }

    // 4) Create RHipoDS from ALL input files; let RDF handle parallelism
    std::cout << "[Events] Creating RHipoDS from " << inputFiles.size() << " input file(s)...\n";
    dataSource = std::make_unique<RHipoDS>(inputFiles);

    auto rdf = ROOT::RDataFrame(std::move(dataSource));
    dfNode_.emplace(rdf);

    std::cout << "[Events] DataFrame initialized successfully.\n";

  } catch (const std::exception& e) {
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