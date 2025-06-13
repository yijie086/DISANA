#include "Events.h"

#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

// Constructor
Events::Events(const std::string& directory, bool fIsReprocessRootFile, const std::string& fInputROOTtreeName, const std::string& fInputROOTfileName)
    : fIsReprocessRootFile(fIsReprocessRootFile), fInputROOTtreeName(fInputROOTtreeName), fInputROOTfileName(fInputROOTfileName) {
  if (fIsReprocessRootFile) {
    std::string inputfile_Root = directory + fInputROOTfileName;
    std::cout << "Reprocessing ROOT files is enabled." << std::endl;

    auto rdf = ROOT::RDataFrame(fInputROOTtreeName, inputfile_Root);
    dfNodePtr = std::make_shared<ROOT::RDF::RNode>(rdf);
  } else {
    std::cout << "Reprocessing ROOT files is disabled." << std::endl;

    inputFiles = GetHipoFilesInPath(directory);
    if (inputFiles.empty()) {
      std::cerr << "No .hipo files found in directory: " << directory << std::endl;
      return;
    }

    std::cout << "Creating RHipoDS from input files..." << std::endl;
    dataSource = std::make_unique<RHipoDS>(inputFiles);

    auto rdf = ROOT::RDataFrame(std::move(dataSource));
    dfNodePtr = std::make_shared<ROOT::RDF::RNode>(rdf);
  }

  dfNode = std::make_optional<ROOT::RDF::RNode>(*dfNodePtr);

  std::cout << "DataFrame initialized with " << inputFiles.size() << " input files." << std::endl;
}

// Helper to get HIPO files in a path
std::vector<std::string> Events::GetHipoFilesInPath(const std::string& directory) {
  std::vector<std::string> files;
  for (const auto& entry : std::filesystem::recursive_directory_iterator(directory)) {
    if (entry.path().extension() == ".hipo") {
      files.push_back(entry.path().string());
    }
  }
  std::cout << "================ " << files.size() << " Files Found ================" << std::endl;
  return files;
}

// Accessor methods
std::optional<ROOT::RDF::RNode> Events::getNode() const { return dfNode; }

size_t Events::getFileCount() const { return inputFiles.size(); }
