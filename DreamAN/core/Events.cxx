#include "Events.h"
#include <iostream>

Events::Events(const std::string& directory) {
    inputFiles = GetHipoFilesInPath(directory);

    if (inputFiles.empty()) {
        std::cerr << "No .hipo files found in directory: " << directory << std::endl;
        return;
    }

    std::cout << "Creating RHipoDS from input files..." << std::endl;
    dataSource = std::make_unique<RHipoDS>(inputFiles);
    
    ROOT::RDF::RNode df = ROOT::RDataFrame(std::move(dataSource));
    dfNodePtr = std::make_shared<ROOT::RDF::RNode>(df);
    dfNode = std::make_optional<ROOT::RDF::RNode>(*dfNodePtr);
    std::cout << "DataFrame initialized with " << inputFiles.size() << " input files." << std::endl;
}

std::vector<std::string> Events::GetHipoFilesInPath(const std::string& directory) {
    std::vector<std::string> files;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(directory)) {
        if (entry.path().extension() == ".hipo") {
            files.push_back(entry.path().string());
            //std::cout << "Found file: " << entry.path() << std::endl;
        }
    }
    std::cout << "================ " << files.size() << " Files Found ================" << std::endl;
    return files;
}

std::optional<ROOT::RDF::RNode> Events::getNode() const {
    return dfNode;
}

size_t Events::getFileCount() const {
    return inputFiles.size();
}
