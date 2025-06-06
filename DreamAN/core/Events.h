#ifndef EVENTS_H
#define EVENTS_H

#include "FilesInPath.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "RHipoDS.hxx"

#include <memory>
#include <optional>
#include <string>
#include <vector>

class Events {
public:
    // Constructor
    explicit Events(const std::string& directory);

    // Accessor for the ROOT RDataFrame node
    std::optional<ROOT::RDF::RNode> getNode() const;

    // Number of input files
    size_t getFileCount() const;

private:
    std::vector<std::string> inputFiles;
    std::unique_ptr<RHipoDS> dataSource;
    std::shared_ptr<ROOT::RDF::RNode> dfNodePtr;
    std::optional<ROOT::RDF::RNode> dfNode;
};

#endif // EVENTS_H
