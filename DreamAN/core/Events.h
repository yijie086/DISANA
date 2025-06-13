#ifndef EVENTS_H
#define EVENTS_H

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "RHipoDS.hxx"

#include <string>
#include <vector>
#include <memory>
#include <optional>

class Events {
public:
  Events(const std::string& directory, bool fIsReprocessRootFile,
         const std::string& fInputROOTtreeName, const std::string& fInputROOTfileName);

  std::optional<ROOT::RDF::RNode> getNode() const;
  size_t getFileCount() const;

private:
  std::vector<std::string> GetHipoFilesInPath(const std::string& directory);

  bool fIsReprocessRootFile;
  std::string fInputROOTtreeName;
  std::string fInputROOTfileName;
  std::vector<std::string> inputFiles;

  std::unique_ptr<RHipoDS> dataSource;
  std::shared_ptr<ROOT::RDF::RNode> dfNodePtr;
  std::optional<ROOT::RDF::RNode> dfNode;
};
#endif // EVENTS_H