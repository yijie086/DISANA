#ifndef EVENTS_H
#define EVENTS_H

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include "hipo4/RHipoDS.hxx"

class Events {
public:
  Events(const std::string& inputDirectory,
         const std::string& outputDirectory,
         bool fIsReprocessRootFile,
         const std::string& fInputROOTtreeName,
         const std::string& fOutputROOTfileName,
         int nfiles,
         int nthreads);

  ~Events();

  std::optional<ROOT::RDF::RNode> getNode() const;
  std::size_t getFileCount() const;
  std::string getFinalInputPath() const;

private:
  std::vector<std::string> GetHipoFilesInPath(const std::string& directory, int nfiles);

  std::string fOutputDir_;

  bool        fIsReprocessRootFile_{false};
  int         fnfiles_{0};
  int         fnthreads_{0};
  std::string fInputROOTtreeName_;
  std::string fOutputROOTfileName_;

  std::optional<ROOT::RDF::RNode> dfNode_;
  std::size_t fileCount_{0};
  std::string finalInputPath_;

  std::vector<std::string> inputFiles;
  std::unique_ptr<RHipoDS> dataSource;
};

#endif // EVENTS_H