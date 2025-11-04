#pragma once
#include <string>
#include <vector>
#include <cstddef>

class HipoToRootConverter {
public:
  static constexpr const char* kSnapshotTreeName = "dst";

  HipoToRootConverter(const std::string& inputDir,
                      const std::string& outputDir,
                      int nFiles,
                      int nThreads);

  // CHANGED: This method now returns a vector of temporary ROOT file paths.
  std::vector<std::string> convert(const std::string& tempFilePrefix);

  std::size_t lastInputCount() const { return lastInputCount_; }
  std::vector<std::string> getHipoFilesInPath(const std::string& directory, int nfiles) const;

private:


  std::string fInputDir_;
  std::string fOutputDir_;
  int fnFiles_{0};
  int fnThreads_{0};
  std::size_t lastInputCount_{0};
};