#ifndef HIPO_TO_ROOT_CONVERTER_H
#define HIPO_TO_ROOT_CONVERTER_H

#include <string>
#include <string_view>
#include <vector>

class HipoToRootConverter {
 public:
  // Tree name written by Snapshot
  static inline constexpr std::string_view kSnapshotTreeName{"clas12"};

  HipoToRootConverter(const std::string& inputDir, const std::string& outputDir, int nFiles, int nThreads);

  // Convert per-slice using RHipoDS -> Snapshot(), then merge temp files.
  // Returns the merged ROOT path (empty on failure/no inputs).
  std::string convertAndMerge(const std::string& outputFileName);

  std::size_t lastInputCount() const { return lastInputCount_; }

 private:
  std::vector<std::string> getHipoFilesInPath(const std::string& directory, int nfiles) const;

  const std::string fInputDir_;
  const std::string fOutputDir_;
  const int fnFiles_;
  const int fnThreads_;
  std::size_t lastInputCount_{0};
};

#endif  // HIPO_TO_ROOT_CONVERTER_H
