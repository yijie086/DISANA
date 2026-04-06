#ifndef ANALYSISTASK_H
#define ANALYSISTASK_H

#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>

#include <ROOT/RDF/RInterface.hxx>
#include <map>
#include <string>

#include "RHipoDS.hxx"

class AnalysisTaskManager;  // forward declare

class AnalysisTask {
 public:
  AnalysisTask();
  virtual ~AnalysisTask();

  virtual void UserCreateOutputObjects() = 0;
  virtual void UserExec(ROOT::RDF::RNode& df) = 0;
  virtual void SaveOutput() = 0;
  virtual void Terminate(Option_t* = nullptr) {}

  void SetTaskManager(AnalysisTaskManager* mgr) { fTaskManager = mgr; }

  // New virtual method to receive output file pointer
  virtual void SetOutputFile(TFile* file) {}
  virtual void SetOutputDir(const std::string& dir) {}
  template <typename Lambda>
  ROOT::RDF::RNode DefineOrRedefine(ROOT::RDF::RNode df, const std::string& name, Lambda&& lambda, const std::vector<std::string>& columns) {
    auto existingCols = df.GetColumnNames();
    if (std::find(existingCols.begin(), existingCols.end(), name) != existingCols.end()) {
      return df.Redefine(name, std::forward<Lambda>(lambda), columns);
    }
    return df.Define(name, std::forward<Lambda>(lambda), columns);
  }

  void SafeSnapshot(ROOT::RDF::RNode df, const std::string& treename, const std::string& filename, const std::vector<std::string>& excludeCols = {"EventCutResult"}) {
    auto allCols = df.GetColumnNames();
    std::vector<std::string> outputCols;

    for (const auto& col : allCols) {
      if (std::find(excludeCols.begin(), excludeCols.end(), col) == excludeCols.end()) {
        outputCols.push_back(col);
      }
    }

    df.Snapshot(treename, filename, outputCols);
  }

  // SelectiveSnapshot: snapshot only the columns in wantedCols that actually
  // exist in df.  Columns not present in df are silently skipped, so this is
  // safe to call on both the full HIPO-backed dataframe and on a re-processing
  // dataframe that was read back from a previously snapshotted .root file.
  void SelectiveSnapshot(ROOT::RDF::RNode df, const std::string& treename,
                         const std::string& filename,
                         const std::vector<std::string>& wantedCols) {
    auto allCols = df.GetColumnNames();
    std::vector<std::string> outputCols;
    outputCols.reserve(wantedCols.size());
    for (const auto& col : wantedCols) {
      if (std::find(allCols.begin(), allCols.end(), col) != allCols.end()) {
        outputCols.push_back(col);
      }
    }
    df.Snapshot(treename, filename, outputCols);
  }

  // ResolveSnapshotColumns: returns the subset of wantedCols that exist in df.
  // Use this to build the column list before calling lazy df.Snapshot() directly,
  // so that Count() and Snapshot() can be booked together and share one event loop.
  // Used by the optimised (fOptimizeColumns=true) snapshot path.
  std::vector<std::string> ResolveSnapshotColumns(ROOT::RDF::RNode df,
                                                   const std::vector<std::string>& wantedCols) {
    auto allCols = df.GetColumnNames();
    std::vector<std::string> result;
    result.reserve(wantedCols.size());
    for (const auto& col : wantedCols) {
      if (std::find(allCols.begin(), allCols.end(), col) != allCols.end())
        result.push_back(col);
    }
    return result;
  }

  // SafeSnapshotColumns: returns all columns present in df, minus any in excludeCols.
  // Used by the full (fOptimizeColumns=false) snapshot path.  The default exclusion
  // list contains "EventCutResult", which is an internal struct ROOT cannot serialise.
  std::vector<std::string> SafeSnapshotColumns(ROOT::RDF::RNode df,
                                                const std::vector<std::string>& excludeCols = {"EventCutResult"}) {
    auto allCols = df.GetColumnNames();
    std::vector<std::string> result;
    result.reserve(allCols.size());
    for (const auto& col : allCols) {
      if (std::find(excludeCols.begin(), excludeCols.end(), col) == excludeCols.end())
        result.push_back(col);
    }
    return result;
  }

 protected:
  AnalysisTaskManager* fTaskManager = nullptr;
};

#endif
