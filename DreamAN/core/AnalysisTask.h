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

 protected:
  AnalysisTaskManager* fTaskManager = nullptr;
};

#endif
