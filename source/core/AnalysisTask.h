#ifndef ANALYSISTASK_H
#define ANALYSISTASK_H

#include <ROOT/RDF/RInterface.hxx>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <map>
#include <string>

class AnalysisTaskManager; // forward declare

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

protected:
    AnalysisTaskManager* fTaskManager = nullptr;
};

#endif
