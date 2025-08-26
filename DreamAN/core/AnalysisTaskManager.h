#ifndef ANALYSISTASKMANAGER_H
#define ANALYSISTASKMANAGER_H

#include <memory>
#include <vector>
#include <map>
#include <string>
#include <ROOT/RDF/RInterface.hxx>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

class AnalysisTask;

class AnalysisTaskManager {
public:
    AnalysisTaskManager();
    ~AnalysisTaskManager();

    void AddTask(std::shared_ptr<AnalysisTask> task);
    void UserCreateOutputObjects();
    void Execute(ROOT::RDF::RNode& df);
    void SaveOutput();
    void SetOututDir(const std::string& Outputdir="./");
private:
    std::vector<std::shared_ptr<AnalysisTask>> tasks;
    std::string outputDir;

};

#endif
