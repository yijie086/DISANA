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

    void AddTask(std::unique_ptr<AnalysisTask> task);
    void UserCreateOutputObjects();
    void Execute(ROOT::RDF::RNode& df);
    void SaveOutput();

    void SetOututDir(const std::string& Outputdir="./",const std::string& filename ="AnalysisResults.root", const std::string& directory = "AnalysisResults");
    void AddHistogram(const std::string& name, TH1* hist);
    void AddTree(const std::string& name, TTree* tree);

    // New: Notify tasks of output file
    void SetOutputFileForTasks();

private:
    std::vector<std::unique_ptr<AnalysisTask>> tasks;
    std::map<std::string, TH1*> histograms;
    std::map<std::string, TTree*> trees;
    std::unique_ptr<TFile> outputFile;
    std::string outputDir;
    std::string outputRootDir;
};

#endif
