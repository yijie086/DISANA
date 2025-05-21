/*#ifndef ANALYSISTASKMANAGER_H
#define ANALYSISTASKMANAGER_H

#include "AnalysisTask.h"
#include <vector>

class AnalysisTaskManager {
    public:
        void AddTask(std::unique_ptr<AnalysisTask> task) {
            tasks.push_back(std::move(task)); // Move ownership
        }
        void ExecuteTasks(Event& event) {
            for (auto& task : tasks) {
                task->UserExec(event);
            }
        }
    
    private:
        std::vector<std::unique_ptr<AnalysisTask>> tasks; // Store smart pointers
};

#endif // ANALYSISTASKMANAGER_H
*/
#ifndef ANALYSISTASKMANAGER_H
#define ANALYSISTASKMANAGER_H

#include "AnalysisTask.h"
#include <vector>
#include <map>
#include <memory>
#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TH1F.h>

class AnalysisTaskManager {
public:
    AnalysisTaskManager() = default;

    // Set output directory for all ROOT files
    void SetOutputDir(const std::string& dir) {
        outputDir = dir;
    }

    // Create all output objects (files, trees, histograms)
    void UserCreateOutputObjects();

    // Create a ROOT output file
    void CreateOutputFile(const std::string& name);

    // Create a TTree inside an output file
    void CreateTree(const std::string& fileName, const std::string& treeName);

    // Create a histogram list container
    void CreateHistogramList(const std::string& fileName, const std::string& listName);

    // Add a histogram to a list
    void AddHistogram(const std::string& listName, TH1F* hist);

    // Add a branch to a tree
    template <typename T>
    void AddTreeBranch(const std::string& treeName, const std::string& branchName, T* address);

       // Add a tree to the list of trees
    void AddTree(const std::string& name, TTree* tree);
    
    // Add an analysis task
    void AddTask(std::unique_ptr<AnalysisTask> task);

    // Execute tasks for an event
    void ExecuteTasks(Event& event);

    // Save all output
    void SaveOutput();

    ~AnalysisTaskManager() = default;

private:
    std::vector<std::unique_ptr<AnalysisTask>> tasks;
    std::map<std::string, std::unique_ptr<TFile>> outputFiles;
    std::map<std::string, std::unique_ptr<TTree>> trees;
    std::map<std::string, std::unique_ptr<TList>> histLists;
    std::string outputDir = "./"; // Default output directory

    // Output variables
    float momentum;
    TH1F* hMomentum;
};

#endif // ANALYSISTASKMANAGER_H

#ifndef ANALYSISTASKMANAGER_H
#define ANALYSISTASKMANAGER_H

#include "AnalysisTask.h"
#include <vector>
#include <map>
#include <memory>
#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TH1F.h>

class AnalysisTaskManager {
public:
    AnalysisTaskManager() = default;

    // Set output directory for all ROOT files
    void SetOutputDir(const std::string& dir) {
        outputDir = dir;
    }

    // Create all output objects (files, trees, histograms)
    void UserCreateOutputObjects();

    // Create a ROOT output file
    void CreateOutputFile(const std::string& name);

    // Create a TTree inside an output file
    void CreateTree(const std::string& fileName, const std::string& treeName);

    // Create a histogram list container
    void CreateHistogramList(const std::string& fileName, const std::string& listName);

    // Add a histogram to a list
    void AddHistogram(const std::string& listName, TH1F* hist);

    // Add a tree to the list of trees
    void AddTree(const std::string& name, TTree* tree);

    // Add an analysis task
    void AddTask(std::unique_ptr<AnalysisTask> task);

    // Execute tasks for an event
    void ExecuteTasks(Event& event);

    // Save all output
    void SaveOutput();

    ~AnalysisTaskManager() = default;

private:
    std::vector<std::unique_ptr<AnalysisTask>> tasks;
    std::map<std::string, std::unique_ptr<TFile>> outputFiles;
    std::map<std::string, std::unique_ptr<TTree>> trees;
    std::map<std::string, std::unique_ptr<TList>> histLists;
    std::string outputDir = "./"; // Default output directory

    // Output variables
    float momentum;
    TH1F* hMomentum;
};

#endif // ANALYSISTASKMANAGER_H

