#ifndef ANALYSISTASK_H
#define ANALYSISTASK_H

#include "Event.h"
#include <memory>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

class AnalysisTask {
public:
    AnalysisTask();
    virtual ~AnalysisTask();

    class AnalysisTaskManager* fTaskManager; // Reference to the task manager for output management
    
    // Function to initialize output objects (e.g., histograms, trees, etc.)
    virtual void UserCreateOutputObjects() = 0;

    // Function to execute the task on a single event
    virtual void UserExec(Event& event) = 0;

    // Save the output to file (to be called after event processing)
    void SaveOutput();

    // Set the task manager for the task (for managing output objects)
    void SetTaskManager(class AnalysisTaskManager* manager) { fTaskManager = manager; }

protected:
    // Add histograms, trees, etc., to the task manager
    void AddHistogram(const std::string& name, TH1F* hist);
    void AddTree(const std::string& name, TTree* tree);

private:
    
};

#endif // ANALYSISTASK_H
