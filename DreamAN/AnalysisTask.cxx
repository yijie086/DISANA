#include "AnalysisTask.h"
#include "AnalysisTaskManager.h"
#include <iostream>

AnalysisTask::AnalysisTask()
    : fTaskManager(nullptr) {}

AnalysisTask::~AnalysisTask() {}

void AnalysisTask::SaveOutput() {
    if (fTaskManager) {
        fTaskManager->SaveOutput();
    } else {
        std::cerr << "Error: Task Manager not set!" << std::endl;
    }
}

void AnalysisTask::AddHistogram(const std::string& name, TH1F* hist) {
    if (fTaskManager) {
        fTaskManager->AddHistogram(name, hist);
    }
}

void AnalysisTask::AddTree(const std::string& name, TTree* tree) {
    if (fTaskManager) {
        fTaskManager->AddTree(name, tree);
    }
}
