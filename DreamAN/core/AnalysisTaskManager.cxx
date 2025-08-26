#include "AnalysisTaskManager.h"
#include "AnalysisTask.h"
#include <TFile.h>

AnalysisTaskManager::AnalysisTaskManager() {}
AnalysisTaskManager::~AnalysisTaskManager() {
}

void AnalysisTaskManager::AddTask(std::shared_ptr<AnalysisTask> task) {
    task->SetTaskManager(this);
    tasks.push_back(std::move(task));
}

void AnalysisTaskManager::UserCreateOutputObjects() {
    for (auto& task : tasks) task->UserCreateOutputObjects();
}

void AnalysisTaskManager::Execute(ROOT::RDF::RNode& df) {
    for (auto& task : tasks) task->UserExec(df);
}

void AnalysisTaskManager::SetOututDir(const std::string& Outputdir) {
    outputDir = Outputdir;
}

void AnalysisTaskManager::SaveOutput() {
    if (outputDir =="./") std::cerr << "[SaveOutput] the default output dir is ./!" << std::endl;
    for (auto& task : tasks) {
        task->SetOutputDir(outputDir);
        task->SaveOutput();
    }
}

