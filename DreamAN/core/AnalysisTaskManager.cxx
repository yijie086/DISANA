#include "AnalysisTaskManager.h"
#include "AnalysisTask.h"
#include <TFile.h>

AnalysisTaskManager::AnalysisTaskManager() {}
AnalysisTaskManager::~AnalysisTaskManager() {
}

void AnalysisTaskManager::AddTask(std::unique_ptr<AnalysisTask> task) {
    task->SetTaskManager(this);
    tasks.push_back(std::move(task));
}

void AnalysisTaskManager::UserCreateOutputObjects() {
    for (auto& task : tasks) task->UserCreateOutputObjects();
}

void AnalysisTaskManager::Execute(ROOT::RDF::RNode& df) {
    for (auto& task : tasks) task->UserExec(df);
}

void AnalysisTaskManager::CreateOutputFile(const std::string& filename, const std::string& directory) {
    outputFile = std::make_unique<TFile>(filename.c_str(), "RECREATE");
    outputDir = directory;
    outputFile->mkdir(directory.c_str());
}

void AnalysisTaskManager::AddHistogram(const std::string& name, TH1* hist) {
    if (!hist) {
        std::cerr << "[AddHistogram] Warning: NULL histogram for " << name << std::endl;
        return;
    }
    histograms[name] = hist;
}

void AnalysisTaskManager::AddTree(const std::string& name, TTree* tree) {
    trees[name] = tree;
}

void AnalysisTaskManager::SetOutputFileForTasks() {
    if (!outputFile) return;
    for (auto& task : tasks) {
        task->SetOutputFile(outputFile.get());
    }
}

void AnalysisTaskManager::SaveOutput() {
    if (!outputFile) {
        std::cerr << "[SaveOutput] No output file!" << std::endl;
        return;
    }
    outputFile->cd(outputDir.c_str());
    SetOutputFileForTasks();
    for (size_t i = 0; i < tasks.size(); ++i) {
        tasks[i]->SaveOutput();
    }
    for (const auto& [name, hist] : histograms) {
        if (hist) {
            hist->Write(name.c_str());
        } else {
            std::cerr << "  Null histogram: " << name << std::endl;
        }
    }

    for (const auto& [name, tree] : trees) {
        if (tree) {
            std::cout << "  Writing tree: " << name << std::endl;
            tree->Write(name.c_str());
        } else {
            std::cerr << "  Null tree: " << name << std::endl;
        }
    }
    outputFile->Close();
    std::cout << "Outuput file saved: " << outputFile->GetName() << std::endl;
}

