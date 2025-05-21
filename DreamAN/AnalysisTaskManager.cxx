#include "AnalysisTaskManager.h"
#include <iostream>

void AnalysisTaskManager::UserCreateOutputObjects() {
    std::cout << "Initializing output objects..." << std::endl;

    CreateOutputFile("analysis_output.root");
    CreateHistogramList("analysis_output.root", "DVCS_Histograms");

    for (auto& task : tasks) {
        task->UserCreateOutputObjects();
    }

    std::cout << "Output objects created!" << std::endl;
}


void AnalysisTaskManager::CreateOutputFile(const std::string& name) {
    std::string filePath = outputDir + "/" + name;
    outputFiles[name] = std::make_unique<TFile>(filePath.c_str(), "RECREATE");
}

void AnalysisTaskManager::CreateTree(const std::string& fileName, const std::string& treeName) {
    if (outputFiles.count(fileName)) {
        outputFiles[fileName]->cd();
        trees[treeName] = std::make_unique<TTree>(treeName.c_str(), "Event Data");
    } else {
        std::cerr << "Error: Output file " << fileName << " not found!" << std::endl;
    }
}

void AnalysisTaskManager::CreateHistogramList(const std::string& fileName, const std::string& listName) {
    if (outputFiles.count(fileName)) {
        outputFiles[fileName]->cd();
        histLists[listName] = std::make_unique<TList>();
        histLists[listName]->SetName(listName.c_str()); // Set TList name
    } else {
        std::cerr << "Error: Output file " << fileName << " not found!" << std::endl;
    }
}


void AnalysisTaskManager::AddHistogram(const std::string& listName, TH1F* hist) {
    if (histLists.count(listName)) {
        histLists[listName]->Add(hist);
    } else {
        std::cerr << "Error: Histogram list " << listName << " not found!" << std::endl;
    }
}

template <typename T>
void AnalysisTaskManager::AddTreeBranch(const std::string& treeName, const std::string& branchName, T* address) {
    if (trees.count(treeName)) {
        trees[treeName]->Branch(branchName.c_str(), address);
    } else {
        std::cerr << "Error: Tree " << treeName << " not found!" << std::endl;
    }
}

void AnalysisTaskManager::AddTask(std::unique_ptr<AnalysisTask> task) {
    tasks.push_back(std::move(task));
}

void AnalysisTaskManager::AddTree(const std::string& name, TTree* tree) {
    // Check if the tree already exists
    if (trees.count(name)) {
        std::cerr << "Error: Tree with name " << name << " already exists!" << std::endl;
    } else {
        // Clone the TTree instead of copying
        if (outputFiles.count("analysis_output.root")) {
            outputFiles["analysis_output.root"]->cd();  // Make sure we are in the correct ROOT file
            trees[name] = std::make_unique<TTree>(tree->GetName(), tree->GetTitle());
            tree->CopyAddresses(trees[name].get());  // Copy branches and structure of the original tree
        } else {
            std::cerr << "Error: Output file not found!" << std::endl;
        }
    }
}

void AnalysisTaskManager::ExecuteTasks(Event& event) {
    for (auto& task : tasks) {
        task->UserExec(event);
    }
}

void AnalysisTaskManager::SaveOutput() {
    for (auto& filePair : outputFiles) {
        filePair.second->cd();
        
        // Write histograms grouped by their lists
        for (auto& listPair : histLists) {
            TDirectory *dir = filePair.second->mkdir(listPair.first.c_str());
            dir->cd();
            listPair.second->Write(); 
        }

        // Write trees
        for (auto& treePair : trees) {
            treePair.second->Write();  
        }

        filePair.second->Close();
    }
    std::cout << "ROOT output saved!" << std::endl;
}
