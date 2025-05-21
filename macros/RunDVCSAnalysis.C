#include "./../DreamAN/DVCSAnalysis.h"
#include "./../DreamAN/EventProcessor.h"
#include "./../DreamAN/AnalysisTaskManager.h"
#include <iostream>

void RunDVCSAnalysis(const std::string& inputFile) {
    bool testFile = false; // just to see if file exists;
    if (inputFile == "") {
        std::cerr << "Usage: " << inputFile << " <input_file>" << std::endl;
        return;
    }

    if(testFile){
    std::string fileName = "/cache/clas12/rg-b/production/recon/spring2019/torus-1/pass2/v0/dst/train/jpsi/jpsi_006334.hipo";
    clas12::clas12reader reader(fileName);
    
    std::cout << "File name: " << reader.getFilename() << std::endl;
    std::cout << "Is file open: " << reader.isOpen() << std::endl;
    
    if (!reader.isOpen()) {
        std::cerr << "Failed to open file." << std::endl;
        return ;
    }
    return ;
    }

    AnalysisTaskManager taskManager;

    TrackCuts *PhotonCuts =TrackCuts::PhotonCuts();
    PhotonCuts->SetChargeCut(0);
    PhotonCuts->SetPIDCut(22,3.0); // set pdg and the chisquared

    TrackCuts *ElectronCuts =TrackCuts::ElectronCuts();
    ElectronCuts->SetChargeCut(-1);
    ElectronCuts->SetPIDCut(11,3.0); // set pdg and the chisquared

    TrackCuts *ProtonCuts =TrackCuts::ProtonCuts();
    ProtonCuts->SetChargeCut(1);
    ProtonCuts->SetPIDCut(2212,3.0); // set pdg and the chisquared


    /// here we call our task and set the track event cut settings
    std::unique_ptr<DVCSAnalysis> DVCStask = std::make_unique<DVCSAnalysis>();

    DVCStask->SetPhotonCuts(PhotonCuts);
    DVCStask->SetElectronCuts(ElectronCuts);
    DVCStask->SetProtonCuts(ProtonCuts);

    /// pass the task to the task manager to do the laborious stuff
    taskManager.AddTask(std::move(DVCStask)); // Use smart pointer to memory stuff
    std::cout << "adding the task in the DVCS analysis" << std::endl;
   
    std::cout << "Going to Event Processor" << std::endl;
    EventProcessor processor(inputFile, taskManager);
    std::cout << "Event Proceesing is in progress " << std::endl;
    processor.ProcessEvents(); // Start event loop
}
