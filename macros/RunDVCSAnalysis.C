#include "./../DreamAN/core/DVCSAnalysis.h"
#include "./../DreamAN/core/EventProcessor.h"
#include "./../DreamAN/core/AnalysisTaskManager.h"
#include <iostream>

void RunDVCSAnalysis(const std::string& inputDir) {
    if (inputDir.empty()) {
        std::cerr << "Usage: RunDVCSAnalysis(<inputDir>)" << std::endl;
        return;
    }

    AnalysisTaskManager mgr;
    mgr.CreateOutputFile("/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/rdf_dvcs_output.root", "DVCS");

    // Track Cuts
    auto* trackCuts = new TrackCut();
    trackCuts->SetECALEdgeCut(9, 100000000);
    trackCuts->SetDCEdgeCut(1, 100000000); 

    // Photon Cuts
    auto* photonCuts = new EventCut();
    photonCuts->SetChargeCut(0);
    photonCuts->SetPIDCountCut(22, 1, 1);

    // Electron Cuts
    auto* electronCuts = new EventCut();
    electronCuts->SetChargeCut(-1);
    electronCuts->SetPIDCountCut(11, 1, 1);

    // Proton Cuts
    auto* protonCuts = new EventCut();
    protonCuts->SetChargeCut(1);
    protonCuts->SetPIDCountCut(2212, 1, 1);

    // Task
    auto dvcsTask = std::make_unique<DVCSAnalysis>();
    dvcsTask->SetTrackCuts(trackCuts);
    dvcsTask->SetPhotonCuts(photonCuts);
    dvcsTask->SetElectronCuts(electronCuts);
    dvcsTask->SetProtonCuts(protonCuts);

    mgr.AddTask(std::move(dvcsTask));

    // Processor
    EventProcessor processor(inputDir, mgr);
    processor.ProcessEvents();
}
