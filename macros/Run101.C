#include <string>
#include <iostream>
#include "RHipoDS.hxx"
#include "TCanvas.h"
#include <optional>
#include <vector>

#include "../source/ParticleInformation/RECParticle.h"
#include "../source/ParticleInformation/RECTraj.h"
#include "../source/Cuts/ElectronCut.h"
#include "../source/DrawHist/DrawAndSave.h"
#include "../source/core/FilesInPath.h"
#include "../source/core/Columns.h"
#include "../source/Cuts/EventCut.h"
#include "../source/Cuts/TrackCut.h"
#include "../source/Math/RECParticleKinematic.h"
#include "../source/Math/ParticleMassTable.h"

void Run101(const float beam_energy, const std::string& FilePath) {
    std::vector<std::string> inputFiles = GetHipoFilesInPath(FilePath);
    if(inputFiles.empty()) {
        std::cout << "No .hipo files found in directory!\n";
    }

    auto ds = std::make_unique<RHipoDS>(inputFiles);
    ROOT::RDF::RNode df = RDataFrame(std::move(ds));
    auto dfN = std::make_shared<ROOT::RDF::RNode>(df);
    auto dfNode = std::make_optional<ROOT::RDF::RNode>(*dfN);
    
    TFile* fout = new TFile("test.root", "RECREATE");

    EventCut Electron_cut;
    Electron_cut.SetChargeCut(-1);
    Electron_cut.SetPIDCountCut(11, 1, 1);
    EventCut Proton_cut;
    Proton_cut.SetChargeCut(1);
    Proton_cut.SetPIDCountCut(2212, 1, 1);

    EventCut Photon_cut;
    Photon_cut.SetChargeCut(0);
    Photon_cut.SetPIDCountCut(22, 1, 1);

    TrackCut trackCut1;
    trackCut1.SetECALEdgeCut(9, 100000000); // 设置边缘距离范围
    trackCut1.SetDCEdgeCut(1, 100000000); // 设置边缘距离范围


    //ROOT::RDF::RNode dfSelected = dfNode.value();
    ROOT::RDF::RNode dfSelected = dfNode.value()
    .Define("REC_Particle_num", [](const std::vector<int>& pid) {
        return static_cast<int>(pid.size()); // 显式转换为 int
    }, {"REC_Particle_pid"});

    dfSelected = dfSelected.Define("REC_Traj_pass", trackCut1.RECTrajPass(), 
                                CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    
    //if want to take a look at the edge of the DC or ECAL, uncomment the following lines
    /*
    dfSelected = dfSelected.Define("REC_Particle_edgeDC6", trackCut1.RECTrajedge(6, 6), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    dfSelected = dfSelected.Define("REC_Particle_edgeDC18", trackCut1.RECTrajedge(6, 18), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    dfSelected = dfSelected.Define("REC_Particle_edgeDC36", trackCut1.RECTrajedge(6, 36), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    dfSelected = dfSelected.Define("REC_Particle_edgeECAL1", trackCut1.RECTrajedge(7, 1), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    dfSelected = dfSelected.Define("REC_Particle_edgeECAL4", trackCut1.RECTrajedge(7, 4), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    dfSelected = dfSelected.Define("REC_Particle_edgeECAL7", trackCut1.RECTrajedge(7, 7), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    */

    dfSelected = dfSelected.Filter(Electron_cut, CombineColumns(RECParticle::All(), std::vector<std::string>{"REC_Traj_pass"}));
    dfSelected = dfSelected.Filter(Proton_cut, CombineColumns(RECParticle::All(), std::vector<std::string>{"REC_Traj_pass"}));
    dfSelected = dfSelected.Filter(Photon_cut, CombineColumns(RECParticle::All(), std::vector<std::string>{"REC_Traj_pass"}));

    dfSelected = dfSelected.Define("REC_Particle_theta", RECParticletheta(), RECParticle::All());
    dfSelected = dfSelected.Define("REC_Particle_phi", RECParticlephi(), RECParticle::All());
    dfSelected = dfSelected.Define("REC_Particle_p", RECParticleP(), RECParticle::All());

    dfSelected = dfSelected.Define("REC_Event_Q2", EventQ2(beam_energy,11,-1), RECParticle::All());
    dfSelected = dfSelected.Define("REC_Event_xB", EventxB(beam_energy,11,-1,getParticleMass(2212)), RECParticle::All());
    dfSelected = dfSelected.Define("REC_Event_Nu", EventNu(beam_energy,11,-1), RECParticle::All());
    dfSelected = dfSelected.Define("REC_Event_W", EventW(beam_energy,11,-1,getParticleMass(2212)), RECParticle::All());
    dfSelected = dfSelected.Define("REC_Event_mt", Eventmt(beam_energy,2212,1,getParticleMass(2212)), RECParticle::All());

    //auto disp = dfSelected.Display({ "REC_Traj_pass","REC_Particle_edgeDC6","REC_Particle_edgeDC18","REC_Particle_edgeDC36"},100,48);
    //disp->Print();

    std::cout << "df_selected count: " << *dfSelected.Count() << std::endl;
    
    dfSelected.Snapshot("dfSelected", "snapshot.root", CombineColumns(RECParticle::Extend(),
                                                                      RECTraj::All(),
                                                                      std::vector<std::string>{"REC_Particle_num"}
                                                                      //std::vector<std::string>{"REC_Particle_edgeDC6"},
                                                                      //std::vector<std::string>{"REC_Particle_edgeDC18"},
                                                                      //std::vector<std::string>{"REC_Particle_edgeDC36"},
                                                                      //std::vector<std::string>{"REC_Particle_edgeECAL1"},
                                                                      //std::vector<std::string>{"REC_Particle_edgeECAL4"},
                                                                      //std::vector<std::string>{"REC_Particle_edgeECAL7"}
                                                                      ));

    //example of how to draw histograms
    DrawAndSaveEventQ2(dfSelected, 11, -1, fout, 
                              500, 0, 5, ""); // Q2 histogram: bins, min, max
    DrawAndSaveEventxB(dfSelected, 11, -1, fout,
                              500, 0, 1.5, ""); // xB histogram: bins, min, max
    DrawAndSaveEventNu(dfSelected, 11, -1, fout,
                              500, 0, 6, ""); // Nu histogram: bins, min, max
    DrawAndSaveEventW(dfSelected, 11, -1, fout,
                              500, 0, 4, ""); // W histogram: bins, min, max
    DrawAndSaveEventmt(dfSelected, 2212, 1, fout,
                              500, 0, 6, ""); // t histogram: bins, min, max
    DrawAndSaveQ2vsxB(dfSelected, 11, -1, fout,
                              500, 0, 5,        // Q2 histogram: bins, min, max
                              500, 0, 1.5, "");     // xB histogram: bins, min, max

    DrawAndSaveParticleHistograms(dfSelected, 11, -1, fout, 
                              500, 0, 1,        // Theta histogram: bins, min, max
                              500, 0, 2 * M_PI, // Phi histogram: bins, min, max
                              500, 0, 9, "");       // Momentum histogram: bins, min, max
    DrawAndSaveParticleHistograms(dfSelected, 2212, 1, fout, 
                              500, 0, 2,        // Theta histogram: bins, min, max
                              500, 0, 2 * M_PI, // Phi histogram: bins, min, max
                              500, 0, 5, "");       // Momentum histogram: bins, min, max
    DrawAndSaveParticleHistograms(dfSelected, 22, 0, fout,
                              500, 0, 1,        // Theta histogram: bins, min, max
                              500, 0, 2 * M_PI, // Phi histogram: bins, min, max
                              500, 0, 9, "");       // Momentum histogram: bins, min, max
    

    fout->Close();
}
