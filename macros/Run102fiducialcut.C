#include <string>
#include <iostream>
#include "RHipoDS.hxx"
#include "TCanvas.h"
#include <optional>
#include <vector>

#include "../DreamAN/ParticleInformation/RECParticle.h"
#include "../DreamAN/ParticleInformation/RECTraj.h"
#include "../DreamAN/ParticleInformation/RECTrack.h"
#include "../DreamAN/Cuts/ElectronCut.h"
#include "../DreamAN/DrawHist/DrawAndSave.h"
#include "../DreamAN/core/FilesInPath.h"
#include "../DreamAN/core/Columns.h"
#include "../DreamAN/Cuts/EventCut.h"
#include "../DreamAN/Cuts/TrackCut.h"
#include "../DreamAN/Math/RECParticleKinematic.h"
#include "../DreamAN/Math/ParticleMassTable.h"

void Run102fiducialcut(const float beam_energy, const std::string& FilePath) {
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

    TrackCut trackCut1;
    //trackCut1.SetSectorCut(1, true); // 设置扇区范围
    //trackCut1.SetECALEdgeCut(9, 100000000);
    //trackCut1.SetDCEdgeCut(1, 100000000);

    int detector_investigate = 6; // 6 for for DC
    int layer_investigate1 = 6; // 6,18,36 for DC
    int layer_investigate2 = 18;
    int layer_investigate3 = 36;


    //ROOT::RDF::RNode dfSelected = dfNode.value();
    ROOT::RDF::RNode dfSelected = dfNode.value()
    .Define("REC_Particle_num", [](const std::vector<int>& pid) {
        return static_cast<int>(pid.size()); // 显式转换为 int
    }, {"REC_Particle_pid"});

    //dfSelected = dfSelected.Define("REC_Traj_pass", trackCut1.RECTrajPass(), 
    //                            CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    dfSelected = dfSelected.Define("REC_Traj_pedgeR1", RECTrajedge(detector_investigate,layer_investigate1),
                                CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    dfSelected = dfSelected.Define("REC_Traj_pedgeR2", RECTrajedge(detector_investigate,layer_investigate2),
                                CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    dfSelected = dfSelected.Define("REC_Traj_pedgeR3", RECTrajedge(detector_investigate,layer_investigate3),
                                CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    dfSelected = dfSelected.Define("REC_Track_chi2perndf", RECTrackchi2perndf(detector_investigate), 
                                CombineColumns(RECTrack::All(), std::vector<std::string>{"REC_Particle_num"}));
    //if want to take a look at the edge of the DC or ECAL, uncomment the following lines
    //dfSelected = dfSelected.Define("REC_Traj_ECAL_PCAL_X", RECTrajXYZ(7,1,1), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    //dfSelected = dfSelected.Define("REC_Traj_ECAL_PCAL_Y", RECTrajXYZ(7,1,2), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    //dfSelected = dfSelected.Define("REC_Traj_ECAL_PCAL_Z", RECTrajXYZ(7,1,3), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    

    //dfSelected = dfSelected.Define("REC_Traj_ECAL_ECIN_X", RECTrajXYZ(7,4,1), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    //dfSelected = dfSelected.Define("REC_Traj_ECAL_ECIN_Y", RECTrajXYZ(7,4,2), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    //dfSelected = dfSelected.Define("REC_Traj_ECAL_ECIN_Z", RECTrajXYZ(7,4,3), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    //dfSelected = dfSelected.Define("REC_Traj_ECAL_ECOUT_X", RECTrajXYZ(7,7,1), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    //dfSelected = dfSelected.Define("REC_Traj_ECAL_ECOUT_Y", RECTrajXYZ(7,7,2), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    //dfSelected = dfSelected.Define("REC_Traj_ECAL_ECOUT_Z", RECTrajXYZ(7,7,3), CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));

    dfSelected = dfSelected.Define("REC_Particle_theta", RECParticletheta(), RECParticle::All());
    dfSelected = dfSelected.Define("REC_Particle_phi", RECParticlephi(), RECParticle::All());
    dfSelected = dfSelected.Define("REC_Particle_p", RECParticleP(), RECParticle::All());
    auto dfSelected0 = dfSelected;

    //auto disp = dfSelected.Display({ "REC_Particle_pid","REC_Track_chi2perndf"},100,48);
    //disp->Print();


    for (int i=0; i<5; i++){
            for (int j=1; j<=6; j++){
            dfSelected = dfSelected0;
            int theta_min = (5+i*5); // 5, 10, 15, 20, 25, 30 degrees
            int theta_max = (10+i*5); // 10, 15, 20, 25, 30, 35 degrees
            if (i==4) {
                theta_max = 40;
            }
            std::string outputname = "theta"+std::to_string(theta_min)+"to"+std::to_string(theta_max)+"_sector"+std::to_string(j);
            std::string outputname1 = "theta"+std::to_string(theta_min)+"to"+std::to_string(theta_max)+"_sector"+std::to_string(j)+"_R1";
            std::string outputname2 = "theta"+std::to_string(theta_min)+"to"+std::to_string(theta_max)+"_sector"+std::to_string(j)+"_R2";
            std::string outputname3 = "theta"+std::to_string(theta_min)+"to"+std::to_string(theta_max)+"_sector"+std::to_string(j)+"_R3";
            

            trackCut1.SetSectorCut(j, true); // 设置扇区范围
            dfSelected = dfSelected.Define("REC_Traj_pass", trackCut1.RECTrajPass(), 
                                CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_num"}));
    
            Electron_cut.SetthetaCut(theta_min*M_PI/180,theta_max*M_PI/180);
            dfSelected = dfSelected.Filter(Electron_cut, CombineColumns(RECParticle::All(), std::vector<std::string>{"REC_Traj_pass"}));
            std::cout << "df_selected count: " << *dfSelected.Count() << std::endl;

            //DrawAndSaveedge(dfSelected, detector_investigate, layer_investigate1, 11, -1, fout, 
            //                        25, 0, 25, {"REC_Traj_pedgeR1"}, outputname1); // Edge histogram: bins, min, max
            DrawAndSavechi2perndfvsedge(dfSelected, detector_investigate, layer_investigate1, 11, -1, fout,
                                    5000, 0, 500,       // Chi2 histogram: bins, min, max
                                    25, 0, 25, {"REC_Traj_pedgeR1"}, outputname1);     // Edge histogram: bins, min, max
            //DrawAndSaveedge(dfSelected, detector_investigate, layer_investigate2, 11, -1, fout,
            //                        25, 0, 25, {"REC_Traj_pedgeR2"}, outputname2); // Edge histogram: bins, min, max
            DrawAndSavechi2perndfvsedge(dfSelected, detector_investigate, layer_investigate2, 11, -1, fout,
                                    5000, 0, 500,       // Chi2 histogram: bins, min, max
                                    25, 0, 25, {"REC_Traj_pedgeR2"}, outputname2);     // Edge histogram: bins, min, max
            //DrawAndSaveedge(dfSelected, detector_investigate, layer_investigate3, 11, -1, fout,
            //                        25, 0, 25, {"REC_Traj_pedgeR3"}, outputname3); // Edge histogram: bins, min, max
            DrawAndSavechi2perndfvsedge(dfSelected, detector_investigate, layer_investigate3, 11, -1, fout,
                                    5000, 0, 500,       // Chi2 histogram: bins, min, max
                                    25, 0, 25, {"REC_Traj_pedgeR3"}, outputname3);     // Edge histogram: bins, min, max
            //DrawAndSavechi2perndf(dfSelected, detector_investigate, layer_investigate1, 11, -1, fout,
            //                        5000, 0, 500, outputname1); // Chi2 histogram: bins, min, max
            //DrawAndSavechi2perndf(dfSelected, detector_investigate, layer_investigate2, 11, -1, fout,
            //                        5000, 0, 500, outputname2); // Chi2 histogram: bins, min, max
            //DrawAndSavechi2perndf(dfSelected, detector_investigate, layer_investigate3, 11, -1, fout,
            //                        5000, 0, 500, outputname3); // Chi2 histogram: bins, min, max

            //DrawAndSaveParticleHistograms(dfSelected, 11, -1, fout, 
            //                        500, 0, 1,        // Theta histogram: bins, min, max
            //                        500, 0, 2 * M_PI, // Phi histogram: bins, min, max
            //                        500, 0, 9,outputname);       // Momentum histogram: bins, min, max
        }
    }

    fout->Close();
}
