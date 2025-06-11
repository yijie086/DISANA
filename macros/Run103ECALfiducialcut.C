#include <string>
#include <iostream>
#include "RHipoDS.hxx"
#include "TCanvas.h"
#include <optional>
#include <vector>

#include "../DreamAN/ParticleInformation/RECParticle.h"
#include "../DreamAN/ParticleInformation/RECTraj.h"
#include "../DreamAN/ParticleInformation/RECTrack.h"
#include "../DreamAN/ParticleInformation/RECCalorimeter.h"
#include "../DreamAN/Cuts/ElectronCut.h"
#include "../DreamAN/DrawHist/DrawAndSave.h"
#include "../DreamAN/core/FilesInPath.h"
#include "../DreamAN/core/Columns.h"
#include "../DreamAN/Cuts/EventCut.h"
#include "../DreamAN/Cuts/TrackCut.h"
#include "../DreamAN/Math/RECParticleKinematic.h"
#include "../DreamAN/Math/ParticleMassTable.h"

void Run103ECALfiducialcut(const float beam_energy, const std::string& FilePath) {
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

    //ROOT::RDF::RNode dfSelected = dfNode.value();
    ROOT::RDF::RNode dfSelected = dfNode.value()
    .Define("REC_Particle_num", [](const std::vector<int>& pid) {
        return static_cast<int>(pid.size()); // 显式转换为 int
    }, {"REC_Particle_pid"});

    dfSelected = dfSelected.Define("REC_Particle_PCALlu", RECCalorimeterluvw(7,1,1),
                                CombineColumns(RECCalorimeter::All(), std::vector<std::string>{"REC_Particle_num"}));
    dfSelected = dfSelected.Define("REC_Particle_PCALlv", RECCalorimeterluvw(7,1,2),
                                CombineColumns(RECCalorimeter::All(), std::vector<std::string>{"REC_Particle_num"}));
    dfSelected = dfSelected.Define("REC_Particle_PCALlw", RECCalorimeterluvw(7,1,3),
                                CombineColumns(RECCalorimeter::All(), std::vector<std::string>{"REC_Particle_num"}));
    
    dfSelected = dfSelected.Define("REC_Particle_ECINlu", RECCalorimeterluvw(7,4,1),
                                CombineColumns(RECCalorimeter::All(), std::vector<std::string>{"REC_Particle_num"}));
    dfSelected = dfSelected.Define("REC_Particle_ECINlv", RECCalorimeterluvw(7,4,2),
                                CombineColumns(RECCalorimeter::All(), std::vector<std::string>{"REC_Particle_num"}));
    dfSelected = dfSelected.Define("REC_Particle_ECINlw", RECCalorimeterluvw(7,4,3),
                                CombineColumns(RECCalorimeter::All(), std::vector<std::string>{"REC_Particle_num"}));
    
    dfSelected = dfSelected.Define("REC_Particle_ECOUTlu", RECCalorimeterluvw(7,7,1),
                                CombineColumns(RECCalorimeter::All(), std::vector<std::string>{"REC_Particle_num"}));
    dfSelected = dfSelected.Define("REC_Particle_ECOUTlv", RECCalorimeterluvw(7,7,2),
                                CombineColumns(RECCalorimeter::All(), std::vector<std::string>{"REC_Particle_num"}));
    dfSelected = dfSelected.Define("REC_Particle_ECOUTlw", RECCalorimeterluvw(7,7,3),
                                CombineColumns(RECCalorimeter::All(), std::vector<std::string>{"REC_Particle_num"}));

    

    dfSelected = dfSelected.Define("REC_Particle_theta", RECParticletheta(), RECParticle::All());
    dfSelected = dfSelected.Define("REC_Particle_phi", RECParticlephi(), RECParticle::All());
    dfSelected = dfSelected.Define("REC_Particle_p", RECParticleP(), RECParticle::All());
    auto dfSelected0 = dfSelected;

    //auto disp = dfSelected.Display({ "REC_Particle_pid","REC_Calorimeter_pass"},100,48);
    //disp->Print();

    for (int i=1; i<=6; i++){
        dfSelected = dfSelected0;
        trackCut1.SetSectorCut(i, 11, 7, true);
        dfSelected = dfSelected.Define("REC_Calorimeter_pass", trackCut1.RECCalorimeterPass(), 
                                CombineColumns(RECCalorimeter::All(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_num"}));
        dfSelected = dfSelected.Define("REC_Traj_pass", trackCut1.RECTrajPass(),
                                CombineColumns(RECTraj::All(), std::vector<std::string>{"REC_Particle_pid"}, std::vector<std::string>{"REC_Particle_num"}));
        dfSelected = dfSelected.Define("REC_Track_pass", Columns::LogicalAND2(), {"REC_Traj_pass", "REC_Calorimeter_pass"});
        dfSelected = dfSelected.Filter(Electron_cut, CombineColumns(RECParticle::All(), std::vector<std::string>{"REC_Track_pass"}));
        std::cout << "df_selected count: " << *dfSelected.Count() << std::endl;

        std::string outputname1 = "sector"+std::to_string(i)+"_PCALlw_REC_PCALlv";
        std::string outputname2 = "sector"+std::to_string(i)+"_ECINlw_REC_ECINlv";
        std::string outputname3 = "sector"+std::to_string(i)+"_ECOUTlu_REC_ECOUTlv";

        DrawAndSavelxvsly(dfSelected, 11, -1, fout, 
                  500, 0, 500, 500, 0, 500, 
                  outputname1, 
                  "REC_Particle_PCALlw", "REC_Particle_PCALlv");
        DrawAndSavelxvsly(dfSelected, 11, -1, fout, 
                  500, 0, 500, 500, 0, 500, 
                  outputname2, 
                  "REC_Particle_ECINlw", "REC_Particle_ECINlv");
        DrawAndSavelxvsly(dfSelected, 11, -1, fout, 
                  500, 0, 500, 500, 0, 500, 
                  outputname3, 
                  "REC_Particle_ECOUTlu", "REC_Particle_ECOUTlv");

    }

    fout->Close();
}
