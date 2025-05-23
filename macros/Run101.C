#include <string>
#include <iostream>
#include "RHipoDS.hxx"
#include "TCanvas.h"
#include <optional>

#include "../source/ParticleInformation/RECParticle.h"
#include "../source/Cuts/ElectronCut.h"
#include "../source/DrawHist/DrawAndSave.h"
#include "../source/core/FilesInPath.h"
#include "../source/Cuts/EventCut.h"
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
    Electron_cut.SetChargeCut(-1,true);
    Electron_cut.SetPIDCountCut(11, 1, 1);

    EventCut Proton_cut;
    Proton_cut.SetChargeCut(1,true);
    Proton_cut.SetPIDCountCut(2212, 1, 1);

    EventCut Photon_cut;
    Photon_cut.SetChargeCut(0,true);
    Photon_cut.SetPIDCountCut(22, 1, 1);

    ROOT::RDF::RNode dfSelected = dfNode.value().Filter(Electron_cut, RECParticle::All());
    dfSelected = dfSelected.Filter(Proton_cut, RECParticle::All());
    dfSelected = dfSelected.Filter(Photon_cut, RECParticle::All());

    dfSelected = dfSelected.Define("REC_Particle_theta", RECParticletheta(), RECParticle::All());
    dfSelected = dfSelected.Define("REC_Particle_phi", RECParticlephi(), RECParticle::All());
    dfSelected = dfSelected.Define("REC_Particle_p", RECParticleP(), RECParticle::All());

    dfSelected = dfSelected.Define("REC_Event_Q2", EventQ2(beam_energy,11,-1), RECParticle::All());
    dfSelected = dfSelected.Define("REC_Event_xB", EventxB(beam_energy,11,-1,getParticleMass(2212)), RECParticle::All());
    dfSelected = dfSelected.Define("REC_Event_Nu", EventNu(beam_energy,11,-1), RECParticle::All());
    dfSelected = dfSelected.Define("REC_Event_W", EventW(beam_energy,11,-1,getParticleMass(2212)), RECParticle::All());
    dfSelected = dfSelected.Define("REC_Event_mt", Eventmt(beam_energy,2212,1,getParticleMass(2212)), RECParticle::All());

    std::cout << "df_selected count: " << *dfSelected.Count() << std::endl;
    
    dfSelected.Snapshot("dfSelected", "snapshot.root", RECParticle::Extend());

    //example of how to draw a 1D histogram

    DrawAndSaveEventQ2(dfSelected, 11, -1, fout, 
                              500, 0, 5); // Q2 histogram: bins, min, max
    DrawAndSaveEventxB(dfSelected, 11, -1, fout,
                              500, 0, 1.5); // xB histogram: bins, min, max
    DrawAndSaveEventNu(dfSelected, 11, -1, fout,
                              500, 0, 6); // Nu histogram: bins, min, max
    DrawAndSaveEventW(dfSelected, 11, -1, fout,
                              500, 0, 4); // W histogram: bins, min, max
    DrawAndSaveEventmt(dfSelected, 2212, 1, fout,
                              500, 0, 6); // t histogram: bins, min, max
    DrawAndSaveQ2vsxB(dfSelected, 11, -1, fout,
                              500, 0, 5,        // Q2 histogram: bins, min, max
                              500, 0, 1.5);     // xB histogram: bins, min, max

    DrawAndSaveParticleHistograms(dfSelected, 11, -1, fout, 
                              500, 0, 1,        // Theta histogram: bins, min, max
                              500, 0, 2 * M_PI, // Phi histogram: bins, min, max
                              500, 0, 9);       // Momentum histogram: bins, min, max
    DrawAndSaveParticleHistograms(dfSelected, 2212, 1, fout, 
                              500, 0, 2,        // Theta histogram: bins, min, max
                              500, 0, 2 * M_PI, // Phi histogram: bins, min, max
                              500, 0, 5);       // Momentum histogram: bins, min, max
    DrawAndSaveParticleHistograms(dfSelected, 22, 0, fout,
                              500, 0, 1,        // Theta histogram: bins, min, max
                              500, 0, 2 * M_PI, // Phi histogram: bins, min, max
                              500, 0, 9);       // Momentum histogram: bins, min, max
    

    fout->Close();
}
