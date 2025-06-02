#include <string>
#include <iostream>
#include <chrono>
#include "RHipoDS.hxx"
#include "TCanvas.h"

#include "../DreamAN/ParticleInformation/RECParticle.h"
#include "../DreamAN/Cuts/ElectronCut.h"
#include "../DreamAN/DrawHist/DrawAndSave.h"
#include "../DreamAN/core/FilesInPath.h"

using namespace ROOT;
using namespace ROOT::RDF;


int main(int argc, char **argv) {
   //ROOT::EnableImplicitMT(2);

    if(argc < 2){
        std::cout << "Please specify a HIPO data file on the command line. (Only one file.) \n";
        return 1;
    }else{
        std::cout << "Opening file " << argv[1] << std::endl;
    }
    std::vector<std::string> inputFiles = GetHipoFilesInPath(argv[1]);
    if(inputFiles.empty()) {
        std::cout << "No .hipo files found in directory!\n";
        return 1;
    }

    auto ds = std::make_unique<RHipoDS>(inputFiles);
    auto df = RDataFrame(std::move(ds));

    TFile* fout = new TFile("test.root", "RECREATE");

    auto df_selected = df.Filter(ElectronCut, RECParticle::All());

    std::cout << "df_selected count: " << *df_selected.Count() << std::endl;

    auto h_electron_px = df_selected
        .Define("electron_px", get_RECParticle_float_var(11,-1), {"REC_Particle_px", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Histo1D({"h_electron_px", "Electron Px", 500, -3, 3}, "electron_px");
    Draw1DHist(h_electron_px.GetPtr(), "h_electron_px", 500, -3, 3);
    Save1DHist(h_electron_px.GetPtr(), "h_electron_px", fout);

    auto h_electron_py = df_selected
        .Define("electron_py", get_RECParticle_float_var(11,-1), {"REC_Particle_py", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Histo1D({"h_electron_py", "Electron Py", 500, -3, 3}, "electron_py");
    Draw1DHist(h_electron_py.GetPtr(), "h_electron_py", 500, -3, 3);
    Save1DHist(h_electron_py.GetPtr(), "h_electron_py", fout);

    auto h_electron_pz = df_selected
        .Define("electron_pz", get_RECParticle_float_var(11,-1), {"REC_Particle_pz", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Histo1D({"h_electron_pz", "Electron Pz", 500, 0, 9}, "electron_pz");
    Draw1DHist(h_electron_pz.GetPtr(), "h_electron_pz", 500, 0, 9);
    Save1DHist(h_electron_pz.GetPtr(), "h_electron_pz", fout);

    fout->Close();

}