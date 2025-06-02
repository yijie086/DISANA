#include <string>
#include <iostream>
#include <chrono>
#include "RHipoDS.hxx"
#include "TCanvas.h"
#include <optional>

#include "../DreamAN/ParticleInformation/RECParticle.h"
#include "../DreamAN/Cuts/ElectronCut.h"
#include "../DreamAN/DrawHist/DrawAndSave.h"
#include "../DreamAN/core/FilesInPath.h"
#include "../DreamAN/Cuts/EventCut.h"

using namespace ROOT;
using namespace ROOT::RDF;

auto func0 = [](const std::vector<float>& var0, const std::vector<float>& var1) {
    std::vector<float> result;
    size_t n = std::min(var0.size(), var1.size());
    result.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        result.push_back(var0[i] + var1[i]);
    }
    return result;
};


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
    ROOT::RDF::RNode df = RDataFrame(std::move(ds));

    auto dfNode = std::make_shared<ROOT::RDF::RNode>(df);
    auto dfSelected = std::make_optional<ROOT::RDF::RNode>(*dfNode);

    // example for using a lambda function to define a new column
    //if (dfSelected) {
    //    dfSelected = dfSelected->Define("new_col", func0, {"REC_Particle_px", "REC_Particle_py"});
    //}
    
    TFile* fout = new TFile("test.root", "RECREATE");

    EventCut Electron_cut;
    Electron_cut.SetChargeCut(-1,true);
    Electron_cut.SetPIDCountCut(11, 1, 1);

    ROOT::RDF::RNode dfRECSelected = dfSelected.value().Filter(Electron_cut, RECParticle::All());

    std::cout << "df_selected count: " << *dfRECSelected.Count() << std::endl;

    //example of how to draw a 1D histogram
    auto h_electron_pz = dfRECSelected
        .Define("electron_pz", get_RECParticle_float_var(11,-1), {"REC_Particle_phi", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Histo1D({"h_electron_pz", "Electron Pz", 500, 0, 9}, "electron_pz");
    Draw1DHist(h_electron_pz.GetPtr(), "h_electron_pz", 500, 0, 9);
    Save1DHist(h_electron_pz.GetPtr(), "h_electron_pz", fout);

    fout->Close();

}