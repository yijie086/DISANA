#ifndef DISANA_PLOTTER_H
#define DISANA_PLOTTER_H

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include <ROOT/RDataFrame.hxx>
#include <memory>
#include <string>
#include <vector>

#include "DISANAMath.h"

class DISANAplotter {
 public:
DISANAplotter(ROOT::RDF::RNode df, double beamEnergy) : rdf(df), beam_energy(beamEnergy) {}

  ~DISANAplotter() {
    disHistos.clear();
    kinematicHistos.clear();
  }
   
ROOT::RDF::RNode GetRDF() { return rdf; }


  void GenerateKinematicHistos(const std::string& type) {
    // std::cout << " type is " <<type << std::endl;
    std::vector<std::string> vars = {"p", "theta", "phi"};
    for (const auto& v : vars) {
      std::string base = "rec" + type + "_" + v;
      //std::cout << " base is " << base << std::endl;
      double histMin, histMax;
      auto it = kinematicAxisRanges.find(base);
      if (it != kinematicAxisRanges.end()) {
        histMin = it->second.first;
        histMax = it->second.second;
      } else {
        auto minVal = rdf.Min(base);
        auto maxVal = rdf.Max(base);
        double lo = *minVal;
        double hi = *maxVal;
        double margin = std::max(1e-3, 0.05 * (hi - lo));
        histMin = lo - margin;
        histMax = hi + margin;
        kinematicAxisRanges[base] = {histMin, histMax};  // store it
      }
    
       // auto h_all = rdf.Histo1D({(base + "_all").c_str(), "", 100, histMin, histMax}, base);
    auto h_acc = rdf.Histo1D({(base).c_str(), "", 100, histMin, histMax}, base);
    kinematicHistos.push_back(h_acc);
    }

  std::map<std::string, std::pair<double, double>> cachedRanges;

  for (const auto& var : disvars) {
    double histMin = axisRanges[var].first;
    double histMax = axisRanges[var].second;
    // auto h_all = rdf.Histo1D({(var + "_all").c_str(), (var + " All").c_str(), 100, histMin, histMax}, var);
    auto h_acc = rdf.Histo1D({var.c_str(), var.c_str(), 100, histMin, histMax}, var);
    disHistos.push_back(h_acc);

    //auto h_acceptance = dynamic_cast<TH1*>(h_acc->Clone((var + "_acceptance").c_str()));
    //h_acceptance->SetTitle(("Acceptance for " + var).c_str());
    //acceptHistos.push_back(std::shared_ptr<TH1>(h_acceptance));
  }
  }

  std::vector<TH1*> GetDISHistograms() {
    std::vector<TH1*> allDIShisto;
    for (auto& h : disHistos) allDIShisto.push_back(h.GetPtr());
    return allDIShisto;
  }

  std::vector<TH1*> GetAllHistograms() {
    std::vector<TH1*> all;
    for (auto& h : kinematicHistos) all.push_back(h.GetPtr());
    for (auto& h : disHistos) all.push_back(h.GetPtr());
    return all;
  }

  std::vector<std::vector<std::vector<TH1D*>>> ComputeDVCS_CrossSection(const BinManager& bins, double luminosity) {
    return kinCalc.ComputeDVCS_CrossSection(rdf, bins, luminosity);
  }

  /// BSA computations // we need to have refined version of this codes
  std::vector<std::vector<std::vector<TH1D*>>> ComputeBSA(const BinManager& bins, double luminosity, double pol = 1.0) {
    // 1) Helicity selection
    auto rdf_pos = rdf.Filter("REC_Event_helicity ==  1");
    auto rdf_neg = rdf.Filter("REC_Event_helicity == -1");
    // 2) One-pass cross-section per helicity

    auto sigma_pos_3d = kinCalc.ComputeDVCS_CrossSection(rdf_pos, bins, luminosity);
    auto sigma_neg_3d = kinCalc.ComputeDVCS_CrossSection(rdf_neg, bins, luminosity);
    // 3) 3-D BSA
    //auto bsa_3d = kinCalc.ComputeBeamSpinAsymmetry(sigma_pos_3d, sigma_neg_3d, pol);
    // 4) Flatten to keep legacy return-type
    return kinCalc.ComputeBeamSpinAsymmetry(sigma_pos_3d, sigma_neg_3d, pol);

  }

   /// BSA computations // we need to have refined version of this codes
  std::vector<std::vector<std::vector<TH1D*>>> ComputePi0Corr(const BinManager& bins) {
    // 4) Flatten to keep legacy return-type
    return kinCalc.computePi0Corr(bins);

  }
  
  void ApplyPi0BkgCorr(THnSparseD *correctionHist ){
      kinCalc.SetApplyCorrPi0BKG(true);
      if (correctionHist) {
        kinCalc.SetCorrHist(correctionHist);
      } else {
        std::cerr << "Error: Correction histogram is not set.\n";
      }
  }

 private:

  std::vector<std::string> disvars = {"Q2", "xB", "t", "W", "phi"};
  std::map<std::string, std::pair<double, double>> axisRanges = {{"Q2", {0.0, 15.0}}, {"xB", {0.0, 1.0}}, {"W", {1.0, 10.0}}, {"t", {0.0, 10.0}}, {"phi", {-180.0, 180.0}}};
  std::map<std::string, std::pair<double, double>> kinematicAxisRanges = {
      {"recel_p", {-0.05, 13.0}},    {"recel_theta", {-0.01, 2.0}}, {"recel_phi", {-0.01, 6.1}}, {"recpho_p", {-0.01, 5.0}},
      {"recpho_theta", {-0.01, 2.0}}, {"recpho_phi", {-0.01, 6.1}},   {"recpro_p", {-0.01, 12.0}},  {"recpro_theta", {-0.01, 2.0}},
      {"recpro_phi", {-0.01, 6.1}}
      // Add more as needed
  };

  DISANAMath kinCalc;
  std::unique_ptr<TFile> predFile;
  std::string predFileName = ".";
  double beam_energy;
  std::string ttreeName;
  std::vector<ROOT::RDF::RResultPtr<TH1>> kinematicHistos, disHistos;
  std::vector<std::shared_ptr<TH1>> acceptHistos;
  ROOT::RDF::RNode rdf;
};

#endif  // DISANA_PLOTTER_H