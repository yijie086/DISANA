// usage:
// root -l -q 'DrawStatBoostEffectsForGEMC.C("boost_map_Q2_t_rgasp19_inb.txt", "/path/to/slim_sp19_MC_exclusive_qadb.root","slim_sp19_MC_exclusive_qadb","Q2","t","boost_effect")'
// eg. root -l -q './../../../../source/rootmacros/DrawStatBoostEffectsForGEMC.C("boost_map_Q2_t_rgasp19_inb.txt", "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_results/boostmaps/slim_sp19_MC_exclusive_qadb.root","slim_sp19_MC_exclusive_qadb","Q2","t","boost_effect")'
//root -l -q './../../../../source/rootmacros/DrawStatBoostEffectsForGEMC.C("boost_map_Q2_t_rgasp19_inb.txt", "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_results/boostmaps/slim_sp19_MC_exclusive_qadb.root","slim_sp19_MC_exclusive_qadb","Q2","t","boost_effect")'
//root -l -q './../../DrawBoostEffectsFromLund.C(  "./../before_boost/lager_run/CLAS-ep-phi-10GeV.ep-phi.4pi.run00001-5000_clean.txt","./../after_boost/lager_run/CLAS-ep-phi-10GeV.ep-phi.4pi.run00001-5000_clean.txt","boost_test_15",10.2,15, 0.05, 8.0,15, -8.0, 0.0)'

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TPaletteAxis.h>

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <random>

// ============================================================
struct BoostBin {
  double q2lo, q2hi;
  double tlo, thi;
  double w;
};

static std::vector<BoostBin> gBins;

// ------------------------------------------------------------
static void LoadBoostMap(const std::string& txt) {
  gBins.clear();
  std::ifstream in(txt);
  if(!in){ std::cerr << "Cannot open " << txt << "\n"; return; }
  std::string line;
  while(std::getline(in, line)){
    if(line.empty() || line[0]=='#') continue;
    for(char& c : line) if(c==',') c=' ';
    BoostBin b;
    std::istringstream ss(line);
    if(ss >> b.q2lo >> b.q2hi >> b.tlo >> b.thi >> b.w)
      gBins.push_back(b);
  }
  std::cout << "[INFO] Loaded " << gBins.size() << " boost bins\n";
}

// ------------------------------------------------------------
// Lookup boost weight for the (Q2,t) point.
// Convention:
//   - If point is OUTSIDE the map coverage: return 1.0 (keep event unboosted)
//   - If point is INSIDE a masked bin: map stores w=0, and we return 0
static double BoostLUT(double Q2, double tUsed){
  for(const auto& b : gBins)
    if(Q2 >= b.q2lo && Q2 < b.q2hi && tUsed >= b.tlo && tUsed < b.thi)
      return b.w;
  return 1.0;
}

// ------------------------------------------------------------
// Build unique sorted bin edges from the map
static void GetEdges(std::vector<double>& q2E, std::vector<double>& tE){
  q2E.clear(); tE.clear();
  for(const auto& b : gBins){
    q2E.push_back(b.q2lo); q2E.push_back(b.q2hi);
    tE.push_back(b.tlo);   tE.push_back(b.thi);
  }
  auto uniq = [](std::vector<double>& v){
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end(),
            [](double a, double b){ return std::fabs(a-b)<1e-12; }), v.end());
  };
  uniq(q2E); uniq(tE);
}

// ============================================================
void DrawStatBoostEffectsForGEMC(const char* boostTxt,
                                  const char* rootFile,
                                  const char* treeName,
                                  const char* xVar      = "Q2",
                                  const char* tVar      = "t",
                                  const char* outPrefix = "boost_effect")
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);

  // ---- load map ----
  LoadBoostMap(boostTxt);
  if(gBins.empty()){ std::cerr << "[ERROR] Boost map empty\n"; return; }

  std::vector<double> q2Edges, tEdges;
  GetEdges(q2Edges, tEdges);
  if(q2Edges.size()<2 || tEdges.size()<2){
    std::cerr << "[ERROR] Not enough edges\n"; return;
  }

  // Auto-detect whether the boost-map t-axis is negative (Mandelstam t) or positive (-t).
  const double tEdgeMin = *std::min_element(tEdges.begin(), tEdges.end());
  const double tEdgeMax = *std::max_element(tEdges.begin(), tEdges.end());
  const bool mapIsNegativeT = (tEdgeMax <= 0.0);
  const char* tLabel = mapIsNegativeT ? "t [GeV^{2}]" : "-t [GeV^{2}]";

  // ---- open ROOT file ----
  TFile* f = TFile::Open(rootFile, "READ");
  if(!f || f->IsZombie()){
    std::cerr << "[ERROR] Cannot open " << rootFile << "\n"; return;
  }
  TTree* tr = (TTree*)f->Get(treeName);
  if(!tr){
    std::cerr << "[ERROR] Cannot find tree " << treeName << "\n";
    f->Close(); return;
  }

  double Q2 = 0, tval = 0;
  tr->SetBranchStatus("*", 0);
  tr->SetBranchStatus(xVar, 1);
  tr->SetBranchStatus(tVar, 1);
  tr->SetBranchAddress(xVar, &Q2);
  tr->SetBranchAddress(tVar, &tval);

  // ---- four histograms on the same grid ----
  int nQ2 = (int)q2Edges.size()-1;
  int nT  = (int)tEdges.size()-1;

  // 1) Before  – raw counts per bin
  TH2D* hBefore = new TH2D("hBefore",
    Form("Before (unweighted);Q^{2} [GeV^{2}];%s", tLabel),
    nQ2, q2Edges.data(), nT, tEdges.data());

  // 2) Weights – the boost map value in each bin (just the LUT visualised)
  TH2D* hWeights = new TH2D("hWeights",
    Form("Boost weight map;Q^{2} [GeV^{2}];%s", tLabel),
    nQ2, q2Edges.data(), nT, tEdges.data());

  // 3) After   – SIMULATED post-duplication counts
  //    Each event contributes floor(w) copies + 1 with prob frac(w)
  //    => mean contribution = w  => hAfter ~ hBefore * w  (bin-by-bin)
  TH2D* hAfter = new TH2D("hAfter",
    Form("After (simulated duplication);Q^{2} [GeV^{2}];%s", tLabel),
    nQ2, q2Edges.data(), nT, tEdges.data());

  // 4) Ratio   – After / Before  (should be ~w everywhere)
  TH2D* hRatio = new TH2D("hRatio",
    Form("Ratio After/Before;Q^{2} [GeV^{2}];%s", tLabel),
    nQ2, q2Edges.data(), nT, tEdges.data());

  // ---- fill weights map directly from LUT (independent of tree) ----
  for(int ix=1; ix<=nQ2; ix++){
    double q2c = 0.5*(q2Edges[ix-1]+q2Edges[ix]);
    for(int iy=1; iy<=nT; iy++){
      double tc = 0.5*(tEdges[iy-1]+tEdges[iy]);
      double ww = BoostLUT(q2c, tc);
      hWeights->SetBinContent(ix, iy, ww);
    }
  }

  // ---- RNG for fractional duplication ----
  std::mt19937_64 rng(42);
  std::uniform_real_distribution<double> U(0.0, 1.0);

  // ---- range limits ----
  const double q2Min = q2Edges.front(), q2Max = q2Edges.back();
  const double tMin  = tEdges.front(),  tMax  = tEdges.back();

  long long nRead=0, nInRange=0;
  double minQ2=1e9, maxQ2=-1e9, minMinusT=1e9, maxMinusT=-1e9;

  const long long N = tr->GetEntries();
  for(long long i=0; i<N; i++){
    tr->GetEntry(i);
    nRead++;

    // Convert the tree's t branch to match the boost-map convention.
    // - If map is negative (Mandelstam t): ensure tUsed <= 0.
    // - If map is positive (-t): ensure tUsed >= 0.
    double tUsed = tval;
    if(mapIsNegativeT){
      if(tUsed > 0) tUsed = -tUsed; // common case: stored +(-t) -> make it t
    } else {
      if(tUsed < 0) tUsed = -tUsed; // stored t -> make it -t
    }

    minQ2    = std::min(minQ2,    Q2);
    maxQ2    = std::max(maxQ2,    Q2);
    minMinusT = std::min(minMinusT, tUsed);
    maxMinusT = std::max(maxMinusT, tUsed);

    const bool in = (Q2 >= q2Min && Q2 < q2Max &&
                     tUsed >= tMin && tUsed < tMax);
    if(!in) continue;
    nInRange++;

    const double w = BoostLUT(Q2, tUsed);

    // hBefore: one count regardless of weight
    hBefore->Fill(Q2, tUsed, 1.0);

    // hAfter: simulate event duplication
    //   integer part: always written
    //   fractional part: written with probability = frac
    int nCopies = (int)w;
    double frac = w - nCopies;
    if(U(rng) < frac) nCopies++;
    hAfter->Fill(Q2, tUsed, (double)nCopies);
  }

  // ---- compute ratio ----
  for(int ix=1; ix<=nQ2; ix++){
    for(int iy=1; iy<=nT; iy++){
      double before = hBefore->GetBinContent(ix, iy);
      double after  = hAfter ->GetBinContent(ix, iy);
      hRatio->SetBinContent(ix, iy, before>0 ? after/before : 0.0);
    }
  }

  // ---- diagnostics ----
  std::cout << "[INFO] Entries read   : " << nRead    << "\n";
  std::cout << "[INFO] In-range       : " << nInRange << "\n";
  std::cout << "[INFO] Q2  min/max    : " << minQ2    << " / " << maxQ2    << "\n";
  std::cout << "[INFO] " << (mapIsNegativeT ? "t" : "-t") << "  min/max    : "
            << minMinusT << " / " << maxMinusT << "\n";
  std::cout << "[INFO] hBefore integral = " << hBefore->Integral() << "\n";
  std::cout << "[INFO] hAfter  integral = " << hAfter ->Integral() << "\n";
  std::cout << "[INFO] Expected ratio   = " << hAfter->Integral()/hBefore->Integral() << "\n";

  // Occupancy report: tells you if binning is appropriate for the MC statistics
  int nEmpty=0, nOne=0, nFew=0, nGood=0;
  double sumOcc=0; int nFilled=0;
  for(int ix=1; ix<=nQ2; ix++){
    for(int iy=1; iy<=nT; iy++){
      double bc = hBefore->GetBinContent(ix,iy);
      double ww = hWeights->GetBinContent(ix,iy);
      if(ww <= 0) continue;  // skip masked bins
      if(bc == 0)      nEmpty++;
      else if(bc < 2)  nOne++;
      else if(bc < 5)  nFew++;
      else             nGood++;
      if(bc > 0){ sumOcc += bc; nFilled++; }
    }
  }
  double meanOcc = nFilled>0 ? sumOcc/nFilled : 0;
  std::cout << "[OCCUPANCY] Mean events/bin (boosted region): " << meanOcc << "\n";
  std::cout << "[OCCUPANCY] Bins with 0 events : " << nEmpty << " (these stay empty after boost!)\n";
  std::cout << "[OCCUPANCY] Bins with 1 event  : " << nOne   << " (boosted but isolated)\n";
  std::cout << "[OCCUPANCY] Bins with 2-4 events: " << nFew  << " (marginal)\n";
  std::cout << "[OCCUPANCY] Bins with 5+ events : " << nGood << " (good for flattening)\n";
  if(meanOcc < 5)
    std::cout << "[WARNING] Mean occupancy < 5. Binning is TOO FINE for your MC statistics.\n"
              << "          Reduce nBins in ExportBoostMap_Q2_minusT until mean >= 5.\n"
              << "          Suggested nBins ~ sqrt(" << (int)hBefore->Integral() << "/5) = "
              << (int)std::sqrt(hBefore->Integral()/5) << " per axis.\n";
  else
    std::cout << "[OK] Binning is appropriate for your MC statistics.\n";

  // ---- colour axis ranges ----
  // find max of Before so After colour scale is comparable
  double beforeMax = hBefore->GetMaximum();
  double weightMax = hWeights->GetMaximum();
  double afterMax  = hAfter ->GetMaximum();
  double ratioMax  = hRatio ->GetMaximum();

  hBefore ->GetZaxis()->SetRangeUser(0, beforeMax);
  hWeights->GetZaxis()->SetRangeUser(0, weightMax);
  hAfter  ->GetZaxis()->SetRangeUser(0, afterMax);
  hRatio  ->GetZaxis()->SetRangeUser(0, std::ceil(ratioMax));

  // ---- draw ----
  TCanvas* c = new TCanvas("c", "Boost Effect", 2200, 550);
  c->Divide(4, 1, 0.002, 0.01);

  auto stylepad = [](TVirtualPad* p){
    p->SetRightMargin(0.15);
    p->SetLeftMargin(0.10);
    p->SetBottomMargin(0.12);
    p->SetTopMargin(0.10);
  };

  c->cd(1); stylepad(gPad); hBefore ->Draw("COLZ");
  c->cd(2); stylepad(gPad); hWeights->Draw("COLZ");
  c->cd(3); stylepad(gPad); hAfter  ->Draw("COLZ");
  c->cd(4); stylepad(gPad); hRatio  ->Draw("COLZ");

  c->Update();

  c->SaveAs(Form("%s.pdf",  outPrefix));
  c->SaveAs(Form("%s.png",  outPrefix));
  std::cout << "[INFO] Saved " << outPrefix << ".pdf/.png\n";

  f->Close();
}