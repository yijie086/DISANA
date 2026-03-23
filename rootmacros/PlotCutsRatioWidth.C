// PlotCutsRatioWidth_GlobalY.C
// Plot ratio: (hi-lo)_file1 / (hi-lo)_file2 versus -t bin center
// Requirement: tEdges of two files MUST be identical (t bin 一致)
// Feature: global Y range shared by ALL pads (scan first)
//
// Usage:
//   root -l -b -q 'PlotCutsRatioWidth_GlobalY.C("cuts1.csv","cuts2.csv")'
//
// Output: cuts_cfg0_widthRatio.pdf/png, cuts_cfg1_widthRatio.pdf/png, cuts_cfg2_widthRatio.pdf/png

#include <TCanvas.h>
#include <TPad.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLatex.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include <limits>

static std::string Trim(const std::string& s) {
  size_t b = 0, e = s.size();
  while (b < e && std::isspace((unsigned char)s[b])) ++b;
  while (e > b && std::isspace((unsigned char)s[e - 1])) --e;
  return s.substr(b, e - b);
}

static std::vector<std::string> SplitCSV(const std::string& line) {
  std::vector<std::string> out;
  std::string tok;
  std::stringstream ss(line);
  while (std::getline(ss, tok, ',')) out.push_back(Trim(tok));
  return out;
}

struct CutBin {
  bool has = false;
  double lo = 0.0;
  double hi = 0.0;
};

struct CutFileData {
  std::vector<double> tEdges;
  // key: var -> cfg -> vector bins
  std::map<std::string, std::map<int, std::vector<CutBin>>> cuts;
};

static bool IsSentinel999(double lo, double hi) {
  // treat lo=-999 and hi=999 as "do not plot"
  return (lo <= -998.5 && hi >= 998.5);
}

static CutFileData ReadCutFile(const std::string& fname) {
  CutFileData data;
  if (fname.empty()) return data;

  std::ifstream in(fname);
  if (!in) throw std::runtime_error("Cannot open file: " + fname);

  std::string line;
  bool gotEdges = false;

  while (std::getline(in, line)) {
    line = Trim(line);
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    // tEdges line
    if (!gotEdges) {
      auto toks = SplitCSV(line);
      if (!toks.empty() && toks[0] == "tEdges") {
        if (toks.size() < 3) throw std::runtime_error("Invalid tEdges line in: " + fname);
        data.tEdges.clear();
        for (size_t i = 1; i < toks.size(); i++) data.tEdges.push_back(std::stod(toks[i]));
        gotEdges = true;
        continue;
      }
    }

    auto toks = SplitCSV(line);
    if (toks.size() != 5) continue;

    if (!gotEdges || data.tEdges.size() < 2) {
      throw std::runtime_error("Missing/invalid tEdges line before cut rows in: " + fname);
    }

    const std::string var = toks[0];
    const int tbin = std::stoi(toks[1]);
    const int cfg  = std::stoi(toks[2]);
    const double lo = std::stod(toks[3]);
    const double hi = std::stod(toks[4]);

    const int nbin = (int)data.tEdges.size() - 1;

    auto& vec = data.cuts[var][cfg];
    if ((int)vec.size() != nbin) vec.assign(nbin, CutBin{});

    if (tbin < 0 || tbin >= nbin) continue;

    if (IsSentinel999(lo, hi)) {
      vec[tbin].has = false;
      vec[tbin].lo  = lo;
      vec[tbin].hi  = hi;
    } else {
      vec[tbin].has = true;
      vec[tbin].lo  = lo;
      vec[tbin].hi  = hi;
    }
  }

  return data;
}

static std::vector<double> BinCenters(const std::vector<double>& edges) {
  std::vector<double> x;
  if (edges.size() < 2) return x;
  x.reserve(edges.size() - 1);
  for (size_t i = 0; i + 1 < edges.size(); i++) x.push_back(0.5 * (edges[i] + edges[i + 1]));
  return x;
}

static TGraph* MakeGraph(const std::vector<double>& x, const std::vector<double>& y,
                         int mstyle, double msize) {
  auto* g = new TGraph((int)x.size());
  for (int i = 0; i < (int)x.size(); i++) g->SetPoint(i, x[i], y[i]);
  g->SetMarkerStyle(mstyle);
  g->SetMarkerSize(msize);
  g->SetLineWidth(1);
  return g;
}

// Build ratio points for bins where BOTH files have valid cuts
static void BuildWidthRatioXY(const std::vector<double>& xcenters,
                              const std::vector<CutBin>* Abins,
                              const std::vector<CutBin>* Bbins,
                              std::vector<double>& xout,
                              std::vector<double>& yout) {
  xout.clear(); yout.clear();
  if (!Abins || !Bbins) return;

  const int n = std::min({(int)xcenters.size(), (int)Abins->size(), (int)Bbins->size()});
  for (int i = 0; i < n; i++) {
    const auto& a = (*Abins)[i];
    const auto& b = (*Bbins)[i];
    if (!a.has) continue;
    if (!b.has) continue;

    const double wa = a.hi - a.lo;
    const double wb = b.hi - b.lo;
    if (!std::isfinite(wa) || !std::isfinite(wb)) continue;
    if (std::fabs(wb) < 1e-12) continue;

    const double r = wa / wb;
    if (!std::isfinite(r)) continue;

    xout.push_back(xcenters[i]);
    yout.push_back(r);
  }
}

// Scan all cfg+var ratio points -> global ymin/ymax
static bool ScanGlobalY(const std::vector<double>& xcenters,
                        const CutFileData& A,
                        const CutFileData& B,
                        const std::vector<std::string>& vars,
                        double& yMinOut,
                        double& yMaxOut) {
  double ymin = std::numeric_limits<double>::infinity();
  double ymax = -std::numeric_limits<double>::infinity();
  bool found = false;

  for (int cfg = 0; cfg <= 2; cfg++) {
    for (const auto& var : vars) {
      const std::vector<CutBin>* Abins = nullptr;
      const std::vector<CutBin>* Bbins = nullptr;

      auto itAvar = A.cuts.find(var);
      if (itAvar != A.cuts.end()) {
        auto itAcfg = itAvar->second.find(cfg);
        if (itAcfg != itAvar->second.end()) Abins = &itAcfg->second;
      }
      auto itBvar = B.cuts.find(var);
      if (itBvar != B.cuts.end()) {
        auto itBcfg = itBvar->second.find(cfg);
        if (itBcfg != itBvar->second.end()) Bbins = &itBcfg->second;
      }

      std::vector<double> xr, yr;
      BuildWidthRatioXY(xcenters, Abins, Bbins, xr, yr);

      for (double v : yr) {
        if (!std::isfinite(v)) continue;
        ymin = std::min(ymin, v);
        ymax = std::max(ymax, v);
        found = true;
      }
    }
  }

  if (!found) return false;

  // add padding once globally
  const double pad = 0.10 * (ymax - ymin + 1e-12);
  yMinOut = ymin - pad;
  yMaxOut = ymax + pad;
  return true;
}

void PlotCutsRatioWidth(const char* file1 = "dvcs_cuts_rec_bin.csv",
                        const char* file2 = "dvcs_cuts_gen_bin.csv") {
  gStyle->SetOptStat(0);

  if (std::string(file2).empty()) {
    std::cerr << "[ERROR] Need TWO csv files to make ratio (hi-lo)_1/(hi-lo)_2.\n";
    return;
  }

  CutFileData A = ReadCutFile(file1);
  CutFileData B = ReadCutFile(file2);

  if (A.tEdges.empty() || B.tEdges.empty()) {
    std::cerr << "[ERROR] Missing tEdges in one of the files.\n";
    return;
  }

  // t bin must be identical
  if (A.tEdges != B.tEdges) {
    std::cerr << "[ERROR] tEdges are NOT identical between two files (t bin 不一致).\n";
    std::cerr << "        file1=" << file1 << "\n";
    std::cerr << "        file2=" << file2 << "\n";
    std::cerr << "        Please rebin/make them consistent before plotting ratio.\n";
    return;
  }

  const auto x = BinCenters(A.tEdges);

  const std::vector<std::string> vars = {
    "Mx2_ep",
    "Emiss",
    "PTmiss",
    "Theta_gamma_gamma",
    "DeltaPhi",
    "Mx2_epg",
    "Mx2_eg"
  };

  // Global Y scan
  double globalYmin = 0.5, globalYmax = 1.5;
  if (!ScanGlobalY(x, A, B, vars, globalYmin, globalYmax)) {
    std::cerr << "[WARN] No valid ratio points found anywhere; using default y-range [0.5,1.5].\n";
  } else {
    std::cout << "[INFO] Global Y range: [" << globalYmin << ", " << globalYmax << "]\n";
  }

  // Global X range
  double xmin = 0.0, xmax = 1.0;
  if (!x.empty()) {
    xmin = *std::min_element(x.begin(), x.end());
    xmax = *std::max_element(x.begin(), x.end());
  }

  const int nx = 4, ny = 2;

  for (int cfg = 0; cfg <= 2; cfg++) {
    TString cname = Form("c_cfg%d", cfg);
    auto* c = new TCanvas(cname, cname, 1800, 900);
    c->Divide(nx, ny, 0.01, 0.01);

    for (int iv = 0; iv < (int)vars.size(); iv++) {
      c->cd(iv + 1);
      gPad->SetGrid(1, 1);

      const std::string& var = vars[iv];

      const std::vector<CutBin>* Abins = nullptr;
      const std::vector<CutBin>* Bbins = nullptr;

      auto itAvar = A.cuts.find(var);
      if (itAvar != A.cuts.end()) {
        auto itAcfg = itAvar->second.find(cfg);
        if (itAcfg != itAvar->second.end()) Abins = &itAcfg->second;
      }
      auto itBvar = B.cuts.find(var);
      if (itBvar != B.cuts.end()) {
        auto itBcfg = itBvar->second.find(cfg);
        if (itBcfg != itBvar->second.end()) Bbins = &itBcfg->second;
      }

      std::vector<double> xr, yr;
      BuildWidthRatioXY(x, Abins, Bbins, xr, yr);

      auto* frame = new TGraph(2);
      frame->SetPoint(0, xmin, globalYmin);
      frame->SetPoint(1, xmax, globalYmax);
      frame->SetTitle("");
      frame->Draw("AP");
      frame->GetXaxis()->SetTitle("-t (bin center) [GeV^{2}]");
      frame->GetYaxis()->SetTitle("(hi-lo)_{file1} / (hi-lo)_{file2}");
      frame->GetXaxis()->SetTitleSize(0.055);
      frame->GetYaxis()->SetTitleSize(0.055);
      frame->GetXaxis()->SetLabelSize(0.050);
      frame->GetYaxis()->SetLabelSize(0.050);

      if (!xr.empty()) {
        auto* g = MakeGraph(xr, yr, 20, 1.0);
        g->Draw("P SAME");
      }

      TLatex lat;
      lat.SetNDC(true);
      lat.SetTextSize(0.060);
      lat.DrawLatex(0.12, 0.92, Form("cfg=%d : %s", cfg, var.c_str()));

      auto* leg = new TLegend(0.12, 0.12, 0.88, 0.28);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.045);
      leg->AddEntry((TObject*)0, Form("file1: %s", file1), "");
      leg->AddEntry((TObject*)0, Form("file2: %s", file2), "");
      leg->Draw();
    }

    TString outPdf = Form("cuts_cfg%d_widthRatio.pdf", cfg);
    TString outPng = Form("cuts_cfg%d_widthRatio.png", cfg);
    c->SaveAs(outPdf);
    c->SaveAs(outPng);
  }

  std::cout << "[OK] Saved cuts_cfg0/1/2_widthRatio.(pdf/png)\n";
}
