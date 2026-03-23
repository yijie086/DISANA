// PlotRecCuts.C
// Usage:
//   root -l -b -q 'PlotRecCuts.C("cuts1.csv","cuts2.csv")'
//   root -l -b -q 'PlotRecCuts.C("cuts1.csv","")'
//
// Output: cuts_cfg0.pdf/png, cuts_cfg1.pdf/png, cuts_cfg2.pdf/png

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
      vec[tbin].has = false; // don't plot
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

static void FillXY(const std::vector<double>& xcenters,
                   const std::vector<CutBin>* bins,
                   bool useLo,
                   std::vector<double>& xout,
                   std::vector<double>& yout) {
  xout.clear(); yout.clear();
  if (!bins) return;
  const int n = std::min((int)xcenters.size(), (int)bins->size());
  for (int i = 0; i < n; i++) {
    if (!(*bins)[i].has) continue;
    xout.push_back(xcenters[i]);
    yout.push_back(useLo ? (*bins)[i].lo : (*bins)[i].hi);
  }
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

static void AppendAll(std::vector<double>& dst, const std::vector<double>& src) {
  dst.insert(dst.end(), src.begin(), src.end());
}

void PlotRecCuts(const char* file1 = "dvcs_cuts_rec_bin.csv",
                 const char* file2 = "dvcs_cuts_rec.csv") {
  gStyle->SetOptStat(0);

  CutFileData A = ReadCutFile(file1);
  CutFileData B = ReadCutFile(file2);

  const bool hasA = !std::string(file1).empty() && !A.tEdges.empty();
  const bool hasB = !std::string(file2).empty() && !B.tEdges.empty();

  const auto xA = BinCenters(A.tEdges);
  const auto xB = BinCenters(B.tEdges);

  const std::vector<std::string> vars = {
    "Mx2_ep",
    "Emiss",
    "PTmiss",
    "Theta_gamma_gamma",
    "DeltaPhi",
    "Mx2_epg",
    "Mx2_eg"
  };

  const int nx = 4, ny = 2; // 8 pads

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

      std::vector<double> xA_lo, yA_lo, xA_hi, yA_hi;
      std::vector<double> xB_lo, yB_lo, xB_hi, yB_hi;

      FillXY(xA, Abins, true,  xA_lo, yA_lo);
      FillXY(xA, Abins, false, xA_hi, yA_hi);

      FillXY(xB, Bbins, true,  xB_lo, yB_lo);
      FillXY(xB, Bbins, false, xB_hi, yB_hi);

      // x-range: union of existing centers (not points) to keep stable framing
      double xmin = 0.0, xmax = 1.0;
      std::vector<double> allX;
      if (hasA) AppendAll(allX, xA);
      if (hasB) AppendAll(allX, xB);
      if (!allX.empty()) {
        xmin = *std::min_element(allX.begin(), allX.end());
        xmax = *std::max_element(allX.begin(), allX.end());
      }

      // y-range: union of actually plotted points
      std::vector<double> allY;
      AppendAll(allY, yA_lo); AppendAll(allY, yA_hi);
      AppendAll(allY, yB_lo); AppendAll(allY, yB_hi);

      double ymin = 0.0, ymax = 1.0;
      if (!allY.empty()) {
        ymin = *std::min_element(allY.begin(), allY.end());
        ymax = *std::max_element(allY.begin(), allY.end());
        double pad = 0.10 * (ymax - ymin + 1e-12);
        ymin -= pad; ymax += pad;
      }

      // frame
      auto* frame = new TGraph(2);
      frame->SetPoint(0, xmin, ymin);
      frame->SetPoint(1, xmax, ymax);
      frame->SetTitle("");
      frame->Draw("AP");
      frame->GetXaxis()->SetTitle("-t (bin center) [GeV^{2}]");
      frame->GetYaxis()->SetTitle(Form("%s cut value", var.c_str()));
      frame->GetXaxis()->SetTitleSize(0.055);
      frame->GetYaxis()->SetTitleSize(0.055);
      frame->GetXaxis()->SetLabelSize(0.050);
      frame->GetYaxis()->SetLabelSize(0.050);

      // draw graphs
      // File1: lo/hi
      if (!xA_lo.empty()) { auto* gAlo = MakeGraph(xA_lo, yA_lo, 20, 1.0); gAlo->Draw("P SAME"); }
      if (!xA_hi.empty()) { auto* gAhi = MakeGraph(xA_hi, yA_hi, 24, 1.0); gAhi->Draw("P SAME"); }

      // File2: lo/hi
      if (!xB_lo.empty()) { auto* gBlo = MakeGraph(xB_lo, yB_lo, 21, 1.0); gBlo->Draw("P SAME"); }
      if (!xB_hi.empty()) { auto* gBhi = MakeGraph(xB_hi, yB_hi, 25, 1.0); gBhi->Draw("P SAME"); }

      // title
      TLatex lat;
      lat.SetNDC(true);
      lat.SetTextSize(0.060);
      lat.DrawLatex(0.12, 0.92, Form("cfg=%d : %s", cfg, var.c_str()));

      // legend
      auto* leg = new TLegend(0.12, 0.12, 0.88, 0.33);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.045);

      if (!xA_lo.empty() || !xA_hi.empty()) {
        leg->AddEntry((TObject*)0, Form("File1: %s", file1), "");
        //if (!xA_lo.empty()) leg->AddEntry((TObject*)0, "  lo: filled circle", "");
        //if (!xA_hi.empty()) leg->AddEntry((TObject*)0, "  hi: open circle", "");
      }
      if (!std::string(file2).empty() && (!xB_lo.empty() || !xB_hi.empty())) {
        leg->AddEntry((TObject*)0, Form("File2: %s", file2), "");
        //if (!xB_lo.empty()) leg->AddEntry((TObject*)0, "  lo: filled square", "");
        //if (!xB_hi.empty()) leg->AddEntry((TObject*)0, "  hi: open square", "");
      }
      leg->Draw();
    }

    TString outPdf = Form("cuts_cfg%d.pdf", cfg);
    TString outPng = Form("cuts_cfg%d.png", cfg);
    c->SaveAs(outPdf);
    c->SaveAs(outPng);
  }

  std::cout << "[OK] Saved cuts_cfg0/1/2.(pdf/png)\n";
}
