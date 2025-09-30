// File: plotDeltap_simple_extended.C
// ROOT macro to plot Δp(θ, p) with two models:
// (1) "baseline" parameterization (yours)
// (2) "fit-from-figure" parameterization using the A_p, B_p, C_p fits (quadratic)
// and including the C_p(θ)/p^2 term.

#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TColor.h"
#include <vector>

// Parameterizations
inline double A_p(double theta) { return  0.0135790  - 0.0005303  * theta; }
inline double B_p(double theta) { return -0.02165929 + 0.00121123 * theta; }
inline double C_p(double /*theta*/) { return 0.0; }

// Δp(θ, p)
inline double DeltaP(double theta, double p) {
  return A_p(theta) + B_p(theta)/p + C_p(theta);
}

// ---------------------- Fit-from-figure parameterization ----------------------
// Read off from the provided figure: Ap, Bp, Cp are quadratic in θ (degrees)
inline double A_p_fit(double theta) {
  return  0.0074727 + (-0.0001729)*theta + (-0.0000007)*theta*theta;
}
inline double B_p_fit(double theta) {
  return  0.0098664 + (-0.0019167)*theta + ( 0.0000719)*theta*theta;
}
inline double C_p_fit(double theta) {
  return  0.0046117 + ( 0.0004211)*theta + (-0.0000184)*theta*theta;
}

inline double DeltaP_fit(double theta, double p) {
  // include Cp term as 1/p^2 (units in the figure are GeV, GeV^2, GeV^3)
  return A_p_fit(theta) + B_p_fit(theta)/p + C_p_fit(theta)/(p*p);
}
void plotDeltap_simple_extended();
void test_energyloss_corrWithTimothy_par() {
  gStyle->SetOptStat(0);

  // Theta range and resolution
  const double theta_min = 0.0;    // degrees
  const double theta_max = 40.0;   // degrees
  const int    npts      = 901;    // step = 0.1 deg

  // Momenta
  std::vector<double> mom = {0.75, 1.10, 1.75, 2.75};

  // Colors
  std::vector<int> colors  = {kBlue+1, kRed+1, kGreen+2, kMagenta+2};

  // Graphs
  std::vector<TGraph*> graphs;
  for (size_t i = 0; i < mom.size(); ++i) {
    auto *g = new TGraph(npts);
    for (int j = 0; j < npts; ++j) {
      double t = theta_min + (theta_max - theta_min) * j / (npts - 1.0);
      g->SetPoint(j, t, DeltaP(t, mom[i]));
    }
    g->SetLineColor(colors[i % colors.size()]);
    g->SetLineWidth(2);
    graphs.push_back(g);
  }

  // Canvas + multigraph
  TCanvas *c = new TCanvas("c", "Delta p vs theta", 900, 600);
  TMultiGraph *mg = new TMultiGraph();
  for (auto *g : graphs) mg->Add(g, "L");

  mg->SetTitle("#Delta p(#theta, p) = A_{p}(#theta) + B_{p}(#theta)/p;#theta [deg];#Delta p");
  mg->Draw("A");
  mg->GetXaxis()->SetLimits(theta_min, theta_max);

  // Legend
  TLegend *leg = new TLegend(0.58, 0.68, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  for (size_t i = 0; i < graphs.size(); ++i) {
    leg->AddEntry(graphs[i], Form("p = %.2f", mom[i]), "l");
  }
  leg->Draw();

  c->SetGrid();
  c->Update();

  // Save
  c->Print("DeltaP_simple.pdf");
  c->Print("DeltaP_simple.png");
  plotDeltap_simple_extended();
}


// -----------------------------------------------------------------------
void plotDeltap_simple_extended() {
  gStyle->SetOptStat(0);

  // θ range
  const double theta_min = 5.0;     // degrees (plots start ~5° in the figure)
  const double theta_max = 40.0;    // degrees
  const int    npts      = 351;     // 0.1° steps

  // Momenta to show
  std::vector<double> mom = {0.75, 1.10, 1.75, 2.75};

  // Color palette for momenta
  std::vector<int> colors = {kBlue+1, kRed+1, kGreen+2, kMagenta+2};

  // Build graphs
  std::vector<TGraph*> g_base, g_fit;

  for (size_t i = 0; i < mom.size(); ++i) {
    auto *gf = new TGraph(npts);
    for (int j = 0; j < npts; ++j) {
      double t = theta_min + (theta_max - theta_min) * j / (npts - 1.0);
      gf->SetPoint(j, t, DeltaP_fit (t, mom[i]));
    }

    gf->SetLineColor(colors[i % colors.size()]);
    gf->SetLineWidth(3);
    gf->SetLineStyle(kDashed); // fit curves dashed for clarity
    g_fit.push_back(gf);
  }

  // Canvas + multigraph
  TCanvas *c = new TCanvas("c", "Delta p vs theta (baseline vs fit)", 1100, 700);
  TMultiGraph *mg = new TMultiGraph();

  //for (auto *g : g_base) mg->Add(g, "L");
  for (auto *g : g_fit ) mg->Add(g, "L");

  mg->SetTitle("#Delta p(#theta, p FD, sp2018 inb);#theta (degrees);#Delta p (GeV)");
  mg->Draw("A");
  mg->GetXaxis()->SetLimits(theta_min, theta_max);

  // Legend
  TLegend *leg = new TLegend(0.52, 0.62, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  for (size_t i = 0; i < mom.size(); ++i) {
    //leg->AddEntry(g_base[i], Form("p = %.2f  (baseline)", mom[i]), "l");
    leg->AddEntry(g_fit [i], Form("p = %.2f  (fit from figure)", mom[i]), "l");
  }
  leg->Draw();

  c->SetGrid();
  c->Update();

  // Save
  c->Print("DeltaP_FD_sp18inb_simple_extended.pdf");
  c->Print("DeltaP_FD_sp18inb_simple_extended.png");
}
