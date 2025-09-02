// DISANA_FitUtils.h
#pragma once
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include <utility>
#include <stdexcept>
#include <ROOT/RDataFrame.hxx>


// Model: A*Gaussian(x; mu, sigma) + B*exp(C*x)
// (simple, robust; if you prefer a one-sided tail, switch to a Crystal-Ball)
inline TF1* FitGaussianExpTail(TH1* h, double mean_guess = 0.4937, double sigma_guess = 0.030) {
  if (!h) throw std::runtime_error("FitGaussianExpTail: null histogram");
  const double xmin = h->GetXaxis()->GetXmin();
  const double xmax = h->GetXaxis()->GetXmax();

  TF1* f = new TF1("fGausExp",
                   "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp([4]*x)",
                   xmin, xmax);
  // [0]=A, [1]=mu, [2]=sigma, [3]=B, [4]=C
  f->SetParNames("A","mu","sigma","B","C");

  // start values & sane limits
  f->SetParameters(h->GetMaximum(), mean_guess, sigma_guess, h->GetMaximum()*0.1, -6.0);
  f->SetParLimits(2, 0.005, 0.120);  // 5–120 MeV width
  f->SetParLimits(4, -50.0, -0.01);  // negative slope for tail

  h->Fit(f, "RQ");   // Quiet, bounded fit
  return f;
}

// Turn the Gaussian core into a window
inline std::pair<double,double> FitWindowFromGaussianCore(TF1* f, double nSigma=3.0) {
  if (!f) throw std::runtime_error("FitWindowFromGaussianCore: null TF1");
  const double mu = f->GetParameter(1);
  const double sg = f->GetParameter(2);
  if (!(sg > 0) || !std::isfinite(mu) || !std::isfinite(sg))
    throw std::runtime_error("Fit failed: invalid mu/sigma");
  return {mu - nSigma*sg, mu + nSigma*sg};
}

// --- Add to DISANAMathFitUtils.h -------------------------------------------
inline TF1* FitGausPlusPoly3(TH1* h,
                             double mu_guess = 0.4937,
                             double sigma_guess = 0.030,
                             double xmin = -1, double xmax = -1) {
  if (!h) throw std::runtime_error("FitGausPlusPoly3: null histogram");
  if (xmin < 0 || xmax < 0) {
    xmin = h->GetXaxis()->GetXmin();
    xmax = h->GetXaxis()->GetXmax();
  }
  // unique name per hist to avoid TF1 name clashes
  std::string fname = Form("fGausP3_%p", (void*)h);
  TF1* f = new TF1(fname.c_str(), "gaus(0)+pol3(3)", xmin, xmax);
  f->SetParNames("A","mu","sigma","p0","p1","p2","p3");
  f->SetParameters(h->GetMaximum(), mu_guess, sigma_guess, 0,0,0,0);
  f->SetParLimits(2, 0.003, 0.120); // 3–120 MeV width
  h->Fit(f, "RQ0");                 // quiet, bounded, don't draw
  return f;
}
// NOTE: FitWindowFromGaussianCore you already have works here too because
// gaus(0) uses [1]=mu and [2]=sigma.
