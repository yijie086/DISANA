#ifndef DRAWANDSAVE_H
#define DRAWANDSAVE_H

#include "TH1.h"
#include <string>
#include "TFile.h"
#include <iostream>
#include "RHipoDS.hxx"
#include "TCanvas.h"

void Draw1DHist(TH1* hist, const std::string& histname, int nbins, double xmin, double xmax);
void Save1DHist(TH1* hist, const std::string& histname, TFile* fout);

void DrawAndSaveParticleHistograms(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                                   int bins_theta, double min_theta, double max_theta,
                                   int bins_phi, double min_phi, double max_phi,
                                   int bins_p, double min_p, double max_p);

#endif