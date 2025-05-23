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

void Draw2DHist(TH2* hist, const std::string& histname, int nbinsx, double xmin, double xmax,
                int nbinsy, double ymin, double ymax);
void Save2DHist(TH2* hist, const std::string& histname, TFile* fout);

void DrawAndSaveEventQ2(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                                   int bins_Q2, double min_Q2, double max_Q2);
void DrawAndSaveEventxB(ROOT::RDF::RNode df, int pid, int charge, TFile* fout,
                                   int bins_xB, double min_xB, double max_xB);
void DrawAndSaveEventNu(ROOT::RDF::RNode df, int pid, int charge, TFile* fout,
                                   int bins_nu, double min_nu, double max_nu);
void DrawAndSaveEventW(ROOT::RDF::RNode df, int pid, int charge, TFile* fout,
                                   int bins_W, double min_W, double max_W);
void DrawAndSaveEventmt(ROOT::RDF::RNode df, int pid, int charge, TFile* fout,
                                   int bins_t, double min_t, double max_t);

void DrawAndSaveQ2vsxB(ROOT::RDF::RNode df, int pid, int charge, TFile* fout,
                        int bins_Q2, double min_Q2, double max_Q2, 
                        int bins_xB, double min_xB, double max_xB);

void DrawAndSaveParticleHistograms(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                                   int bins_theta, double min_theta, double max_theta,
                                   int bins_phi, double min_phi, double max_phi,
                                   int bins_p, double min_p, double max_p);

#endif