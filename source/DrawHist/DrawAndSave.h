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

void DrawTProfile(TProfile* prof, const std::string& profname, int nbinsx, double xmin, double xmax);
void SaveTProfile(TProfile* prof, const std::string& profname, TFile* fout);

void DrawAndSaveEventQ2(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                                   int bins_Q2, double min_Q2, double max_Q2, std::string output_name);
void DrawAndSaveEventxB(ROOT::RDF::RNode df, int pid, int charge, TFile* fout,
                                   int bins_xB, double min_xB, double max_xB, std::string output_name);
void DrawAndSaveEventNu(ROOT::RDF::RNode df, int pid, int charge, TFile* fout,
                                   int bins_nu, double min_nu, double max_nu, std::string output_name);
void DrawAndSaveEventW(ROOT::RDF::RNode df, int pid, int charge, TFile* fout,
                                   int bins_W, double min_W, double max_W, std::string output_name);
void DrawAndSaveEventmt(ROOT::RDF::RNode df, int pid, int charge, TFile* fout,
                                   int bins_t, double min_t, double max_t, std::string output_name);

void DrawAndSaveQ2vsxB(ROOT::RDF::RNode df, int pid, int charge, TFile* fout,
                        int bins_Q2, double min_Q2, double max_Q2, 
                        int bins_xB, double min_xB, double max_xB, std::string output_name);

void DrawAndSavechi2perndfvsedge(ROOT::RDF::RNode df, int detector, int layer, int pid, int charge, TFile* fout, 
                           int bins_chi2, double min_chi2, double max_chi2,
                           int bins_edge, double min_edge, double max_edge, std::vector<std::string> input_name_edge, std::string output_name);
void DrawAndSaveedge(ROOT::RDF::RNode df, int detector, int layer, int pid, int charge, TFile* fout, 
                           int bins_edge, double min_edge, double max_edge, std::vector<std::string> input_name_edge, std::string output_name);
void DrawAndSavechi2perndf(ROOT::RDF::RNode df, int detector, int layer, int pid, int charge, TFile* fout,
                           int bins_chi2, double min_chi2, double max_chi2, std::string output_name);

void DrawAndSaveParticleHistograms(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                                   int bins_theta, double min_theta, double max_theta,
                                   int bins_phi, double min_phi, double max_phi,
                                   int bins_p, double min_p, double max_p, std::string output_name);

#endif