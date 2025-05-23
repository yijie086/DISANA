#include "DrawAndSave.h"
#include "TH1.h"
#include <string>
#include "TFile.h"
#include <iostream>
#include "RHipoDS.hxx"
#include "TCanvas.h"

#include "../ParticleInformation/RECParticle.h"

void Draw1DHist(TH1* hist, const std::string& histname, int nbins, double xmin, double xmax) {
    hist->SetBins(nbins, xmin, xmax);
    TCanvas* c = new TCanvas(histname.c_str(), histname.c_str(), 2000, 1500);
    hist->Draw();
    c->SaveAs((histname + ".png").c_str());
    delete c;
}
void Save1DHist(TH1* hist, const std::string& histname, TFile* fout) {
    hist->Write(histname.c_str(), TObject::kOverwrite);
}

void Draw2DHist(TH2* hist, const std::string& histname, int nbinsx, double xmin, double xmax,
                int nbinsy, double ymin, double ymax) {
    hist->SetBins(nbinsx, xmin, xmax, nbinsy, ymin, ymax);
    TCanvas* c = new TCanvas(histname.c_str(), histname.c_str(), 2000, 1500);
    hist->Draw("COLZ");
    c->SetLogz();
    hist->SetStats(0); // Disable the statistics box
    c->SaveAs((histname + ".png").c_str());
    delete c;
}
void Save2DHist(TH2* hist, const std::string& histname, TFile* fout) {
    hist->Write(histname.c_str(), TObject::kOverwrite);
}

void DrawAndSaveEventQ2(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                                   int bins_Q2, double min_Q2, double max_Q2) {
    auto h_Q2 = df
        .Define("Q2", get_RECParticle_float_var(pid,charge), {"REC_Event_Q2", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Histo1D({"Q2", "Q2", bins_Q2, min_Q2, max_Q2}, "Q2");
    Draw1DHist(h_Q2.GetPtr(), "h_Q2", bins_Q2, min_Q2, max_Q2);
    Save1DHist(h_Q2.GetPtr(), "Q2", fout);
}
void DrawAndSaveEventxB(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                                   int bins_xB, double min_xB, double max_xB) {
    auto h_xB = df
        .Define("xB", get_RECParticle_float_var(pid,charge), {"REC_Event_xB", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Histo1D({"xB", "xB", bins_xB, min_xB, max_xB}, "xB");
    Draw1DHist(h_xB.GetPtr(), "h_xB", bins_xB, min_xB, max_xB);
    Save1DHist(h_xB.GetPtr(), "xB", fout);
}
void DrawAndSaveEventNu(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                                   int bins_nu, double min_nu, double max_nu) {
    auto h_nu = df
        .Define("Nu", get_RECParticle_float_var(pid,charge), {"REC_Event_Nu", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Histo1D({"Nu", "Nu", bins_nu, min_nu, max_nu}, "Nu");
    Draw1DHist(h_nu.GetPtr(), "h_Nu", bins_nu, min_nu, max_nu);
    Save1DHist(h_nu.GetPtr(), "Nu", fout);
}

void DrawAndSaveEventW(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                                   int bins_W, double min_W, double max_W) {
    auto h_W = df
        .Define("W", get_RECParticle_float_var(pid,charge), {"REC_Event_W", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Histo1D({"W", "W", bins_W, min_W, max_W}, "W");
    Draw1DHist(h_W.GetPtr(), "h_W", bins_W, min_W, max_W);
    Save1DHist(h_W.GetPtr(), "W", fout);
}

void DrawAndSaveEventmt(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                                   int bins_t, double min_t, double max_t) {
    auto h_t = df
        .Define("mt", get_RECParticle_float_var(pid,charge), {"REC_Event_mt", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Histo1D({"mt", "-t", bins_t, min_t, max_t}, "mt");
    Draw1DHist(h_t.GetPtr(), "h_mt", bins_t, min_t, max_t);
    Save1DHist(h_t.GetPtr(), "mt", fout);
}

void DrawAndSaveQ2vsxB(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                       int bins_Q2, double min_Q2, double max_Q2, 
                       int bins_xB, double min_xB, double max_xB) {
    auto h_Q2vsxB = df
        .Define("Q2", get_RECParticle_float_var(pid, charge), {"REC_Event_Q2", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Define("xB", get_RECParticle_float_var(pid, charge), {"REC_Event_xB", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Histo2D({"Q2_vs_xB", "Q2 vs xB", bins_xB, min_xB, max_xB, bins_Q2, min_Q2, max_Q2}, "xB", "Q2");
    Draw2DHist(h_Q2vsxB.GetPtr(), "h_Q2vsxB", bins_xB, min_xB, max_xB, bins_Q2, min_Q2, max_Q2);
    Save2DHist(h_Q2vsxB.GetPtr(), "Q2_vs_xB", fout);
}

void DrawAndSaveParticleHistograms(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                                   int bins_theta, double min_theta, double max_theta,
                                   int bins_phi, double min_phi, double max_phi,
                                   int bins_p, double min_p, double max_p) {
    // 动态生成名称和标题
    std::string particle_name;
    if (pid == 11) {
        particle_name = "Electron";
    } else if (pid == 22) {
        particle_name = "Photon";
    } else if (pid == 2212) {
        particle_name = "Proton";
    } else {
        particle_name = "PID_" + std::to_string(pid);
    }

    // Theta histogram
    std::string theta_name = "h_" + particle_name + "_theta";
    std::string theta_title = particle_name + "_theta";
    auto h_theta = df
        .Define("theta", get_RECParticle_float_var(pid,charge), {"REC_Particle_theta", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Histo1D({theta_title.c_str(), theta_title.c_str(), bins_theta, min_theta, max_theta}, "theta");
    Draw1DHist(h_theta.GetPtr(), theta_name.c_str(), bins_theta, min_theta, max_theta);
    Save1DHist(h_theta.GetPtr(), theta_title.c_str(), fout);

    // Phi histogram
    std::string phi_name = "h_" + particle_name + "_phi";
    std::string phi_title = particle_name + "_phi";
    auto h_phi = df
        .Define("phi", get_RECParticle_float_var(pid,charge), {"REC_Particle_phi", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Histo1D({phi_title.c_str(), phi_title.c_str(), bins_phi, min_phi, max_phi}, "phi");
    Draw1DHist(h_phi.GetPtr(), phi_name.c_str(), bins_phi, min_phi, max_phi);
    Save1DHist(h_phi.GetPtr(), phi_title.c_str(), fout);

    // Momentum histogram
    std::string p_name = "h_" + particle_name + "_p";
    std::string p_title = particle_name + "_p";
    auto h_p = df
        .Define("p", get_RECParticle_float_var(pid,charge), {"REC_Particle_p", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Histo1D({p_title.c_str(), p_title.c_str(), bins_p, min_p, max_p}, "p");
    Draw1DHist(h_p.GetPtr(), p_name.c_str(), bins_p, min_p, max_p);
    Save1DHist(h_p.GetPtr(), p_title.c_str(), fout);

    std::string charg_name = "h_" + particle_name + "_charge";
    std::string charg_title = particle_name + "_charge";
    auto h_charge = df
        .Define("charge", get_RECParticle_int_var(pid,charge), {"REC_Particle_charge", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Histo1D({charg_title.c_str(), charg_title.c_str(), 500, -2, 2}, "charge");
    Draw1DHist(h_charge.GetPtr(), charg_name.c_str(), 500, -2, 2);
    Save1DHist(h_charge.GetPtr(), charg_title.c_str(), fout);
}