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
    std::string theta_title = particle_name + " Theta";
    auto h_theta = df
        .Define("theta", get_RECParticle_float_var(pid,charge), {"REC_Particle_theta", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Histo1D({theta_name.c_str(), theta_title.c_str(), bins_theta, min_theta, max_theta}, "theta");
    Draw1DHist(h_theta.GetPtr(), theta_name.c_str(), bins_theta, min_theta, max_theta);
    Save1DHist(h_theta.GetPtr(), theta_name.c_str(), fout);

    // Phi histogram
    std::string phi_name = "h_" + particle_name + "_phi";
    std::string phi_title = particle_name + " Phi";
    auto h_phi = df
        .Define("phi", get_RECParticle_float_var(pid,charge), {"REC_Particle_phi", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Histo1D({phi_name.c_str(), phi_title.c_str(), bins_phi, min_phi, max_phi}, "phi");
    Draw1DHist(h_phi.GetPtr(), phi_name.c_str(), bins_phi, min_phi, max_phi);
    Save1DHist(h_phi.GetPtr(), phi_name.c_str(), fout);

    // Momentum histogram
    std::string p_name = "h_" + particle_name + "_p";
    std::string p_title = particle_name + " Momentum";
    auto h_p = df
        .Define("p", get_RECParticle_float_var(pid,charge), {"REC_Particle_p", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Histo1D({p_name.c_str(), p_title.c_str(), bins_p, min_p, max_p}, "p");
    Draw1DHist(h_p.GetPtr(), p_name.c_str(), bins_p, min_p, max_p);
    Save1DHist(h_p.GetPtr(), p_name.c_str(), fout);

    std::string histname = "h_" + particle_name + "_charge";
    std::string histtitle = particle_name + " Charge";
    auto h_charge = df
        .Define("charge", get_RECParticle_int_var(pid,charge), {"REC_Particle_charge", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p"})
        .Histo1D({histname.c_str(), histtitle.c_str(), 500, -2, 2}, "charge");
    Draw1DHist(h_charge.GetPtr(), histname.c_str(), 500, -2, 2);
    Save1DHist(h_charge.GetPtr(), histname.c_str(), fout);
}