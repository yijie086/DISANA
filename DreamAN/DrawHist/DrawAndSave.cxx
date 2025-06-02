#include "DrawAndSave.h"
#include "TH1.h"
#include <string>
#include "TFile.h"
#include <iostream>
#include "RHipoDS.hxx"
#include "TCanvas.h"

#include "../ParticleInformation/RECParticle.h"
#include "../ParticleInformation/RECTraj.h"
#include "../core/Columns.h"

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

void DrawTProfile(TProfile* prof, const std::string& profname, int nbinsx, double xmin, double xmax) {
    prof->SetBins(nbinsx, xmin, xmax);
    TCanvas* c = new TCanvas(profname.c_str(), profname.c_str(), 2000, 1500);
    prof->Draw();
    c->SaveAs((profname + ".png").c_str());
    delete c;
}
void SaveTProfile(TProfile* prof, const std::string& profname, TFile* fout) {
    prof->Write(profname.c_str(), TObject::kOverwrite);
}

void DrawAndSaveEventQ2(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                        int bins_Q2, double min_Q2, double max_Q2, std::string output_name) {
    std::string Q2_name = "h_" + output_name + "Q2";
    std::string Q2_title = output_name + "Q2";
    auto h_Q2 = df
        .Define("Q2", get_RECParticle_float_var(pid, charge), {"REC_Event_Q2", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p", "REC_Traj_pass"})
        .Histo1D({Q2_title.c_str(), Q2_title.c_str(), bins_Q2, min_Q2, max_Q2}, "Q2");
    Draw1DHist(h_Q2.GetPtr(), Q2_name.c_str(), bins_Q2, min_Q2, max_Q2);
    Save1DHist(h_Q2.GetPtr(), Q2_title.c_str(), fout);
}

void DrawAndSaveEventxB(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                        int bins_xB, double min_xB, double max_xB, std::string output_name) {
    std::string xB_name = "h_" + output_name + "xB";
    std::string xB_title = output_name + "xB";
    auto h_xB = df
        .Define("xB", get_RECParticle_float_var(pid, charge), {"REC_Event_xB", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p", "REC_Traj_pass"})
        .Histo1D({xB_title.c_str(), xB_title.c_str(), bins_xB, min_xB, max_xB}, "xB");
    Draw1DHist(h_xB.GetPtr(), xB_name.c_str(), bins_xB, min_xB, max_xB);
    Save1DHist(h_xB.GetPtr(), xB_title.c_str(), fout);
}

void DrawAndSaveEventNu(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                        int bins_nu, double min_nu, double max_nu, std::string output_name) {
    std::string nu_name = "h_" + output_name + "Nu";
    std::string nu_title = output_name + "Nu";
    auto h_nu = df
        .Define("Nu", get_RECParticle_float_var(pid, charge), {"REC_Event_Nu", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p", "REC_Traj_pass"})
        .Histo1D({nu_title.c_str(), nu_title.c_str(), bins_nu, min_nu, max_nu}, "Nu");
    Draw1DHist(h_nu.GetPtr(), nu_name.c_str(), bins_nu, min_nu, max_nu);
    Save1DHist(h_nu.GetPtr(), nu_title.c_str(), fout);
}

void DrawAndSaveEventW(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                       int bins_W, double min_W, double max_W, std::string output_name) {
    std::string W_name = "h_" + output_name + "W";
    std::string W_title = output_name + "W";
    auto h_W = df
        .Define("W", get_RECParticle_float_var(pid, charge), {"REC_Event_W", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p", "REC_Traj_pass"})
        .Histo1D({W_title.c_str(), W_title.c_str(), bins_W, min_W, max_W}, "W");
    Draw1DHist(h_W.GetPtr(), W_name.c_str(), bins_W, min_W, max_W);
    Save1DHist(h_W.GetPtr(), W_title.c_str(), fout);
}

void DrawAndSaveEventmt(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                        int bins_t, double min_t, double max_t, std::string output_name) {
    std::string mt_name = "h_" + output_name + "mt";
    std::string mt_title = output_name + "mt";
    auto h_t = df
        .Define("mt", get_RECParticle_float_var(pid, charge), {"REC_Event_mt", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p", "REC_Traj_pass"})
        .Histo1D({mt_title.c_str(), mt_title.c_str(), bins_t, min_t, max_t}, "mt");
    Draw1DHist(h_t.GetPtr(), mt_name.c_str(), bins_t, min_t, max_t);
    Save1DHist(h_t.GetPtr(), mt_title.c_str(), fout);
}

void DrawAndSaveQ2vsxB(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                       int bins_Q2, double min_Q2, double max_Q2, 
                       int bins_xB, double min_xB, double max_xB, std::string output_name) {
    std::string Q2vsxB_name = "h_" + output_name + "Q2vsxB";
    std::string Q2vsxB_title = output_name + "Q2_vs_xB";
    auto h_Q2vsxB = df
        .Define("Q2", get_RECParticle_float_var(pid, charge), {"REC_Event_Q2", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p", "REC_Traj_pass"})
        .Define("xB", get_RECParticle_float_var(pid, charge), {"REC_Event_xB", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p", "REC_Traj_pass"})
        .Histo2D({Q2vsxB_title.c_str(), Q2vsxB_title.c_str(), bins_xB, min_xB, max_xB, bins_Q2, min_Q2, max_Q2}, "xB", "Q2");
    Draw2DHist(h_Q2vsxB.GetPtr(), Q2vsxB_name.c_str(), bins_xB, min_xB, max_xB, bins_Q2, min_Q2, max_Q2);
    Save2DHist(h_Q2vsxB.GetPtr(), Q2vsxB_title.c_str(), fout);
}

void DrawAndSavechi2perndfvsedge(ROOT::RDF::RNode df, int detector, int layer, int pid, int charge, TFile* fout, 
                          int bins_chi2, double min_chi2, double max_chi2,
                          int bins_edge, double min_edge, double max_edge, std::vector<std::string> input_name_edge, std::string output_name) {
    std::string char_name = "h_" + output_name +"_chi2perndf_vs_edge";
    std::string char_title = output_name + "_chi2perndf_vs_edge";
    std::string profile_name = "h_" + char_title + "_profile";
    std::string profile_title = char_title + "_profile";
    auto h_chi2vsedge = df
        .Define("chi2", get_RECTraj_float_var(detector, layer, pid, charge), 
                {"REC_Track_chi2perndf", "REC_Traj_detector", "REC_Traj_layer", "REC_Traj_pindex", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p", "REC_Traj_pass"})
        .Define("edge", get_RECTraj_float_var(detector, layer, pid, charge), 
                CombineColumns(input_name_edge, std::vector<std::string>{"REC_Traj_detector", "REC_Traj_layer", "REC_Traj_pindex", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p", "REC_Traj_pass"}))
                //{"REC_Traj_pedge", "REC_Traj_detector", "REC_Traj_layer", "REC_Traj_pindex", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p", "REC_Traj_pass"})
        .Histo2D({char_title.c_str(), char_title.c_str(), bins_edge, min_edge, max_edge, bins_chi2, min_chi2, max_chi2}, "edge", "chi2");
    Draw2DHist(h_chi2vsedge.GetPtr(), char_name.c_str(), bins_edge, min_edge, max_edge, bins_chi2, min_chi2, max_chi2);
    Save2DHist(h_chi2vsedge.GetPtr(), char_title.c_str(), fout);
    DrawTProfile(h_chi2vsedge.GetPtr()->ProfileX(profile_title.c_str(), 1, -1, "s"), profile_name.c_str(), 
                  bins_edge, min_edge, max_edge);
    SaveTProfile(h_chi2vsedge.GetPtr()->ProfileX(profile_title.c_str(), 1, -1, "s"), profile_title.c_str(), fout);
}
void DrawAndSaveedge(ROOT::RDF::RNode df, int detector, int layer, int pid, int charge, TFile* fout, 
                                   int bins_edge, double min_edge, double max_edge, std::vector<std::string> input_name_edge, std::string output_name) {
    std::string edge_name = "h_" + output_name + "_edge";
    std::string edge_title = output_name + "_edge";
    auto h_edge = df
        .Define("edge", get_RECTraj_float_var(detector, layer, pid, charge), 
                CombineColumns(input_name_edge, std::vector<std::string>{"REC_Traj_detector", "REC_Traj_layer", "REC_Traj_pindex", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p", "REC_Traj_pass"}))
                //{"REC_Traj_pedge", "REC_Traj_detector", "REC_Traj_layer", "REC_Traj_pindex", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p", "REC_Traj_pass"})
        .Histo1D({edge_title.c_str(), edge_title.c_str(), bins_edge, min_edge, max_edge}, "edge");
    Draw1DHist(h_edge.GetPtr(), edge_name.c_str(), bins_edge, min_edge, max_edge);
    Save1DHist(h_edge.GetPtr(), edge_title.c_str(), fout);
}
void DrawAndSavechi2perndf(ROOT::RDF::RNode df, int detector, int layer, int pid, int charge, TFile* fout, 
                          int bins_chi2, double min_chi2, double max_chi2, std::string output_name) {
    std::string chi2_name = "h_" + output_name + "_chi2perndf";
    std::string chi2_title = output_name + "_chi2perndf";
    auto h_chi2 = df
        .Define("chi2", get_RECTraj_float_var(detector, layer, pid, charge), {"REC_Track_chi2perndf", "REC_Traj_detector", "REC_Traj_layer", "REC_Traj_pindex", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p", "REC_Traj_pass"})
        .Histo1D({chi2_title.c_str(), chi2_title.c_str(), bins_chi2, min_chi2, max_chi2}, "chi2");
    Draw1DHist(h_chi2.GetPtr(), chi2_name.c_str(), bins_chi2, min_chi2, max_chi2);
    Save1DHist(h_chi2.GetPtr(), chi2_title.c_str(), fout);
}

void DrawAndSaveParticleHistograms(ROOT::RDF::RNode df, int pid, int charge, TFile* fout, 
                                   int bins_theta, double min_theta, double max_theta,
                                   int bins_phi, double min_phi, double max_phi,
                                   int bins_p, double min_p, double max_p, std::string output_name) {
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
    particle_name = particle_name + "_" + output_name;

    // Theta histogram
    std::string theta_name = "h_" + particle_name + "_theta";
    std::string theta_title = particle_name + "_theta";
    auto h_theta = df
        .Define("theta", get_RECParticle_float_var(pid,charge), {"REC_Particle_theta", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p", "REC_Traj_pass"})
        .Histo1D({theta_title.c_str(), theta_title.c_str(), bins_theta, min_theta, max_theta}, "theta");
    Draw1DHist(h_theta.GetPtr(), theta_name.c_str(), bins_theta, min_theta, max_theta);
    Save1DHist(h_theta.GetPtr(), theta_title.c_str(), fout);

    // Phi histogram
    std::string phi_name = "h_" + particle_name + "_phi";
    std::string phi_title = particle_name + "_phi";
    auto h_phi = df
        .Define("phi", get_RECParticle_float_var(pid,charge), {"REC_Particle_phi", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p", "REC_Traj_pass"})
        .Histo1D({phi_title.c_str(), phi_title.c_str(), bins_phi, min_phi, max_phi}, "phi");
    Draw1DHist(h_phi.GetPtr(), phi_name.c_str(), bins_phi, min_phi, max_phi);
    Save1DHist(h_phi.GetPtr(), phi_title.c_str(), fout);

    // Momentum histogram
    std::string p_name = "h_" + particle_name + "_p";
    std::string p_title = particle_name + "_p";
    auto h_p = df
        .Define("p", get_RECParticle_float_var(pid,charge), {"REC_Particle_p", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p", "REC_Traj_pass"})
        .Histo1D({p_title.c_str(), p_title.c_str(), bins_p, min_p, max_p}, "p");
    Draw1DHist(h_p.GetPtr(), p_name.c_str(), bins_p, min_p, max_p);
    Save1DHist(h_p.GetPtr(), p_title.c_str(), fout);

    std::string charg_name = "h_" + particle_name + "_charge";
    std::string charg_title = particle_name + "_charge";
    auto h_charge = df
        .Define("charge", get_RECParticle_int_var(pid,charge), {"REC_Particle_charge", "REC_Particle_pid", "REC_Particle_charge", "REC_Particle_p", "REC_Traj_pass"})
        .Histo1D({charg_title.c_str(), charg_title.c_str(), 500, -2, 2}, "charge");
    Draw1DHist(h_charge.GetPtr(), charg_name.c_str(), 500, -2, 2);
    Save1DHist(h_charge.GetPtr(), charg_title.c_str(), fout);
}