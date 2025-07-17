#include <THnSparse.h>

#include "../DreamAN/DrawHist/DISANAcomparer.h"
#include "../DreamAN/DrawHist/DrawStyle.h"
#include "../DreamAN/DrawHist/DISANAMath.h"

ROOT::RDF::RNode RejectPi0TwoPhoton(ROOT::RDF::RNode df_);
ROOT::RDF::RNode SelectPi0Event(ROOT::RDF::RNode df);

ROOT::RDF::RNode ApplyFinalDVCSSelections(ROOT::RDF::RNode df, bool inbending);

ROOT::RDF::RNode InitKinematics(const std::string& filename_ = "", const std::string& treename_ = "", float beam_energy = 0);

static double MomentumFunc(float px, float py, float pz) { return std::sqrt(px * px + py * py + pz * pz); }
static double ThetaFunc(float px, float py, float pz) { return std::acos(pz / std::sqrt(px * px + py * py + pz * pz)); }
static double PhiFunc(float px, float py) {
  double phi = std::atan2(py, px);
  return phi < 0 ? phi + 2 * M_PI : phi;
}

/// styling plots
// double double titleSize = 0.05, double labelSize = 0.04,double xTitleOffset = 1.1, double yTitleOffset = 1.6, int font = 42, int maxDigits = 5, int nDivisions = 510, double
// leftMargin = 0.16, double rightMargin = 0.07, double bottomMargin = 0.13, double topMargin = 0.06
DrawStyle KinStyle(0.07, 0.06, 0.9, 1.2);                                       // For Kin plots
DrawStyle dvcsStyle(0.06, 0.06, 1.2, 1.4, 42, 5, 510, 0.14, 0.07, 0.13, 0.06);  // For DVCS plots
DrawStyle csStyle(0.05, 0.04, 1.0, 1.3);                                        // For Cross-Sections
DrawStyle bsaStyle(0.06, 0.045, 1.0, 1.2);                                      // For BSA

// for exclusivity plots
std::vector<std::pair<std::string, std::string>> detCuts = {
    {"pho_det_region == 0 && pro_det_region == 2", "FT-CD"},
    {"pho_det_region == 1 && pro_det_region == 2", "FD-CD"},
    {"pho_det_region == 1 && pro_det_region == 1", "FD-FD"},
};

template <typename Method>
ROOT::RDF::RNode define_DISCAT(ROOT::RDF::RNode node, const std::string& name, const Method method, float beam_energy) {
  return node.Define(name,
                     [method, beam_energy](double recel_p, double recel_theta, double recel_phi, double recpro_p, double recpro_theta, double recpro_phi, double recpho_p,
                                           double recpho_theta, double recpho_phi) {
                       return (DISANAMath(beam_energy, recel_p, recel_theta, recel_phi, recpro_p, recpro_theta, recpro_phi, recpho_p, recpho_theta, recpho_phi).*method)();
                     },
                     {"recel_p", "recel_theta", "recel_phi", "recpro_p", "recpro_theta", "recpro_phi", "recpho_p", "recpho_theta", "recpho_phi"});
}

void DISANA_XplotterOut() {
  bool ComputeBgk_core = false;  // Set to true if you want to compute background
  bool DoBkgCorr = true;       // Set to true if you want to apply background correction

  ROOT::EnableImplicitMT();
 
  std::string input_path_from_analysisRun_rgasp18outb_data = "./../build/rgasp18outdatanoSF/";
  std::string input_path_from_analysisRun_rgasp18outb_data_mc = "./../build/";
  std::string input_path_from_analysisRun_rgasp18outb_MC = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/RGA_spring2018_Analysis/pi0Sims/Inb/";

  std::string filename_afterFid_rgasp18outb_data = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_rgasp18outb_data.c_str());
  std::string filename_afterFid_rgasp18outb_data_mc = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_rgasp18outb_data_mc.c_str());
  std::string filename_afterFid_rgasp18outb_MC = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_rgasp18outb_MC.c_str());
  float beam_energy = 10.6;

  ROOT::RDF::RNode df_afterFid_rgasp18outb_data = InitKinematics(filename_afterFid_rgasp18outb_data, "dfSelected_afterFid", beam_energy);
  ROOT::RDF::RNode df_afterFid_rgasp18outb_data_mc = InitKinematics(filename_afterFid_rgasp18outb_data_mc, "dfSelected_afterFid", beam_energy);
  ROOT::RDF::RNode df_afterFid_rgasp18outb_MC = InitKinematics(filename_afterFid_rgasp18outb_MC, "dfSelected_afterFid", beam_energy);



  DrawStyle fitStyle(0.06, 0.05, 1.0, 1.3);  // You can tweak this
  // DVCS event selection cuts

  // Apply final DVCS cuts
  auto df_final_dvcs_rgasp18outb_data = ApplyFinalDVCSSelections(df_afterFid_rgasp18outb_data, true);
  auto df_final_dvcsPi_rejected_rgasp18outb_data = RejectPi0TwoPhoton(df_final_dvcs_rgasp18outb_data);
  auto df_final_dvcs_rgasp18outb_data_mc = ApplyFinalDVCSSelections(df_afterFid_rgasp18outb_data_mc, false);
  auto df_final_dvcsPi_rejected_rgasp18outb_data_mc = RejectPi0TwoPhoton(df_final_dvcs_rgasp18outb_data_mc);  
  auto df_final_dvcs_rgasp18outb_MC = ApplyFinalDVCSSelections(df_afterFid_rgasp18outb_MC, false);
  auto df_final_dvcsPi_rejected_rgasp18outb_MC = RejectPi0TwoPhoton(df_final_dvcs_rgasp18outb_MC);

  auto df_final_OnlPi0_rgasp18outb_data = SelectPi0Event(df_final_dvcs_rgasp18outb_data);
  auto df_final_OnlPi0_rgasp18outb_data_mc = SelectPi0Event(df_final_dvcs_rgasp18outb_data_mc);
  auto df_final_OnlPi0_rgasp18outb_MC = SelectPi0Event(df_final_dvcs_rgasp18outb_MC);

  DISANAcomparer comparer;
  comparer.SetOutputDir("./");
  comparer.SetKinStyle(KinStyle);
  comparer.SetDVCSStyle(dvcsStyle);
  comparer.SetCrossSectionStyle(csStyle);
  comparer.SetBSAStyle(bsaStyle);

  comparer.PlotIndividual(false);
  /// bins for cross-section plots
  BinManager xBins;
  // xBins.SetQ2Bins({.11,1.3,1.6,2.1,2.8,3.6,8.0});
  // xBins.SetTBins({0.0, 1.2});
  // xBins.SetXBBins({0.0, 0.08,.1,.14,.18,.23,.3,.39,.50});
  xBins.SetQ2Bins({1.0, 10.0});
  //xBins.SetTBins({0.0, 1.0});
  xBins.SetXBBins({0.0, 1.0});
  //xBins.SetQ2Bins({1.0, 1.2, 1.456, 1.912, 2.51, 3.295, 4.326, 5.761});
  xBins.SetTBins({0.11, 1.0});
  //xBins.SetXBBins({0.062, 0.09, 0.118, 0.155, 0.204, 0.268, 0.357, 0.446, 0.581});
  //xBins.SetQ2Bins({1.0, 5.0});
  //xBins.SetTBins({0.0, 1.0});
  //xBins.SetXBBins({0.0, 1.0});
  comparer.SetXBinsRanges(xBins);

  comparer.AddModelwithPi0Corr(df_final_dvcsPi_rejected_rgasp18outb_data_mc,df_final_OnlPi0_rgasp18outb_data_mc,df_final_dvcsPi_rejected_rgasp18outb_MC,df_final_OnlPi0_rgasp18outb_MC, "RGA Sp18 Outb mc ", beam_energy, false);
  comparer.AddModelwithPi0Corr(df_final_dvcsPi_rejected_rgasp18outb_data,df_final_OnlPi0_rgasp18outb_data,df_final_dvcsPi_rejected_rgasp18outb_MC,df_final_OnlPi0_rgasp18outb_MC, "RGA Sp18 Outb", beam_energy, false);
  
  double luminosity = 1.0;  // Set your desired luminosity here
  double polarisation = 0.85;  // Set your desired polarisation here

  comparer.PlotKinematicComparison();
  //comparer.PlotDVCSKinematicsComparison();
  //comparer.PlotDISCrossSectionComparison(luminosity);  // argument is Luminosity, polarisation
  //comparer.PlotDIS_BSA_Comparison(luminosity, polarisation);         // argument is Luminosity
  //comparer.PlotDIS_Pi0CorrComparison();
  //comparer.PlotExclusivityComparisonByDetectorCases(detCuts);

  gApplication->Terminate(0);
}

ROOT::RDF::RNode InitKinematics(const std::string& filename_, const std::string& treename_, float beam_energy) {
  ROOT::RDataFrame rdf(treename_, filename_);
  auto df_ = std::make_unique<ROOT::RDF::RNode>(rdf);
  *df_ = df_->Define("ele_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("ele_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("ele_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("pho_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i] && maxEpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass", "REC_Photon_MaxE"})
             .Define("pho_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i] && maxEpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass", "REC_Photon_MaxE"})
             .Define("pho_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i] && maxEpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass", "REC_Photon_MaxE"})
             .Define("recpho_beta",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& beta, const ROOT::VecOps::RVec<bool>& trackpass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i] && maxEpass[i]) return beta[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_beta", "REC_Particle_pass", "REC_Photon_MaxE"})
             .Define("pro_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("pro_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("pro_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Filter([](float ex, float gx, float px) { return ex != -999 && gx != -999 && px != -999; }, {"ele_px", "pho_px", "pro_px"})
             .Define("recel_p", MomentumFunc, {"ele_px", "ele_py", "ele_pz"})
             .Define("recel_theta", ThetaFunc, {"ele_px", "ele_py", "ele_pz"})
             .Define("recel_phi", PhiFunc, {"ele_px", "ele_py"})
             .Define("recpho_p", MomentumFunc, {"pho_px", "pho_py", "pho_pz"})
             .Define("recpho_theta", ThetaFunc, {"pho_px", "pho_py", "pho_pz"})
             .Define("recpho_phi", PhiFunc, {"pho_px", "pho_py"})
             .Define("recpro_p", MomentumFunc, {"pro_px", "pro_py", "pro_pz"})
             .Define("recpro_theta", ThetaFunc, {"pro_px", "pro_py", "pro_pz"})
             .Define("recpro_phi", PhiFunc, {"pro_px", "pro_py"})
             .Define("pho_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 22 && pass[i] && maxEpass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;  // Unknown/Other
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass", "REC_Photon_MaxE"})

             .Define("pro_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass) {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 2212 && pass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT (probably rare for protons)
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass"})
             .Define("ele_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass) {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 11 && pass[i]) {
                           int abs_status = std::abs(status[i]);
                           if (abs_status >= 1000 && abs_status < 2000)
                             return 0;  // FT (probably rare for protons)
                           else if (abs_status >= 2000 && abs_status < 3000)
                             return 1;  // FD
                           else if (abs_status >= 4000 && abs_status < 5000)
                             return 2;  // CD
                           else
                             return -1;
                         }
                       }
                       return -1;
                     },
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass"});

  *df_ = define_DISCAT(*df_, "Q2", &DISANAMath::GetQ2, beam_energy);
  *df_ = define_DISCAT(*df_, "xB", &DISANAMath::GetxB, beam_energy);
  *df_ = define_DISCAT(*df_, "t", &DISANAMath::GetT, beam_energy);
  *df_ = define_DISCAT(*df_, "phi", &DISANAMath::GetPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "W", &DISANAMath::GetW, beam_energy);
  *df_ = define_DISCAT(*df_, "nu", &DISANAMath::GetNu, beam_energy);
  *df_ = define_DISCAT(*df_, "y", &DISANAMath::Gety, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_ep", &DISANAMath::GetMx2_ep, beam_energy);
  *df_ = define_DISCAT(*df_, "Emiss", &DISANAMath::GetEmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "PTmiss", &DISANAMath::GetPTmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epg", &DISANAMath::GetMx2_epg, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_eg", &DISANAMath::GetMx2_egamma, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_e_gamma", &DISANAMath::GetTheta_e_gamma, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaE", &DISANAMath::GetDeltaE, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaPhi", &DISANAMath::GetDeltaPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_gamma_gamma", &DISANAMath::GetTheta_gamma_gamma, beam_energy);

  return *df_;
}

//
ROOT::RDF::RNode RejectPi0TwoPhoton(ROOT::RDF::RNode df_) {
  return df_.Filter(
      [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<bool>& pass, const ROOT::VecOps::RVec<bool>& daughterPass) {
        int e = 0, g = 0, p = 0;
        for (size_t i = 0; i < pid.size(); ++i) {
          if (!pass[i]) continue;
          if (daughterPass[i]) return false;  // reject if any daughter particle is a pi0
          if (pid[i] == 11)
            e++;
          else if (pid[i] == 22)
            g++;  // photon must NOT be from pi0
          else if (pid[i] == 2212)
            p++;
        }
        return (e == 1 && g == 1 && p == 1);
      },
      {"REC_Particle_pid", "REC_Particle_pass", "REC_DaughterParticle_pass"}, "Cut: one good e, γ (not π⁰-like), p");
}
// pi-0 event selection cuts for single photon contaminations
ROOT::RDF::RNode SelectPi0Event(ROOT::RDF::RNode df_) {
  return df_.Filter(
      [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<bool>& pass, const ROOT::VecOps::RVec<bool>& daughterPass) {
        int e = 0, g = 0, p = 0;
        for (size_t i = 0; i < pid.size(); ++i) {
          if (!pass[i]) continue;
          if (pid[i] == 11)
            e++;
          else if (pid[i] == 22 && daughterPass[i])
            g++;  // photon must NOT be from pi0
          else if (pid[i] == 2212)
            p++;
        }
        return (e == 1 && g == 2 && p == 1);
      },
      {"REC_Particle_pid", "REC_Particle_pass", "REC_DaughterParticle_pass"}, "Cut: one good e, γ (not π⁰-like), p");
}
// exclusivity cuts
ROOT::RDF::RNode ApplyFinalDVCSSelections(ROOT::RDF::RNode df, bool inbending) {
  return df
      // 4. Q2 > 1
      .Filter("Q2 > 1.0", "Cut: Q2 > 1 GeV^2")
      .Filter("t < 1.0", "Cut: t > 1 GeV^2")
      //.Filter("recel_p > 6.0", "Cut: recel_p > 0.6")

      // 5. W > 2
      .Filter("W > 2.0", "Cut: W > 1.8 GeV")
      //.Filter("phi > 100.0 && phi < 300 ", "Cut: phi")


  // 6. Electron and photon in different sectors
  //.Filter("ele_sector != pho_sector", "Cut: e and gamma in different sectors")

  // 7. Proton and photon in different sectors if ECAL hit
  //.Filter(
  //    [](int p_sec, int g_sec, bool has_ecal_hit) {
  //     return (p_sec != g_sec) || !has_ecal_hit;
  //   },
  //  {"pro_sector", "pho_sector", "pro_has_ECAL_hit"},
  //  "Cut: p and gamma different sector if ECAL hit")
  //
  // 9. 3σ exclusivity cuts
  //.Filter("Mx2_ep > -1.5 && Mx2_ep < 1.5", "Cut: MM^2(ep) in 3sigma")
  .Filter("Emiss < 1.0", "Cut: Missing energy")
  .Filter("PTmiss < 0.2", "Cut: Transverse missing momentum")
  .Filter("Theta_e_gamma > 5 ", "Cut: Theta_e_gamma")
  .Filter("Theta_gamma_gamma < 2.0", "Cut: photon-missing angle");
  //.Filter("DeltaPhi < 25.0", "Cut: Coplanarity");
  //.Filter("(pho_det_region==0&&pro_det_region==2)||(pho_det_region==1&&pro_det_region==1)||(pho_det_region==1&&pro_det_region==2)", "Cut: three config")
  //.Filter("(pho_det_region==0&&pro_det_region==2&&Mx2_ep<0.25&&Mx2_ep>-0.23)||(pho_det_region==1&&pro_det_region==1&&Mx2_ep<0.33&&Mx2_ep>-0.33)||(pho_det_region==1&&pro_det_region==2&&Mx2_ep<0.25&&Mx2_ep>-0.23)", "Cut: Mx2_ep in 3sigma")
  //.Filter("(pho_det_region==0&&pro_det_region==2&&Emiss<0.37&&Emiss>-0.29)||(pho_det_region==1&&pro_det_region==1&&Emiss<0.64&&Emiss>-0.56)||(pho_det_region==1&&pro_det_region==2&&Emiss<0.72&&Emiss>-0.6)", "Cut: Emiss in 3sigma")
  //.Filter("(pho_det_region==0&&pro_det_region==2&&PTmiss<0.08&&PTmiss>-0.04)||(pho_det_region==1&&pro_det_region==1&&PTmiss<0.18&&PTmiss>-0.06)||(pho_det_region==1&&pro_det_region==2&&PTmiss<0.13&&PTmiss>-0.05)", "Cut: PTmiss in 3sigma")
  //.Filter("(pho_det_region==0&&pro_det_region==2&&Theta_gamma_gamma<1.16&&Theta_gamma_gamma>-0.58)||(pho_det_region==1&&pro_det_region==1&&Theta_gamma_gamma<1.81&&Theta_gamma_gamma>-0.71)||(pho_det_region==1&&pro_det_region==2&&Theta_gamma_gamma<1.26&&Theta_gamma_gamma>-0.6)", "Cut: Theta_gamma_gamma in 3sigma")
  //.Filter("(pho_det_region==0&&pro_det_region==2&&DeltaPhi<3.93&&DeltaPhi>-2.55)||(pho_det_region==1&&pro_det_region==1&&DeltaPhi<11.75&&DeltaPhi>-6.73)||(pho_det_region==1&&pro_det_region==2&&DeltaPhi<7.42&&DeltaPhi>-4.76)", "Cut: DeltaPhi in 3sigma")
  //.Filter("(pho_det_region==0&&pro_det_region==2&&Mx2_epg<0.03&&Mx2_epg>-0.03)||(pho_det_region==1&&pro_det_region==1&&Mx2_epg<0.03&&Mx2_epg>-0.03)||(pho_det_region==1&&pro_det_region==2&&Mx2_epg<0.03&&Mx2_epg>-0.03)", "Cut: Mx2_epg in 3sigma")
  //.Filter("(pho_det_region==0&&pro_det_region==2&&Mx2_eg<1.55&&Mx2_eg>0.35)||(pho_det_region==1&&pro_det_region==1&&Mx2_eg<1.68&&Mx2_eg>0.12)||(pho_det_region==1&&pro_det_region==2&&Mx2_eg<1.99&&Mx2_eg>-0.05)", "Cut: Mx2_eg in 3sigma")
  //.Filter("(pho_det_region==0&&pro_det_region==2&&Theta_e_gamma<27.92&&Theta_e_gamma>5.42)||(pho_det_region==1&&pro_det_region==1&&Theta_e_gamma<44.68&&Theta_e_gamma>29.32)||(pho_det_region==1&&pro_det_region==2&&Theta_e_gamma<36.36&&Theta_e_gamma>10.36)", "Cut: Theta_e_gamma in 3sigma");

  // 10. Quality Assurance Cut
  //.Filter("REC_Event_pass == true", "Cut: QA pass");
}

/*
 *df_ = define_DISCAT(*df_, "Mx2_ep", &DISANAMath::GetMx2_ep, beam_energy);
  *df_ = define_DISCAT(*df_, "Emiss", &DISANAMath::GetEmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "PTmiss", &DISANAMath::GetPTmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epg", &DISANAMath::GetMx2_epg, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_eg", &DISANAMath::GetMx2_egamma, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_e_gamma", &DISANAMath::GetTheta_e_gamma, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaE", &DISANAMath::GetDeltaE, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaPhi", &DISANAMath::GetDeltaPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_gamma_gamma", &DISANAMath::GetTheta_gamma_gamma, beam_energy);*/