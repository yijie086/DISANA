#include <THnSparse.h>

#include "../DreamAN/DrawHist/DISANAcomparer.h"
#include "../DreamAN/DrawHist/DrawStyle.h"

ROOT::RDF::RNode RejectPi0TwoPhoton(ROOT::RDF::RNode df_);
ROOT::RDF::RNode SelectPi0Event(ROOT::RDF::RNode df);
void CreateCorrectionHistogram4D(ROOT::RDF::RNode df_dvcs_mc, ROOT::RDF::RNode df_pi0_mc, ROOT::RDF::RNode df_dvcs_data, ROOT::RDF::RNode df_pi0_data,
                                 const std::string& out_file_name);

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

void DISANA_Xplotter() {
  bool ComputeBgk_core = false;  // Set to true if you want to compute background
  bool DoBkgCorr = true;       // Set to true if you want to apply background correction

  ROOT::EnableImplicitMT();
  // std::string input_path_from_analysisRun = "/work/clas12/singh/CrossSectionAN/RGA_spring2018_Analysis/fromDVCS_wagon/Inb/";
  // test case
  std::string input_path_from_analysisRun_inb_data = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/testing_outupt/afterFiducialCuts/DVCS_wagon/inb/";
  std::string input_path_from_analysisRun_inb_MC = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/RGA_spring2018_Analysis/pi0Sims/Inb/";
  std::string input_path_from_analysisRun_out_data = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/RGA_spring2018_Analysis/fromDVCS_wagon/Outb/";
  std::string input_path_from_analysisRun_out_MC = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/RGA_spring2018_Analysis/pi0Sims/Outb/";

  // std::string input_path_from_analysisRun = "./../build";
  std::string filename_afterFid_inb_data = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_inb_data.c_str());
  std::string filename_afterFid_inb_MC = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_inb_MC.c_str());
  std::string filename_afterFid_outb_data = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_out_data.c_str());
  std::string filename_afterFid_outb_MC = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_out_MC.c_str());
  // std::string filename_afterFid_afterCorr = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun.c_str());
  float beam_energy = 10.6;

  ROOT::RDF::RNode df_afterFid_inb_data = InitKinematics(filename_afterFid_inb_data, "dfSelected_afterFid", beam_energy);
  ROOT::RDF::RNode df_afterFid_inb_MC = InitKinematics(filename_afterFid_inb_MC, "dfSelected_afterFid", beam_energy);

  ROOT::RDF::RNode df_afterFid_outb_data = InitKinematics(filename_afterFid_outb_data, "dfSelected_afterFid", beam_energy);
  ROOT::RDF::RNode df_afterFid_outb_MC = InitKinematics(filename_afterFid_outb_MC, "dfSelected_afterFid", beam_energy);

  // ROOT::RDF::RNode df_afterFid_afterCorr = InitKinematics(filename_afterFid_afterCorr, "dfSelected_afterFid_afterCorr");
  // input files for the data

  gSystem->Exec("mkdir -p ExclusivityFits");

  DrawStyle fitStyle(0.06, 0.05, 1.0, 1.3);  // You can tweak this
  // DVCS event selection cuts

  // Apply final DVCS cuts
  // inb
  auto df_final_dvcs_inb_data = ApplyFinalDVCSSelections(df_afterFid_inb_data, true);
  auto df_final_dvcsPi_rejected_inb_data = RejectPi0TwoPhoton(df_final_dvcs_inb_data);
  auto df_final_dvcs_inb_MC = ApplyFinalDVCSSelections(df_afterFid_inb_MC, true);
  auto df_final_dvcsPi_rejected_inb_MC = RejectPi0TwoPhoton(df_final_dvcs_inb_MC);

  // outb
  auto df_final_dvcs_outb_data = ApplyFinalDVCSSelections(df_afterFid_outb_data, false);
  auto df_final_dvcsPi_rejected_outb_data = RejectPi0TwoPhoton(df_final_dvcs_outb_data);
  auto df_final_dvcs_outb_MC = ApplyFinalDVCSSelections(df_afterFid_outb_MC, false);
  auto df_final_dvcsPi_rejected_outb_MC = RejectPi0TwoPhoton(df_final_dvcs_outb_MC);

  // pi0 event selection cuts
  // inb
  auto df_final_OnlPi0_inb_data = SelectPi0Event(df_final_dvcs_inb_data);
  auto df_final_OnlPi0_inb_MC = SelectPi0Event(df_final_dvcs_inb_MC);
  // outb
  auto df_final_OnlPi0_outb_data = SelectPi0Event(df_final_dvcs_outb_data);
  auto df_final_OnlPi0_outb_MC = SelectPi0Event(df_final_dvcs_outb_MC);

  // final single photon from pi0 correction factors here
  if (ComputeBgk_core) {
    CreateCorrectionHistogram4D(df_final_dvcsPi_rejected_inb_MC, df_final_OnlPi0_inb_MC, df_final_dvcsPi_rejected_inb_data, df_final_OnlPi0_inb_data, "correction_factorsInb.root");
    CreateCorrectionHistogram4D(df_final_dvcsPi_rejected_outb_MC, df_final_OnlPi0_outb_MC, df_final_dvcsPi_rejected_outb_data, df_final_OnlPi0_outb_data, "correction_factorsOutb.root");
  }

  // for inbending data

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
  //xBins.SetQ2Bins({1.4, 2.0});
  //xBins.SetTBins({0.4, 0.6});
  //xBins.SetXBBins({0.15, 0.25});
  xBins.SetQ2Bins({1.0, 1.5, 2.0, 3.0, 5.0});
  xBins.SetTBins({0.1, 1.0});
  xBins.SetXBBins({0.06, 0.1, 0.16, 0.24, 0.36, 0.48, 0.6});
  comparer.SetXBinsRanges(xBins);

  // comparer.AddModel(df_afterFid_afterCorr, "after Correction", beam_energy);
  // comparer.AddModel(df_afterFid, "Before Exclusivity cuts", beam_energy);

  //comparer.AddModel(df_final_dvcsPi_rejected_inb_data, "Sp18 Inb C", beam_energy, false, "./../build/correction_factorsInb.root");
  comparer.AddModel(df_final_dvcsPi_rejected_inb_data, "Sp18 Inb", beam_energy);
  //comparer.AddModel(df_final_dvcsPi_rejected_outb_data, "Sp18 OutB C", beam_energy, false, "./../build/correction_factorsOutb.root");
  comparer.AddModel(df_final_dvcsPi_rejected_outb_data, "Sp18 OutB", beam_energy);
  //comparer.AddModel(df_final_dvcsPi_rejected_outb_data, "Sp18 OutB", beam_energy);
  //comparer.PlotKinematicComparison();
  //comparer.PlotDVCSKinematicsComparison();
  double luminosity = 1.0;  // Set your desired luminosity here
  double polarisation = 0.85;  // Set your desired polarisation here

  comparer.PlotDISCrossSectionComparison(luminosity);  // argument is Luminosity, polarisation
  comparer.PlotDIS_BSA_Comparison(luminosity, polarisation);         // argument is Luminosity
  // comparer.PlotExclusivityComparisonByDetectorCases(detCuts);

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
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("pho_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("pho_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("recpho_beta",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& beta, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i]) return beta[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_beta", "REC_Particle_pass"})
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
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass) {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 22 && pass[i]) {
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
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass"})

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
          if (pid[i] == 11)
            e++;
          else if (pid[i] == 22 && !daughterPass[i])
            g++;  // photon must NOT be from pi0
          else if (pid[i] == 2212)
            p++;
        }
        return (e >= 1 && g >= 1 && p >= 1);
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
        return (e >= 1 && g >= 1 && p >= 1);
      },
      {"REC_Particle_pid", "REC_Particle_pass", "REC_DaughterParticle_pass"}, "Cut: one good e, γ (not π⁰-like), p");
}
// exclusivity cuts
ROOT::RDF::RNode ApplyFinalDVCSSelections(ROOT::RDF::RNode df, bool inbending) {
  return df
      // 4. Q2 > 1
      .Filter("Q2 > 1.0", "Cut: Q2 > 1 GeV^2")
      //.Filter("t > 1.0", "Cut: t > 1 GeV^2");

      // 5. W > 2
      .Filter("W > 1.8", "Cut: W > 1.8 GeV");

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
  //.Filter("Emiss > 0 && Emiss < 2.0", "Cut: Missing energy")
  //.Filter("PTmiss < 0.35", "Cut: Transverse missing momentum");
  //.Filter("Theta_gamma_gamma < 2.0", "Cut: photon-missing angle")
  //.Filter("DeltaPhi < 25.0", "Cut: Coplanarity");

  // 10. Quality Assurance Cut
  //.Filter("REC_Event_pass == true", "Cut: QA pass");
}
void CreateCorrectionHistogram4D(ROOT::RDF::RNode df_dvcs_mc, ROOT::RDF::RNode df_pi0_mc, ROOT::RDF::RNode df_dvcs_data, ROOT::RDF::RNode df_pi0_data,
                                 const std::string& out_file_name) {
  std::cout << "⏳ Creating 4D correction histograms...\n";

  // Define binning
  const int nQ2Bins = 20, nTBins = 20, nXBBins = 20, nPhiBins = 36;
  std::vector<double> Q2Bins(nQ2Bins + 1), tBins(nTBins + 1), xBBins(nXBBins + 1), phiBins(nPhiBins + 1);

  // Compute min and max from data
  double Q2_min = df_dvcs_data.Min("Q2").GetValue();
  double Q2_max = df_dvcs_data.Max("Q2").GetValue();

  double t_min = df_dvcs_data.Min("t").GetValue();
  double t_max = df_dvcs_data.Max("t").GetValue();

  double xB_min = df_dvcs_data.Min("xB").GetValue();
  double xB_max = df_dvcs_data.Max("xB").GetValue();

  double phi_min = df_dvcs_data.Min("phi").GetValue();
  double phi_max = df_dvcs_data.Max("phi").GetValue();

  // Debug print
  std::cout << "Variable Ranges (from df_dvcs_data):\n";
  std::cout << "  Q²:  " << Q2_min << " to " << Q2_max << "\n";
  std::cout << "  t:   " << t_min << " to " << t_max << "\n";
  std::cout << "  xB:  " << xB_min << " to " << xB_max << "\n";
  std::cout << "  φ:   " << phi_min << "° to " << phi_max << "°\n";

  // Create bin edges based on these ranges
  for (int i = 0; i <= nQ2Bins; ++i) Q2Bins[i] = Q2_min + i * (Q2_max - Q2_min) / nQ2Bins;

  for (int i = 0; i <= nTBins; ++i) tBins[i] = t_min + i * (t_max - t_min) / nTBins;

  for (int i = 0; i <= nXBBins; ++i) xBBins[i] = xB_min + i * (xB_max - xB_min) / nXBBins;

  for (int i = 0; i <= nPhiBins; ++i) phiBins[i] = phi_min + i * (phi_max - phi_min) / nPhiBins;

  int nbins[4] = {nQ2Bins, nTBins, nXBBins, nPhiBins};
  double xmin[4] = {Q2Bins.front(), tBins.front(), xBBins.front(), phiBins.front()};
  double xmax[4] = {Q2Bins.back(), tBins.back(), xBBins.back(), phiBins.back()};

  auto h_dvcs_data = new THnSparseD("h_dvcs_data", "DVCS Data;Q2;t;xB;phi", 4, nbins, xmin, xmax);
  auto h_pi0_data = new THnSparseD("h_pi0_data", "Pi0 Data;Q2;t;xB;phi", 4, nbins, xmin, xmax);
  auto h_dvcs_mc = new THnSparseD("h_dvcs_mc", "DVCS MC;Q2;t;xB;phi", 4, nbins, xmin, xmax);
  auto h_pi0_mc = new THnSparseD("h_pi0_mc", "Pi0 MC;Q2;t;xB;phi", 4, nbins, xmin, xmax);
  auto h_corr = new THnSparseD("h_correction", "Correction;Q2;t;xB;phi", 4, nbins, xmin, xmax);

  // Set variable bin edges
  for (auto* h : {h_dvcs_data, h_pi0_data, h_dvcs_mc, h_pi0_mc, h_corr}) {
    h->SetBinEdges(0, Q2Bins.data());
    h->SetBinEdges(1, tBins.data());
    h->SetBinEdges(2, xBBins.data());
    h->SetBinEdges(3, phiBins.data());
  }

  // Fill helper
  auto fill4D = [](THnSparseD* h, float Q2, float t, float xB, float phi) {
    double vals[4] = {Q2, t, xB, phi * 180.0 / M_PI};  // Convert phi to degrees
    h->Fill(vals);
  };

  auto fillFromRDF = [&](ROOT::RDF::RNode df, THnSparseD* h) {
    auto Q2 = df.Take<double>("Q2");
    auto t = df.Take<double>("t");
    auto xB = df.Take<double>("xB");
    auto phi = df.Take<double>("phi");
    for (size_t i = 0; i < Q2->size(); ++i) fill4D(h, (*Q2)[i], (*t)[i], (*xB)[i], (*phi)[i]);
  };

  fillFromRDF(df_dvcs_data, h_dvcs_data);
  fillFromRDF(df_pi0_data, h_pi0_data);
  fillFromRDF(df_dvcs_mc, h_dvcs_mc);
  fillFromRDF(df_pi0_mc, h_pi0_mc);

  // Correction factor: 1 - (N_pi0_MC / N_dvcs_MC) * (N_pi0_data / N_dvcs_data)
  int nQ2BinsEff = h_corr->GetAxis(0)->GetNbins();
  int nTBinsEff = h_corr->GetAxis(1)->GetNbins();
  int nXBBinsEff = h_corr->GetAxis(2)->GetNbins();
  int nPhiBinsEff = h_corr->GetAxis(3)->GetNbins();

  std::cout << "Entries in df_dvcs_data: " << df_dvcs_data.Count().GetValue() << std::endl;
  std::cout << "Entries in df_pi0_data: " << df_pi0_data.Count().GetValue() << std::endl;
  std::cout << "Entries in df_dvcs_mc:   " << df_dvcs_mc.Count().GetValue() << std::endl;
  std::cout << "Entries in df_pi0_mc:    " << df_pi0_mc.Count().GetValue() << std::endl;

  for (int iq2 = 1; iq2 <= nQ2BinsEff; ++iq2) {
    for (int it = 1; it <= nTBinsEff; ++it) {
      for (int ixb = 1; ixb <= nXBBinsEff; ++ixb) {
        for (int iphi = 1; iphi <= nPhiBinsEff; ++iphi) {
          int indices[4] = {iq2, it, ixb, iphi};

          double N_dvcs_data = h_dvcs_data->GetBinContent(indices);
          double N_pi0_data = h_pi0_data->GetBinContent(indices);
          double N_dvcs_mc = h_dvcs_mc->GetBinContent(indices);
          double N_pi0_mc = h_pi0_mc->GetBinContent(indices);

          double factor = 1.0;
          if (N_dvcs_data > 0 && N_pi0_mc > 0&& N_dvcs_mc > 0 && N_pi0_data > 0) {
            factor = 1.0 - (N_dvcs_mc / N_pi0_mc) * (N_pi0_data / N_dvcs_data);
           // std::cout<<"N_dvcs_data: "<<N_dvcs_data<<" N_pi0_data: "<<N_pi0_data<<" N_dvcs_mc: "<<N_dvcs_mc<<" N_pi0_mc: "<<N_pi0_mc<<std::endl;
            //std::cout << "Q2: " << h_corr->GetAxis(0)->GetBinLowEdge(iq2) << " - " << h_corr->GetAxis(0)->GetBinUpEdge(iq2) << ", t: " << h_corr->GetAxis(1)->GetBinLowEdge(it)
            ///          << " - " << h_corr->GetAxis(1)->GetBinUpEdge(it) << ", xB: " << h_corr->GetAxis(2)->GetBinLowEdge(ixb) << " - " << h_corr->GetAxis(2)->GetBinUpEdge(ixb)
             //         << ", phi: " << h_corr->GetAxis(3)->GetBinLowEdge(iphi) << " - " << h_corr->GetAxis(3)->GetBinUpEdge(iphi) << " => Factor: " << factor << "\n";
          }
          if (factor < 0) {
            std::cout << "Warning: Negative correction factor at Q2 bin " << iq2 << ", t bin " << it << ", xB bin " << ixb << ", phi bin " << iphi
                      << ". Setting to 1.\n";
            factor = 1.0;
          }
          h_corr->SetBinContent(indices, factor);
        }
      }
    }
  }

  // Save output
  TFile fout(out_file_name.c_str(), "RECREATE");
  h_dvcs_mc->Write();
  h_pi0_mc->Write();
  h_dvcs_data->Write();
  h_pi0_data->Write();
  h_corr->Write();
  fout.Close();

  std::cout << "✅ Correction histogram written to: " << out_file_name << "\n";
}