#include <THnSparse.h>

#include "../DreamAN/DrawHist/DISANAMath.h"
#include "../DreamAN/DrawHist/DISANAcomparer.h"
#include "../DreamAN/DrawHist/DrawStyle.h"

// ROOT::RDF::RNode RejectPi0TwoPhoton(ROOT::RDF::RNode df_);
ROOT::RDF::RNode SelectPhiEvent(ROOT::RDF::RNode df);

ROOT::RDF::RNode ApplyFinalPhiSelections(ROOT::RDF::RNode df, bool inbending);

ROOT::RDF::RNode InitKinematics(const std::string& filename_ = "", const std::string& treename_ = "", float beam_energy = 0);

ROOT::RDF::RNode ApplyFinalPhiSelections_NoMass(ROOT::RDF::RNode df, bool inbending);

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
DrawStyle csStyle(0.05, 0.05, .95, 1.1, 42, 5, 510, 0.12, 0.03, 0.12, 0.02);                                        // For Cross-Sections
DrawStyle bsaStyle(0.06, 0.045, .8, .8,42, 5, 510, 0.15, 0.07, 0.16, 0.06);                                      // For BSA

// for exclusivity plots
std::vector<std::pair<std::string, std::string>> detCuts = {
    {"pro_det_region == 2", "CD"},
    {"pro_det_region == 1", "FD"},
};

template <typename Method>
ROOT::RDF::RNode define_DISCAT(ROOT::RDF::RNode node, const std::string& name, const Method method, float beam_energy) {
  return node.Define(name,
                     [method, beam_energy](double recel_p, double recel_theta, double recel_phi, double recpro_p, double recpro_theta, double recpro_phi, double reckMinus_p,
                                           double reckMinus_theta, double reckMinus_phi, double reckPlus_p, double reckPlus_theta, double reckPlus_phi) {
                       return (DISANAMath(beam_energy, recel_p, recel_theta, recel_phi, recpro_p, recpro_theta, recpro_phi, reckMinus_p, reckMinus_theta, reckMinus_phi, reckPlus_p,
                                          reckPlus_theta, reckPlus_phi).*
                               method)();
                     },
                     {"recel_p", "recel_theta", "recel_phi", "recpro_p", "recpro_theta", "recpro_phi", "reckMinus_p", "reckMinus_theta", "reckMinus_phi", "reckPlus_p",
                      "reckPlus_theta", "reckPlus_phi"});
}
//// Phi Analysis Plotter
void DISANA_PhiAnalysisPlotter() {
  bool ComputeBgk_core = false;  // Set to true if you want to compute background
  bool DoBkgCorr = true;         // Set to true if you want to apply background correction

  ROOT::EnableImplicitMT();


  std::string input_path_from_analysisRun_SP18inb_data = "./../data_processed/spring2018/inb/DVKpKm_wagon/";
  std::string input_path_from_analysisRun_SP18outb_data = "./../data_processed/spring2018/outb/DVKpKm_wagon/";

  std::string input_path_from_analysisRun_Fall18inb_data = "./../data_processed/fall2018/inb/DVKpKm_wagon/";
  std::string input_path_from_analysisRun_Fall18outb_data = "./../data_processed/fall2018/outb/DVKpKm_wagon/";

  std::string input_path_from_analysisRun_SP19inb_data = "./../data_processed/spring2019/inb/DVKpKm_wagon/";

 

  std::string filename_afterFid_SP18inb_data = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_SP18inb_data.c_str());
  std::string filename_afterFid_SP18outb_data = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_SP18outb_data.c_str());

  std::string filename_afterFid_Fall18inb_data = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_Fall18inb_data.c_str());
  std::string filename_afterFid_Fall18outb_data = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_Fall18outb_data.c_str());

  std::string filename_afterFid_SP19inb_data = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_SP19inb_data.c_str());

  // std::string filename_afterFid_7546_MC = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_7546_MC.c_str());
  float beam_energy_sp2018 =10.5940;
  float beam_energy_fall2018 =10.6000;
  float beam_energy_sp2019 =10.1998;
 
  ROOT::RDF::RNode df_afterFid_sp18inb_data = InitKinematics(filename_afterFid_SP18inb_data, "dfSelected_afterFid", beam_energy_sp2018);
  ROOT::RDF::RNode df_afterFid_sp18outb_data = InitKinematics(filename_afterFid_SP18outb_data, "dfSelected_afterFid", beam_energy_sp2018);
  ROOT::RDF::RNode df_afterFid_fall18inb_data = InitKinematics(filename_afterFid_Fall18inb_data, "dfSelected_afterFid", beam_energy_fall2018);
  ROOT::RDF::RNode df_afterFid_fall18outb_data = InitKinematics(filename_afterFid_Fall18outb_data, "dfSelected_afterFid", beam_energy_fall2018);
  ROOT::RDF::RNode df_afterFid_sp19inb_data = InitKinematics(filename_afterFid_SP19inb_data, "dfSelected_afterFid", beam_energy_sp2019);

  // ROOT::RDF::RNode df_afterFid_7546_data_mc = InitKinematics(filename_afterFid_7546_data_mc, "dfSelected_afterFid_afterCorr", beam_energy);
  // ROOT::RDF::RNode df_afterFid_7546_MC = InitKinematics(filename_afterFid_7546_MC, "dfSelected_afterFid_afterCorr", beam_energy);

  std::cout << "Applying further cuts and plottings" << std::endl;

  DrawStyle fitStyle(0.06, 0.05, 1.0, 1.3);  // You can tweak this
  // DVCS event selection cuts

  // Apply final DVCS cuts
  auto df_afterFid_sp18inb_with_all_data = ApplyFinalPhiSelections(df_afterFid_sp18inb_data, true);
  auto df_afterFid_sp18outb_with_all_data = ApplyFinalPhiSelections(df_afterFid_sp18outb_data, true);
  auto df_afterFid_sp18inb_with_phi_data = SelectPhiEvent(df_afterFid_sp18inb_with_all_data);
  auto df_afterFid_sp18outb_with_phi_data = SelectPhiEvent(df_afterFid_sp18outb_with_all_data);

  auto df_afterFid_fall18inb_with_all_data = ApplyFinalPhiSelections(df_afterFid_fall18inb_data, true);
  auto df_afterFid_fall18outb_with_all_data = ApplyFinalPhiSelections(df_afterFid_fall18outb_data, true);
  auto df_afterFid_fall18inb_with_phi_data = SelectPhiEvent(df_afterFid_fall18inb_with_all_data);
  auto df_afterFid_fall18outb_with_phi_data = SelectPhiEvent(df_afterFid_fall18outb_with_all_data);

  auto df_afterFid_sp19inb_with_all_data = ApplyFinalPhiSelections(df_afterFid_sp19inb_data, true);
  auto df_afterFid_sp19inb_with_phi_data = SelectPhiEvent(df_afterFid_sp19inb_with_all_data);


  // auto df_final_dvcsPi_rejected_7546_data = RejectPi0TwoPhoton(df_final_dvcs_7546_data);

  DISANAcomparer comparer;
  comparer.SetOutputDir("./");
  comparer.SetKinStyle(KinStyle);
  comparer.SetDVCSStyle(dvcsStyle);
  comparer.SetCrossSectionStyle(csStyle);
  comparer.SetBSAStyle(bsaStyle);

  comparer.PlotIndividual(false);
  /// bins for cross-section plots
  BinManager xBins;
  xBins.SetQ2Bins({0.4, 1.3, 1.86, 2.71, 8.35});
  //xBins.SetQ2Bins({0.4, 8.35});
  xBins.SetXBBins({0, 0.99});
  xBins.SetTBins({0.2, .5,  0.8,  1.2, 8.0});// 5.4, 5.8, 6.2, 6.6, 7.0, 7.4, 7.8, 8.2, 8.6, 9.0, 9.4, 9.8});
 
  comparer.SetXBinsRanges(xBins);
  comparer.UseFittedPhiYields(true);
   
   comparer.AddModelPhi(df_afterFid_fall18outb_with_phi_data, "Fall18 outb", beam_energy_fall2018);
   comparer.AddModelPhi(df_afterFid_fall18inb_with_phi_data, "Fall18 inb", beam_energy_fall2018);

   comparer.AddModelPhi(df_afterFid_sp18outb_with_phi_data, "Sp18 outb", beam_energy_sp2018);
   comparer.AddModelPhi(df_afterFid_sp18inb_with_phi_data, "Sp18 inb", beam_energy_sp2018);

   comparer.AddModelPhi(df_afterFid_sp19inb_with_phi_data, "Sp19 inb", beam_energy_sp2019);

  double luminosity = 24.3065 * pow(10, 6);  // Set your desired luminosity here nb^-1
  double luminosity_rga_fall18 = 1.0;        // 5.47*pow(10,40);  // Set your desired luminosity here nb^-1
  double polarisation = 0.85;                // Set your desired polarisation here
  double branching = 0.49;                // Set your desired polarisation here

  
  comparer.PlotPhiElectroProKinematicsComparison();
  comparer.PlotKinematicComparison_phiAna();
  comparer.PlotPhiAnaExclusivityComparisonByDetectorCases(detCuts);
  comparer.PlotPhiInvMassPerBin_AllModels("PhiInvMassFits", /*nBins*/30, /*mMin*/0.9974, /*mMax*/1.2, /*constrainSigma*/true, luminosity_rga_fall18,branching);
  comparer.PlotPhiDSigmaDt_FromCache();

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
             .Define("kMinus_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("kMinus_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("kMinus_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == -321 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})
             .Define("kPlus_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && trackpass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass"})
             .Define("kPlus_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && trackpass[i]) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass"})
             .Define("kPlus_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 321 && trackpass[i]) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass"})

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
             .Filter([](float ex, float kMinusx, float kPlusx, float px) { return ex != -999 && kMinusx != -999 && kPlusx != -999 && px != -999; },
                     {"ele_px", "kMinus_px", "kPlus_px", "pro_px"})
             .Define("recel_p", MomentumFunc, {"ele_px", "ele_py", "ele_pz"})
             .Define("recel_theta", ThetaFunc, {"ele_px", "ele_py", "ele_pz"})
             .Define("recel_phi", PhiFunc, {"ele_px", "ele_py"})
             .Define("reckMinus_p", MomentumFunc, {"kMinus_px", "kMinus_py", "kMinus_pz"})
             .Define("reckMinus_theta", ThetaFunc, {"kMinus_px", "kMinus_py", "kMinus_pz"})
             .Define("reckMinus_phi", PhiFunc, {"kMinus_px", "kMinus_py"})
             .Define("reckPlus_p", MomentumFunc, {"kPlus_px", "kPlus_py", "kPlus_pz"})
             .Define("reckPlus_theta", ThetaFunc, {"kPlus_px", "kPlus_py", "kPlus_pz"})
             .Define("reckPlus_phi", PhiFunc, {"kPlus_px", "kPlus_py"})
             .Define("recpro_p", MomentumFunc, {"pro_px", "pro_py", "pro_pz"})
             .Define("recpro_theta", ThetaFunc, {"pro_px", "pro_py", "pro_pz"})
             .Define("recpro_phi", PhiFunc, {"pro_px", "pro_py"})
             .Define("kMinus_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass) {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == -321 && pass[i]) {
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
             .Define("kPlus_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<short>& status, const ROOT::VecOps::RVec<bool>& pass) {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 321 && pass[i]) {
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
                     {"REC_Particle_pid", "REC_Particle_status", "REC_Particle_pass"})
             .Define("invMass_KpKm",
                     [](float px1, float py1, float pz1, float px2, float py2, float pz2) {
                       constexpr float mK = 0.493677;  // Kaon mass in GeV/c²
                       float E1 = std::sqrt(px1 * px1 + py1 * py1 + pz1 * pz1 + mK * mK);
                       float E2 = std::sqrt(px2 * px2 + py2 * py2 + pz2 * pz2 + mK * mK);
                       float px = px1 + px2;
                       float py = py1 + py2;
                       float pz = pz1 + pz2;
                       float E = E1 + E2;
                       return std::sqrt(E * E - (px * px + py * py + pz * pz));
                     },
                     {"kPlus_px", "kPlus_py", "kPlus_pz", "kMinus_px", "kMinus_py", "kMinus_pz"});

         * df_ = define_DISCAT(*df_, "Q2", &DISANAMath::GetQ2, beam_energy);
  *df_ = define_DISCAT(*df_, "xB", &DISANAMath::GetxB, beam_energy);
  *df_ = define_DISCAT(*df_, "t", &DISANAMath::GetT, beam_energy);
  *df_ = define_DISCAT(*df_, "phi", &DISANAMath::GetPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "W", &DISANAMath::GetW, beam_energy);
  *df_ = define_DISCAT(*df_, "nu", &DISANAMath::GetNu, beam_energy);
  *df_ = define_DISCAT(*df_, "y", &DISANAMath::Gety, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_ep", &DISANAMath::GetMx2_ep, beam_energy);
  *df_ = define_DISCAT(*df_, "Emiss", &DISANAMath::GetEmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "PTmiss", &DISANAMath::GetPTmiss, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKpKm", &DISANAMath::GetMx2_epKpKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_eKpKm", &DISANAMath::GetMx2_eKpKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKp", &DISANAMath::GetMx2_epKp, beam_energy);
  *df_ = define_DISCAT(*df_, "Mx2_epKm", &DISANAMath::GetMx2_epKm, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_e_phimeson", &DISANAMath::GetTheta_e_phimeson, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaE", &DISANAMath::GetDeltaE, beam_energy);
  *df_ = define_DISCAT(*df_, "DeltaPhi", &DISANAMath::GetDeltaPhi, beam_energy);
  *df_ = define_DISCAT(*df_, "Theta_g_phimeson", &DISANAMath::GetTheta_g_phimeson, beam_energy);

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

ROOT::RDF::RNode ApplyFinalPhiSelections_NoMass(ROOT::RDF::RNode df, bool inbending) {
  return df
    .Filter("reckPlus_p<3.5",  "Cut: reckPlus_p < 3.5")
    .Filter("reckMinus_p<3.5", "Cut: reckMinus_p < 3.5")
    .Filter("W > 2.0",         "Cut: W > 2 GeV")
    .Filter("Emiss < 10",     "Cut: Missing energy");
}

// pi-0 event selection cuts for single photon contaminations
ROOT::RDF::RNode SelectPhiEvent(ROOT::RDF::RNode df_) {
  return df_.Filter(
      [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<bool>& pass, const ROOT::VecOps::RVec<bool>& daughterPass) {
        int e = 0, km = 0, kp = 0, p = 0;
        bool hasPhiDaughter = true;

        for (size_t i = 0; i < pid.size(); ++i) {
          if (!pass[i]) continue;

          if (pid[i] == 11) {
            e++;
          } else if (pid[i] == 321) {
            kp++;
            hasPhiDaughter = hasPhiDaughter /*|| daughterPass[i]/*/;  // check if kaon is from phi
          } else if (pid[i] == -321) {
            km++;
            hasPhiDaughter = hasPhiDaughter /*|| daughterPass[i]*/;  // check if kaon is from phi
          } else if (pid[i] == 2212) {
            p++;
          }
        }

        return (e == 1 && kp >= 1 && km >= 1 && p >= 1 && hasPhiDaughter);
      },
      {"REC_Particle_pid", "REC_Particle_pass", "REC_DaughterParticle_pass"}, "Cut: 1 e⁻, ≥1 K⁺, ≥1 K⁻, 1 proton, ≥1 kaon from phi");
}

// exclusivity cuts
ROOT::RDF::RNode ApplyFinalPhiSelections(ROOT::RDF::RNode df, bool inbending) {
  return df
      // 4. Q2 > 1
      //.Filter("Q2 > 1.0", "Cut: Q2 > 1 GeV^2")
      //.Filter("t < 1.0", "Cut: t > 1 GeV^2")
      .Filter("reckPlus_p<3.5", "Cut: reckPlus_p < 0.6")
      .Filter("reckMinus_p<3.5", "Cut: reckMinus_p < 0.6")

      // 5. W > 2
      .Filter("W > 2.0", "Cut: W > 1.8 GeV")
      //.Filter("invMass_KpKm > 1.02-0.0138 && invMass_KpKm < 1.02+0.0138", "Cut: Invariant mass of K⁺K⁻ around ϕ") //systemtaics +-0.3
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
      //.Filter("Mx2_ep > 0.8 && Mx2_ep < 1.5", "Cut: MM^2(ep) in 3sigma");
      .Filter("Emiss < 0.200", "Cut: Missing energy");
}
