#include <THnSparse.h>

#include "../DreamAN/DrawHist/DISANAcomparer.h"
#include "../DreamAN/DrawHist/DrawStyle.h"
#include "../DreamAN/DrawHist/DISANAMath.h"

ROOT::RDF::RNode RejectPi0TwoPhoton(ROOT::RDF::RNode df_);
ROOT::RDF::RNode SelectPi0Event(ROOT::RDF::RNode df);

ROOT::RDF::RNode ApplyFinalDVCSSelections(ROOT::RDF::RNode df, bool inbending);
ROOT::RDF::RNode ApplyFinalGenDVCSSelections(ROOT::RDF::RNode df, bool inbending);

ROOT::RDF::RNode InitKinematics(const std::string& filename_ = "", const std::string& treename_ = "", float beam_energy = 0);
ROOT::RDF::RNode InitGenKinematics(const std::string& filename_ = "", const std::string& treename_ = "", float beam_energy = 0);

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

void DISANA_Xplotter2() {
  bool ComputeBgk_core = false;  // Set to true if you want to compute background
  bool DoBkgCorr = true;       // Set to true if you want to apply background correction

  ROOT::EnableImplicitMT();
 
  //std::string input_path_from_analysisRun_7546_data = "/w/hallb-scshelf2102/clas12/yijie/clas12ana/analysis401/DISANA/build/rgk7546dataCorr/";
  std::string input_path_from_analysisRun_7546_data = "./../build";
  //std::string input_path_from_analysisRun_7546_data_mc = "./../build/rgk7546mcSFCorr";
  std::string input_path_from_analysisRun_7546_pi0MC = "/w/hallb-scshelf2102/clas12/yijie/clas12ana/analysis401/DISANA/build/rgk7546dvpiomcCorr/";
  
  std::string input_path_from_analysisRun_7546_dvcsmc_gen = "/w/hallb-scshelf2102/clas12/yijie/clas12ana/analysis401/DISANA/build/rgk7546dvcsmcAll2000/";
  std::string input_path_from_analysisRun_7546_dvcsmc_rec = "/w/hallb-scshelf2102/clas12/yijie/clas12ana/analysis401/DISANA/build/rgk7546dvcsmcSel2000/";

  std::string input_path_from_analysisRun_7546_dvcsmc_bkg = "/w/hallb-scshelf2102/clas12/yijie/clas12ana/analysis501/DISANA/build/rgkdvcs7546bkg/";
  std::string input_path_from_analysisRun_7546_dvcsmc_nobkg = "/w/hallb-scshelf2102/clas12/yijie/clas12ana/analysis501/DISANA/build/rgkdvcs7546nobkg/";

  std::string input_path_from_analysisRun_7546_dvcsmc_rad = "/w/hallb-scshelf2102/clas12/yijie/clas12ana/analysis501/DISANA/build/radcorr/";
  std::string input_path_from_analysisRun_7546_dvcsmc_norad = "/w/hallb-scshelf2102/clas12/yijie/clas12ana/analysis501/DISANA/build/radcorr/";

  std::string filename_afterFid_7546_data = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_7546_data.c_str());
  //std::string filename_afterFid_7546_data_mc = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_7546_data_mc.c_str());
  std::string filename_afterFid_7546_pi0MC = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_7546_pi0MC.c_str());

  std::string filename_afterFid_7546_dvcsmc_gen = Form("%s/dfSelected.root", input_path_from_analysisRun_7546_dvcsmc_gen.c_str());
  std::string filename_afterFid_7546_dvcsmc_rec = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_7546_dvcsmc_rec.c_str());

  std::string filename_afterFid_7546_dvcsmc_bkg = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_7546_dvcsmc_bkg.c_str());
  std::string filename_afterFid_7546_dvcsmc_nobkg = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_7546_dvcsmc_nobkg.c_str());

  std::string filename_afterFid_7546_dvcsmc_rad = Form("%s/rad100M.root", input_path_from_analysisRun_7546_dvcsmc_rad.c_str());
  std::string filename_afterFid_7546_dvcsmc_norad = Form("%s/norad100M.root", input_path_from_analysisRun_7546_dvcsmc_norad.c_str());
  //filename_afterFid_7546_dvcsmc_rad = "../build/radtest/dfSelected.root";
  //filename_afterFid_7546_dvcsmc_norad = "../build/nortest/dfSelected.root";
  filename_afterFid_7546_dvcsmc_rad = "/w/hallb-scshelf2102/clas12/yijie/clas12ana/analysis701/DISANA/rootmacros/rad.root";
  filename_afterFid_7546_dvcsmc_norad = "/w/hallb-scshelf2102/clas12/yijie/clas12ana/analysis701/DISANA/rootmacros/nor.root";
  

  float beam_energy = 7.546;

  ROOT::RDF::RNode df_afterFid_7546_data = InitKinematics(filename_afterFid_7546_data, "dfSelected_afterFid_afterCorr", beam_energy);
  //ROOT::RDF::RNode df_afterFid_7546_data_mc = InitKinematics(filename_afterFid_7546_data_mc, "dfSelected_afterFid_afterCorr", beam_energy);
  ROOT::RDF::RNode df_afterFid_7546_pi0MC = InitKinematics(filename_afterFid_7546_pi0MC, "dfSelected_afterFid_afterCorr", beam_energy);

  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_gen = InitGenKinematics(filename_afterFid_7546_dvcsmc_gen, "dfSelected", beam_energy);
  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_rec = InitKinematics(filename_afterFid_7546_dvcsmc_rec, "dfSelected_afterFid_afterCorr", beam_energy);

  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_bkg = InitKinematics(filename_afterFid_7546_dvcsmc_bkg, "dfSelected_afterFid_afterCorr", beam_energy);
  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_nobkg = InitKinematics(filename_afterFid_7546_dvcsmc_nobkg, "dfSelected_afterFid_afterCorr", beam_energy);

  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_rad = InitGenKinematics(filename_afterFid_7546_dvcsmc_rad, "MC", beam_energy);
  ROOT::RDF::RNode df_afterFid_7546_dvcsmc_norad = InitGenKinematics(filename_afterFid_7546_dvcsmc_norad, "MC", beam_energy);
  //auto df_afterFid_7546_dvcsmc_rad_sel=ApplyFinalGenDVCSSelections(df_afterFid_7546_dvcsmc_rad, true);
  //auto df_afterFid_7546_dvcsmc_norad_sel=ApplyFinalGenDVCSSelections(df_afterFid_7546_dvcsmc_norad, true);

  DrawStyle fitStyle(0.06, 0.05, 1.0, 1.3);  // You can tweak this
  // DVCS event selection cuts

  // Apply final DVCS cuts
  auto df_final_dvcs_7546_data = ApplyFinalDVCSSelections(df_afterFid_7546_data, true);
  auto df_final_dvcsPi_rejected_7546_data = RejectPi0TwoPhoton(df_final_dvcs_7546_data);
  //auto df_final_dvcs_7546_data_mc = ApplyFinalDVCSSelections(df_afterFid_7546_data_mc, false);
  //auto df_final_dvcsPi_rejected_7546_data_mc = RejectPi0TwoPhoton(df_final_dvcs_7546_data_mc);  
  auto df_final_dvcs_7546_pi0MC = ApplyFinalDVCSSelections(df_afterFid_7546_pi0MC, true);
  auto df_final_dvcsPi_rejected_7546_pi0MC = RejectPi0TwoPhoton(df_final_dvcs_7546_pi0MC);

  auto df_final_OnlPi0_7546_data = SelectPi0Event(df_final_dvcs_7546_data);
  //auto df_final_OnlPi0_7546_data_mc = SelectPi0Event(df_final_dvcs_7546_data_mc);
  auto df_final_OnlPi0_7546_pi0MC = SelectPi0Event(df_final_dvcs_7546_pi0MC);

  auto df_final_dvcs_7546_dvcsmc_rec = ApplyFinalDVCSSelections(df_afterFid_7546_dvcsmc_rec, true);
  auto df_final_dvcsPi_rejected_7546_dvcsmc_rec = RejectPi0TwoPhoton(df_final_dvcs_7546_dvcsmc_rec);

  auto df_final_dvcs_7546_dvcsmc_bkg = ApplyFinalDVCSSelections(df_afterFid_7546_dvcsmc_bkg, true);
  auto df_final_dvcsPi_rejected_7546_dvcsmc_bkg = RejectPi0TwoPhoton(df_final_dvcs_7546_dvcsmc_bkg);

  auto df_final_dvcs_7546_dvcsmc_nobkg = ApplyFinalDVCSSelections(df_afterFid_7546_dvcsmc_nobkg, true);
  auto df_final_dvcsPi_rejected_7546_dvcsmc_nobkg = RejectPi0TwoPhoton(df_final_dvcs_7546_dvcsmc_nobkg);


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
  //xBins.SetTBins({0.0, 10.0});
  // xBins.SetXBBins({0.0, 0.08,.1,.14,.18,.23,.3,.39,.50});
  //xBins.SetQ2Bins({1.0,1.5});
  //xBins.SetTBins({0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0});
  //xBins.SetXBBins({0.0,10.0});
  xBins.SetQ2Bins({1.0, 1.5, 2.0});
  xBins.SetTBins({0.2, 0.3, 0.4, 0.5});
  xBins.SetXBBins({0.15, 0.2, 0.25, 0.3});
  //xBins.SetQ2Bins({0.0, 10.0});
  //xBins.SetTBins({0.2, 0.3});
  //xBins.SetXBBins({0.0, 1.0});
  comparer.SetXBinsRanges(xBins);

  df_final_dvcsPi_rejected_7546_data.Count();
  df_final_OnlPi0_7546_data.Count();
  df_final_dvcsPi_rejected_7546_pi0MC.Count();
  df_final_OnlPi0_7546_pi0MC.Count();
  df_afterFid_7546_dvcsmc_gen.Count();
  df_final_dvcsPi_rejected_7546_dvcsmc_rec.Count();
  df_final_dvcsPi_rejected_7546_dvcsmc_bkg.Count();
  df_final_dvcsPi_rejected_7546_dvcsmc_nobkg.Count();
  //std::cout << df_afterFid_7546_dvcsmc_rad.Count().GetValue() << std::endl;
  //std::cout << df_afterFid_7546_dvcsmc_norad.Count().GetValue() << std::endl;


  comparer.AddModelwithPi0Corr(df_final_dvcsPi_rejected_7546_data,
                              //df_afterFid_7546_dvcsmc_gen,
                              df_final_OnlPi0_7546_data,
                              df_final_dvcsPi_rejected_7546_pi0MC,
                              df_final_OnlPi0_7546_pi0MC,
                              df_afterFid_7546_dvcsmc_gen,
                              df_final_dvcsPi_rejected_7546_dvcsmc_rec,
                              df_final_dvcsPi_rejected_7546_dvcsmc_bkg,
                              df_final_dvcsPi_rejected_7546_dvcsmc_nobkg,
                              df_afterFid_7546_dvcsmc_rad,
                              df_afterFid_7546_dvcsmc_norad,
                              "RGK 7.5GeV", beam_energy, true, true, true, true);

  //comparer.AddModelwithPi0Corr(df_final_dvcsPi_rejected_7546_data,
  //                              //df_afterFid_7546_dvcsmc_norad,
  //                              df_final_OnlPi0_7546_data,
  //                              df_final_dvcsPi_rejected_7546_pi0MC,
  //                              df_final_OnlPi0_7546_pi0MC,
  //                              df_afterFid_7546_dvcsmc_gen,
  //                              df_final_dvcsPi_rejected_7546_dvcsmc_rec,
  //                              df_final_dvcsPi_rejected_7546_dvcsmc_bkg,
  //                              df_final_dvcsPi_rejected_7546_dvcsmc_nobkg,
  //                              df_afterFid_7546_dvcsmc_rad,
  //                              df_afterFid_7546_dvcsmc_norad,
  //                              "born", beam_energy, false, false, false, false);
  

  //comparer.AddModelwithPi0Corr(df_final_dvcsPi_rejected_7546_dvcsmc_rec,df_final_OnlPi0_7546_data,df_final_dvcsPi_rejected_7546_pi0MC,df_final_OnlPi0_7546_pi0MC, df_afterFid_7546_dvcsmc_gen, df_final_dvcs_7546_dvcsmc_rec, "RGK 7.5GeV rec", beam_energy, false, false);
  //comparer.AddModelwithPi0Corr(df_afterFid_7546_dvcsmc_gen,df_final_OnlPi0_7546_data,df_final_dvcsPi_rejected_7546_pi0MC,df_final_OnlPi0_7546_pi0MC, df_afterFid_7546_dvcsmc_gen, df_final_dvcs_7546_dvcsmc_rec, "RGK 7.5GeV gen", beam_energy, false, false);

  //comparer.AddModelwithPi0Corr(df_afterFid_7546_dvcsmc_gen,df_final_OnlPi0_7546_data,df_final_dvcsPi_rejected_7546_pi0MC,df_final_OnlPi0_7546_pi0MC, "RGK 7.5GeV Gen", beam_energy, false);
  //comparer.AddModelwithPi0Corr(df_final_dvcsPi_rejected_7546_data,df_final_OnlPi0_7546_data,df_final_dvcsPi_rejected_7546_MC,df_final_OnlPi0_7546_MC, "RGK 7.5GeV C", beam_energy, true);
  //comparer.AddModelwithPi0Corr(df_afterFid_7546_MC,df_final_OnlPi0_7546_data,df_final_dvcsPi_rejected_7546_MC,df_final_OnlPi0_7546_MC, "RGK 7.5GeV", beam_energy, false);
  //comparer.AddModelwithPi0Corr(df_final_dvcsPi_rejected_7546_data_mc,df_final_OnlPi0_7546_data_mc,df_final_dvcsPi_rejected_7546_MC,df_final_OnlPi0_7546_MC, "RGK 7.5GeV mc ", beam_energy, false);
  
  //double luminosity = 18.8425*pow(10,6);  // Set your desired luminosity here nb^-1
  //double luminosity = (5.516893390230349+2.2356409321666653+2.296226768624658)*1.33*pow(10,6);  // Set your desired luminosity here nb^-1
  double luminosity = (1.515618+3.039676)*1.33*pow(10,6);  // Set your desired luminosity here nb^-1
  double polarisation = 0.85;  // Set your desired polarisation here

  //comparer.PlotKinematicComparison();
  //comparer.PlotxBQ2tBin();
  //comparer.PlotDVCSKinematicsComparison();
  comparer.PlotDIS_BSA_Cross_Section_AndCorr_Comparison(luminosity, polarisation, true, true, true, true, true, true, false);   
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
             .Define("recel_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& vz, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && pass[i]) return vz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})
             .Define("pho_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& trackpass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i] /*&& maxEpass[i]*/) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_px", "REC_Particle_pass", "REC_Photon_MaxE"})
             .Define("pho_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<bool>& trackpass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i] /*&& maxEpass[i]*/) return py[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_py", "REC_Particle_pass", "REC_Photon_MaxE"})
             .Define("pho_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<bool>& trackpass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i] /*&& maxEpass[i]*/) return pz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_pz", "REC_Particle_pass", "REC_Photon_MaxE"})
             .Define("recpho_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& vz, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && pass[i]) return vz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})
             .Define("recpho_beta",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& beta, const ROOT::VecOps::RVec<bool>& trackpass, const ROOT::VecOps::RVec<bool>& maxEpass) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && trackpass[i] /*&& maxEpass[i]*/) return beta[i];
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
             .Define("recpro_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& vz, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && pass[i]) return vz[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})
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
                         if (pid[i] == 22 && pass[i] /*&& maxEpass[i]*/) {
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
        bool hasPi0Daughter = false;
        bool result = false;
        for (size_t i = 0; i < pid.size(); ++i) {
          if (!pass[i]) continue;
          if (pid[i] == 11){
            e++;
          }
          else if (pid[i] == 22 /*&& daughterPass[i]*/){
            g++;  // photon must NOT be from pi0
            hasPi0Daughter = hasPi0Daughter || daughterPass[i];  // check if this photon is a daughter of pi0
            //std::cout << "pid[i] = " << pid[i] << " pass[i] = " << pass[i] << " daughterPass[i] = " << daughterPass[i] << std::endl;
          }
          else if (pid[i] == 2212){
            p++;
          }
        }
        //if(g==2) std::cout << "e = " << e << " g = " << g << " p = " << p << std::endl;
        result = (e == 1 && g >= 1 && p == 1 && hasPi0Daughter);  // at least one photon, no pi0 daughter
        hasPi0Daughter = false;  // reset for next event
        return result;  // at least one photon, no pi0 daughter
      },
      {"REC_Particle_pid", "REC_Particle_pass", "REC_DaughterParticle_pass"}, "Cut: one good e, γ (not π⁰-like), p");
}
// exclusivity cuts
ROOT::RDF::RNode ApplyFinalGenDVCSSelections(ROOT::RDF::RNode df, bool inbending) {
  return df
      .Filter("Q2 > 1.0", "Cut: Q2 > 1 GeV^2")
      .Filter("y < 0.95", "Cut: y < 0.95")
      .Filter("W > 2.0", "Cut: W > 1.8 GeV");
}

ROOT::RDF::RNode ApplyFinalDVCSSelections(ROOT::RDF::RNode df, bool inbending) {
  return df
      //.Filter("(RUN_config_run >= 5689 && RUN_config_run <=5761)||(RUN_config_run >= 5764 && RUN_config_run <=5799)||(RUN_config_run >= 5801 && RUN_config_run <=5840)||RUN_config_run ==11", "Cut: Golden Run")
      .Filter("(RUN_config_run >= 5694 && RUN_config_run <=5705)||(RUN_config_run >= 5707 && RUN_config_run <=5739)||RUN_config_run ==11", "Cut: Golden Run")
      // 4. Q2 > 1
      .Filter("Q2 > 1.0", "Cut: Q2 > 1 GeV^2")
      .Filter("t < 1.0", "Cut: t < 1 GeV^2")
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
  .Filter("Theta_gamma_gamma < 2.0", "Cut: photon-missing angle")
  //.Filter("DeltaPhi < 25.0", "Cut: Coplanarity");
  .Filter("(pho_det_region==0&&pro_det_region==2)||(pho_det_region==1&&pro_det_region==1)||(pho_det_region==1&&pro_det_region==2)", "Cut: three config")
  .Filter("(pho_det_region==0&&pro_det_region==2&&Mx2_ep<0.25&&Mx2_ep>-0.23)||(pho_det_region==1&&pro_det_region==1&&Mx2_ep<0.35&&Mx2_ep>-0.31)||(pho_det_region==1&&pro_det_region==2&&Mx2_ep<0.25&&Mx2_ep>-0.23)", "Cut: Mx2_ep in 3sigma")
  .Filter("(pho_det_region==0&&pro_det_region==2&&Emiss<0.37&&Emiss>-0.29)||(pho_det_region==1&&pro_det_region==1&&Emiss<0.62&&Emiss>-0.58)||(pho_det_region==1&&pro_det_region==2&&Emiss<0.72&&Emiss>-0.6)", "Cut: Emiss in 3sigma")
  .Filter("(pho_det_region==0&&pro_det_region==2&&PTmiss<0.08&&PTmiss>-0.04)||(pho_det_region==1&&pro_det_region==1&&PTmiss<0.18&&PTmiss>-0.06)||(pho_det_region==1&&pro_det_region==2&&PTmiss<0.13&&PTmiss>-0.05)", "Cut: PTmiss in 3sigma")
  .Filter("(pho_det_region==0&&pro_det_region==2&&Theta_gamma_gamma<1.16&&Theta_gamma_gamma>-0.58)||(pho_det_region==1&&pro_det_region==1&&Theta_gamma_gamma<1.73&&Theta_gamma_gamma>-0.79)||(pho_det_region==1&&pro_det_region==2&&Theta_gamma_gamma<1.26&&Theta_gamma_gamma>-0.6)", "Cut: Theta_gamma_gamma in 3sigma")
  .Filter("(pho_det_region==0&&pro_det_region==2&&DeltaPhi<3.93&&DeltaPhi>-2.55)||(pho_det_region==1&&pro_det_region==1&&DeltaPhi<11.87&&DeltaPhi>-6.85)||(pho_det_region==1&&pro_det_region==2&&DeltaPhi<7.42&&DeltaPhi>-4.76)", "Cut: DeltaPhi in 3sigma")
  .Filter("(pho_det_region==0&&pro_det_region==2&&Mx2_epg<0.03&&Mx2_epg>-0.03)||(pho_det_region==1&&pro_det_region==1&&Mx2_epg<0.03&&Mx2_epg>-0.03)||(pho_det_region==1&&pro_det_region==2&&Mx2_epg<0.03&&Mx2_epg>-0.03)", "Cut: Mx2_epg in 3sigma")
  .Filter("(pho_det_region==0&&pro_det_region==2&&Mx2_eg<1.55&&Mx2_eg>0.35)||(pho_det_region==1&&pro_det_region==1&&Mx2_eg<1.68&&Mx2_eg>0.12)||(pho_det_region==1&&pro_det_region==2&&Mx2_eg<1.99&&Mx2_eg>-0.05)", "Cut: Mx2_eg in 3sigma")
  .Filter("(pho_det_region==0&&pro_det_region==2&&Theta_e_gamma<27.92&&Theta_e_gamma>5.42)||(pho_det_region==1&&pro_det_region==1&&Theta_e_gamma<44.71&&Theta_e_gamma>29.17)||(pho_det_region==1&&pro_det_region==2&&Theta_e_gamma<36.36&&Theta_e_gamma>10.36)", "Cut: Theta_e_gamma in 3sigma");

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


ROOT::RDF::RNode InitGenKinematics(const std::string& filename_, const std::string& treename_, float beam_energy) {
  ROOT::RDataFrame rdf(treename_, filename_);
  auto df_ = std::make_unique<ROOT::RDF::RNode>(rdf);
  *df_ = df_->Define("ele_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11) return px[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_px"})
             .Define("ele_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11) return py[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_py"})
             .Define("ele_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11) return pz[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_pz"})
             .Define("pho_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<float>& pz) {
                      for (size_t i = 0; i < pid.size(); ++i){
                         if (pid[i]==22) return px[i];
                      }
                         //if (pid[i] == 22 /*&& maxEpass[i]*/) return px[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_px","MC_Particle_py","MC_Particle_pz"})
             .Define("pho_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<float>& pz) {
                      for (size_t i = 0; i < pid.size(); ++i){
                         if (pid[i]==22) return py[i];
                       }
                         //if (pid[i] == 22 /*&& maxEpass[i]*/) return px[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_px","MC_Particle_py","MC_Particle_pz"})
             .Define("pho_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<float>& pz) {
                      float temp_p = 0.0; 
                      float result = -999.0f;
                      for (size_t i = 0; i < pid.size(); ++i){
                         if (pid[i]==22) return pz[i];
                       }
                         //if (pid[i] == 22 /*&& maxEpass[i]*/) return px[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_px","MC_Particle_py","MC_Particle_pz"})
             .Define("pro_px",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 ) return px[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_px"})
             .Define("pro_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 ) return py[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_py"})
             .Define("pro_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 ) return pz[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_pz"})
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
             .Define("recel_vz", [](){ return 0.0; })
             .Define("recpho_vz", [](){ return 0.0; })
             .Define("recpro_vz", [](){ return 0.0; })
             .Define("pho_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid) {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 22  /*&& maxEpass[i]*/) {
                           return 1;
                         }
                       }
                       return -1;
                     },
                     {"MC_Particle_pid"})

             .Define("pro_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid) {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 2212 ) {
                           return 1;
                         }
                       }
                       return -1;
                     },
                     {"MC_Particle_pid"})
             .Define("ele_det_region",
                     [](const ROOT::VecOps::RVec<int>& pid) {
                       for (size_t i = 0; i < pid.size(); ++i) {
                         if (pid[i] == 11 ) {
                           return 1;
                         }
                       }
                       return -1;
                     },
                     {"MC_Particle_pid"});

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