#include"../DreamAN/DrawHist/DrawStyle.h"
#include "../DreamAN/DrawHist/DISANAcomparer.h"

ROOT::RDF::RNode InitKinematics(const std::string& filename_ = "", const std::string& treename_ = "");
 void PlotDVCSKinematicsComparison(ROOT::RDF::RNode& rdf);
static double MomentumFunc(float px, float py, float pz) { return std::sqrt(px * px + py * py + pz * pz); }
static double ThetaFunc(float px, float py, float pz) { return std::acos(pz / std::sqrt(px * px + py * py + pz * pz)); }
static double PhiFunc(float px, float py) {
  double phi = std::atan2(py, px);
  return phi < 0 ? phi + 2 * M_PI : phi;
}

template <typename Method>
ROOT::RDF::RNode define_DISCAT(ROOT::RDF::RNode node, const std::string& name, const Method method, float beam_energy) {
  return node.Define(name,
                     [method, beam_energy](double recel_p, double recel_theta, double recel_phi, double recpho_p, double recpho_theta, double recpho_phi, double recpro_p,
                                           double recpro_theta, double recpro_phi) {
                       return (DISANAMath(beam_energy, recel_p, recel_theta, recel_phi, recpho_p, recpho_theta, recpho_phi, recpro_p, recpro_theta, recpro_phi).*method)();
                     },
                     {"recel_p", "recel_theta", "recel_phi", "recpho_p", "recpho_theta", "recpho_phi", "recpro_p", "recpro_theta", "recpro_phi"});
}

void DISANA_Xplotter() {
  ROOT::EnableImplicitMT();
  //std::string input_path_from_analysisRun = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/testing_outupt/afterFiducialCuts/DC_fiducialcuts/";
  std::string input_path_from_analysisRun = "./../build";
  std::string filename_afterFid = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun.c_str());
  std::string filename_afterFid_afterCorr = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun.c_str());
  float beam_energy = 7.546;
  
  //// this is where you can plot the comparisions
  ROOT::RDF::RNode df_afterFid = InitKinematics(filename_afterFid, "dfSelected_afterFid");
  ROOT::RDF::RNode df_afterFid_afterCorr = InitKinematics(filename_afterFid_afterCorr, "dfSelected_afterFid_afterCorr");
  // input files for the data


  // compute the dvcs variables
  df_afterFid = define_DISCAT(df_afterFid, "Q2", &DISANAMath::GetQ2, beam_energy);
  df_afterFid = define_DISCAT(df_afterFid, "xB", &DISANAMath::GetxB, beam_energy);
  df_afterFid = define_DISCAT(df_afterFid, "t", &DISANAMath::GetT, beam_energy);
  df_afterFid = define_DISCAT(df_afterFid, "phi", &DISANAMath::GetPhi, beam_energy);
  df_afterFid = define_DISCAT(df_afterFid, "W", &DISANAMath::GetW, beam_energy);
  df_afterFid = define_DISCAT(df_afterFid, "nu", &DISANAMath::GetNu, beam_energy);
  df_afterFid = define_DISCAT(df_afterFid, "y", &DISANAMath::Gety, beam_energy);

  df_afterFid_afterCorr = define_DISCAT(df_afterFid_afterCorr, "Q2", &DISANAMath::GetQ2, beam_energy);
  df_afterFid_afterCorr = define_DISCAT(df_afterFid_afterCorr, "xB", &DISANAMath::GetxB, beam_energy);
  df_afterFid_afterCorr = define_DISCAT(df_afterFid_afterCorr, "t", &DISANAMath::GetT, beam_energy);
  df_afterFid_afterCorr = define_DISCAT(df_afterFid_afterCorr, "phi", &DISANAMath::GetPhi, beam_energy);
  df_afterFid_afterCorr = define_DISCAT(df_afterFid_afterCorr, "W", &DISANAMath::GetW, beam_energy);
  df_afterFid_afterCorr = define_DISCAT(df_afterFid_afterCorr, "nu", &DISANAMath::GetNu, beam_energy);
  df_afterFid_afterCorr = define_DISCAT(df_afterFid_afterCorr, "y", &DISANAMath::Gety, beam_energy);
  

  //// styling plots 
  //double double titleSize = 0.05, double labelSize = 0.04,double xTitleOffset = 1.1, double yTitleOffset = 1.6, int font = 42, int maxDigits = 5, int nDivisions = 510, double leftMargin = 0.16, double rightMargin = 0.07, double bottomMargin = 0.13, double topMargin = 0.06
  DrawStyle KinStyle(0.07, 0.06, 0.9, 1.2);     // For Kin plots
  DrawStyle dvcsStyle(0.06, 0.06,1.2,1.4,42,5,510,0.14,0.07,0.13,0.06);     // For DVCS plots
  DrawStyle csStyle(0.05, 0.04, 1.0, 1.3);        // For Cross-Sections
  DrawStyle bsaStyle(0.06, 0.045, 1.0, 1.2);      // For BSA

  DISANAcomparer comparer;
  comparer.SetOutputDir("./");
  comparer.SetKinStyle(KinStyle);
  comparer.SetDVCSStyle(dvcsStyle);
  comparer.SetCrossSectionStyle(csStyle);
  comparer.SetBSAStyle(bsaStyle);
  comparer.PlotIndividual(false);
  /// bins for cross-section plots
  BinManager xBins;
   //xBins.SetQ2Bins({.11,1.3,1.6,2.1,2.8,3.6,8.0});
   //xBins.SetTBins({0.0, 1.2});
   //xBins.SetXBBins({0.0, 0.08,.1,.14,.18,.23,.3,.39,.50});
  xBins.SetQ2Bins({.1,3.0, 8.0});
  xBins.SetTBins({0.0,0.5, 1.2});
  xBins.SetXBBins({0.0,0.5, 1});
  comparer.SetXBinsRanges(xBins);

  
  //comparer.AddModel(df_afterFid_afterCorr, "after Correction", beam_energy);
  comparer.AddModel(df_afterFid, "before Correction", beam_energy);
  comparer.PlotKinematicComparison();
  comparer.PlotDVCSKinematicsComparison();
  comparer.PlotDISCrossSectionComparison(1);  // argument is Luminosity
  comparer.PlotDIS_BSA_Comparison(1);         // argument is Luminosity
  gApplication->Terminate(0);
}

ROOT::RDF::RNode InitKinematics(const std::string& filename_, const std::string& treename_) {
  ROOT::RDataFrame rdf(treename_, filename_);

  auto df_ = std::make_unique<ROOT::RDF::RNode>(rdf.Define("ele_px",
                                                           [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& trackpass) {
                                                             for (size_t i = 0; i < pid.size(); ++i)
                                                               if (pid[i] == 11 && trackpass[i] == 1) return px[i];
                                                             return -999.0f;
                                                           },
                                                           {"REC_Particle_pid", "REC_Particle_px" ,"REC_Track_pass_fid"})
                                                    .Define("ele_py",
                                                            [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<int>& trackpass) {
                                                              for (size_t i = 0; i < pid.size(); ++i)
                                                                if (pid[i] == 11 && trackpass[i] == 1) return py[i];
                                                              return -999.0f;
                                                            },
                                                            {"REC_Particle_pid", "REC_Particle_py" ,"REC_Track_pass_fid"})
                                                    .Define("ele_pz",
                                                            [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<int>& trackpass) {
                                                              for (size_t i = 0; i < pid.size(); ++i)
                                                                if (pid[i] == 11 && trackpass[i] == 1) return pz[i];
                                                              return -999.0f;
                                                            },
                                                            {"REC_Particle_pid", "REC_Particle_pz" ,"REC_Track_pass_fid"})
                                                    .Define("pho_px",
                                                            [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& trackpass) {
                                                              for (size_t i = 0; i < pid.size(); ++i)
                                                                if (pid[i] == 22 && trackpass[i] == 1) return px[i];
                                                              return -999.0f;
                                                            },
                                                            {"REC_Particle_pid", "REC_Particle_px" ,"REC_Track_pass_fid"})
                                                    .Define("pho_py",
                                                            [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<int>& trackpass) {
                                                              for (size_t i = 0; i < pid.size(); ++i)
                                                                if (pid[i] == 22 && trackpass[i] == 1) return py[i];
                                                              return -999.0f;
                                                            },
                                                            {"REC_Particle_pid", "REC_Particle_py" ,"REC_Track_pass_fid"})
                                                    .Define("pho_pz",
                                                            [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<int>& trackpass) {
                                                              for (size_t i = 0; i < pid.size(); ++i)
                                                                if (pid[i] == 22 && trackpass[i] == 1) return pz[i];
                                                              return -999.0f;
                                                            },
                                                            {"REC_Particle_pid", "REC_Particle_pz" ,"REC_Track_pass_fid"})
                                                    .Define("pro_px",
                                                            [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<int>& trackpass) {
                                                              for (size_t i = 0; i < pid.size(); ++i)
                                                                if (pid[i] == 2212 && trackpass[i] == 1) return px[i];
                                                              return -999.0f;
                                                            },
                                                            {"REC_Particle_pid", "REC_Particle_px" ,"REC_Track_pass_fid"})
                                                    .Define("pro_py",
                                                            [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py, const ROOT::VecOps::RVec<int>& trackpass) {
                                                              for (size_t i = 0; i < pid.size(); ++i)
                                                                if (pid[i] == 2212 && trackpass[i] == 1) return py[i];
                                                              return -999.0f;
                                                            },
                                                            {"REC_Particle_pid", "REC_Particle_py" ,"REC_Track_pass_fid"})
                                                    .Define("pro_pz",
                                                            [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz, const ROOT::VecOps::RVec<int>& trackpass) {
                                                              for (size_t i = 0; i < pid.size(); ++i)
                                                                if (pid[i] == 2212 && trackpass[i] == 1) return pz[i];
                                                              return -999.0f;
                                                            },
                                                            {"REC_Particle_pid", "REC_Particle_pz" ,"REC_Track_pass_fid"})
                                                    .Filter([](float ex, float gx, float px) { return ex != -999 && gx != -999 && px != -999; }, {"ele_px", "pho_px", "pro_px"})
                                                    .Define("recel_p", MomentumFunc, {"ele_px", "ele_py", "ele_pz"})
                                                    .Define("recel_theta", ThetaFunc, {"ele_px", "ele_py", "ele_pz"})
                                                    .Define("recel_phi", PhiFunc, {"ele_px", "ele_py"})
                                                    .Define("recpho_p", MomentumFunc, {"pho_px", "pho_py", "pho_pz"})
                                                    .Define("recpho_theta", ThetaFunc, {"pho_px", "pho_py", "pho_pz"})
                                                    .Define("recpho_phi", PhiFunc, {"pho_px", "pho_py"})
                                                    .Define("recpro_p", MomentumFunc, {"pro_px", "pro_py", "pro_pz"})
                                                    .Define("recpro_theta", ThetaFunc, {"pro_px", "pro_py", "pro_pz"})
                                                    .Define("recpro_phi", PhiFunc, {"pro_px", "pro_py"}));

  return *df_;
}