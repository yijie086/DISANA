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
  std::string input_path_from_analysisRun = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/NewAnalysisFrameWork/testing_outupt/afterFiducialCuts/DC_fiducialcuts/";
  //std::string input_path_from_analysisRun = "./../build";
  std::string filename_after = Form("%s/dfSelected_after_fiducialCuts.root", input_path_from_analysisRun.c_str());
  float beam_energy = 10.6;
  
  //// this is where you can plot the comparisions
  ROOT::RDF::RNode df = InitKinematics(filename_after, "dfSelected_after");
  // input files for the data


  // compute the dvcs variables
  df = define_DISCAT(df, "Q2", &DISANAMath::GetQ2, beam_energy);
  df = define_DISCAT(df, "xB", &DISANAMath::GetxB, beam_energy);
  df = define_DISCAT(df, "t", &DISANAMath::GetT, beam_energy);
  df = define_DISCAT(df, "phi", &DISANAMath::GetPhi, beam_energy);
  df = define_DISCAT(df, "W", &DISANAMath::GetW, beam_energy);
  df = define_DISCAT(df, "nu", &DISANAMath::GetNu, beam_energy);
  df = define_DISCAT(df, "y", &DISANAMath::Gety, beam_energy);
  

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

  comparer.AddModel(df, "RGA-A inb", 10.6);
  comparer.PlotKinematicComparison();
  comparer.PlotDVCSKinematicsComparison();
  comparer.PlotDISCrossSectionComparison(1);  // argument is Luminosity
  comparer.PlotDIS_BSA_Comparison(1);         // argument is Luminosity
}

ROOT::RDF::RNode InitKinematics(const std::string& filename_, const std::string& treename_) {
  ROOT::RDataFrame rdf(treename_, filename_);

  auto df_ = std::make_unique<ROOT::RDF::RNode>(rdf.Define("ele_px",
                                                           [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px) {
                                                             for (size_t i = 0; i < pid.size(); ++i)
                                                               if (pid[i] == 11) return px[i];
                                                             return -999.0f;
                                                           },
                                                           {"REC_Particle_pid", "REC_Particle_px"})
                                                    .Define("ele_py",
                                                            [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py) {
                                                              for (size_t i = 0; i < pid.size(); ++i)
                                                                if (pid[i] == 11) return py[i];
                                                              return -999.0f;
                                                            },
                                                            {"REC_Particle_pid", "REC_Particle_py"})
                                                    .Define("ele_pz",
                                                            [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz) {
                                                              for (size_t i = 0; i < pid.size(); ++i)
                                                                if (pid[i] == 11) return pz[i];
                                                              return -999.0f;
                                                            },
                                                            {"REC_Particle_pid", "REC_Particle_pz"})
                                                    .Define("pho_px",
                                                            [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px) {
                                                              for (size_t i = 0; i < pid.size(); ++i)
                                                                if (pid[i] == 22) return px[i];
                                                              return -999.0f;
                                                            },
                                                            {"REC_Particle_pid", "REC_Particle_px"})
                                                    .Define("pho_py",
                                                            [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py) {
                                                              for (size_t i = 0; i < pid.size(); ++i)
                                                                if (pid[i] == 22) return py[i];
                                                              return -999.0f;
                                                            },
                                                            {"REC_Particle_pid", "REC_Particle_py"})
                                                    .Define("pho_pz",
                                                            [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz) {
                                                              for (size_t i = 0; i < pid.size(); ++i)
                                                                if (pid[i] == 22) return pz[i];
                                                              return -999.0f;
                                                            },
                                                            {"REC_Particle_pid", "REC_Particle_pz"})
                                                    .Define("pro_px",
                                                            [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px) {
                                                              for (size_t i = 0; i < pid.size(); ++i)
                                                                if (pid[i] == 2212) return px[i];
                                                              return -999.0f;
                                                            },
                                                            {"REC_Particle_pid", "REC_Particle_px"})
                                                    .Define("pro_py",
                                                            [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py) {
                                                              for (size_t i = 0; i < pid.size(); ++i)
                                                                if (pid[i] == 2212) return py[i];
                                                              return -999.0f;
                                                            },
                                                            {"REC_Particle_pid", "REC_Particle_py"})
                                                    .Define("pro_pz",
                                                            [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz) {
                                                              for (size_t i = 0; i < pid.size(); ++i)
                                                                if (pid[i] == 2212) return pz[i];
                                                              return -999.0f;
                                                            },
                                                            {"REC_Particle_pid", "REC_Particle_pz"})
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