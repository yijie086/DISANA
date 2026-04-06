#include <TApplication.h>
#include <THnSparse.h>
#include <TROOT.h>
#include <TString.h>

#include <ROOT/RDataFrame.hxx>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

#include "../DreamAN/DrawHist/DISANAMMUtils.h"
#include "../DreamAN/DrawHist/DISANAMath.h"
#include "../DreamAN/DrawHist/DISANAMathFitUtils.h"
#include "../DreamAN/DrawHist/DISANA_PhiMassUtils.h"
#include "../DreamAN/DrawHist/DISANAcomparer.h"
#include "../DreamAN/DrawHist/DrawStyle.h"

ROOT::RDF::RNode ApplyFinalDVEPSelections(ROOT::RDF::RNode df);
ROOT::RDF::RNode ApplyFinalGenDVEPSelections(ROOT::RDF::RNode df);
static void ExportBinningCSV(const std::vector<double>& q2Edges, const std::vector<double>& tpEdges, double beamMom, double W2lo, double W2hi, const std::string& outPath);
std::vector<std::pair<std::string, std::string>> detCuts = {{"pro_det_region == 2", "CD"}, {"pro_det_region == 1", "FD"}};
void PrintEdges(const std::string& label, const std::vector<double>& edges);
// Plot styling (unchanged)
DrawStyle KinStyle(0.07, 0.06, 0.9, 1.2);                                       // For Kin plots
DrawStyle dvcsStyle(0.06, 0.06, 1.0, 1.2, 42, 5, 510, 0.14, 0.07, 0.14, 0.06);  // For DVCS plots
DrawStyle csStyle(0.05, 0.05, .95, 1.1, 42, 5, 510, 0.12, 0.03, 0.12, 0.02);    // For Cross-Sections
DrawStyle bsaStyle(0.06, 0.045, .8, .8, 42, 5, 510, 0.15, 0.07, 0.16, 0.06);    // For BSA
double computeLuminosity(double Q_coulombs);

// Main plotter with toggles
void DISANA_PhiAnalysisPlotter()  // subset toggle inside missing-mass
{
  bool runExclusive = true;
  bool runMissingMass = false;
  bool plotIntInvMass = true;  // toggle for plotting invariant mass (phi peak) distributions
  bool GenOnly = false;
  bool RecGEMCOnly = false;
  bool RecDataOnly = true;
  bool TryHStyle = false;  // toggle to try HStyle for plotting (instead of custom DrawStyle)
  ROOT::EnableImplicitMT(40);

  std::string outputDir = "phi_analysis_plots";

  // Spring 2018 INB (10.2 GeV)
  // -----------------------------
  // Input locations
  // Exclusive reconstruction K+K-
  std::string input_result_folder_excl = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/DVKpKm/";
  //std::string input_result_folder_excl_test = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/DSTs/script/";
  std::string input_result_folder_excl_test = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/DSTs/standalone/Column_optimisation_test/";
  std::string input_result_folder_excl_missing = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/nSIDIS/hipo2root/Final_reprocessed/";
  /*std::string input_path_from_analysisRun_SP18inb_data =
   * "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/spring2018/inb/DVKpKm_wagon/after_fids/SF_momentum_corr";
   */
  // std::string input_path_from_analysisRun_SP18outb_data
  // ="/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/spring2018/outb/DVKpKm_wagon/after_fids/SF_momentum_corr";
  std::string input_path_from_analysisRun_SP18outb_data = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_skiming/build/resutls/skims/sp2019_inb/n_sidis_missingKm/";
  /*std::string input_path_from_analysisRun_Fall18inb_data =
  "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/fall2018/inb/DVKpKm_wagon/after_fids/SF_momentum_corr"; std::string
  input_path_from_analysisRun_Fall18outb_data = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/fall2018/outb/DVKpKm_wagon/after_fids/SF_momentum_corr";
  std::string input_path_from_analysisRun_SP19inb_data =
  "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/spring2019/inb/DVKpKm_wagon/after_fids/SF_momentum_corr";*/
  /// MC path for exclusivity fits

  // MC RECONSTRUCTED GEMC
  std::string input_path_from_analysisRun_SP18inb_MC = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/sp2018_inb/50nA/";
  std::string input_path_from_analysisRun_SP19inb_MC = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/sp2019/50nA/";
  std::string input_path_from_analysisRun_SP19inb_MC_clasdis = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/clasdis/sp2019_50nA_inb/";
  std::string input_path_from_analysisRun_SP19inb_MC_Harut_p39 = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/Harut_MC/p39_10p2/";

  std::string input_path_from_analysisRun_Fall18inb_MC = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/fall2018_inb/55nA/";
  std::string input_path_from_analysisRun_Fall18outb_MC = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/fall2018_outb/50nA/";
  std::string input_path_from_analysisRun_SP18outb_MC = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/sp2018_outb/45nA/";

  // MC Reconstructed for efficiencies
  std::string input_path_from_analysisRun_SP18inb_MC_eff_nobkg =
      "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/Efficiency/sp2018_inb/no_bkg/";
  std::string input_path_from_analysisRun_SP18inb_MC_eff_50nA = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/Efficiency/sp2018_inb/50nA/";
  std::string input_path_from_analysisRun_SP18outb_MC_eff_45nA =
      "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/Efficiency/sp2018_outb/45nA/";
  std::string input_path_from_analysisRun_SP18outb_MC_eff_nobkg =
      "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/Efficiency/sp2018_outb/no_bkg/";

  std::string input_path_from_analysisRun_SP19inb_MC_eff_50nA = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/Efficiency/sp2019_inb/50nA/";
  std::string input_path_from_analysisRun_SP19inb_MC_eff_nobkg =
      "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/Efficiency/sp2019_inb/no_bkg/";

  // MC GEN
  std::string input_path_from_analysisRun_SP18inb_MCgen = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/sp2018_inb/50nA/";
  std::string input_path_from_analysisRun_SP19inb_MCgen = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/sp2019/50nA/";
  std::string input_path_from_analysisRun_SP19inb_MCgen_Harut_p39 = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/Harut_MC/p39_10p2/";
  std::string input_path_from_analysisRun_Fall18inb_MCgen = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/fall2018_inb/55nA/";
  std::string input_path_from_analysisRun_Fall18outb_MCgen = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/fall2018_outb/50nA/";
  std::string input_path_from_analysisRun_SP18outb_MCgen = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/sims/lager/sp2018_outb/45nA/";

  // MC GEN test paths Harut
  std::string input_path_from_analysisRun_SP18inb_data_qadb = Form("%s/rgasp18_inb", input_result_folder_excl.c_str());
  std::string input_path_from_analysisRun_SP18outb_data_qadb = Form("%s/rgasp18_outb", input_result_folder_excl.c_str());
  std::string input_path_from_analysisRun_Fall18inb_data_qadb = Form("%s/rgafall18_inb", input_result_folder_excl.c_str());
  std::string input_path_from_analysisRun_Fall18outb_data_qadb = Form("%s/rgafall18_outb", input_result_folder_excl.c_str());
  std::string input_path_from_analysisRun_SP19inb_data_qadb = Form("%s", input_result_folder_excl_test.c_str());

  // std::string input_path_from_analysisRun_SP19inb_data_qadb = Form("%s", "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/spring2019/inb/DVKpKm_wagon/");
  //  Reprocessed missing–mass
  std::string input_path_from_analysisRun_SP19inb_data_missingKm = "./../data_processed/spring2019/inb/nsidis_wagon/missing_Km_output/";
  std::string input_path_from_analysisRun_Fall18inb_data_missingKm = "./../data_processed/fall2018/inb/nsidis_wagon/missing_Km_output/";
  std::string input_path_from_analysisRun_SP18inb_data_missingKm = "./../data_processed/spring2018/inb/nsidis_wagon/missing_Km_output/";
  std::string input_path_from_analysisRun_SP18outb_data_missingKp = "./../data_processed/spring2018/outb/nsidis_wagon/missing_Kp_output/";
  /// std::string input_path_from_analysisRun_Fall18outb_data_missingKp =
  /// "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/spring2019/inb/nsidis_wagon/missing_Km_output/";
  std::string input_path_from_analysisRun_Sp2019_data_missingKp = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/Phi_data_processed/n_sidis/sp2019_inb_missingKm/";

  /*std::string input_path_from_analysisRun_SP18inb_data_missingKm = "/w/hallb-scshelf2102/clas12/singh/data_repo/phi_analysis/skim_from_nsidis/sp2018/inb/missing_Km_output/";
  std::string input_path_from_analysisRun_SP18outb_data_missingKp = "/w/hallb-scshelf2102/clas12/singh/data_repo/phi_analysis/skim_from_nsidis/sp2018/outb/missing_Kp_output/";
  std::string input_path_from_analysisRun_Fall18inb_data_missingKm = "/w/hallb-scshelf2102/clas12/singh/data_repo/phi_analysis/skim_from_nsidis/fall2018/inb/missing_km_output/";
  std::string input_path_from_analysisRun_Fall18outb_data_missingKp = "/w/hallb-scshelf2102/clas12/singh/data_repo/phi_analysis/skim_from_nsidis/fall2018/outb/missing_Kp_output/";
  std::string input_path_from_analysisRun_SP19inb_data_qadb_missingKm =
  "/w/hallb-scshelf2102/clas12/singh/data_repo/phi_analysis/skim_from_nsidis/sp2019/sp2019_inb/missing_km_output/";
*/
  std::string input_path_from_analysisRun_SP18inb_data_qadb_missingKm = Form("%s/rgasp18_inb", input_result_folder_excl_missing.c_str());
  std::string input_path_from_analysisRun_SP18outb_data_qadb_missingKp = Form("%s/rgasp18_outb", input_result_folder_excl_missing.c_str());
  std::string input_path_from_analysisRun_Fall18inb_data_qadb_missingKm = Form("%s/rgafall18_inb", input_result_folder_excl_missing.c_str());
  std::string input_path_from_analysisRun_Fall18outb_data_qadb_missingKp = Form("%s/rgafall18_outb", input_result_folder_excl_missing.c_str());
  std::string input_path_from_analysisRun_SP19inb_data_qadb_missingKm = Form("%s/rgasp19_inb", input_result_folder_excl_missing.c_str());

  // File names
  std::string filename_afterFid_SP18inb_data = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_SP18inb_data_qadb.c_str());
  std::string filename_afterFid_SP18outb_data = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_SP18outb_data_qadb.c_str());
  std::string filename_afterFid_Fall18inb_data = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_Fall18inb_data_qadb.c_str());
  std::string filename_afterFid_Fall18outb_data = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_Fall18outb_data_qadb.c_str());
  std::string filename_afterFid_SP19inb_data = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_SP19inb_data_qadb.c_str());

  // no qadb files for exclusive analysis
  std::string filename_afterFid_SP18outb_data_no_qadb = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_SP18outb_data.c_str());

  std::string filename_afterFid_SP19inb_missingKm_data_org = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_Sp2019_data_missingKp.c_str());
  std::string filename_afterFid_SP19inb_missingKm_data = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_SP19inb_data_qadb_missingKm.c_str());
  std::string filename_afterFid_Fall18inb_missingKm_data = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_Fall18inb_data_qadb_missingKm.c_str());
  std::string filename_afterFid_SP18inb_missingKm_data = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_SP18inb_data_qadb_missingKm.c_str());
  std::string filename_afterFid_SP18outb_missingKp_data = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_SP18outb_data_qadb_missingKp.c_str());
  std::string filename_afterFid_Fall18outb_missingKp_data = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_Fall18outb_data_qadb_missingKp.c_str());

  // MC rec files for exclusivity fits
  std::string filename_afterFid_SP18inb_MC = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_SP18inb_MC.c_str());
  std::string filename_afterFid_SP19inb_MC = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_SP19inb_MC.c_str());
  std::string filename_afterFid_SP19inb_MC_clasdis = Form("%s/dfSelected.root", input_path_from_analysisRun_SP19inb_MC_clasdis.c_str());
  std::string filename_afterFid_SP19inb_MC_Harut_p39 = Form("%s/dfSelected.root", input_path_from_analysisRun_SP19inb_MC_Harut_p39.c_str());
  std::string filename_afterFid_Fall18inb_MC = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_Fall18inb_MC.c_str());
  std::string filename_afterFid_Fall18outb_MC = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_Fall18outb_MC.c_str());
  std::string filename_afterFid_SP18outb_MC = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_SP18outb_MC.c_str());

  // MC rec for efficiencies
  std::string filename_afterFid_SP18inb_MC_eff_nobkg = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_SP18inb_MC_eff_nobkg.c_str());
  std::string filename_afterFid_SP18outb_MC_eff_nobkg = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_SP18outb_MC_eff_nobkg.c_str());
  std::string filename_afterFid_SP19inb_MC_eff_nobkg = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_SP19inb_MC_eff_nobkg.c_str());
  std::string filename_afterFid_SP18inb_MC_eff_50nA = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_SP18inb_MC_eff_50nA.c_str());
  std::string filename_afterFid_SP19inb_MC_eff_50nA = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_SP19inb_MC_eff_50nA.c_str());
  std::string filename_afterFid_SP18outb_MC_eff_45nA = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_SP18outb_MC_eff_45nA.c_str());

  // MC GEN files for exclusivity fits
  std::string filename_afterFid_SP18inb_MCgen = Form("%s/dfSelectedMC.root", input_path_from_analysisRun_SP18inb_MCgen.c_str());
  std::string filename_afterFid_SP19inb_MCgen = Form("%s/dfSelectedMC.root", input_path_from_analysisRun_SP19inb_MCgen.c_str());
  std::string filename_afterFid_SP19inb_MCgen_Harut_p39 = Form("%s/dfSelectedMC.root", input_path_from_analysisRun_SP19inb_MCgen_Harut_p39.c_str());
  std::string filename_afterFid_Fall18inb_MCgen = Form("%s/dfSelectedMC.root", input_path_from_analysisRun_Fall18inb_MCgen.c_str());
  std::string filename_afterFid_Fall18outb_MCgen = Form("%s/dfSelectedMC.root", input_path_from_analysisRun_Fall18outb_MCgen.c_str());
  std ::string filename_afterFid_SP18outb_MCgen = Form("%s/dfSelectedMC.root", input_path_from_analysisRun_SP18outb_MCgen.c_str());
  // Beam energies
  float beam_energy_sp2018 = 10.5940f;
  float beam_energy_fall2018 = 10.6000f;
  float beam_energy_sp2019 = 10.1998f;
  double polarization_sp18 = 0.85;

  // -----------------------------
  // Comparer setup (we’ll add models conditionally below)
  DISANAcomparer comparer;
  comparer.SetOutputDir("./");
  comparer.SetKinStyle(KinStyle);
  comparer.SetDVCSStyle(dvcsStyle);
  comparer.SetCrossSectionStyle(csStyle);
  comparer.SetBSAStyle(bsaStyle);
  comparer.PlotIndividual(false);

  // Binning (unchanged)
  BinManager xBins;
  xBins.SetXBBins({0.0, 0.7});
  xBins.SetWBins({1.8, 4.0});
  xBins.SetQ2Bins({0.0000, 2.2667, 2.9333, 3.7333, 6.0000});
  xBins.SetTprimeBins({0.000, 0.2250, 0.3750, 0.5250, 0.6750, 0.8250, 1.0500, 1.3500, 1.8750, 4.5000});
  xBins.SetCosThetaKKBins({-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0});
  comparer.SetXBinsRanges(xBins);
  comparer.UseFittedPhiYields(true);

  // Some global constants you had
  double luminosity_rga_sp18_outb = computeLuminosity(16.465615589217398);  // in nb^-1
  // double luminosity_rga_sp18_outb = computeLuminosity(1.7447073374328613E7); // in nb^-1
  double luminosity_rga_sp18_inb = computeLuminosity(44.183881173812);        // in nb^-1
  double luminosity_rga_fall18_outb = computeLuminosity(44.478816492680613);  // in nb^-1
  double luminosity_rga_fall18_inb = computeLuminosity(35.023407601823784);   // in nb^-1
  double luminosity_rga_sp19_inb = computeLuminosity(44.183881173812);        // in nb^-1
  double polarisation = 0.85;
  double branching = 0.49;

  cout << "Using luminosity (nb^-1): " << luminosity_rga_sp18_outb << " for RGA SP18 OUTB" << endl;
  cout << "Using luminosity (nb^-1): " << luminosity_rga_sp18_inb << " for RGA SP18 INB" << endl;
  cout << "Using luminosity (nb^-1): " << luminosity_rga_fall18_outb << " for RGA FALL18 OUTB" << endl;
  cout << "Using luminosity (nb^-1): " << luminosity_rga_fall18_inb << " for RGA FALL18 INB" << endl;
  cout << "Using luminosity (nb^-1): " << luminosity_rga_sp19_inb << " for RGA SP19 INB" << endl;

  //
  double charge = 16.465615589217398;  //(mC)//4.815525219658029+(8.88177914805192-0.2128897513862203)*0.5; // mC (5681-5757, 5757-5870(trigger prescale 2), 5758removed)
  std::cout << "Total effective charge (mC): " << charge << std::endl;

  double luminosity = charge * 1.33 * pow(10, 6);
  ;  // Set your desired luminosity here nb^-1
  std::cout << "Computed luminosity (nb^-1): " << computeLuminosity(charge) << std::endl;
  std::cout << "Computed luminosity (nb^-1): " << luminosity << std::endl;
  // -----------------------------can
  // -----------------------------
  // 1) Exclusive measured K+K- analysis (toggle)
  if (runExclusive) {
    // Initialize RDataFrames up front (shared by both modes)
    ROOT::RDF::RNode df_afterFid_sp18inb_data_init = InitKinematics(filename_afterFid_SP18inb_data, "dfSelected_afterFid_afterCorr", beam_energy_sp2018);
    ROOT::RDF::RNode df_afterFid_sp18outb_data_init = InitKinematics(filename_afterFid_SP18outb_data, "dfSelected_afterFid_afterCorr", beam_energy_sp2018);
    // ROOT::RDF::RNode df_afterFid_sp18outb_data_init_noqadb = InitKinematics(filename_afterFid_SP18outb_data_no_qadb, "dfSelected_afterFid_afterCorr", beam_energy_sp2018);
    ROOT::RDF::RNode df_afterFid_fall18inb_data_init = InitKinematics(filename_afterFid_Fall18inb_data, "dfSelected_afterFid_afterCorr", beam_energy_fall2018);
    ROOT::RDF::RNode df_afterFid_fall18outb_data_init = InitKinematics(filename_afterFid_Fall18outb_data, "dfSelected_afterFid_afterCorr", beam_energy_fall2018);
    ROOT::RDF::RNode df_afterFid_sp19inb_data_init = InitKinematics(filename_afterFid_SP19inb_data, "dfSelected_afterFid_afterCorr", beam_energy_sp2019);

    // MC RECONSTRUCTED GEMC
    // ROOT::RDF::RNode df_afterFid_sp18inb_MC_init = InitKinematics(filename_afterFid_SP18inb_MC, "dfSelected_afterFid_afterCorr", beam_energy_sp2018);
    // auto df_afterFid_sp18inb_MC = GetSlim_exclusive(df_afterFid_sp18inb_MC_init, "slim_sp18_MC_exclusive_qadb.root", "slim_sp18_MC_exclusive_qadb");
    ROOT::RDF::RNode df_afterFid_sp18inb_MC_init = InitKinematics(filename_afterFid_SP18inb_MC, "dfSelected_afterFid_afterCorr", beam_energy_sp2018);
    auto df_afterFid_sp18inb_MC = GetSlim_exclusive(df_afterFid_sp18inb_MC_init, "slim_sp18_MC_exclusive_qadb.root", "slim_sp18_MC_exclusive_qadb");
    ROOT::RDF::RNode df_afterFid_sp19inb_MC_init = InitKinematics(filename_afterFid_SP19inb_MC, "dfSelected_afterFid_afterCorr", beam_energy_sp2019);

    auto df_afterFid_sp19inb_MC = GetSlim_exclusive(df_afterFid_sp19inb_MC_init, "slim_sp19_MC_exclusive_qadb.root", "slim_sp19_MC_exclusive_qadb");
    ROOT::RDF::RNode df_afterFid_sp19inb_MC_clasdis_init = InitKinematics(filename_afterFid_SP19inb_MC_clasdis, "dfSelected", beam_energy_sp2019);
    auto df_afterFid_sp19inb_MC_clasdis = GetSlim_exclusive(df_afterFid_sp19inb_MC_clasdis_init, "slim_sp19_MC_clasdis_exclusive_qadb.root", "slim_sp19_MC_clasdis_exclusive_qadb");

    ROOT::RDF::RNode df_afterFid_sp19inb_MC_Harut_p39_init = InitKinematics(filename_afterFid_SP19inb_MC_Harut_p39, "dfSelected", beam_energy_sp2019);
    auto df_afterFid_sp19inb_MC_Harut_p39 =
        GetSlim_exclusive(df_afterFid_sp19inb_MC_Harut_p39_init, "slim_sp19_MC_Harut_p39_exclusive_qadb.root", "slim_sp19_MC_Harut_p39_exclusive_qadb");

    auto df_afterFid_fall18inb_MC_init = InitKinematics(filename_afterFid_Fall18inb_MC, "dfSelected_afterFid_afterCorr", beam_energy_fall2018);
    auto df_afterFid_fall18inb_MC = GetSlim_exclusive(df_afterFid_fall18inb_MC_init, "slim_fall18inb_MC_exclusive_qadb.root", "slim_fall18inb_MC_exclusive_qadb");
    auto df_afterFid_fall18outb_MC_init = InitKinematics(filename_afterFid_Fall18outb_MC, "dfSelected_afterFid_afterCorr", beam_energy_fall2018);
    auto df_afterFid_fall18outb_MC = GetSlim_exclusive(df_afterFid_fall18outb_MC_init, "slim_fall18outb_MC_exclusive_qadb.root", "slim_fall18outb_MC_exclusive_qadb");
    auto df_afterFid_sp18outb_MC_init = InitKinematics(filename_afterFid_SP18outb_MC, "dfSelected_afterFid_afterCorr", beam_energy_sp2018);
    auto df_afterFid_sp18outb_MC = GetSlim_exclusive(df_afterFid_sp18outb_MC_init, "slim_sp18outb_MC_exclusive_qadb.root", "slim_sp18outb_MC_exclusive_qadb");

    // MC RECONSTRUCTED GEMC for efficiencies
    auto df_afterFid_SP18inb_MC_eff_nobkg_init = InitKinematics(filename_afterFid_SP18inb_MC_eff_nobkg, "dfSelected_afterFid_afterCorr", beam_energy_sp2018);
    auto df_afterFid_SP18inb_MC_eff_nobkg = GetSlim_exclusive(df_afterFid_SP18inb_MC_eff_nobkg_init, "slim_sp18inb_MC_eff_nobkg.root", "slim_sp18inb_MC_eff_nobkg");
    auto df_afterFid_SP18inb_MC_eff_50nA_init = InitKinematics(filename_afterFid_SP18inb_MC_eff_50nA, "dfSelected_afterFid_afterCorr", beam_energy_sp2018);
    auto df_afterFid_SP18inb_MC_eff_50nA = GetSlim_exclusive(df_afterFid_SP18inb_MC_eff_50nA_init, "slim_sp18inb_MC_eff_50nA.root", "slim_sp18inb_MC_eff_50nA");
    auto df_afterFid_SP18outb_MC_eff_45nA_init = InitKinematics(filename_afterFid_SP18outb_MC_eff_45nA, "dfSelected_afterFid_afterCorr", beam_energy_sp2018);
    auto df_afterFid_SP18outb_MC_eff_45nA = GetSlim_exclusive(df_afterFid_SP18outb_MC_eff_45nA_init, "slim_sp18outb_MC_eff_45nA.root", "slim_sp18outb_MC_eff_45nA");
    auto df_afterFid_SP18outb_MC_eff_nobkg_init = InitKinematics(filename_afterFid_SP18outb_MC_eff_nobkg, "dfSelected_afterFid_afterCorr", beam_energy_sp2018);
    auto df_afterFid_SP18outb_MC_eff_nobkg = GetSlim_exclusive(df_afterFid_SP18outb_MC_eff_nobkg_init, "slim_sp18outb_MC_eff_nobkg.root", "slim_sp18outb_MC_eff_nobkg");
    auto df_afterFid_SP19inb_MC_eff_50nA_init = InitKinematics(filename_afterFid_SP19inb_MC_eff_50nA, "dfSelected_afterFid_afterCorr", beam_energy_sp2019);
    auto df_afterFid_SP19inb_MC_eff_50nA = GetSlim_exclusive(df_afterFid_SP19inb_MC_eff_50nA_init, "slim_sp19inb_MC_eff_50nA.root", "slim_sp19inb_MC_eff_50nA");
    auto df_afterFid_SP19inb_MC_eff_nobkg_init = InitKinematics(filename_afterFid_SP19inb_MC_eff_nobkg, "dfSelected_afterFid_afterCorr", beam_energy_sp2019);
    auto df_afterFid_SP19inb_MC_eff_nobkg = GetSlim_exclusive(df_afterFid_SP19inb_MC_eff_nobkg_init, "slim_sp19inb_MC_eff_nobkg.root", "slim_sp19inb_MC_eff_nobkg");

    // MC GEN
    ROOT::RDF::RNode df_afterFid_sp18inb_MCgen = InitGenKinematics(filename_afterFid_SP18inb_MCgen, "dfSelectedMC", beam_energy_sp2018);
    ROOT::RDF::RNode df_afterFid_sp19inb_MCgen = InitGenKinematics(filename_afterFid_SP19inb_MCgen, "dfSelectedMC", beam_energy_sp2019);
    ROOT::RDF::RNode df_afterFid_sp19inb_MCgen_Harut_p39 = InitGenKinematics(filename_afterFid_SP19inb_MCgen_Harut_p39, "dfSelectedMC", beam_energy_sp2019);
    ROOT::RDF::RNode df_afterFid_fall18inb_MCgen = InitGenKinematics(filename_afterFid_Fall18inb_MCgen, "dfSelectedMC", beam_energy_fall2018);
    ROOT::RDF::RNode df_afterFid_fall18outb_MCgen = InitGenKinematics(filename_afterFid_Fall18outb_MCgen, "dfSelectedMC", beam_energy_fall2018);
    ROOT::RDF::RNode df_afterFid_sp18outb_MCgen = InitGenKinematics(filename_afterFid_SP18outb_MCgen, "dfSelectedMC", beam_energy_sp2018);

    auto df_afterFid_sp19inb_data = GetSlim_exclusive(df_afterFid_sp19inb_data_init, "slim_sp19_exlcusive_qadb.root", "slim_sp19_exlcusive_qadb");
    auto df_afterFid_sp18inb_data = GetSlim_exclusive(df_afterFid_sp18inb_data_init, "slim_sp18inb_exlcusive_qadb.root", "slim_sp18inb_exlcusive_qadb");
    auto df_afterFid_sp18outb_data = GetSlim_exclusive(df_afterFid_sp18outb_data_init, "slim_sp18outb_exlcusive_qadb.root", "slim_sp18outb_exlcusive_qadb");
    // auto df_afterFid_sp18outb_data_noqadb = GetSlim_exclusive(df_afterFid_sp18outb_data_init_noqadb, "slim_sp18outb_exlcusive.root", "slim_sp18outb_exlcusive");
    auto df_afterFid_fall18inb_data = GetSlim_exclusive(df_afterFid_fall18inb_data_init, "slim_fall18inb_exlcusive_qadb.root", "slim_fall18inb_exlcusive_qadb");
    auto df_afterFid_fall18outb_data = GetSlim_exclusive(df_afterFid_fall18outb_data_init, "slim_fall18outb_exlcusive_qadb.root", "slim_fall18outb_exlcusive_qadb");

    // Apply your “final” DVEP selections and then pick exclusive phi events
    auto df_sp18inb_all = SelectExclusivePhiEvent(df_afterFid_sp18inb_data);
    auto df_sp18outb_all = SelectExclusivePhiEvent(df_afterFid_sp18outb_data);
    // auto df_sp18outb_all_noqadb = SelectExclusivePhiEvent(df_afterFid_sp18outb_data_noqadb);
    auto df_fall18inb_all = SelectExclusivePhiEvent(df_afterFid_fall18inb_data);
    auto df_fall18outb_all = SelectExclusivePhiEvent(df_afterFid_fall18outb_data);
    auto df_sp19inb_all = SelectExclusivePhiEvent(df_afterFid_sp19inb_data);

    // Apply final selections MC
    auto df_sp18inb_all_MC = SelectExclusivePhiEvent(df_afterFid_sp18inb_MC);
    auto df_sp19inb_all_MC = SelectExclusivePhiEvent(df_afterFid_sp19inb_MC);
    auto df_sp19inb_all_MC_clasdis = SelectExclusivePhiEvent(df_afterFid_sp19inb_MC_clasdis);
    auto df_sp19inb_all_MC_Harut_p39 = SelectExclusivePhiEvent(df_afterFid_sp19inb_MC_Harut_p39);
    auto df_fall18inb_all_MC = SelectExclusivePhiEvent(df_afterFid_fall18inb_MC);
    auto df_fall18outb_all_MC = SelectExclusivePhiEvent(df_afterFid_fall18outb_MC);
    auto df_sp18outb_all_MC = SelectExclusivePhiEvent(df_afterFid_sp18outb_MC);

    // Apply final selections MC rec for efficiencies
    auto df_afterFid_SP18inb_MC_eff_nobkg_final = SelectExclusivePhiEvent(df_afterFid_SP18inb_MC_eff_nobkg);
    auto df_afterFid_SP18inb_MC_eff_50nA_final = SelectExclusivePhiEvent(df_afterFid_SP18inb_MC_eff_50nA);
    auto df_afterFid_SP18outb_MC_eff_45nA_final = SelectExclusivePhiEvent(df_afterFid_SP18outb_MC_eff_45nA);
    auto df_afterFid_SP18outb_MC_eff_nobkg_final = SelectExclusivePhiEvent(df_afterFid_SP18outb_MC_eff_nobkg);
    auto df_afterFid_SP19inb_MC_eff_50nA_final = SelectExclusivePhiEvent(df_afterFid_SP19inb_MC_eff_50nA);
    auto df_afterFid_SP19inb_MC_eff_nobkg_final = SelectExclusivePhiEvent(df_afterFid_SP19inb_MC_eff_nobkg);

    // MC gen

    auto df_sp18inb_phi = ApplyFinalDVEPSelections(df_sp18inb_all);
    // auto df_sp18outb_phi_noqadb = ApplyFinalDVEPSelections(df_sp18outb_all_noqadb);
    auto df_sp18outb_phi = ApplyFinalDVEPSelections(df_sp18outb_all);
    auto df_fall18inb_phi = ApplyFinalDVEPSelections(df_fall18inb_all);
    auto df_fall18outb_phi = ApplyFinalDVEPSelections(df_fall18outb_all);
    auto df_sp19inb_phi = ApplyFinalDVEPSelections(df_sp19inb_all);

    /// MC
    auto df_sp18inb_phi_MC = ApplyFinalDVEPSelections(df_sp18inb_all_MC);
    auto df_sp19inb_phi_MC = ApplyFinalDVEPSelections(df_sp19inb_all_MC);
    auto df_sp19inb_phi_MC_Harut_p39 = ApplyFinalDVEPSelections(df_sp19inb_all_MC_Harut_p39);
    auto df_fall18inb_phi_MC = ApplyFinalDVEPSelections(df_fall18inb_all_MC);
    auto df_fall18outb_phi_MC = ApplyFinalDVEPSelections(df_fall18outb_all_MC);
    auto df_sp18outb_phi_MC = ApplyFinalDVEPSelections(df_sp18outb_all_MC);

    // MC rec for efficiencies
    auto df_sp18inb_phi_MC_eff_nobkg = ApplyFinalDVEPSelections(df_afterFid_SP18inb_MC_eff_nobkg_final);
    auto df_sp18inb_phi_MC_eff_50nA = ApplyFinalDVEPSelections(df_afterFid_SP18inb_MC_eff_50nA_final);
    auto df_sp18outb_phi_MC_eff_45nA = ApplyFinalDVEPSelections(df_afterFid_SP18outb_MC_eff_45nA_final);
    auto df_sp18outb_phi_MC_eff_nobkg = ApplyFinalDVEPSelections(df_afterFid_SP18outb_MC_eff_nobkg_final);
    auto df_sp19inb_phi_MC_eff_50nA = ApplyFinalDVEPSelections(df_afterFid_SP19inb_MC_eff_50nA_final);
    auto df_sp19inb_phi_MC_eff_nobkg = ApplyFinalDVEPSelections(df_afterFid_SP19inb_MC_eff_nobkg_final);
    // auto df_sp19inb_phi_MC_clasdis = df_afterFid_sp19inb_MC_clasdis;

    // Gen DVEP selections
    auto df_afterFid_sp19inb_MCgen_Harut_p39final = ApplyFinalGenDVEPSelections(df_afterFid_sp19inb_MCgen_Harut_p39);
    auto df_afterFid_sp19inb_MCgenfinal = ApplyFinalGenDVEPSelections(df_afterFid_sp19inb_MCgen);
    auto df_afterFid_sp18inb_MCgenfinal = ApplyFinalGenDVEPSelections(df_afterFid_sp18inb_MCgen);
    auto df_afterFid_fall18inb_MCgenfinal = ApplyFinalGenDVEPSelections(df_afterFid_fall18inb_MCgen);
    auto df_afterFid_fall18outb_MCgenfinal = ApplyFinalGenDVEPSelections(df_afterFid_fall18outb_MCgen);
    auto df_afterFid_sp18outb_MCgenfinal = ApplyFinalGenDVEPSelections(df_afterFid_sp18outb_MCgen);

    if (plotIntInvMass) {
      /* code */
      // PlotEventOverview(df_sp19inb_phi, "PhiMassPlots/spring2019/inb/DVKpKm", "exclusiveKpKm");
      DISANA::PhiMass::DrawPhiMass_Measured(df_sp19inb_phi, "PhiMassPlots/spring2019/inb/DVKpKm/exclusiveKpKm/", "Sp19 INB DVKpKm", /*nBins*/ 200, /*mMin*/ 0.8, /*mMax*/ 1.6,
                                            /*mPhiLo*/ 0.987, /*mPhiHi*/ 1.2, /*nSigma*/ 8.0);
      DISANA::PhiMass::DrawPhiMass_Measured(df_fall18inb_phi, "PhiMassPlots/fall2018/inb/DVKpKm/exclusiveKpKm/", "Fall18 INB DVKpKm", /*nBins*/ 200, /*mMin*/ 0.8, /*mMax*/ 1.6,
                                            /*mPhiLo*/ 0.987, /*mPhiHi*/ 1.2, /*nSigma*/ 8.0);
      DISANA::PhiMass::DrawPhiMass_Measured(df_sp18inb_phi, "PhiMassPlots/spring2018/inb/DVKpKm/exclusiveKpKm/", "Sp18 INB DVKpKm", /*nBins*/ 200, /*mMin*/ 0.8, /*mMax*/ 1.6,
                                            /*mPhiLo*/ 0.9874, /*mPhiHi*/ 1.2, /*nSigma*/ 8.0);
      DISANA::PhiMass::DrawPhiMass_Measured(df_sp18outb_phi, "PhiMassPlots/spring2018/outb/DVKpKm/exclusiveKpKm/", "Sp18 OUTB DVKpKm ", /*nBins*/ 200, /*mMin*/ 0.8, /*mMax*/ 1.6,
                                            /*mPhiLo*/ 0.9874, /*mPhiHi*/ 1.2, /*nSigma*/ 8.0);
      // DISANA::PhiMass::DrawPhiMass_Measured(df_sp18outb_phi_noqadb, "PhiMassPlots/spring2018/outb_noqadb/DVKpKm/exclusiveKpKm/", "Sp18 OUTB DVKpKm no qa", /*nBins*/ 200,
      //                                   /*mMin*/ 0.8, /*mMax*/ 1.6, /*mPhiLo*/ 0.9874, /*mPhiHi*/ 1.2, /*nSigma*/ 8.0);
      DISANA::PhiMass::DrawPhiMass_Measured(df_fall18outb_phi, "PhiMassPlots/fall2018/outb/DVKpKm/exclusiveKpKm/", "Fall18 OUTB DVKpKm", /*nBins*/ 200, /*mMin*/ 0.8, /*mMax*/ 1.6,
                                            /*mPhiLo*/ 0.9874, /*mPhiHi*/ 1.2, /*nSigma*/ 8.0);

      // MC RECONSTRUCTED
      DISANA::PhiMass::DrawPhiMass_Measured(df_sp18inb_phi_MC, "PhiMassPlots/spring2018/inb/MC_reconstructed/DVKpKm/", "Sp18 INB DVKpKm MC Reconstructed", /*nBins*/ 200,
                                            /*mMin*/ 0.8, /*mMax*/ 1.8,
                                            /*mPhiLo*/ 0.987, /*mPhiHi*/ 1.2, /*nSigma*/ 8.0);
      DISANA::PhiMass::DrawPhiMass_Measured(df_sp19inb_phi_MC, "PhiMassPlots/spring2019/inb/MC_reconstructed/DVKpKm/", "Sp19 INB DVKpKm MC Reconstructed", /*nBins*/ 200,
                                            /*mMin*/ 0.8, /*mMax*/ 1.8,
                                            /*mPhiLo*/ 0.987, /*mPhiHi*/ 1.2, /*nSigma*/ 8.0);
      DISANA::PhiMass::DrawPhiMass_Measured(df_sp19inb_phi_MC_Harut_p39, "PhiMassPlots/spring2019/inb/MC_reconstructed/DVKpKm/Harut_p39/",
                                            "Sp19 INB DVKpKm MC Reconstructed Harut p39", /*nBins*/ 200,
                                            /*mMin*/ 0.8, /*mMax*/ 1.8,
                                            /*mPhiLo*/ 0.987, /*mPhiHi*/ 1.2, /*nSigma*/ 8.0);

      // DISANA::PhiMass::DrawPhiMass_Measured(df_sp19inb_phi_MC_clasdis, "PhiMassPlots/spring2019/inb/MC_reconstructed/DVKpKm/clasdis/", "Sp19 INB DVKpKm MC Reconstructed
      // clasdis", /*nBins*/ 200, /*mMin*/ 0.8, /*mMax*/ 1.8,
      //                                    /*mPhiLo*/ 0.987, /*mPhiHi*/ 1.2, /*nSigma*/ 8.0);
      DISANA::PhiMass::DrawPhiMass_Measured(df_fall18inb_phi_MC, "PhiMassPlots/fall2018/inb/MC_reconstructed/DVKpKm/", "Fall18 INB DVKpKm MC Reconstructed", /*nBins*/ 200,
                                            /*mMin*/ 0.8, /*mMax*/ 1.8,
                                            /*mPhiLo*/ 0.987, /*mPhiHi*/ 1.2, /*nSigma*/ 8.0);
      DISANA::PhiMass::DrawPhiMass_Measured(df_fall18outb_phi_MC, "PhiMassPlots/fall2018/outb/MC_reconstructed/DVKpKm/", "Fall18 OUTB DVKpKm MC Reconstructed", /*nBins*/ 200,
                                            /*mMin*/ 0.8, /*mMax*/ 1.8,
                                            /*mPhiLo*/ 0.987, /*mPhiHi*/ 1.2, /*nSigma*/ 8.0);
      DISANA::PhiMass::DrawPhiMass_Measured(df_sp18outb_phi_MC, "PhiMassPlots/spring2018/outb/MC_reconstructed/DVKpKm/", "Sp18 OUTB DVKpKm MC Reconstructed", /*nBins*/ 200,
                                            /*mMin*/ 0.8, /*mMax*/ 1.8,
                                            /*mPhiLo*/ 0.987, /*mPhiHi*/ 1.2, /*nSigma*/ 8.0);
    }

    // Add any datasets you want to compare
    // Binning setup for phi analysis
    auto hQ2t_fall18outb = df_fall18outb_phi.Histo2D({"hQ2t_fall18outb", "t' vs Q^{2};t' [GeV^{2}];Q^{2} [GeV^{2}]", 100, 0.01, 4.0,  // x:  t' range
                                                      100, 0.8, 6.0},                                                                 // y:  Q^2 range
                                                     "mtprime", "Q2"                                                                 // x var, y var
    );
    // MC Reconstructed
    // auto hQ2t_sp18inb_MC = df_sp18inb_phi_MC.Histo2D({"hQ2t_sp18inb_MC", "t' vs Q^{2};t' [GeV^{2}];Q^{2} [GeV^{2}]", 60, 0.01, 4.5,  // x:  t' range
    //                                                  60, 0.1, 8.0},                                                                // y:  Q^2 range
    //                                                  "mtprime", "Q2"                                                                // x var, y var
    // );
    auto* hptr = hQ2t_fall18outb.GetPtr();
    std::cout << "[DEBUG] hQ2t_fall18outb entries = " << hptr->GetEntries() << " integral = " << hptr->Integral() << std::endl;

    int nQ2Bins = 5;       // n Q^2 bins
    int nTprimeBins = 10;  // m t' bins

    EqualStatBinningResult eqBins = xBins.MakeEqualStatBinning(hQ2t_fall18outb.GetPtr(), nQ2Bins, nTprimeBins);
    std::cout << "Q2 edges size = " << eqBins.q2Edges.size() << std::endl;
    std::cout << "t' edges size = " << eqBins.tprimeEdges.size() << std::endl;

    PrintEdges("Q^{2}", eqBins.q2Edges);
    PrintEdges("t'", eqBins.tprimeEdges);

    xBins.DrawQ2TprimeWithGrid(hQ2t_fall18outb.GetPtr(), eqBins.q2Edges, {0.0100, 0.329, 0.401, 0.52, 0.64, 0.76, 0.9250, 1.1080, 1.38, 1.75, 2.2, 2.7, 4.0},
                               "cQ2t_equalStat", "PhiMassPlots/Q2_vs_tprime_equalStatBins.pdf");
    // Fix t bin bin
    xBins.SetQ2Bins(eqBins.q2Edges);
    xBins.SetTprimeBins({0.0100, 0.329, 0.401, 0.52, 0.64, 0.76, 0.9250, 1.1080, 1.38, 1.75, 2.2, 2.7, 4.0});
    // keep your xB, W binning as-is
    xBins.SetWBins({1.8, 4.0});
    //xBins.SetTprimeBins(eqBins.tprimeEdges);
    {
      std::filesystem::create_directories(outputDir);

      // Compute mean Q^2 and mean x_B per Q^2 bin from each data sample.
      // Both quantities are booked in a single lazy RDataFrame event loop.
      const auto means_sp19 = ComputeMeanKinPerQ2Bin(df_sp19inb_phi, eqBins.q2Edges);
      const auto means_sp18 = ComputeMeanKinPerQ2Bin(df_sp18inb_phi, eqBins.q2Edges);

      // Spring 2019 inb  (10.2 GeV)
      PlotVmCut(eqBins.q2Edges, means_sp19, static_cast<double>(beam_energy_sp2019), m_phi, outputDir + "/vm_cut_sp19inb_10p2GeV.pdf");

      // Spring 2018 inb  (10.6 GeV)
      PlotVmCut(eqBins.q2Edges, means_sp18, static_cast<double>(beam_energy_sp2018), m_phi, outputDir + "/vm_cut_sp18inb_10p6GeV.pdf");

      // ── Egamma* and Emiss distributions per Q^2 bin ────────────────────
      PlotRadiativeKinematics(df_sp19inb_phi, eqBins.q2Edges, means_sp19, outputDir + "/rad_kinematics_sp19inb_10p2GeV.pdf");

      PlotRadiativeKinematics(df_sp18inb_phi, eqBins.q2Edges, means_sp18, outputDir + "/rad_kinematics_sp18inb_10p6GeV.pdf");
      // 2D heat map — integrated over all Q²
      PlotInvMassVsMissingMass(df_sp19inb_phi, outputDir + "/invmass_vs_missmass.pdf");

      // Side-by-side 1D projections from the same 2D distribution
      PlotInvMassVsMissingMassProjections(df_sp19inb_phi, outputDir + "/invmass_missmass_projections.pdf",
                                          /*nBinsX=*/100, /*nBinsY=*/80,
                                          /*xLo=*/0.98, /*xHi=*/1.10,
                                          /*yLo=*/0.90, /*yHi=*/1.15,
                                          /*sigWin=*/0.010);  // ← widen to 0.015 if stats are low
    }

    {
      double W2lo = 3.24;       // W > 1.8 GeV
      double W2hi = 4.0 * 4.0;  // W < 3.5 GeV
      std::filesystem::create_directories(outputDir);

      // 10.2 GeV CSV
      std::string csv102 = outputDir + "/diffrad_binning_10p2GeV.csv";
      ExportBinningCSV(eqBins.q2Edges, {0.0100, 0.3760, 0.4675, 0.5590, 0.6505, 0.8335, 0.9250, 1.1080, 1.3825, 1.8400, 2.2, 2.7, 4.0}, static_cast<double>(beam_energy_sp2019),
                       W2lo, W2hi, csv102);

      // 10.6 GeV CSV
      std::string csv106 = outputDir + "/diffrad_binning_10p6GeV.csv";
      ExportBinningCSV(eqBins.q2Edges, {0.0100, 0.3760, 0.4675, 0.5590, 0.6505, 0.8335, 0.9250, 1.1080, 1.3825, 1.8400, 2.2, 2.7, 4.0}, static_cast<double>(beam_energy_fall2018),
                       W2lo, W2hi, csv106);
    }

    // ── Load rad-corr RDF trees from RunMDiffradNew output ROOT files ─────────
    std::string kRadRootFile_10p2 = "/w/hallb-scshelf2102/clas12/singh/Softwares/Generators/PhiEventGen/DIFFRAD/diffradstn/outputs/10p2GeV/mdiffrad_output.root";
    std::string kRadRootFile_10p6 = "/w/hallb-scshelf2102/clas12/singh/Softwares/Generators/PhiEventGen/DIFFRAD/diffradstn/outputs/10p6GeV/mdiffrad_output.root";
    if (TryHStyle) {
      kRadRootFile_10p2 = "/w/hallb-scshelf2102/clas12/singh/Softwares/Generators/PhiEventGen/DIFFRAD/diffradstn/outputs/HS_rad/10p2GeV/mdiffrad_output.root";
      kRadRootFile_10p6 = "/w/hallb-scshelf2102/clas12/singh/Softwares/Generators/PhiEventGen/DIFFRAD/diffradstn/outputs/HS_rad/10p6GeV/mdiffrad_output.root";
    } else {
      kRadRootFile_10p2 = "/w/hallb-scshelf2102/clas12/singh/Softwares/Generators/PhiEventGen/DIFFRAD/diffradstn/outputs/10p2GeV/mdiffrad_output.root";
      kRadRootFile_10p6 = "/w/hallb-scshelf2102/clas12/singh/Softwares/Generators/PhiEventGen/DIFFRAD/diffradstn/outputs/10p6GeV/mdiffrad_output.root";
    }
    // Helper: load radcorr_rdf tree, return RNode. Falls back to a 1-row dummy
    // (rad_corr=1.0) if the file is not yet available so the rest of the analysis
    // still compiles and runs without corrections applied.
    auto LoadRadRDF = [](const std::string& rootFile) -> ROOT::RDF::RNode {
      if (std::filesystem::exists(rootFile)) {
        std::cout << "[RadCorr] Loading radcorr_rdf from " << rootFile << "\n";
        return ROOT::RDataFrame("radcorr_rdf", rootFile);
      }
      std::cerr << "[RadCorr] WARNING: " << rootFile << " not found — using rad_corr=1.0 (no correction).\n"
                << "  Run RunMDiffradNew and set the output path above.\n";
      // Return a minimal dummy RDF with the required columns so downstream
      // code does not crash.  MakePhiRadiativeCorrection3D will fall back to 1.0
      // for all bins with fewer than minEntries (default 10) entries.
      return ROOT::RDataFrame(0)
          .Define("Q2", [] { return 1.0; })
          .Define("W", [] { return 3.0; })
          .Define("mtprime", [] { return 0.1; })
          .Define("rad_corr", [] { return 1.0; })
          .Define("beam_mom", [] { return 10.6; });
    };

    ROOT::RDF::RNode df_rad_10p2 = LoadRadRDF(kRadRootFile_10p2);
    ROOT::RDF::RNode df_rad_10p6 = LoadRadRDF(kRadRootFile_10p6);

    // Period → rad RDF mapping
    ROOT::RDF::RNode& df_rad_sp19inb = df_rad_10p2;  // 10.2 GeV
    ROOT::RDF::RNode& df_rad_sp18inb = df_rad_10p6;  // 10.6 GeV
    ROOT::RDF::RNode& df_rad_sp18outb = df_rad_10p6;
    ROOT::RDF::RNode& df_rad_fall18inb = df_rad_10p6;
    ROOT::RDF::RNode& df_rad_fall18outb = df_rad_10p6;

    comparer.SetXBinsRanges(xBins);
    // toggle to skip reconstructed and only compare gen-level distributions to models (e.g. for acceptance/efficiency studies)
    // ── AddModelPhi — pass rad RDF + enable rad correction per dataset ────────
    if (RecDataOnly) {
      comparer.AddModelPhi(df_sp19inb_phi, "Sp19 inb", beam_energy_sp2019, luminosity_rga_sp19_inb, df_afterFid_sp19inb_MCgenfinal, df_sp19inb_phi_MC, df_sp19inb_phi_MC_eff_50nA,
                           df_sp19inb_phi_MC_eff_nobkg, df_rad_sp19inb,
                           /*doAcc=*/false, /*doEff=*/false, /*doRadCorr=*/false);
      // comparer.AddModelPhi(df_sp18inb_phi, "Sp18 inb", beam_energy_sp2018, luminosity_rga_sp18_inb, df_afterFid_sp18inb_MCgenfinal, df_sp18inb_phi_MC,
      // df_sp18inb_phi_MC_eff_nobkg, df_sp18inb_phi_MC_eff_50nA, df_rad_sp18inb,
      //                     /*doAcc=*/false, /*doEff=*/false, /*doRadCorr=*/false);
      // comparer.AddModelPhi(df_sp18outb_phi, "Sp18 outb", beam_energy_sp2018, luminosity_rga_sp18_outb, df_afterFid_sp18outb_MCgenfinal, df_sp18outb_phi_MC,
      // df_sp18outb_phi_MC_eff_45nA, df_sp18outb_phi_MC_eff_nobkg, df_rad_sp18outb,
      //                    /*doAcc=*/false, /*doEff=*/false, /*doRadCorr=*/false);

    } else if (GenOnly) {
      xBins.SetQ2Bins({0.8, 8.0});
      comparer.AddModelPhi(df_afterFid_sp19inb_MCgen_Harut_p39final, "JLab MC", beam_energy_sp2019, luminosity_rga_sp19_inb);
      comparer.AddModelPhi(df_afterFid_sp19inb_MCgenfinal, "Phi Eventgen", beam_energy_sp2019, luminosity_rga_sp19_inb);
      /*comparer.AddModelPhi(df_afterFid_sp19inb_MCgenfinal, "Sp18 inb gen", beam_energy_sp2018, luminosity_rga_sp18_inb);
      comparer.AddModelPhi(df_afterFid_sp19inb_MCgenfinal, "Sp18 outb gen", beam_energy_sp2018, luminosity_rga_sp18_outb);
      comparer.AddModelPhi(df_afterFid_sp19inb_MCgenfinal, "Fall18 inb gen", beam_energy_fall2018, luminosity_rga_fall18_inb);
      comparer.AddModelPhi(df_afterFid_sp19inb_MCgenfinal, "Fall18 outb gen", beam_energy_fall2018, luminosity_rga_fall18_outb);*/
    } else if (RecGEMCOnly) {
      comparer.AddModelPhi(df_sp19inb_phi_MC_Harut_p39, "JLab MC", beam_energy_sp2019, luminosity_rga_sp19_inb);
      comparer.AddModelPhi(df_sp19inb_phi_MC, "Phi Eventgen", beam_energy_sp2019, luminosity_rga_sp19_inb);
      comparer.AddModelPhi(df_sp19inb_phi, "Sp19 inb (data)", beam_energy_sp2019, luminosity_rga_sp19_inb);
      /*comparer.AddModelPhi(df_afterFid_sp19inb_MCgenfinal, "Sp18 inb gen", beam_energy_sp2018, luminosity_rga_sp18_inb);
      comparer.AddModelPhi(df_afterFid_sp19inb_MCgenfinal, "Sp18 outb gen", beam_energy_sp2018, luminosity_rga_sp18_outb);
      comparer.AddModelPhi(df_afterFid_sp19inb_MCgenfinal, "Fall18 inb gen", beam_energy_fall2018, luminosity_rga_fall18_inb);
      comparer.AddModelPhi(df_afterFid_sp19inb_MCgenfinal, "Fall18 outb gen", beam_energy_fall2018, luminosity_rga_fall18_outb);*/
    }

    // outb (no MC gen/rec available — rad corr only)
    // comparer.AddModelPhi(df_fall18outb_phi, "Fall18 outb", beam_energy_fall2018, luminosity_rga_fall18_outb);
    // comparer.AddModelPhi(df_sp18outb_phi, "Sp18 outb", beam_energy_sp2018, luminosity_rga_sp18_outb);
    // inb
    // comparer.AddModelPhi(df_fall18inb_phi, "Fall18 inb", beam_energy_fall2018, luminosity_rga_fall18_inb);

    // MC RECONSTRUCTED
    // comparer.AddModelPhi(df_sp18inb_phi_MC, "Sp18 inb MC Recon", beam_energy_sp2018, luminosity_rga_sp18_inb);
    // comparer.AddModelPhi(df_sp19inb_phi_MC, "Sp19inb MC Rec", beam_energy_sp2019, luminosity_rga_sp19_inb);
    // comparer.AddModelPhi(df_sp19inb_phi_MC_clasdis, "Sp19inb MC Rec clasdis", beam_energy_sp2019, luminosity_rga_sp19_inb);
    // DumpExclusiveTxt(df_sp18inb_phi_MC, "exclusive_events_dump_sp2018_inb_MC_reconstructed.txt");
    // DumpExclusiveTxt(df_sp19inb_phi_MC, "exclusive_events_dump_sp2019_inb_MC_reconstructed.txt");
  }

  if (runMissingMass) {
    // bining is done same as the exlcusive case,I am considering only sp2019 inb
    if (!runExclusive) {
      xBins.SetQ2Bins({0.0000, 1.2000, 1.7333, 2.5333, 8.0000});
      xBins.SetTprimeBins({0.0000, 0.3750, 0.4500, 0.6000, 0.6750, 0.9000, 1.0500, 1.3500, 1.7250, 4.5000});
      comparer.SetXBinsRanges(xBins);
    }
    // ROOT::RDF::RNode df_afterFid_sp19inb_missingKm_data_org = InitKinematics_MissingKm(filename_afterFid_SP19inb_missingKm_data_org, "dfSelected", beam_energy_sp2019);

    ROOT::RDF::RNode df_afterFid_sp18inb_missingKm_data_init =
        InitKinematics_MissingKm(filename_afterFid_SP18inb_missingKm_data, "dfSelected_afterFid_afterCorr", beam_energy_sp2018);
    std::cout << "[Debug] Initialized RDataFrame for SP18 INB missing Km data" << std::endl;
    ROOT::RDF::RNode df_afterFid_fall18inb_missingKm_data_init =
        InitKinematics_MissingKm(filename_afterFid_Fall18inb_missingKm_data, "dfSelected_afterFid_afterCorr", beam_energy_fall2018);
    std::cout << "[Debug] Initialized RDataFrame for Fall18 INB missing Km data" << std::endl;
    ROOT::RDF::RNode df_afterFid_sp19inb_missingKm_data_init =
        InitKinematics_MissingKm(filename_afterFid_SP19inb_missingKm_data, "dfSelected_afterFid_afterCorr", beam_energy_sp2019);
    std::cout << "[Debug] Initialized RDataFrame for SP19 INB missing Km data" << std::endl;
    ROOT::RDF::RNode df_afterFid_sp18outb_missingKp_data_init =
        InitKinematics_MissingKp(filename_afterFid_SP18outb_missingKp_data, "dfSelected_afterFid_afterCorr", beam_energy_sp2018);
    std::cout << "[Debug] Initialized RDataFrame for SP18 OUTB missing Kp data" << std::endl;
    ROOT::RDF::RNode df_afterFid_fall18outb_missingKp_data_init = InitKinematics_MissingKp(filename_afterFid_Fall18outb_missingKp_data, "dfSelected_afterFid_afterCorr",
                                                                                           beam_energy_fall2018);  /// check this number as fall 18 is still misssing
    std::cout << "[Debug] Initialized RDataFrame for Fall18 OUTB missing Kp data" << std::endl;

    auto df_afterFid_sp19inb_missingKm_data = GetSlim_missingKm(df_afterFid_sp19inb_missingKm_data_init, "slim_sp19_missingKm_all.root", "slim_sp19_missingKm_all");
    std::cout << "[Debug] Obtained slim dataframe for SP19 INB missing Km data" << std::endl;
    auto df_afterFid_sp18inb_missingKm_data = GetSlim_missingKm(df_afterFid_sp18inb_missingKm_data_init, "slim_sp18_missingKm_all.root", "slim_sp18_missingKm_all");
    std::cout << "[Debug] Obtained slim dataframe for SP18 INB missing Km data" << std::endl;
    auto df_afterFid_sp18outb_missingKp_data = GetSlim_missingKp(df_afterFid_sp18outb_missingKp_data_init, "slim_sp18_missingKp_all.root", "slim_sp18_missingKp_all");
    std::cout << "[Debug] Obtained slim dataframe for SP18 OUTB missing Kp data" << std::endl;
    auto df_afterFid_fall18inb_missingKm_data = GetSlim_missingKm(df_afterFid_fall18inb_missingKm_data_init, "slim_fall18_missingKm_all.root", "slim_fall18_missingKm_all");
    std::cout << "[Debug] Obtained slim dataframe for Fall18 INB missing Km data" << std::endl;
    auto df_afterFid_fall18outb_missingKp_data = GetSlim_missingKp(df_afterFid_fall18outb_missingKp_data_init, "slim_fall18_missingKp_all.root", "slim_fall18_missingKp_all");

    // auto [df_cut_missing_Km_sp19_org, winKm2019_org] =
    //  DISANA::PhiMass::WireKaonMMCut_FromMx2(df_afterFid_sp19inb_missingKm_data_org, "Mx2_epKp", "PhiMassPlots/spring2019/inb/nsidis_org", "Sp19 INB", "K^{+}", 220, 0.25,
    //  0.75, 3.0);
    auto [df_cut_missing_Km_sp19, winKm2019] =
        DISANA::PhiMass::WireKaonMMCut_FromMx2(df_afterFid_sp19inb_missingKm_data, "Mx2_epKp", "PhiMassPlots/spring2019/inb/nsidis", "Sp19 INB", "K^{+}", 220, 0.25, 0.75, 3.0);
    auto [df_cut_missing_Km_fall18, winKmFall18] =
        DISANA::PhiMass::WireKaonMMCut_FromMx2(df_afterFid_fall18inb_missingKm_data, "Mx2_epKp", "PhiMassPlots/fall2018/inb/nsidis", "Fall18 INB", "K^{+}", 220, 0.25, 0.75, 3.0);
    auto [df_cut_missing_Km_sp18, winKmSp18] =
        DISANA::PhiMass::WireKaonMMCut_FromMx2(df_afterFid_sp18inb_missingKm_data, "Mx2_epKp", "PhiMassPlots/spring2018/inb/nsidis", "Sp18 INB", "K^{+}", 220, 0.25, 0.75, 3.0);
    auto [df_cut_missing_Kp_sp18, winKpSp18] =
        DISANA::PhiMass::WireKaonMMCut_FromMx2(df_afterFid_sp18outb_missingKp_data, "Mx2_epKm", "PhiMassPlots/spring2018/outb/nsidis", "Sp18 OUTB", "K^{-}", 220, 0.30, 0.80, 3.0);
    auto [df_cut_missing_Kp_fall18, winKpFall18] = DISANA::PhiMass::WireKaonMMCut_FromMx2(df_afterFid_fall18outb_missingKp_data, "Mx2_epKm", "PhiMassPlots/fall2018/outb/nsidis",
                                                                                          "Fall18 OUTB", "K^{-}", 220, 0.25, 0.75, 3.0);
    std::cout << "[Debug] Applied missing mass cuts and obtained cut dataframes" << std::endl;
    // Draw phi mass with Km missing
    // DISANA::PhiMass::DrawPhiMass_Kp_plus_KmMissing(df_cut_missing_Km_sp19_org, "PhiMassPlots/spring2019/inb/nsidis_org", "Sp19 INB", 200, 0.8, 1.8, 0.987, 1.2, 3.0);
    DISANA::PhiMass::DrawPhiMass_Kp_plus_KmMissing(df_cut_missing_Km_sp19, "PhiMassPlots/spring2019/inb/nsidis", "Sp19 INB", 200, 0.8, 1.8, 0.987, 1.2, 3.0);
    DISANA::PhiMass::DrawPhiMass_Kp_plus_KmMissing(df_cut_missing_Km_fall18, "PhiMassPlots/fall2018/inb/nsidis", "Fall18 INB", 200, 0.8, 1.6, 0.9874, 1.12, 3.0);
    DISANA::PhiMass::DrawPhiMass_Kp_plus_KmMissing(df_cut_missing_Km_sp18, "PhiMassPlots/spring2018/inb/nsidis", "Sp18 INB", 200, 0.8, 1.6, 0.9874, 1.12, 3.0);
    DISANA::PhiMass::DrawPhiMass_KpMissing_plus_Km(df_cut_missing_Kp_sp18, "PhiMassPlots/spring2018/outb/nsidis", "Sp18 OUTB", 200, 0.8, 1.6, 0.9874, 1.12, 3.0);
    DISANA::PhiMass::DrawPhiMass_KpMissing_plus_Km(df_cut_missing_Kp_fall18, "PhiMassPlots/fall2018/outb/nsidis", "Fall18 OUTB", 200, 0.8, 1.6, 0.9874, 1.12, 3.0);
    std::cout << "[Debug] Plotted phi mass with missing kaon for all datasets" << std::endl;
    // Pass through your final DVEP selections & dedicated phi picker for missing-Km streams
    // auto df_sp19_missingKm_phi_org = SelectPhiEvent_MissingKm(df_cut_missing_Km_sp19_org);
    auto df_sp19_missingKm_phi = SelectPhiEvent_MissingKm(df_cut_missing_Km_sp19);
    auto df_fall18_missingKm_phi = SelectPhiEvent_MissingKm(df_cut_missing_Km_fall18);
    auto df_sp18_missingKm_phi = SelectPhiEvent_MissingKm(df_cut_missing_Km_sp18);
    auto df_sp18_missingKp_phi = SelectPhiEvent_MissingKp(df_cut_missing_Kp_sp18);
    auto df_fall18_missingKp_phi = SelectPhiEvent_MissingKp(df_cut_missing_Kp_fall18);
    std::cout << "[Debug] Selected phi events with missing kaon for all datasets" << std::endl;
    // auto df_sp19_missingKm_all_org = ApplyFinalDVEPSelections(df_sp19_missingKm_phi_org);
    auto df_sp19_missingKm_all = ApplyFinalDVEPSelections(df_sp19_missingKm_phi);
    auto df_fall18_missingKm_all = ApplyFinalDVEPSelections(df_fall18_missingKm_phi);
    auto df_sp18_missingKm_all = ApplyFinalDVEPSelections(df_sp18_missingKm_phi);
    auto df_sp18_missingKp_all = ApplyFinalDVEPSelections(df_sp18_missingKp_phi);
    auto df_fall18_missingKp_all = ApplyFinalDVEPSelections(df_fall18_missingKp_phi);
    std::cout << "[Debug] Applied final DVEP selections to missing kaon phi dataframes" << std::endl;
    // Add to comparer
    // DumpExclusiveTxt(df_sp19_missingKm_phi_org, "missingKminus_events_dump_sp2019_inb.txt");
    // comparer.AddModelPhi(df_sp19_missingKm_phi_org, "Sp19 inb", beam_energy_sp2019, luminosity_rga_sp19_inb);
    // comparer.AddModelPhi(df_sp19_missingKm_all_org, "Sp19 inb (Missing K-)", beam_energy_sp2019, luminosity_rga_sp19_inb);
    comparer.AddModelPhi(df_sp19_missingKm_all, "Sp19 inb (Missing K-)", beam_energy_sp2019, luminosity_rga_sp19_inb);
    comparer.AddModelPhi(df_fall18_missingKm_all, "Fall18 inb (Missing K-)", beam_energy_fall2018, luminosity_rga_fall18_inb);
    comparer.AddModelPhi(df_sp18_missingKm_all, "Sp18 inb (Missing K-)", beam_energy_sp2018, luminosity_rga_sp18_inb);
    comparer.AddModelPhi(df_sp18_missingKp_all, "Sp18 outb (Missing K+)", beam_energy_sp2018, luminosity_rga_sp18_outb);
    comparer.AddModelPhi(df_fall18_missingKp_all, "Fall18 outb (Missing K+)", beam_energy_fall2018, luminosity_rga_fall18_outb);
  }
  // -----------------------------
  // Shared summary plots
  std::cout << "Applying further cuts and plotting…" << std::endl;

  // comparer.AddModelPhi(df_sp19inb_phi, "Sp19 inb", beam_energy_sp2019);
  // comparer.AddModelPhi(df_sp19_missingKm_all, "Sp19 inb (Missing K-)", beam_energy_sp2019);
  /*comparer.AddModelPhi(df_fall18outb_phi, "Fall18 outb", beam_energy_fall2018);
  comparer.AddModelPhi(df_fall18inb_phi,  "Fall18 inb",  beam_energy_fall2018);
  comparer.AddModelPhi(df_sp18outb_phi,   "Sp18 outb",   beam_energy_sp2018);
  comparer.AddModelPhi(df_sp18inb_phi,    "Sp18 inb",    beam_energy_sp2018);*/
  // Examples you had

  comparer.PlotPhiDVEPKinematicsPlots();
  comparer.PlotKinematicComparison_phiAna();
  comparer.PlotPhiAnaExclusivityComparisonByDetectorCases(detCuts);

  comparer.PlotPhiInvMassPerBin_AllModels("PhiInvMassFits", 40, 0.988, 1.15, true, 0.004, 0.25,  // sigmaRef, sigmaFrac
                                          branching,                                             // branching (0.492 etc)
                                          true,                                                  // doAcceptanceCorr
                                          true,                                                  // efficiencyCorr
                                          false                                                  // doRadCorr

  );
  comparer.PlotPhiDSigmaDt_FromCache(true, true, true, true, true);
  return;
  // ----------------------------------------------------------------
  //  R = sigma_L / sigma_T  from cos(theta_H) angular distribution
  //  Step 1: run per-(Q2,W,t',cos θ) mass fits and fill R / r04_00 cache
  //  Step 2: plot R vs t' and write summary CSV
  // ----------------------------------------------------------------
  comparer.PlotPhiRLTPerBin_AllModels("PhiRLTFits",  // output subdirectory
                                      40,            // nMassBins for K+K- invariant mass fit
                                      0.988,         // mMin  [GeV]
                                      1.15,          // mMax  [GeV]
                                      true,          // constrainSigma
                                      0.004,         // sigmaRef [GeV]
                                      0.25,          // sigmaFrac
                                      10.6           // beam energy [GeV]
  );
  comparer.PlotPhiRLT_FromCache("PhiRLT");

  return;
  xBins.SetQ2Bins({0.9, 8.35});
  comparer.SetXBinsRanges(xBins);
  xBins.SetTrentoPhiBins({0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360});
  comparer.SetXBinsRanges(xBins);
  // comparer.PlotPhiBSATrentoPhiPerBin_AllModels("PhiBSATrentoPhiFits", 40, 0.988, 1.15,true, 0.004, 0.25, polarisation);
  // comparer.PlotPhiBSATrentoPhi_FromCache();
  // comparer.PlotPhiALUCosThetaPerBin_AllModels("PhiALUCosThetaFits", 40, 0.988, 1.15, true, 0.004, 0.25, polarisation);
  // comparer.PlotPhiALUCosTheta_FromCache();

  // BSA as function of zPhi
  // xBins.SetZPhiBins({0.0, 0.9, 1.0});
  // comparer.SetXBinsRanges(xBins);

  // comparer.PlotPhiALUZPhiPerBin_AllModels("PhiALUZPhiFits",
  //                                      40, 0.988, 1.15,
  //                                      true, 0.004, 0.25, polarisation);
  // comparer.PlotPhiALUZPhi_FromCache();
  double mWindowLo = 1.01;  // your lower invMass_KpKm cut
  double mWindowHi = 1.03;  // your upper invMass_KpKm cut
  // Build A_LU^{sin(phi)}(cos(theta_KK)) using the moment method
  // comparer.PlotPhiALUCosThetaPerBin_AllModels_SinPhiMoment("PhiALUCosTheta_SinPhiMoment",    mWindowLo,    mWindowHi,    polarisation );
  // Then draw using the existing "FromCache" function
  // comparer.PlotPhiALUCosTheta_FromCache("PhiALUCosTheta_sin_phimoments");

  // Build A_LU(cosθ_KK) using the sin(phi)/(1 + b cos(phi)) fit
  comparer.PlotPhiALUCosThetaPerBin_AllModels_SinOver1PlusbCosFit("PhiALUCosTheta_SinOver1PlusbCosFit", mWindowLo, mWindowHi,
                                                                  polarisation  // beam polarization
  );

  // Then plot A_LU vs cosθ_KK using your existing cache-plotter
  comparer.PlotPhiALUCosTheta_FromCache("PhiALUCosTheta_with_SinOver1PlusbCosFit");
}

// -----------------------------------------------------------------------------
// Final DVEP selections (kept, just fixed labels to match values)
ROOT::RDF::RNode ApplyFinalDVEPSelections(ROOT::RDF::RNode df) {
  return df
      // 4. Q2 > 1
      .Filter("Q2 > 0.8 && Q2 < 6.0", "Cut: Q2 > 0.8 GeV^2")
      .Filter("W > 1.8", "Cut: W > 1.8 && W < 2.5 GeV")
      .Filter("recel_p  > 1.5", "Cut: recel_p < 1.5 GeV")
      //.Filter("bestEle_idx > 0", "Cut: reckPlus_p < 3.5 GeV")

      .Filter("reckPlus_p  <7.5", "Cut: reckPlus_p < 3.5 GeV")
      .Filter("reckMinus_p < 7.5", "Cut: reckMinus_p < 3.5 GeV")
      .Filter("ele_det_region == 1", "Cut: ele_det_region == 1")
      .Filter("Mx2_eKpKm > 0.8*0.8 && Mx2_eKpKm < 1.08*1.08", "Cut: Proton Missing Mass Squared in [0.8,1.08] GeV^2")
      .Filter("Mx2_epKm > .08 && Mx2_epKm < 0.48", "Cut: Kaon Missing Mass Squared in [0.08,.48] GeV^2")
      .Filter("Mx2_epKp > .08 && Mx2_epKp < 0.48", "Cut: Kaon Missing Mass Squared in [0.08,.48] GeV^2")
      .Filter("Cone_Kp < 6.0", "Cut: cone Angle  K+ < 6°")
      .Filter("Cone_Km < 6.0", "Cut: cone Angle  K- < 6°")
      .Filter("Cone_p < 6.0", "Cut:  cone between p < 6°")
      .Filter("Coplanarity_had_normals_deg <15", "Cut: Coplanarity angle 15°")
      .Filter("PTmiss < 0.120", "Cut: Total Missing PTmiss < .120 GeV")
      .Filter("Mx2_epKpKm < 0.0075", "Cut: Total Missing Mass squared < 0.0074 GeV")
      //.Filter("z_phi > 0.", "Cut: zPhi between 0.2 and 0.9")
      .Filter("Mx2_ep > 1.11-3*0.11&&Mx2_ep < 1.11+3*0.1", "Cut: Total ep Missing Mass squared < mean +3 sigma GeV")
      .Filter("Emiss <0.32 && Emiss> -0.175", "Cut: Missing energy Emiss < 0.2 GeV");
  //.Filter("invMass_pKminus < 1.312 || (invMass_pKminus > 1.5 && invMass_pKminus < 1.7) || invMass_pKminus > 1.9","Exclude invMass_pKminus in [1.312,1.5] and [1.7,1.9]");
  //.Filter("invMass_pKminus > 1.9","Exclude invMass_pKminus in [1.312,1.5] and [1.7,1.9]")
  //.Filter("invMass_KpKm > 0.9874 && invMass_KpKm < 1.12", "Cut: invMass_KpKm in [0.9874,1.12] GeV");
}

ROOT::RDF::RNode ApplyFinalGenDVEPSelections(ROOT::RDF::RNode df) {
  return df
      // 4. Q2 > 1
      .Filter("Q2 > 0.8", "Cut: Q2 > 0.8 GeV^2")
      .Filter("W > 1.8", "Cut: W > 1.8 && W < 2.5 GeV");
}

/////// mathematica functions
void PrintEdges(const std::string& label, const std::vector<double>& edges) {
  std::cout << "=== " << label << " bin edges (" << edges.size() - 1 << " bins) ===" << std::endl;
  std::cout << std::fixed << std::setprecision(4);
  for (std::size_t i = 0; i < edges.size(); ++i) {
    std::cout << "  edge[" << i << "] = " << edges[i] << std::endl;
  }
  std::cout << std::endl;
}

static void ExportBinningCSV(const std::vector<double>& q2Edges, const std::vector<double>& tpEdges, double beamMom, double W2lo, double W2hi, const std::string& outPath) {
  std::ofstream csv(outPath);
  if (!csv.is_open()) {
    std::cerr << "[ExportBinningCSV] ERROR: cannot open " << outPath << std::endl;
    return;
  }
  csv << std::fixed << std::setprecision(8);

  csv << "# Binning CSV for RunMDiffradNew — exported by DISANA_PhiAnalysisPlotter\n"
      << "# RunMDiffradNew computes tmin/|t| internally; only edges are needed here.\n"
      << "# Format: key,value0,value1,...\n";

  // Line 1: metadata
  csv << "beam_momentum," << beamMom << "\n";
  csv << "W2_range," << W2lo << "," << W2hi << "\n";

  // Line 2: Q² edges
  csv << "Q2_edges";
  for (auto e : q2Edges) csv << "," << e;
  csv << "\n";

  // Line 3: t' edges
  csv << "tprime_edges";
  for (auto e : tpEdges) csv << "," << e;
  csv << "\n";

  csv.close();

  // Terminal summary
  int nQ2 = (int)q2Edges.size() - 1;
  int nTp = (int)tpEdges.size() - 1;
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "[ExportBinningCSV] " << outPath << "  (" << nQ2 << " Q2 bins × " << nTp << " t' bins"
            << ", E_beam=" << beamMom << " GeV)\n"
            << "  Q2 edges:";
  for (auto e : q2Edges) std::cout << " " << e;
  std::cout << "\n  t' edges:";
  for (auto e : tpEdges) std::cout << " " << e;
  std::cout << "\n";
}

double computeLuminosity(double Q_mC) {
  const double e = 1.602e-19;                               // C
  const double NA = 6.022e23;                               // /mol
  const double rho = 0.07151;                               // g/cm³
  const double ell = 5.0;                                   // cm
  const double Aw = 1.00794;                                // g/mol
  const double target_areal_density = NA * rho * ell / Aw;  // ~2.136e22 cm⁻²
  // 1 cm⁻² = 1e33 nb⁻¹
  double n_electrons = (Q_mC * 1e-3) / e;
  return n_electrons * target_areal_density * 1e-33;  // nb⁻¹
}

/*double computeLuminosity(double Q_coulombs) {
  // Physical constants
  const double e = 1.602e-19;                    // electron charge (C)
  const double target_areal_density = 2.136e23;  // NA * rho * ell / (A * H) in cm^-2

  // Unit conversion: 1 nb = 1e-33 cm^2. Therefore, 1 cm^-2 = 1e33 nb^-1.
  // We must multiply the luminosity in cm^-2 by 1e33 to get nb^-1.
  const double cm2_to_nb_inv_factor = 1e33;

  // Integrated luminosity:
  // L = (Q/e) * (Target Areal Density) * (Conversion Factor)
  // Units: (unitless) * (cm^-2) * (nb^-1 / cm^-2) -> nb^-1
  double luminosity_nb_inv = (Q_coulombs / 1000 / e) * target_areal_density / cm2_to_nb_inv_factor;
  // double luminosity_nb_inv = Q_coulombs * 1.316875; // Simplified constant factor
  return luminosity_nb_inv;  // Unit is nb^-1
}*/