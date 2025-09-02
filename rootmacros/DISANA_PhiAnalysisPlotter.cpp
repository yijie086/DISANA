#include <THnSparse.h>

#include "../DreamAN/DrawHist/DISANAMMUtils.h"
#include "../DreamAN/DrawHist/DISANAMath.h"
#include "../DreamAN/DrawHist/DISANAMathFitUtils.h"
#include "../DreamAN/DrawHist/DISANA_PhiMassUtils.h"
#include "../DreamAN/DrawHist/DISANAcomparer.h"
#include "../DreamAN/DrawHist/DrawStyle.h"

// define analysis cuts here
ROOT::RDF::RNode ApplyFinalPhiSelections(ROOT::RDF::RNode df, bool inbending);
/// styling plots
// double double titleSize = 0.05, double labelSize = 0.04,double xTitleOffset = 1.1, double yTitleOffset = 1.6, int font = 42, int maxDigits = 5, int nDivisions = 510, double
// leftMargin = 0.16, double rightMargin = 0.07, double bottomMargin = 0.13, double topMargin = 0.06
DrawStyle KinStyle(0.07, 0.06, 0.9, 1.2);                                       // For Kin plots
DrawStyle dvcsStyle(0.06, 0.06, 1.2, 1.4, 42, 5, 510, 0.14, 0.07, 0.13, 0.06);  // For DVCS plots
DrawStyle csStyle(0.05, 0.05, .95, 1.1, 42, 5, 510, 0.12, 0.03, 0.12, 0.02);    // For Cross-Sections
DrawStyle bsaStyle(0.06, 0.045, .8, .8, 42, 5, 510, 0.15, 0.07, 0.16, 0.06);    // For BSA

// for exclusivity plots
std::vector<std::pair<std::string, std::string>> detCuts = {
    {"pro_det_region == 2", "CD"},
    {"pro_det_region == 1", "FD"},
};

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
  std::string input_path_from_analysisRun_SP19inb_data_missingKm = "./../data_processed/spring2019/inb/nsidis_wagon/missing_Km_output/";

  std::string filename_afterFid_SP18inb_data = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_SP18inb_data.c_str());
  std::string filename_afterFid_SP18outb_data = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_SP18outb_data.c_str());

  std::string filename_afterFid_Fall18inb_data = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_Fall18inb_data.c_str());
  std::string filename_afterFid_Fall18outb_data = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_Fall18outb_data.c_str());

  std::string filename_afterFid_SP19inb_data = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_SP19inb_data.c_str());
  std::string filename_afterFid_SP19inb_missingKm_data = Form("%s/dfSelected_afterFid_reprocessed.root", input_path_from_analysisRun_SP19inb_data_missingKm.c_str());

  // std::string filename_afterFid_7546_MC = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_7546_MC.c_str());
  float beam_energy_sp2018 = 10.5940;
  float beam_energy_fall2018 = 10.6000;
  float beam_energy_sp2019 = 10.1998;

  ROOT::RDF::RNode df_afterFid_sp18inb_data = InitKinematics(filename_afterFid_SP18inb_data, "dfSelected_afterFid", beam_energy_sp2018);
  ROOT::RDF::RNode df_afterFid_sp18outb_data = InitKinematics(filename_afterFid_SP18outb_data, "dfSelected_afterFid", beam_energy_sp2018);
  ROOT::RDF::RNode df_afterFid_fall18inb_data = InitKinematics(filename_afterFid_Fall18inb_data, "dfSelected_afterFid", beam_energy_fall2018);
  ROOT::RDF::RNode df_afterFid_fall18outb_data = InitKinematics(filename_afterFid_Fall18outb_data, "dfSelected_afterFid", beam_energy_fall2018);
  ROOT::RDF::RNode df_afterFid_sp19inb_data = InitKinematics(filename_afterFid_SP19inb_data, "dfSelected_afterFid", beam_energy_sp2019);
  ROOT::RDF::RNode df_afterFid_sp19inb_data_mKm = InitKinematics_MissingKm(filename_afterFid_SP19inb_data, "dfSelected_afterFid", beam_energy_sp2019);
  ROOT::RDF::RNode df_afterFid_sp19inb_data_missingKm_data =
      InitKinematics_MissingKm(filename_afterFid_SP19inb_missingKm_data, "dfSelected_afterFid_reprocessed", beam_energy_sp2019);



auto [df_cut, winKm] = DISANA::PhiMass::WireKaonMMCut_FromMx2(df_afterFid_sp19inb_data_missingKm_data, "Mx2_epKm","PhiMassPlots/spring2019/inb/nsidis" ,"Sp19 INB", "Km_miss", /*nBins*/220, /*xMin*/0.30, /*xMax*/0.70, /*nSigma*/3.0);      
auto [df_cut_exclusive, win_dvkpkm] = DISANA::PhiMass::WireKaonMMCut_FromMx2(df_afterFid_sp19inb_data_mKm,"Mx2_epKm", "PhiMassPlots/spring2019/inb/DVKpKm" ,"Sp19 INB DVKPKM", "Km_miss", /*nBins*/220, /*xMin*/0.30, /*xMax*/0.70, /*nSigma*/3.0);      
//auto [df_cut_exclusive, win_dvkpkm] = DISANA::PhiMass::WireKaonMMCut_FromMx2(df_afterFid_sp19inb_data_mKm,"Mx2_epKm", "PhiMassPlots/spring2019/inb/DVKpKm/exclusiveKpKm/" ,"Sp19 INB exlcusive", "Km_miss", /*nBins*/220, /*xMin*/0.30, /*xMax*/0.70, /*nSigma*/3.0);      


DISANA::PhiMass::DrawPhiMass_Kp_plus_KmMissing(df_cut,"PhiMassPlots/spring2019/inb/nsidis","Sp19 INB",200,0.8,1.6,0.9874,1.12,3.0);
DISANA::PhiMass::DrawPhiMass_Kp_plus_KmMissing(df_cut_exclusive, "PhiMassPlots/spring2019/inb/DVKpKm", "Sp19 INB DVKPKM",200,0.8,1.6,0.9874,1.12,3.0);
DISANA::PhiMass::DrawPhiMass_Measured(df_afterFid_sp19inb_data, "PhiMassPlots/spring2019/inb/PhiMassPlots/spring2019/inb/DVKpKm/exclusiveKpKm/", "Sp19 INB DVKPKM",200,0.8,1.6,0.9874,1.12,3.0);

std::cout << "Applying further cuts and plottings" << std::endl;

  DrawStyle fitStyle(0.06, 0.05, 1.0, 1.3);  
  
  // You can tweak this
  // DVCS event selection cuts

  // Apply final DVCS cuts
  auto df_afterFid_sp18inb_with_all_data = ApplyFinalPhiSelections(df_afterFid_sp18inb_data, true);
  auto df_afterFid_sp18outb_with_all_data = ApplyFinalPhiSelections(df_afterFid_sp18outb_data, true);
  auto df_afterFid_sp18inb_with_phi_data = SelectExclusivePhiEvent(df_afterFid_sp18inb_with_all_data);
  auto df_afterFid_sp18outb_with_phi_data = SelectExclusivePhiEvent(df_afterFid_sp18outb_with_all_data);

  auto df_afterFid_fall18inb_with_all_data = ApplyFinalPhiSelections(df_afterFid_fall18inb_data, true);
  auto df_afterFid_fall18outb_with_all_data = ApplyFinalPhiSelections(df_afterFid_fall18outb_data, true);
  auto df_afterFid_fall18inb_with_phi_data = SelectExclusivePhiEvent(df_afterFid_fall18inb_with_all_data);
  auto df_afterFid_fall18outb_with_phi_data = SelectExclusivePhiEvent(df_afterFid_fall18outb_with_all_data);

  auto df_afterFid_sp19inb_with_all_data = ApplyFinalPhiSelections(df_afterFid_sp19inb_data, true);
  auto df_afterFid_sp19inb_with_all_data_km = ApplyFinalPhiSelections(df_cut_exclusive, true);
  auto df_afterFid_sp19inb_with_phi_data = SelectExclusivePhiEvent(df_afterFid_sp19inb_with_all_data);
  auto df_afterFid_sp19inb_with_phi_data_Km = SelectExclusivePhiEvent(df_afterFid_sp19inb_with_all_data_km);
  auto df_sp19_missingKm_all = ApplyFinalPhiSelections(df_cut, true);
  auto df_sp19_missingKm_phi = SelectPhiEvent_MissingKm(df_sp19_missingKm_all);


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
  // xBins.SetQ2Bins({0.4, 8.35});
  xBins.SetXBBins({0, 0.99});
  xBins.SetTBins({0.2, .5, 0.8, 1.2, 8.0});  // 5.4, 5.8, 6.2, 6.6, 7.0, 7.4, 7.8, 8.2, 8.6, 9.0, 9.4, 9.8});

  comparer.SetXBinsRanges(xBins);
  comparer.UseFittedPhiYields(true);

  //comparer.AddModelPhi(df_afterFid_fall18outb_with_phi_data, "Fall18 outb", beam_energy_fall2018);
  //comparer.AddModelPhi(df_afterFid_fall18inb_with_phi_data, "Fall18 inb", beam_energy_fall2018);

  //comparer.AddModelPhi(df_afterFid_sp18outb_with_phi_data, "Sp18 outb", beam_energy_sp2018);
  //comparer.AddModelPhi(df_afterFid_sp18inb_with_phi_data, "Sp18 inb", beam_energy_sp2018);

  comparer.AddModelPhi(df_afterFid_sp19inb_with_phi_data, "Sp19 inb", beam_energy_sp2019);
  comparer.AddModelPhi(df_afterFid_sp19inb_with_phi_data_Km, "Sp19 inb DVKPKM Wa(Missing K-)", beam_energy_sp2019);
  comparer.AddModelPhi(df_sp19_missingKm_phi, "Sp19 inb (Missing K-)", beam_energy_sp2019);

  //double luminosity = 24.3065 * pow(10, 6);  // Set your desired luminosity here nb^-1
  double luminosity_rga_fall18 = 1.0;        // 5.47*pow(10,40);  // Set your desired luminosity here nb^-1
  double polarisation = 0.85;                // Set your desired polarisation here
  double branching = 0.49;                   // Set your desired polarisation here

  comparer.PlotPhiElectroProKinematicsComparison();
  comparer.PlotKinematicComparison_phiAna();
  //comparer.PlotPhiAnaExclusivityComparisonByDetectorCases(detCuts);
  // comparer.PlotPhiInvMassPerBin_AllModels("PhiInvMassFits", /*nBins*/30, /*mMin*/0.9974, /*mMax*/1.2, /*constrainSigma*/true, luminosity_rga_fall18,branching);
  // comparer.PlotPhiDSigmaDt_FromCache();

  gApplication->Terminate(0);
}

// exclusivity cuts
ROOT::RDF::RNode ApplyFinalPhiSelections(ROOT::RDF::RNode df, bool inbending) {
  return df
      // 4. Q2 > 1
      .Filter("Q2 > 1.0", "Cut: Q2 > 1 GeV^2")
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
