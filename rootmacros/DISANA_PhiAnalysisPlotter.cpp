#include <THnSparse.h>
#include <ROOT/RDataFrame.hxx>
#include <TROOT.h>
#include <TApplication.h>
#include <TString.h>

#include "../DreamAN/DrawHist/DISANAMMUtils.h"
#include "../DreamAN/DrawHist/DISANAMath.h"
#include "../DreamAN/DrawHist/DISANAMathFitUtils.h"
#include "../DreamAN/DrawHist/DISANA_PhiMassUtils.h"
#include "../DreamAN/DrawHist/DISANAcomparer.h"
#include "../DreamAN/DrawHist/DrawStyle.h"

// -----------------------------------------------------------------------------
// Forward declare your final selections (implemented below)
ROOT::RDF::RNode ApplyFinalDVEPSelections(ROOT::RDF::RNode df);

// -----------------------------------------------------------------------------
// Plot styling (unchanged)
DrawStyle KinStyle(0.07, 0.06, 0.9, 1.2);                                       // For Kin plots
DrawStyle dvcsStyle(0.06, 0.06, 1.2, 1.4, 42, 5, 510, 0.14, 0.07, 0.13, 0.06);  // For DVCS plots
DrawStyle csStyle(0.05, 0.05, .95, 1.1, 42, 5, 510, 0.12, 0.03, 0.12, 0.02);    // For Cross-Sections
DrawStyle bsaStyle(0.06, 0.045, .8, .8, 42, 5, 510, 0.15, 0.07, 0.16, 0.06);    // For BSA

// For exclusivity plots
std::vector<std::pair<std::string, std::string>> detCuts = {
    {"pro_det_region == 2", "CD"},
    {"pro_det_region == 1", "FD"},
};

// -----------------------------------------------------------------------------
// Main plotter with toggles
void DISANA_PhiAnalysisPlotter()   // subset toggle inside missing-mass
{
  bool runExclusive   = false;
  bool runMissingMass = true;
  ROOT::EnableImplicitMT();

  // -----------------------------
  // Input locations
  // Exclusive reconstruction K+K-
  std::string input_path_from_analysisRun_SP18inb_data   = "./../data_processed/spring2018/inb/DVKpKm_wagon/";
  std::string input_path_from_analysisRun_SP18outb_data  = "./../data_processed/spring2018/outb/DVKpKm_wagon/";
  std::string input_path_from_analysisRun_Fall18inb_data = "./../data_processed/fall2018/inb/DVKpKm_wagon/";
  std::string input_path_from_analysisRun_Fall18outb_data= "./../data_processed/fall2018/outb/DVKpKm_wagon/";
  std::string input_path_from_analysisRun_SP19inb_data   = "./../data_processed/spring2019/inb/DVKpKm_wagon/";

  // Reprocessed missing–mass
  std::string input_path_from_analysisRun_SP19inb_data_missingKm   = "./../data_processed/spring2019/inb/nsidis_wagon/missing_Km_output/";
  std::string input_path_from_analysisRun_Fall18inb_data_missingKm  = "./../data_processed/fall2018/inb/nsidis_wagon/missing_Km_output/";
  std::string input_path_from_analysisRun_SP18inb_data_missingKm    = "./../data_processed/spring2018/inb/nsidis_wagon/missing_Km_output/";
  std::string input_path_from_analysisRun_SP18outb_data_missingKp   = "./../data_processed/spring2018/outb/nsidis_wagon/missing_Kp_output/";
  std::string input_path_from_analysisRun_Fall18outb_data_missingKp = "./../data_processed/fall2018/outb/nsidis_wagon/missing_Kp_output/";

  // File names
  std::string filename_afterFid_SP18inb_data   = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_SP18inb_data.c_str());
  std::string filename_afterFid_SP18outb_data  = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_SP18outb_data.c_str());
  std::string filename_afterFid_Fall18inb_data = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_Fall18inb_data.c_str());
  std::string filename_afterFid_Fall18outb_data= Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_Fall18outb_data.c_str());
  std::string filename_afterFid_SP19inb_data   = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_SP19inb_data.c_str());

  std::string filename_afterFid_SP19inb_missingKm_data   = Form("%s/dfSelected_afterFid_reprocessed.root", input_path_from_analysisRun_SP19inb_data_missingKm.c_str());
  std::string filename_afterFid_Fall18inb_missingKm_data = Form("%s/dfSelected_afterFid_reprocessed.root", input_path_from_analysisRun_Fall18inb_data_missingKm.c_str());
  std::string filename_afterFid_SP18inb_missingKm_data   = Form("%s/dfSelected_afterFid_reprocessed.root", input_path_from_analysisRun_SP18inb_data_missingKm.c_str());
  std::string filename_afterFid_SP18outb_missingKp_data  = Form("%s/dfSelected_afterFid_reprocessed.root", input_path_from_analysisRun_SP18outb_data_missingKp.c_str());
  std::string filename_afterFid_Fall18outb_missingKp_data= Form("%s/dfSelected_afterFid_reprocessed.root", input_path_from_analysisRun_Fall18outb_data_missingKp.c_str());

  // Beam energies
  float beam_energy_sp2018   = 10.5940f;
  float beam_energy_fall2018 = 10.6000f;
  float beam_energy_sp2019   = 10.1998f;

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
  xBins.SetQ2Bins({0.4, 1.3, 1.86, 2.71, 8.35});
  xBins.SetXBBins({0, 0.99});
  xBins.SetTBins({0.2, .5, 0.8, 1.2, 8.0});
  comparer.SetXBinsRanges(xBins);
  comparer.UseFittedPhiYields(true);

  // Some global constants you had
  double luminosity_rga_fall18 = 1.0;
  double polarisation          = 0.85;
  double branching             = 0.49;

  // -----------------------------
  // Initialize RDataFrames up front (shared by both modes)
  ROOT::RDF::RNode df_afterFid_sp18inb_data   = InitKinematics(filename_afterFid_SP18inb_data,   "dfSelected_afterFid", beam_energy_sp2018);
  ROOT::RDF::RNode df_afterFid_sp18outb_data  = InitKinematics(filename_afterFid_SP18outb_data,  "dfSelected_afterFid", beam_energy_sp2018);
  ROOT::RDF::RNode df_afterFid_fall18inb_data = InitKinematics(filename_afterFid_Fall18inb_data, "dfSelected_afterFid", beam_energy_fall2018);
  ROOT::RDF::RNode df_afterFid_fall18outb_data= InitKinematics(filename_afterFid_Fall18outb_data,"dfSelected_afterFid", beam_energy_fall2018);
  ROOT::RDF::RNode df_afterFid_sp19inb_data   = InitKinematics(filename_afterFid_SP19inb_data,   "dfSelected_afterFid", beam_energy_sp2019);

  ROOT::RDF::RNode df_afterFid_sp18inb_missingKm_data   = InitKinematics_MissingKm(filename_afterFid_SP18inb_missingKm_data,   "dfSelected_afterFid_reprocessed", beam_energy_sp2018);
  ROOT::RDF::RNode df_afterFid_fall18inb_missingKm_data = InitKinematics_MissingKm(filename_afterFid_Fall18inb_missingKm_data, "dfSelected_afterFid_reprocessed", beam_energy_fall2018);
  ROOT::RDF::RNode df_afterFid_sp19inb_missingKm_data   = InitKinematics_MissingKm(filename_afterFid_SP19inb_missingKm_data,   "dfSelected_afterFid_reprocessed", beam_energy_sp2019);
  //ROOT::RDF::RNode df_afterFid_sp18outb_missingKp_data  = InitKinematics_MissingKp(filename_afterFid_SP18outb_missingKp_data,  "dfSelected_afterFid_reprocessed", beam_energy_sp2018);
  ROOT::RDF::RNode df_afterFid_fall18outb_missingKp_data= InitKinematics_MissingKp(filename_afterFid_Fall18outb_missingKp_data,"dfSelected_afterFid_reprocessed", beam_energy_fall2018);

  // -----------------------------
  // 1) Exclusive measured K+K- analysis (toggle)
  if (runExclusive)
  {
    // Apply your “final” DVEP selections and then pick exclusive phi events
    auto df_sp18inb_all    = SelectExclusivePhiEvent(df_afterFid_sp18inb_data);
    auto df_sp18outb_all   = SelectExclusivePhiEvent(df_afterFid_sp18outb_data);
    auto df_fall18inb_all  = SelectExclusivePhiEvent(df_afterFid_fall18inb_data);
    auto df_fall18outb_all = SelectExclusivePhiEvent(df_afterFid_fall18outb_data);
    auto df_sp19inb_all    = SelectExclusivePhiEvent(df_afterFid_sp19inb_data);

    auto df_sp18inb_phi    = ApplyFinalDVEPSelections(df_sp18inb_all);
    auto df_sp18outb_phi   = ApplyFinalDVEPSelections(df_sp18outb_all);
    auto df_fall18inb_phi  = ApplyFinalDVEPSelections(df_fall18inb_all);
    auto df_fall18outb_phi = ApplyFinalDVEPSelections(df_fall18outb_all);
    auto df_sp19inb_phi    = ApplyFinalDVEPSelections(df_sp19inb_all);

    DISANA::PhiMass::DrawPhiMass_Measured(df_sp19inb_phi,"PhiMassPlots/spring2019/inb/DVKpKm/exclusiveKpKm/","Sp19 INB DVKpKm",/*nBins*/200, /*mMin*/0.8, /*mMax*/1.6, /*mPhiLo*/0.9874, /*mPhiHi*/1.12, /*nSigma*/3.0);
    DISANA::PhiMass::DrawPhiMass_Measured(df_fall18inb_phi,"PhiMassPlots/fall2018/inb/DVKpKm/exclusiveKpKm/","Fall18 INB DVKpKm",/*nBins*/200, /*mMin*/0.8, /*mMax*/1.6, /*mPhiLo*/0.9874, /*mPhiHi*/1.12, /*nSigma*/3.0);  
    DISANA::PhiMass::DrawPhiMass_Measured(df_sp18inb_phi,"PhiMassPlots/spring2018/inb/DVKpKm/exclusiveKpKm/","Sp18 INB DVKpKm",/*nBins*/200, /*mMin*/0.8, /*mMax*/1.6, /*mPhiLo*/0.9874, /*mPhiHi*/1.12, /*nSigma*/3.0);
    DISANA::PhiMass::DrawPhiMass_Measured(df_sp18outb_phi,"PhiMassPlots/spring2018/outb/DVKpKm/exclusiveKpKm/","Sp18 OUTB DVKpKm",/*nBins*/200, /*mMin*/0.8, /*mMax*/1.6, /*mPhiLo*/0.9874, /*mPhiHi*/1.12, /*nSigma*/3.0);
    DISANA::PhiMass::DrawPhiMass_Measured(df_fall18outb_phi,"PhiMassPlots/fall2018/outb/DVKpKm/exclusiveKpKm/","Fall18 OUTB DVKpKm",/*nBins*/200,/*mMin*/0.8, /*mMax*/1.6, /*mPhiLo*/0.9874, /*mPhiHi*/1.12, /*nSigma*/3.0);
    
    // Add any datasets you want to compare
    comparer.AddModelPhi(df_sp19inb_phi, "Sp19 inb", beam_energy_sp2019);
    comparer.AddModelPhi(df_fall18outb_phi, "Fall18 outb", beam_energy_fall2018);
    comparer.AddModelPhi(df_fall18inb_phi,  "Fall18 inb",  beam_energy_fall2018);
    comparer.AddModelPhi(df_sp18outb_phi,   "Sp18 outb",   beam_energy_sp2018);
    comparer.AddModelPhi(df_sp18inb_phi,    "Sp18 inb",    beam_energy_sp2018);
  }
  if (runMissingMass)
  {
      auto [df_cut_missing_Km_sp19,  winKm2019]  = DISANA::PhiMass::WireKaonMMCut_FromMx2(df_afterFid_sp19inb_missingKm_data,   "Mx2_epKp", "PhiMassPlots/spring2019/inb/nsidis",  "Sp19 INB",  "K^{+}", 220, 0.25, 0.75, 3.0);
      auto [df_cut_missing_Km_fall18,winKmFall18]= DISANA::PhiMass::WireKaonMMCut_FromMx2(df_afterFid_fall18inb_missingKm_data, "Mx2_epKp", "PhiMassPlots/fall2018/inb/nsidis",   "Fall18 INB","K^{+}", 220, 0.25, 0.75, 3.0);
      auto [df_cut_missing_Km_sp18,  winKmSp18]  = DISANA::PhiMass::WireKaonMMCut_FromMx2(df_afterFid_sp18inb_missingKm_data,   "Mx2_epKp", "PhiMassPlots/spring2018/inb/nsidis", "Sp18 INB",  "K^{+}", 220, 0.25, 0.75, 3.0);
      //auto [df_cut_missing_Kp_sp18,  winKpSp18]   = DISANA::PhiMass::WireKaonMMCut_FromMx2(df_afterFid_sp18outb_missingKp_data,  "Mx2_epKm", "PhiMassPlots/spring2018/outb/nsidis", "Sp18 OUTB", "K^{-}", 220, 0.30, 0.80, 3.0);
      auto [df_cut_missing_Kp_fall18,winKpFall18] = DISANA::PhiMass::WireKaonMMCut_FromMx2(df_afterFid_fall18outb_missingKp_data,"Mx2_epKm", "PhiMassPlots/fall2018/outb/nsidis",   "Fall18 OUTB","K^{-}", 220, 0.25, 0.75, 3.0);

      // Draw phi mass with Km missing
      DISANA::PhiMass::DrawPhiMass_Kp_plus_KmMissing(df_cut_missing_Km_sp19,  "PhiMassPlots/spring2019/inb/nsidis",  "Sp19 INB",  200, 0.8, 1.6, 0.9874, 1.12, 3.0);
      //DISANA::PhiMass::DrawPhiMass_Kp_plus_KmMissing(df_cut_missing_Km_fall18,"PhiMassPlots/fall2018/inb/nsidis",    "Fall18 INB",200, 0.8, 1.6, 0.9874, 1.12, 3.0);
      //DISANA::PhiMass::DrawPhiMass_Kp_plus_KmMissing(df_cut_missing_Km_sp18,  "PhiMassPlots/spring2018/inb/nsidis",  "Sp18 INB",  200, 0.8, 1.6, 0.9874, 1.12, 3.0);
      //DISANA::PhiMass::DrawPhiMass_Km_plus_KpMissing(df_cut_missing_Kp_sp18,  "PhiMassPlots/spring2018/outb/nsidis", "Sp18 OUTB", 200, 0.8, 1.6, 0.9874, 1.12, 3.0);
      //DISANA::PhiMass::DrawPhiMass_KpMissing_plus_Km(df_cut_missing_Kp_fall18,"PhiMassPlots/fall2018/outb/nsidis",   "Fall18 OUTB",200, 0.8, 1.6, 0.9874, 1.12, 3.0);

      // Pass through your final DVEP selections & dedicated phi picker for missing-Km streams
      auto df_sp19_missingKm_phi  = SelectPhiEvent_MissingKm(df_cut_missing_Km_sp19);
      auto df_fall18_missingKm_phi= SelectPhiEvent_MissingKm(df_cut_missing_Km_fall18);
      auto df_sp18_missingKm_phi  = SelectPhiEvent_MissingKm(df_cut_missing_Km_sp18);
      //auto df_sp18_missingKp_phi  = SelectPhiEvent_MissingKp(df_cut_missing_Kp_sp18);
      auto df_fall18_missingKp_phi= SelectPhiEvent_MissingKp(df_cut_missing_Kp_fall18);
      
      auto df_sp19_missingKm_all  = ApplyFinalDVEPSelections(df_sp19_missingKm_phi);
      auto df_fall18_missingKm_all= ApplyFinalDVEPSelections(df_fall18_missingKm_phi);
      auto df_sp18_missingKm_all  = ApplyFinalDVEPSelections(df_sp18_missingKm_phi);
      //auto df_sp18_missingKp_all  = ApplyFinalDVEPSelections(df_sp18_missingKp_phi);
      auto df_fall18_missingKp_all= ApplyFinalDVEPSelections(df_fall18_missingKp_phi);  
      // Add to comparer
      comparer.AddModelPhi(df_sp19_missingKm_all, "Sp19 inb (Missing K-)", beam_energy_sp2019);
      comparer.AddModelPhi(df_fall18_missingKm_all, "Fall18 inb (Missing K-)", beam_energy_fall2018); 
      comparer.AddModelPhi(df_sp18_missingKm_all, "Sp18 inb (Missing K-)", beam_energy_sp2018);
      //comparer.AddModelPhi(df_sp18_missingKp_all, "Sp18 outb (Missing K+)", beam_energy_sp2018);
      comparer.AddModelPhi(df_fall18_missingKp_all, "Fall18 outb (Missing K+)", beam_energy_fall2018);
  }

  // -----------------------------
  // Shared summary plots
  std::cout << "Applying further cuts and plotting…" << std::endl;

  // Examples you had:
  comparer.PlotPhiElectroProKinematicsComparison();
  comparer.PlotKinematicComparison_phiAna();
  comparer.PlotPhiAnaExclusivityComparisonByDetectorCases(detCuts);
  comparer.PlotPhiInvMassPerBin_AllModels("PhiInvMassFits", 30, 0.9974, 1.2, true, luminosity_rga_fall18, branching);
  comparer.PlotPhiDSigmaDt_FromCache();

  gApplication->Terminate(0);
}

// -----------------------------------------------------------------------------
// Final DVEP selections (kept, just fixed labels to match values)
ROOT::RDF::RNode ApplyFinalDVEPSelections(ROOT::RDF::RNode df) {
  return df
      // 4. Q2 > 1
      .Filter("Q2 > 1.0", "Cut: Q2 > 1 GeV^2")
      // Momentum caps (keep your 3.5 value; corrected title text)
      .Filter("reckPlus_p  < 3.5", "Cut: reckPlus_p < 3.5 GeV")
      .Filter("reckMinus_p < 3.5", "Cut: reckMinus_p < 3.5 GeV")
      // 5. W > 2 (your title said 1.8; using 2.0 as in expression)
      .Filter("W > 2.0", "Cut: W > 2.0 GeV")
      // 9. Missing energy / exclusivity
      .Filter("Emiss < 0.200", "Cut: Missing energy Emiss < 0.2 GeV");
}
