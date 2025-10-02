#include <THnSparse.h>

#include "../DreamAN/DrawHist/DISANAcomparer.h"
#include "../DreamAN/DrawHist/DrawStyle.h"
#include "../DreamAN/DrawHist/DISANAMath.h"

ROOT::RDF::RNode RejectPi0TwoPhoton(ROOT::RDF::RNode df_);
ROOT::RDF::RNode SelectPi0Event(ROOT::RDF::RNode df);
void CreateCorrectionHistogram4D(ROOT::RDF::RNode df_dvcs_mc, ROOT::RDF::RNode df_pi0_mc, ROOT::RDF::RNode df_dvcs_data, ROOT::RDF::RNode df_pi0_data,
                                 const std::string& out_file_name);

ROOT::RDF::RNode ApplyFinalDVCSSelections(ROOT::RDF::RNode df, bool inbending);

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

void DISANA_Xplotter_RGA() {
  bool ComputeBgk_core = false;  // Set to true if you want to compute background
  bool DoBkgCorr = true;       // Set to true if you want to apply background correction
  float beam_energy = 10.594;
  ROOT::EnableImplicitMT();


  /// proton momentum correction applied
  
  //data path
  std::string input_path_from_analysisRun_inb_data = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/spring2018/inb/DVCS_wagon/";
  std::string input_path_from_analysisRun_out_data = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/spring2018/outb/DVCS_wagon/";
  
  //pi-0 MC path
  std::string input_path_from_analysisRun_inb_pi0MC  ="/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/sims/aaogen/inb/";
  std::string input_path_from_analysisRun_outb_pi0MC  ="/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/sims/aaogen/inb/";
  
  //DVCS MC path
  std::string input_path_from_analysisRun_inb_DVCSMC_rec  ="/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/sims/DVCSgen/inb/rec/";
  std::string input_path_from_analysisRun_outb_DVCSMC_rec  ="/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/sims/DVCSgen/outb/rec/";
  std::string input_path_from_analysisRun_inb_DVCSMC_gen  ="/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/sims/DVCSgen/inb/accept_all/";
  std::string input_path_from_analysisRun_outb_DVCSMC_gen  ="/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/sims/DVCSgen/outb/accept_all/";

  std::string input_path_from_analysisRun_inb_dvcsmc_bkg = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/sims/DVCSgen/outb/45na_bkg/reprocessed_from/";
  std::string input_path_from_analysisRun_outb_dvcsmc_bkg = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/sims/DVCSgen/outb/45na_bkg/";
  std::string input_path_from_analysisRun_inb_dvcsmc_nobkg = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/sims/DVCSgen/inb/no_bgk/";
  std::string input_path_from_analysisRun_outb_dvcsmc_nobkg = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/sims/DVCSgen/outb/no_bkg/";

  std::string input_path_from_analysisRun_dvcsmc_rad = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/sims/aaogen/rad_gen/";
  std::string input_path_from_analysisRun_dvcsmc_norad = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/sims/aaogen/no_rad_gen/";
  
  // File names
  std::string filename_afterFid_inb_data_corr = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_inb_data.c_str());
  std::string filename_afterFid_outb_data_corr = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_out_data.c_str());
  
  std::string filename_afterFid_inb_pi0MC = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_inb_pi0MC.c_str());
  std::string filename_afterFid_outb_pi0MC = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_outb_pi0MC.c_str());
  
  std::string filename_afterFid_dvcsmc_inb_gen = Form("%s/dfSelected.root", input_path_from_analysisRun_inb_DVCSMC_gen.c_str());
  std::string filename_afterFid_dvcsmc_outb_gen = Form("%s/dfSelected.root", input_path_from_analysisRun_outb_DVCSMC_gen.c_str());
  std::string filename_afterFid_dvcsmc_inb_rec = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_inb_DVCSMC_rec.c_str());
  std::string filename_afterFid_dvcsmc_outb_rec = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_outb_DVCSMC_rec.c_str());

  std::string filename_afterFid_dvcsmc_inb_bkg = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_inb_dvcsmc_bkg.c_str());
  std::string filename_afterFid_dvcsmc_inb_nobkg = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_inb_dvcsmc_nobkg.c_str());

  std::string filename_afterFid_dvcsmc_outb_bkg = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_outb_dvcsmc_bkg.c_str());
  std::string filename_afterFid_dvcsmc_outb_nobkg = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun_outb_dvcsmc_nobkg.c_str());


  std::string filename_afterFid_dvcsmc_rad = Form("%s/rad100M.root", input_path_from_analysisRun_dvcsmc_rad.c_str());
  std::string filename_afterFid_dvcsmc_norad = Form("%s/no_rad100M.root", input_path_from_analysisRun_dvcsmc_norad.c_str());

  

  /// Initialize dataframes
  ROOT::RDF::RNode df_afterFid_inb_data_corr = InitKinematics(filename_afterFid_inb_data_corr, "dfSelected_afterFid_afterCorr", beam_energy);
  ROOT::RDF::RNode df_afterFid_outb_data_corr = InitKinematics(filename_afterFid_outb_data_corr, "dfSelected_afterFid_afterCorr", beam_energy);
  
  ROOT::RDF::RNode df_afterFid_inb_pi0MC = InitKinematics(filename_afterFid_inb_pi0MC, "dfSelected_afterFid_afterCorr", beam_energy);
  ROOT::RDF::RNode df_afterFid_outb_pi0MC = InitKinematics(filename_afterFid_outb_pi0MC, "dfSelected_afterFid_afterCorr", beam_energy);

  ROOT::RDF::RNode df_afterFid_dvcsmc_inb_gen = InitGenKinematics(filename_afterFid_dvcsmc_inb_gen, "dfSelected", beam_energy);
  ROOT::RDF::RNode df_afterFid_dvcsmc_outb_gen = InitGenKinematics(filename_afterFid_dvcsmc_outb_gen, "dfSelected", beam_energy);
  ROOT::RDF::RNode df_afterFid_dvcsmc_inb_rec = InitKinematics(filename_afterFid_dvcsmc_inb_rec, "dfSelected_afterFid_afterCorr", beam_energy);
  ROOT::RDF::RNode df_afterFid_dvcsmc_outb_rec = InitKinematics(filename_afterFid_dvcsmc_outb_rec, "dfSelected_afterFid_afterCorr", beam_energy);

  ROOT::RDF::RNode df_afterFid_dvcsmc_inb_bkg = InitKinematics(filename_afterFid_dvcsmc_inb_bkg, "dfSelected_afterFid_afterCorr", beam_energy);
  ROOT::RDF::RNode df_afterFid_dvcsmc_outb_bkg = InitKinematics(filename_afterFid_dvcsmc_outb_bkg, "dfSelected_afterFid_afterCorr", beam_energy);
  ROOT::RDF::RNode df_afterFid_dvcsmc_inb_nobkg = InitKinematics(filename_afterFid_dvcsmc_inb_nobkg, "dfSelected_afterFid_afterCorr", beam_energy);
  ROOT::RDF::RNode df_afterFid_dvcsmc_outb_nobkg = InitKinematics(filename_afterFid_dvcsmc_outb_nobkg, "dfSelected_afterFid_afterCorr", beam_energy);
 

  ROOT::RDF::RNode df_afterFid_dvcsmc_rad = InitGenKinematics(filename_afterFid_dvcsmc_rad, "MC", beam_energy);
  ROOT::RDF::RNode df_afterFid_dvcsmc_norad = InitGenKinematics(filename_afterFid_dvcsmc_norad, "MC", beam_energy);
 // input files for the data
  gSystem->Exec("mkdir -p ExclusivityFits");

  DrawStyle fitStyle(0.06, 0.05, 1.0, 1.3);  // You can tweak this
  // DVCS event selection cuts

  // Apply final DVCS cuts---------------  
  auto df_final_dvcs_inb_data_corr = ApplyFinalDVCSSelections(df_afterFid_inb_data_corr, true);
  auto df_final_dvcs_outb_data_corr =  ApplyFinalDVCSSelections(df_afterFid_outb_data_corr, false);
  
  auto df_final_dvcs_inb_pi0MC = ApplyFinalDVCSSelections(df_afterFid_inb_pi0MC, true);
  auto df_final_dvcs_outb_pi0MC = ApplyFinalDVCSSelections(df_afterFid_outb_pi0MC, false);

  auto df_final_dvcs_inb_DVCSMC_rec = ApplyFinalDVCSSelections(df_afterFid_dvcsmc_inb_rec, true);
  auto df_final_dvcs_outb_DVCSMC_rec = ApplyFinalDVCSSelections(df_afterFid_dvcsmc_outb_rec, false);

  auto df_final_dvcs_inb_DVCSMC_bkg = ApplyFinalDVCSSelections(df_afterFid_dvcsmc_inb_bkg, true);
  auto df_final_dvcs_outb_DVCSMC_bkg = ApplyFinalDVCSSelections(df_afterFid_dvcsmc_outb_bkg, false);

  auto df_final_dvcs_inb_DVCSMC_nobkg = ApplyFinalDVCSSelections(df_afterFid_dvcsmc_inb_nobkg, true);
  auto df_final_dvcs_outb_DVCSMC_nobkg = ApplyFinalDVCSSelections(df_afterFid_dvcsmc_outb_nobkg, false);

  // reject pi0 events from dvcs sample ----------- 
  auto df_final_dvcsPi_rejected_inb_data_corr = RejectPi0TwoPhoton(df_final_dvcs_inb_data_corr);
  auto df_final_dvcsPi_rejected_outb_data_corr = RejectPi0TwoPhoton(df_final_dvcs_outb_data_corr);
  
  auto df_final_dvcsPi_rejected_inb_pi0MC = RejectPi0TwoPhoton(df_final_dvcs_inb_pi0MC);
  auto df_final_dvcsPi_rejected_outb_pi0MC = RejectPi0TwoPhoton(df_final_dvcs_outb_pi0MC);

  auto df_final_dvcsPi_rejected_inb_DVCSMC_rec = RejectPi0TwoPhoton(df_final_dvcs_inb_DVCSMC_rec);
  auto df_final_dvcsPi_rejected_outb_DVCSMC_rec = RejectPi0TwoPhoton(df_final_dvcs_outb_DVCSMC_rec);
 
  auto df_final_dvcsPi_rejected_inb_DVCSMC_bkg = RejectPi0TwoPhoton(df_final_dvcs_inb_DVCSMC_bkg);
  auto df_final_dvcsPi_rejected_outb_DVCSMC_bkg = RejectPi0TwoPhoton(df_final_dvcs_outb_DVCSMC_bkg);

  auto df_final_dvcsPi_rejected_inb_DVCSMC_nobkg = RejectPi0TwoPhoton(df_final_dvcs_inb_DVCSMC_nobkg);
  auto df_final_dvcsPi_rejected_outb_DVCSMC_nobkg = RejectPi0TwoPhoton(df_final_dvcs_outb_DVCSMC_nobkg);

  
  // pi0 event selection cuts ------------
  auto df_final_OnlPi0_inb_data_corr = SelectPi0Event(df_final_dvcs_inb_data_corr);
  auto df_final_OnlPi0_outb_data_corr = SelectPi0Event(df_final_dvcs_outb_data_corr);

  auto df_final_OnlPi0_inb_pi0MC = SelectPi0Event(df_final_dvcs_inb_pi0MC);
  auto df_final_OnlPi0_outb_pi0MC = SelectPi0Event(df_final_dvcs_outb_pi0MC);

  // final single photon from pi0 correction factors here

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
  xBins.SetQ2Bins({1.0, 1.2, 1.456, 1.912, 2.51, 3.295, 4.326, 5.761});
  xBins.SetTBins({0.11, 1.0});
  xBins.SetXBBins({0.062, 0.09, 0.118, 0.155, 0.204, 0.268, 0.357, 0.446, 0.581});
  comparer.SetXBinsRanges(xBins);

  comparer.AddModelwithPi0Corr(df_final_dvcsPi_rejected_outb_data_corr,df_final_OnlPi0_outb_data_corr,df_final_dvcsPi_rejected_outb_pi0MC,df_final_OnlPi0_outb_pi0MC,/*mcdvcs*/df_afterFid_dvcsmc_outb_gen,/*mcdvcs_acceptance*/df_final_dvcsPi_rejected_outb_DVCSMC_rec, /*mc dvcsbkg*/df_final_dvcsPi_rejected_outb_DVCSMC_bkg, /*dvcs mcnobkg*/df_final_dvcsPi_rejected_outb_DVCSMC_nobkg, /*dvscrad*/df_afterFid_dvcsmc_rad, /*dvcsmc_norad*/df_afterFid_dvcsmc_norad, "Sp18 outb after", beam_energy, true, true,true, true);
  comparer.AddModelwithPi0Corr(df_final_dvcsPi_rejected_inb_data_corr,df_final_OnlPi0_inb_data_corr,df_final_dvcsPi_rejected_inb_pi0MC,df_final_OnlPi0_inb_pi0MC,/*mcdvcs*/df_afterFid_dvcsmc_inb_gen,/*mcdvcs_acceptance*/df_final_dvcsPi_rejected_inb_DVCSMC_rec, /*mc dvcsbkg*/df_final_dvcsPi_rejected_inb_DVCSMC_bkg, /*dvcs mcnobkg*/df_final_dvcsPi_rejected_inb_DVCSMC_nobkg, /*dvscrad*/df_afterFid_dvcsmc_rad, /*dvcsmc_norad*/df_afterFid_dvcsmc_norad, "Sp18 Inb after", beam_energy, true, true,true,true);
  //comparer.AddModelwithPi0Corr(df_final_dvcsPi_rejected_inb_data,df_final_OnlPi0_inb_data,df_final_dvcsPi_rejected_inb_MC,df_final_OnlPi0_inb_MC,/*mcdvcs*/df_final_dvcsPi_rejected_inb_MC,/*mcdvcs_acceptance*/df_final_OnlPi0_inb_MC, /*mc dvcsbkg*/df_final_OnlPi0_inb_MC, /*dvcs mcnobkg*/df_final_OnlPi0_inb_MC, /*dvscrad*/df_final_OnlPi0_inb_MC, /*dvcsmc_norad*/df_final_OnlPi0_inb_MC, "Sp18 Inb no corr", beam_energy, false, false,false,false);
  //comparer.AddModelwithPi0Corr(df_final_dvcsPi_rejected_inb_data_corr,df_final_OnlPi0_inb_data_corr,df_final_dvcsPi_rejected_inb_pi0MC,df_final_OnlPi0_inb_pi0MC,/*mcdvcs*/df_afterFid_dvcsmc_inb_gen,/*mcdvcs_acceptance*/df_final_dvcsPi_rejected_inb_DVCSMC_rec, /*mc dvcsbkg*/df_final_OnlPi0_inb_pi0MC, /*dvcs mcnobkg*/df_final_OnlPi0_inb_pi0MC, /*dvscrad*/df_final_OnlPi0_inb_pi0MC, /*dvcsmc_norad*/df_final_OnlPi0_inb_pi0MC, "Sp18 Inb Before", beam_energy, false, false,false,false);
  //comparer.AddModelwithPi0Corr(df_final_dvcsPi_rejected_inb_data_corr,df_final_OnlPi0_inb_data_corr,df_final_dvcsPi_rejected_inb_pi0MC,df_final_OnlPi0_inb_pi0MC,/*mcdvcs*/df_afterFid_dvcsmc_inb_gen,/*mcdvcs_acceptance*/df_final_dvcsPi_rejected_inb_DVCSMC_rec, /*mc dvcsbkg*/df_final_OnlPi0_inb_pi0MC, /*dvcs mcnobkg*/df_final_OnlPi0_inb_pi0MC, /*dvscrad*/df_final_OnlPi0_inb_pi0MC, /*dvcsmc_norad*/df_final_OnlPi0_inb_pi0MC, "Sp18 Inb after", beam_energy, true, true,false,false);

  //comparer.AddModelwithPi0Corr(df_final_dvcsPi_rejected_outb_data_corr,df_final_OnlPi0_outb_data_corr,df_final_dvcsPi_rejected_outb_pi0MC,df_final_OnlPi0_outb_pi0MC,/*mcdvcs*/df_afterFid_dvcsmc_outb_gen,/*mcdvcs_acceptance*/df_final_dvcsPi_rejected_outb_DVCSMC_rec, /*mc dvcsbkg*/df_final_OnlPi0_outb_pi0MC, /*dvcs mcnobkg*/df_final_OnlPi0_outb_pi0MC, /*dvscrad*/df_afterFid_dvcsmc_rad, /*dvcsmc_norad*/df_afterFid_dvcsmc_norad, "Sp18 outb Before", beam_energy, false, false,false,false);
  
  
  //comparer.AddModelwithPi0Corr(df_final_dvcsPi_rejected_outb_data_corr,df_final_OnlPi0_outb_data,df_final_dvcsPi_rejected_outb_MC,df_final_OnlPi0_outb_MC, "Sp18 outb corr", beam_energy, false);


  double luminosity = 1.0;  // Set your desired luminosity here
  double polarisation = 0.85;  // Set your desired polarisation here

  //comparer.PlotDISCrossSectionComparison(luminosity);  // argument is Luminosity, polarisation
  //comparer.PlotDIS_BSA_Comparison(luminosity, polarisation);         // argument is Luminosity
  comparer.PlotDIS_BSA_Cross_Section_AndCorr_Comparison(luminosity, polarisation, true, true, true, true, true, true);
  //comparer.PlotExclusivityComparisonByDetectorCases(detCuts);
  //comparer.PlotKinematicComparison();
  //comparer.PlotDVCSKinematicsComparison();
  gApplication->Terminate(0);
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
ROOT::RDF::RNode ApplyFinalDVCSSelections(ROOT::RDF::RNode df, bool inbending) {
   auto df1= df
      // 4. Q2 > 1
      .Filter("Q2 > 1.0", "Cut: Q2 > 1 GeV^2")
      .Filter("t < 1.0", "Cut: t > 1 GeV^2")
      //.Filter("recel_p > 6.0", "Cut: recel_p > 0.6")
      .Filter("recpro_p < 1.2", "Cut: recpro_p < 1.2")
      // 5. W > 2
      .Filter("W > 2.0", "Cut: W > 2.0 GeV")
     // 9. 3σ exclusivity cuts intitial loose cuts
  .Filter("Emiss < 1.0", "Cut: Missing energy")
  .Filter("PTmiss < 0.2", "Cut: Transverse missing momentum")
  .Filter("Theta_e_gamma > 5 ", "Cut: Theta_e_gamma")
  .Filter("Theta_gamma_gamma < 2.0", "Cut: photon-missing angle")
  .Filter("(pho_det_region==0&&pro_det_region==2)||(pho_det_region==1&&pro_det_region==1)||(pho_det_region==1&&pro_det_region==2)", "Cut: three config");

  if (inbending){
  return df1 = df1
  .Filter("(pho_det_region==0&&pro_det_region==2&&Mx2_ep<(0.10+3*0.20)&&Mx2_ep>(0.10-3*0.20))||(pho_det_region==1&&pro_det_region==1&&Mx2_ep<(0.06+3*0.18)&&Mx2_ep>(0.06-3*0.18))||(pho_det_region==1&&pro_det_region==2&&Mx2_ep<(0.08+3*0.18)&&Mx2_ep>(0.08-3*0.18))", "Cut: Mx2_ep in 3sigma")
  .Filter("(pho_det_region==0&&pro_det_region==2&&Emiss<(.20+3*0.24)&&Emiss>(.20-3*0.24))||(pho_det_region==1&&pro_det_region==1&&Emiss<(.27+3*0.27)&&Emiss>.27-3*0.27)||(pho_det_region==1&&pro_det_region==2&&Emiss<(0.29+3*0.34)&&Emiss>(0.29-3*0.34))", "Cut: Emiss in 3sigma")
  .Filter("(pho_det_region==0&&pro_det_region==2&&PTmiss<(.05+3*0.04)&&PTmiss>(.05-3*0.04))||(pho_det_region==1&&pro_det_region==1&&PTmiss<(.09+3*0.05)&&PTmiss>(.09-3*0.05))||(pho_det_region==1&&pro_det_region==2&&PTmiss<(.07+3*0.04)&&PTmiss>(.07-3*0.04))", "Cut: PTmiss in 3sigma")
  .Filter("(pho_det_region==0&&pro_det_region==2&&Theta_gamma_gamma<(0.41+3*0.33)&&Theta_gamma_gamma>(0.41-3*0.33))||(pho_det_region==1&&pro_det_region==1&&Theta_gamma_gamma<(.81+3*0.50)&&Theta_gamma_gamma>(.81-3*0.50))||(pho_det_region==1&&pro_det_region==2&&Theta_gamma_gamma<(.63+3*0.47)&&Theta_gamma_gamma>(.63-3*0.47))", "Cut: Theta_gamma_gamma in 3sigma")
  .Filter("(pho_det_region==0&&pro_det_region==2&&DeltaPhi<(2.60+3*3.01)&&DeltaPhi>(2.60-3*3.01))||(pho_det_region==1&&pro_det_region==1&&DeltaPhi<(5.55+3*4.70)&&DeltaPhi>(5.55-3*4.70))||(pho_det_region==1&&pro_det_region==2&&DeltaPhi<(3.92+3*3.90)&&DeltaPhi>(3.92-3*3.90))", "Cut: DeltaPhi in 3sigma")
  .Filter("(pho_det_region==0&&pro_det_region==2&&Mx2_epg<(0.0+3*0.01)&&Mx2_epg>(0.0-3*0.01))||(pho_det_region==1&&pro_det_region==1&&Mx2_epg<(0.0+3*0.01)&&Mx2_epg>(0.0-3*0.01))||(pho_det_region==1&&pro_det_region==2&&Mx2_epg<(0.0+3*0.01)&&Mx2_epg>(0.0-3*0.01))", "Cut: Mx2_epg in 3sigma")
  .Filter("(pho_det_region==0&&pro_det_region==2&&Mx2_eg<(1.17+3*0.35)&&Mx2_eg>(1.17-3*0.35))||(pho_det_region==1&&pro_det_region==1&&Mx2_eg<(1.16+3*0.31)&&Mx2_eg>(1.16-3*0.31))||(pho_det_region==1&&pro_det_region==2&&Mx2_eg<(1.27+3*0.49)&&Mx2_eg>(1.27-3*0.49))", "Cut: Mx2_eg in 3sigma")
  .Filter("(pho_det_region==0&&pro_det_region==2&&Theta_e_gamma<(19.11+3*3.80)&&Theta_e_gamma>(19.11-3*3.80))||(pho_det_region==1&&pro_det_region==1&&Theta_e_gamma<(33.82+3*3.34)&&Theta_e_gamma>(33.82-3*3.34))||(pho_det_region==1&&pro_det_region==2&&Theta_e_gamma<(23.82+3*3.50)&&Theta_e_gamma>(23.82-3*3.50))", "Cut: Theta_e_gamma in 3sigma");
}else{
   return df1 = df1
    .Filter("(pho_det_region==0&&pro_det_region==2&&Mx2_ep<(0.11+3*0.19)&&Mx2_ep>(0.11-3*0.19))||(pho_det_region==1&&pro_det_region==1&&Mx2_ep<(0.07+3*0.17)&&Mx2_ep>(0.07-3*0.17))||(pho_det_region==1&&pro_det_region==2&&Mx2_ep<(0.07+3*0.16)&&Mx2_ep>(0.07-3*0.16))", "Cut: Mx2_ep in 3sigma")
    .Filter("(pho_det_region==0&&pro_det_region==2&&Emiss<(.20+3*0.24)&&Emiss>(.20-3*0.24))||(pho_det_region==1&&pro_det_region==1&&Emiss<(.27+3*0.27)&&Emiss>.27-3*0.27)||(pho_det_region==1&&pro_det_region==2&&Emiss<(0.29+3*0.34)&&Emiss>(0.29-3*0.34))", "Cut: Emiss in 3sigma")
    .Filter("(pho_det_region==0&&pro_det_region==2&&PTmiss<(.05+3*0.04)&&PTmiss>(.05-3*0.04))||(pho_det_region==1&&pro_det_region==1&&PTmiss<(.09+3*0.05)&&PTmiss>(.09-3*0.05))||(pho_det_region==1&&pro_det_region==2&&PTmiss<(.07+3*0.04)&&PTmiss>(.07-3*0.04))", "Cut: PTmiss in 3sigma")
    .Filter("(pho_det_region==0&&pro_det_region==2&&Theta_gamma_gamma<(0.49+3*0.40)&&Theta_gamma_gamma>(0.5-3*0.40))||(pho_det_region==1&&pro_det_region==1&&Theta_gamma_gamma<(.91+3*0.51)&&Theta_gamma_gamma>(.91-3*0.51))||(pho_det_region==1&&pro_det_region==2&&Theta_gamma_gamma<(.71+3*0.50)&&Theta_gamma_gamma>(.71-3*0.50))", "Cut: Theta_gamma_gamma in 3sigma")
    .Filter("(pho_det_region==0&&pro_det_region==2&&DeltaPhi<(2.41+3*2.92)&&DeltaPhi>(2.41-3*2.92))||(pho_det_region==1&&pro_det_region==1&&DeltaPhi<(4.35+3*4.02)&&DeltaPhi>(4.35-3*4.02))||(pho_det_region==1&&pro_det_region==2&&DeltaPhi<(3.71+3*3.67)&&DeltaPhi>(3.71-3*3.67))", "Cut: DeltaPhi in 3sigma")
    .Filter("(pho_det_region==0&&pro_det_region==2&&Mx2_epg<(0.0+3*0.01)&&Mx2_epg>(0.0-3*0.01))||(pho_det_region==1&&pro_det_region==1&&Mx2_epg<(0.0+3*0.01)&&Mx2_epg>(0.0-3*0.01))||(pho_det_region==1&&pro_det_region==2&&Mx2_epg<(0.0+3*0.01)&&Mx2_epg>(0.0-3*0.01))", "Cut: Mx2_epg in 3sigma")
    .Filter("(pho_det_region==0&&pro_det_region==2&&Mx2_eg<(1.19+3*0.37)&&Mx2_eg>(1.19-3*0.37))||(pho_det_region==1&&pro_det_region==1&&Mx2_eg<(1.18+3*0.30)&&Mx2_eg>(1.18-3*0.30))||(pho_det_region==1&&pro_det_region==2&&Mx2_eg<(1.29+3*0.49)&&Mx2_eg>(1.29-3*0.49))", "Cut: Mx2_eg in 3sigma")
    .Filter("(pho_det_region==0&&pro_det_region==2&&Theta_e_gamma<(15.72+3*4.69)&&Theta_e_gamma>(15.72-3*4.69))||(pho_det_region==1&&pro_det_region==1&&Theta_e_gamma<(33.56+3*2.85)&&Theta_e_gamma>(33.56-3*2.85))||(pho_det_region==1&&pro_det_region==2&&Theta_e_gamma<(21.36+3*4.22)&&Theta_e_gamma>(21.36-3*4.22))", "Cut: Theta_e_gamma in 3sigma");
  }
  // 10. Quality Assurance Cut
  //.Filter("REC_Event_pass == true", "Cut: QA pass");
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
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 11 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})
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
              .Define("recpro_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 2212 && pass[i]) return px[i];
                       return -999.0f;
                     },
                     {"REC_Particle_pid", "REC_Particle_vz", "REC_Particle_pass"})
              .Define("recpho_vz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px, const ROOT::VecOps::RVec<bool>& pass) -> float {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 && pass[i]) return px[i];
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
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& px) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 /*&& maxEpass[i]*/) return px[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_px"})
             .Define("pho_py",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& py) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22  /*&& maxEpass[i]*/) return py[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_py"})
             .Define("pho_pz",
                     [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<float>& pz) {
                       for (size_t i = 0; i < pid.size(); ++i)
                         if (pid[i] == 22 /*&& maxEpass[i]*/) return pz[i];
                       return -999.0f;
                     },
                     {"MC_Particle_pid", "MC_Particle_pz"})
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


