#include "../DreamAN/DrawHist/DISANAcomparer.h"
#include "../DreamAN/DrawHist/DrawStyle.h"

template <typename T>
void DrawAllExclusivityVarsOnCanvas(T&& df, const DrawStyle& style, const std::string& outname_prefix) {
  std::vector<std::tuple<std::string, std::string, std::string, double, double>> vars = {
      {"Mx2_ep", "Missing Mass Squared (ep)", "MM^{2}(ep) [GeV^{2}]", -2.0, 2.0},
      {"Emiss", "Missing Energy", "E_{miss} [GeV]", -0.5, 3.0},
      {"PTmiss", "Transverse Missing Momentum", "P_{T}^{miss} [GeV/c]", 0, 1.0},
      {"Theta_gamma_gamma", "#theta(#gamma, #vec{q})", "#theta_{#gamma#gamma'} [deg]", 0, 5},
      {"DeltaPhi", "Coplanarity Angle", "#Delta#phi [deg]", 0, 50},
      {"Mx2_epg", "Missing Mass Squared (ep#gamma)", "MM^{2}(ep#gamma) [GeV^{2}]", -1.0, 1.0},
      {"Mx2_eg", "Invariant Mass (e#gamma)", "M^{2}(e#gamma) [GeV^{2}]", 0.0, 2.5},
      {"Theta_e_gamma", "Angle: e-#gamma", "#theta(e, #gamma) [deg]", 0.0, 180.0},
      {"DeltaE", "Energy Balance", "#DeltaE [GeV]", -2.0, 2.0},
  };

  const int cols = 3;
  const int rows = (vars.size() + cols - 1) / cols;
  TCanvas* c = new TCanvas((outname_prefix + "_canvas").c_str(), "DVCS Exclusivity Variables", 1600, 1000);
  c->Divide(cols, rows);

  std::vector<std::unique_ptr<TH1D>> histos;

  for (size_t i = 0; i < vars.size(); ++i) {
    const auto& [var, title, xlabel, xmin, xmax] = vars[i];

    auto histR = df.Histo1D({("h_" + var).c_str(), (title + ";" + xlabel + ";Counts").c_str(), 100, xmin, xmax}, var);
    histR.GetValue();  // force evaluation

    auto* h_original = (TH1D*)histR.GetPtr();
    if (!h_original) {
      std::cerr << "Warning: null histogram pointer for variable " << var << std::endl;
      histos.emplace_back(nullptr);
      continue;
    }

    // Clone to take ownership and avoid RResultPtr lifetimes
    auto* h_clone = (TH1D*)h_original->Clone(("clone_" + var).c_str());
    h_clone->SetDirectory(0);  // detach from any file
    histos.emplace_back(h_clone);
  }

  for (size_t i = 0; i < histos.size(); ++i) {
    c->cd(i + 1);
    auto h = histos[i].get();
    const auto& [var, title, xlabel, xmin, xmax] = vars[i];

    if (!h || !gPad) {
      std::cerr << "Skipping invalid histogram or gPad for " << var << std::endl;
      continue;
    }

    gPad->SetTicks();
    style.StylePad((TPad*)gPad);
    style.StyleTH1(h, title.c_str());

    h->SetLineWidth(2);
    h->SetLineColor(kRed + 1);
    h->Draw("HISTE");

    double mu = h->GetMean();
    double sigma = h->GetStdDev();
    double ymax = h->GetMaximum();

    TLine* l1 = new TLine(mu - 3 * sigma, 0, mu - 3 * sigma, ymax * 0.6);
    TLine* l2 = new TLine(mu + 3 * sigma, 0, mu + 3 * sigma, ymax * 0.6);
    l1->SetLineColor(kMagenta + 2);
    l2->SetLineColor(kMagenta + 2);
    l1->SetLineStyle(2);
    l2->SetLineStyle(2);
    l1->Draw("SAME");
    l2->Draw("SAME");

    TLatex text;
    text.SetNDC();
    text.SetTextFont(42);
    text.SetTextSize(0.045);
    text.DrawLatex(0.6, 0.83, Form("Mean = %.3f", mu));
    text.DrawLatex(0.6, 0.75, Form("#sigma = %.3f", sigma));

    std::cout << var << ": Entries = " << h->GetEntries() << ", Mean = " << mu << ", RMS = " << h->GetRMS() << std::endl;
  }

  c->SaveAs((outname_prefix + "_ExclusivitySummary.png").c_str());
  delete c;
}

void FitAndDrawExclusivityVar(ROOT::RDF::RNode& df, const std::string& var, const std::string& title, const std::string& xlabel, double xmin, double xmax,
                              const std::string& outname, const DrawStyle& style) {
  auto hist = df.Histo1D({("h_" + var).c_str(), (title + ";" + xlabel + ";Counts").c_str(), 100, xmin, xmax}, var);

  TCanvas* c = new TCanvas(("c_" + var).c_str(), title.c_str(), 1200, 1000);
  style.ApplyGlobalStyle();
  style.StylePad(c);
  style.StyleTH1(hist.GetPtr(), title.c_str());
  style.NormalizeHistogram(hist.GetPtr());
  hist->Draw("h");

  // Get mean and standard deviation
  double mu = hist->GetMean();
  double sigma = hist->GetStdDev();

  // ±3σ lines
  TLine* l1 = new TLine(mu - 3 * sigma, 0, mu - 3 * sigma, hist->GetMaximum() * 0.6);
  TLine* l2 = new TLine(mu + 3 * sigma, 0, mu + 3 * sigma, hist->GetMaximum() * 0.6);
  l1->SetLineColor(kMagenta + 2);
  l2->SetLineColor(kMagenta + 2);
  l1->SetLineStyle(2);
  l2->SetLineStyle(2);
  if ((mu - 3 * sigma) > 0) l1->Draw("SAME");
  l2->Draw("SAME");

  // Draw text
  TLatex text;
  text.SetNDC();
  text.SetTextSize(0.045);
  text.SetTextFont(42);

  text.DrawLatex(0.7, 0.85, Form("Mean = %.3f", mu));
  text.DrawLatex(0.7, 0.80, Form("#sigma = %.3f", sigma));

  // Save
  c->SaveAs((outname + ".png").c_str());
  // c->SaveAs((outname + ".pdf").c_str());
  delete c;
}

ROOT::RDF::RNode SelectOneTriplet(ROOT::RDF::RNode df_) {
  return df_.Filter(
      [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<bool>& pass) {
        int e = 0, g = 0, p = 0;
        for (size_t i = 0; i < pid.size(); ++i) {
          if (!pass[i]) continue;
          if (pid[i] == 11)
            e++;
          else if (pid[i] == 22)
            g++;
          else if (pid[i] == 2212)
            p++;
        }
        return (e == 1 && g == 1 && p == 1);
      },
      {"REC_Particle_pid", "REC_Particle_pass"});
}

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

ROOT::RDF::RNode InitKinematics(const std::string& filename_ = "", const std::string& treename_ = "", float beam_energy = 0);
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
                     [method, beam_energy](double recel_p, double recel_theta, double recel_phi, double recpro_p, double recpro_theta, double recpro_phi, double recpho_p,
                                           double recpho_theta, double recpho_phi) {
                       return (DISANAMath(beam_energy, recel_p, recel_theta, recel_phi, recpro_p, recpro_theta, recpro_phi, recpho_p, recpho_theta, recpho_phi).*method)();
                     },
                     {"recel_p", "recel_theta", "recel_phi", "recpro_p", "recpro_theta", "recpro_phi", "recpho_p", "recpho_theta", "recpho_phi"});
}

void DISANA_Xplotter() {
  ROOT::EnableImplicitMT();
  // std::string input_path_from_analysisRun = "/work/clas12/singh/CrossSectionAN/RGA_spring2018_Analysis/fromDVCS_wagon/Inb/";
  // test case
  std::string input_path_from_analysisRun_inb = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/RGA_spring2018_Analysis/fromDVCS_wagon/Inb/";
  // std::string input_path_from_analysisRun_inb = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/RGA_sims/test/";
  std::string input_path_from_analysisRun_out = "/w/hallb-scshelf2102/clas12/singh/CrossSectionAN/RGA_spring2018_Analysis/fromDVCS_wagon/Outb/";
  // std::string input_path_from_analysisRun = "./../build";
  std::string filename_afterFid_inb = Form("%s/dfSelected_afterFid_reprocessed.root", input_path_from_analysisRun_inb.c_str());
  std::string filename_afterFid_outb = Form("%s/dfSelected_afterFid.root", input_path_from_analysisRun_out.c_str());
  // std::string filename_afterFid_afterCorr = Form("%s/dfSelected_afterFid_afterCorr.root", input_path_from_analysisRun.c_str());
  float beam_energy = 10.6;

  //// this is where you can plot the comparisions
  ROOT::RDF::RNode df_afterFid_inb = InitKinematics(filename_afterFid_inb, "dfSelected_afterFid_reprocessed", beam_energy);
  ROOT::RDF::RNode df_afterFid_all_inb = InitKinematics(filename_afterFid_inb, "dfSelected_afterFid_reprocessed", beam_energy);

  ROOT::RDF::RNode df_afterFid_outb = InitKinematics(filename_afterFid_outb, "dfSelected_afterFid", beam_energy);
  ROOT::RDF::RNode df_afterFid_all_outb = InitKinematics(filename_afterFid_outb, "dfSelected_afterFid", beam_energy);

  // ROOT::RDF::RNode df_afterFid_afterCorr = InitKinematics(filename_afterFid_afterCorr, "dfSelected_afterFid_afterCorr");
  // input files for the data

  gSystem->Exec("mkdir -p ExclusivityFits");

  DrawStyle fitStyle(0.06, 0.05, 1.0, 1.3);  // You can tweak this
  // DVCS event selection cuts
  bool inbending = true;

  // Apply final DVCS cuts
  auto df_final_all_all_inb = ApplyFinalDVCSSelections(df_afterFid_all_inb, inbending);
  auto df_final_all_inb = SelectOneTriplet(df_final_all_all_inb.Filter(
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
        return (e == 1 && g == 1 && p == 1);
      },
      {"REC_Particle_pid", "REC_Particle_pass", "REC_DaughterParticle_pass"}, "Cut: one good e, γ (not π⁰-like), p"));

  auto df_final_all_all_outb = ApplyFinalDVCSSelections(df_afterFid_all_outb, false);
  auto df_final_all_outb = SelectOneTriplet(df_final_all_all_outb.Filter(
      [](const ROOT::VecOps::RVec<int>& pid, const ROOT::VecOps::RVec<bool>& pass, const ROOT::VecOps::RVec<bool>& daughterPass) {
        int e = 0, g = 0, p = 0;
        // float pho_pMax= 0.0;
        for (size_t i = 0; i < pid.size(); ++i) {
          if (!pass[i]) continue;
          if (pid[i] == 11)
            e++;
          else if (pid[i] == 22 && !daughterPass[i])
            // if(pho_p >= pho_pMax) pho_pMax = pho_p;
            g++;  // photon must NOT be from pi0
          else if (pid[i] == 2212)
            p++;
        }
        return (e <= 20 && g <= 10 && p <= 10);
      },
      {"REC_Particle_pid", "REC_Particle_pass", "REC_DaughterParticle_pass"}, "Cut: one good e, γ (not π⁰-like), p"));

  auto df_FD_CD = df_final_all_inb.Filter("pho_det_region == 1 && pro_det_region == 2", "Photon in FD, Proton in CD");
  auto df_FT_CD = df_final_all_inb.Filter("pho_det_region == 0 && pro_det_region == 2", "Photon in FT, Proton in CD");
  auto df_FD_FD = df_final_all_inb.Filter("pho_det_region == 1 && pro_det_region == 1", "Photon in FD, Proton in FD");

  DrawAllExclusivityVarsOnCanvas(df_FT_CD, fitStyle, "ExclusivityFits/Inbending_FT_CD");
  DrawAllExclusivityVarsOnCanvas(df_FD_CD, fitStyle, "ExclusivityFits/Inbending_FD_CD");
  DrawAllExclusivityVarsOnCanvas(df_FD_FD, fitStyle, "ExclusivityFits/Inbending_FD_FD");

  //// styling plots
  // double double titleSize = 0.05, double labelSize = 0.04,double xTitleOffset = 1.1, double yTitleOffset = 1.6, int font = 42, int maxDigits = 5, int nDivisions = 510, double
  // leftMargin = 0.16, double rightMargin = 0.07, double bottomMargin = 0.13, double topMargin = 0.06
  DrawStyle KinStyle(0.07, 0.06, 0.9, 1.2);                                       // For Kin plots
  DrawStyle dvcsStyle(0.06, 0.06, 1.2, 1.4, 42, 5, 510, 0.14, 0.07, 0.13, 0.06);  // For DVCS plots
  DrawStyle csStyle(0.05, 0.04, 1.0, 1.3);                                        // For Cross-Sections
  DrawStyle bsaStyle(0.06, 0.045, 1.0, 1.2);                                      // For BSA

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
  xBins.SetQ2Bins({.1, 3.0, 8.0});
  xBins.SetTBins({0.0, 0.5, 1.2});
  xBins.SetXBBins({0.0, 0.5, 1});
  comparer.SetXBinsRanges(xBins);

  // comparer.AddModel(df_afterFid_afterCorr, "after Correction", beam_energy);
  // comparer.AddModel(df_afterFid, "Before Exclusivity cuts", beam_energy);
  comparer.AddModel(df_final_all_inb, "Sp18 Inb", beam_energy);
  comparer.AddModel(df_final_all_outb, "Sp18 OutB", beam_energy);
  comparer.PlotKinematicComparison();
  comparer.PlotDVCSKinematicsComparison();
  comparer.PlotDISCrossSectionComparison(1);  // argument is Luminosity
  //comparer.PlotDIS_BSA_Comparison(1);         // argument is Luminosity
  comparer.PlotDIS_BSA_Comparison_optimized(1);  // argument is Luminosity
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