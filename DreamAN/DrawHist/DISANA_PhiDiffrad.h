#ifndef DISANA_PHI_DIFFRAD_H
#define DISANA_PHI_DIFFRAD_H

// Utilities to compute and cache Phi radiative-correction ratios using diffrad.
//
// ===== Sign conventions =====
//
//  Column   Sign        Definition
//  ------   ----        ----------
//  t        POSITIVE    std::abs((proton_in - proton_out).Mag2())   → +|t|
//  tmin     NEGATIVE    tmin_phi_from_Q2_xB(Q2, xB)                → physical tmin < 0
//  mtprime  POSITIVE    |t| - |tmin|  = t_col + tmin_col  ≥ 0
//
// ===== DIFFRAD inmdi.dat format (matches RunMDiffrad_backup.C exactly) =====
//
//  9 scalar header lines (one value per line, inline ! comment OK)
//  6 data rows (all npoi values on ONE line, NO inline comments):
//    row 1: -W2_hi  per point  (more negative)
//    row 2: -W2_lo  per point  (less negative)
//    row 3: -Q2_hi  per point
//    row 4: -Q2_lo  per point   ← must never be 0.0 (causes SIGFPE in Fortran)
//    row 5:  t_lo   per point  (physical t, more negative = larger |t|)
//    row 6:  t_hi   per point  (physical t, less negative = smaller |t|)
//
// ===== t ↔ t' conversion =====
//
//   t_phys  = tmin_used - mtprime   (both negative)
//   mtprime = tmin_used - t_phys    (inverse)
//
//   For mtprime bin [mtp_lo, mtp_hi]:
//     t_phys_hi = tmin_used - mtp_lo  (less negative, smaller |t|)
//     t_phys_lo = tmin_used - mtp_hi  (more negative, larger |t|)
//
// ===== Cache tree columns (radcorr) =====
//
//   Q2          GeV^2    bin-centre Q^2
//   W           GeV      bin-centre W from data Mean
//   t_centre    GeV^2    physical t bin centre from DIFFRAD (NEGATIVE)
//   tmin_used   GeV^2    tmin from data Mean("tmin") for this (Q2,W) bin (NEGATIVE)
//   mtprime     GeV^2    t' = tmin_used - t_phys  (POSITIVE, primary quantity)
//   rad_corr             radiative correction factor
//   rad_corr_err         sigma over nev samples
//   iq, iw, it           bin indices

#include <ROOT/RDataFrame.hxx>
#include "DISANAMath.h"   // for m_p, tmin_phi_from_Q2_xB

#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

class BinManager;  // must provide GetQ2Bins(), GetWBins(), GetTprimeBins()

namespace DISANA {
namespace PhiDiffrad {

namespace fs = std::filesystem;

// ============================================================
// SentBin — everything needed to map a DIFFRAD output entry
// back to its analysis bin and reconstruct t' consistently.
// ============================================================
struct SentBin {
  size_t iq, iw, it;
  double tmin_used;     // Mean("tmin") from data for this (Q2,W) bin (NEGATIVE)
  double mtp_lo, mtp_hi;  // t' bin edges (POSITIVE)
  double Q2c;           // Q^2 bin centre
  double Wc;            // W bin centre (from data, used for display only)
};

// ============================================================
// DiffradRunParams
// ============================================================
struct DiffradRunParams {
  double bmom       = 10.60;   // beam lepton momentum [GeV]
  double tmom       = 0.0;
  int    lepton     = 1;       // 1=electron, 2=muon
  int    ivec       = 3;       // 1=rho,2=omega,3=phi,4=J/Psi
  double ann1       = 1e6;
  double ann2       = 1e6;
  double ann3       = 5e5;
  double vcut       = 0.0;
  int    nev        = 3;
  long   seed       = 333522;
  // W window half-width sent to DIFFRAD around the data W-mean [GeV].
  // Backup uses ~0.5 GeV (W range 2.5–3.5, ΔW=1). Sending the full
  // analysis bin (e.g. 2–10 GeV) makes DIFFRAD average over a huge
  // unphysical phase space.
  double wHalfWidth = 0.5;
  // Minimum Q2 lower edge sent to DIFFRAD [GeV^2].
  // DIFFRAD Fortran crashes (SIGFPE at mdiffrad.f:178) if Q2_lo = 0.0
  // because it computes log(Q2) or 1/Q2 internally.
  double q2LoMin    = 0.5;
  // Half-width of the t window sent to DIFFRAD around each t-prime bin centre [GeV^2].
  //
  // CRITICAL: DIFFRAD must receive a NARROW t window centred on the bin midpoint,
  // NOT the full bin edge-to-edge width. The radiative integral in sig_rad_
  // (mdiffrad.f:575) diverges when the t range is wide.
  //
  // Sending full bin widths (up to 0.53 GeV^2) produced diverging RC values:
  //   RC(it=1)=1.33, RC(it=2)=1.77, RC(it=3)=2.77, RC(it=4)=5.38,
  //   RC(it=5)=12.5, RC(it=6)=47.4, RC(it=7)=336 -> SIGFPE in sig_rad_
  //
  // The backup uses t in [-0.25,-0.15] -> half-width = 0.05 GeV^2.
  // For each t-prime bin: t_centre = tmin_used - mtp_mid
  //   t_lo = t_centre - tHalfWidth
  //   t_hi = t_centre + tHalfWidth
  double tHalfWidth = 0.05;
};

// ============================================================
// DiffradCacheConfig
// ============================================================
struct DiffradCacheConfig {
  std::string outDir;
  std::string cacheRoot       = "phi_radcorr.root";
  std::string cacheTree       = "radcorr";
  std::string inmdiName       = "inmdi.dat";
  std::string diffradOutRoot  = "mdiffrad_output.root";
  std::string diffradOutTree  = "rcsum";
  std::string runMacroPath    = "/mnt/data/RunMDiffrad.C";
  std::string diffradWorkDir  = "";
  bool        storeBinIndices = true;
  double      maxMtprimeForDiffrad = 2.0;  // GeV^2
  DiffradRunParams runParams;
};

// ============================================================
// Path helpers
// ============================================================
inline static std::string JoinPath(const std::string& a, const std::string& b) {
  if (a.empty()) return b;
  return (a.back() == '/') ? a + b : a + "/" + b;
}
inline static std::string Dirname(const std::string& p) {
  try { fs::path fp(p); return fp.has_parent_path() ? fp.parent_path().string() : "."; }
  catch (...) { return "."; }
}
inline static std::string AbsPath(const std::string& p) {
  try { return fs::absolute(fs::path(p)).string(); } catch (...) { return p; }
}

// ============================================================
// MakePeriodConfig
// ============================================================
inline DiffradCacheConfig MakePeriodConfig(const std::string& baseDir,
                                           const std::string& tag) {
  DiffradCacheConfig cfg;
  cfg.outDir    = JoinPath(baseDir, tag);
  cfg.cacheRoot = "phi_radcorr_" + tag + ".root";
  return cfg;
}

// ============================================================
// Forward declarations
// ============================================================
inline static std::vector<std::vector<double>>
ComputeTminRef_QW(ROOT::RDF::RNode df_data, const BinManager& bins);

inline static void
WriteInmdiDat(const std::string& path,
              const DiffradRunParams& p,
              const std::vector<double>& W2min_v, const std::vector<double>& W2max_v,
              const std::vector<double>& Q2min_v, const std::vector<double>& Q2max_v,
              const std::vector<double>& tmin_v,  const std::vector<double>& tmax_v);

inline static void
BuildInmdiForBins(const BinManager& bins,
                  const std::vector<std::vector<double>>& tminRef,
                  const std::vector<std::vector<double>>& WmeanRef,
                  const std::string& inmdiPath,
                  const DiffradRunParams& runParams,
                  double maxMtprime,
                  std::vector<SentBin>& sentBins);

inline static void RunDiffrad(const DiffradCacheConfig& cfg);

inline static void
ConvertDiffradToCache(const DiffradCacheConfig& cfg,
                      const std::vector<SentBin>& sentBins);

// ============================================================
// Public entry point
// ============================================================
inline ROOT::RDF::RNode
GetOrBuildPhiRadCorrRDF(ROOT::RDF::RNode          df_data,
                        const BinManager&         bins,
                        const DiffradCacheConfig& cfg,
                        bool forceRecompute = false)
{
  fs::create_directories(cfg.outDir);
  const std::string cachePath = JoinPath(cfg.outDir, cfg.cacheRoot);

  if (!forceRecompute && fs::exists(cachePath)) {
    std::cout << "[PhiDiffrad] Cache found, loading: " << cachePath << "\n";
    ROOT::RDataFrame rdf(cfg.cacheTree, cachePath);
    return rdf;
  }

  std::cout << "[PhiDiffrad] Building cache: " << cfg.outDir << "\n";

  // Compute Mean("tmin") and Mean("W") per (Q2,W) bin from data
  // tminRef[iq][iw] < 0  — used directly as tmin_used (no kinematic override)
  // WmeanRef[iq][iw] > 0 — used as Wc for the narrow W window
  auto tminRef = ComputeTminRef_QW(df_data, bins);

  // Also compute mean W per bin so the W window is centred on data W-mean
  const auto& q2E = bins.GetQ2Bins();
  const auto& wE  = bins.GetWBins();
  const bool  hasW = !wE.empty();
  const size_t nQ  = (q2E.size() > 1) ? q2E.size() - 1 : 0;
  const size_t nW  = hasW ? wE.size() - 1 : 1;
  std::vector<std::vector<double>> WmeanRef(nQ, std::vector<double>(nW, 3.0));

  for (size_t iq = 0; iq < nQ; ++iq) {
    const double qLo = q2E[iq], qHi = q2E[iq + 1];
    auto df_q = df_data.Filter(
        [=](double Q2){ return Q2 > qLo && Q2 <= qHi; }, {"Q2"});
    for (size_t iw = 0; iw < nW; ++iw) {
      const double wLo = hasW ? wE[iw]     : 0.0;
      const double wHi = hasW ? wE[iw + 1] : 1e9;
      auto df_qw = df_q.Filter(
          [=](double W){ return !hasW || (W > wLo && W <= wHi); }, {"W"});
      const auto cnt = *df_qw.Count();
      WmeanRef[iq][iw] = (cnt > 0) ? *df_qw.Mean("W") : 0.5*(wLo + wHi);
      std::cout << "[PhiDiffrad] WmeanRef[" << iq << "][" << iw << "] = "
                << WmeanRef[iq][iw] << "  (N=" << cnt << ")\n";
    }
  }

  const std::string inmdiPath = JoinPath(cfg.outDir, cfg.inmdiName);
  std::vector<SentBin> sentBins;
  BuildInmdiForBins(bins, tminRef, WmeanRef, inmdiPath,
                    cfg.runParams, cfg.maxMtprimeForDiffrad, sentBins);

  RunDiffrad(cfg);
  ConvertDiffradToCache(cfg, sentBins);

  ROOT::RDataFrame rdf(cfg.cacheTree, cachePath);
  return rdf;
}

// ============================================================
// ComputeTminRef_QW
// Returns Mean("tmin") per (Q2,W) bin. tminRef[iq][iw] < 0.
// This is used directly as tmin_used — no kinematic override.
// ============================================================
inline static std::vector<std::vector<double>>
ComputeTminRef_QW(ROOT::RDF::RNode df_data, const BinManager& bins)
{
  const auto& q2E  = bins.GetQ2Bins();
  const auto& wE   = bins.GetWBins();
  const bool  hasW = !wE.empty();
  const size_t nQ  = (q2E.size() > 1) ? q2E.size() - 1 : 0;
  const size_t nW  = hasW ? wE.size() - 1 : 1;

  std::vector<std::vector<double>> tminRef(nQ, std::vector<double>(nW, 0.0));

  for (size_t iq = 0; iq < nQ; ++iq) {
    const double qLo = q2E[iq], qHi = q2E[iq + 1];
    auto df_q = df_data.Filter(
        [=](double Q2){ return Q2 > qLo && Q2 <= qHi; }, {"Q2"});
    for (size_t iw = 0; iw < nW; ++iw) {
      const double wLo = hasW ? wE[iw]     : 0.0;
      const double wHi = hasW ? wE[iw + 1] : 1e9;
      auto df_qw = df_q.Filter(
          [=](double W){ return !hasW || (W > wLo && W <= wHi); }, {"W"});
      const auto cnt   = *df_qw.Count();
      tminRef[iq][iw]  = (cnt > 0) ? *df_qw.Mean("tmin") : 0.0;
      std::cout << "[PhiDiffrad] tminRef[" << iq << "][" << iw << "] = "
                << tminRef[iq][iw] << "  (N=" << cnt << ")\n";
      if (cnt > 0 && tminRef[iq][iw] >= 0.0)
        std::cerr << "[PhiDiffrad] ERROR: tminRef[" << iq << "][" << iw
                  << "] >= 0! Check 'tmin' column sign.\n";
    }
  }
  return tminRef;
}

// ============================================================
// WriteInmdiDat
//
// Matches RunMDiffrad_backup.C WriteInmdiDat EXACTLY:
//   - 9 scalar header lines (inline ! comments OK)
//   - 6 data rows: bare numbers, NO inline ! comments
//     (Fortran READ is newline/whitespace-sensitive on array rows)
//
// Row order:
//   row 1: -W2_hi  (W2min_v[i], more negative)
//   row 2: -W2_lo  (W2max_v[i], less negative)
//   row 3: -Q2_hi  (Q2min_v[i])
//   row 4: -Q2_lo  (Q2max_v[i])  ← must never be 0.0
//   row 5:  t_lo   (tmin_v[i], more negative)
//   row 6:  t_hi   (tmax_v[i], less negative)
// ============================================================
inline static void
WriteInmdiDat(const std::string& path,
              const DiffradRunParams& p,
              const std::vector<double>& W2min_v, const std::vector<double>& W2max_v,
              const std::vector<double>& Q2min_v, const std::vector<double>& Q2max_v,
              const std::vector<double>& tmin_v,  const std::vector<double>& tmax_v)
{
  const int npoi = (int)W2min_v.size();
  if (npoi == 0)
    throw std::runtime_error("WriteInmdiDat: no kinematic points");
  if ((int)W2max_v.size() != npoi || (int)Q2min_v.size() != npoi ||
      (int)Q2max_v.size() != npoi || (int)tmin_v.size()  != npoi ||
      (int)tmax_v.size()  != npoi)
    throw std::runtime_error("WriteInmdiDat: vector size mismatch");

  std::ofstream f(path);
  if (!f) throw std::runtime_error("WriteInmdiDat: cannot open " + path);
  f << std::fixed << std::setprecision(10);

  // 9 scalar header lines
  f << p.bmom   << "     !  bmom\n";
  f << p.tmom   << "     !  tmom\n";
  f << p.lepton << "         !  lepton\n";
  f << p.ivec   << "         !  ivec\n";
  f << p.ann1 << " " << p.ann2 << " " << p.ann3 << "   ! ann1 ann2 ann3\n";
  f << p.vcut   << "        !  vcut\n";
  f << p.nev    << "         !  nev\n";
  f << p.seed   << "    !  seed\n";
  f << npoi     << "         !  npoi\n";

  // 6 data rows: NO inline comments (matches backup writeRow exactly)
  auto writeRow = [&](const std::vector<double>& v) {
    for (int i = 0; i < npoi; ++i)
      f << v[i] << (i + 1 < npoi ? " " : "");
    f << "\n";
  };
  writeRow(W2min_v);
  writeRow(W2max_v);
  writeRow(Q2min_v);
  writeRow(Q2max_v);
  writeRow(tmin_v);
  writeRow(tmax_v);

  // ── Diagnostic table ────────────────────────────────────────────────────
  std::cout << "\n[WriteInmdiDat] Written: " << path
            << "  (npoi=" << npoi << ", nev=" << p.nev << ")\n\n";
  std::cout << std::fixed << std::setprecision(5);
  std::cout
    << "  idx |    W2 range [GeV^2]      |    Q2 range [GeV^2]     "
    << "|      t range [GeV^2]      | Wc[GeV] Qc[GeV] tc[GeV^2]\n"
    << "  ----|--------------------------|-------------------------"
    << "|---------------------------|---------------------------\n";

  for (int i = 0; i < npoi; ++i) {
    double W2_lo = -W2max_v[i],  W2_hi = -W2min_v[i];
    double Q2_lo = -Q2max_v[i],  Q2_hi = -Q2min_v[i];
    double Wc = std::sqrt(std::max(0.0, 0.5*(W2_lo + W2_hi)));
    double Qc = std::sqrt(std::max(0.0, 0.5*(Q2_lo + Q2_hi)));
    double tc = 0.5*(tmin_v[i] + tmax_v[i]);
    std::cout
      << "  " << std::setw(3) << i
      << " | [" << std::setw(7) << W2_lo << ", " << std::setw(7) << W2_hi << "]"
      << " | [" << std::setw(7) << Q2_lo << ", " << std::setw(7) << Q2_hi << "]"
      << " | [" << std::setw(9) << tmin_v[i] << ", " << std::setw(9) << tmax_v[i] << "]"
      << " | Wc=" << std::setw(5) << Wc
      << " Qc=" << std::setw(5) << Qc
      << " tc=" << std::setw(8) << tc
      << "\n";
  }
  std::cout << std::defaultfloat << "\n";
}

// ============================================================
// BuildInmdiForBins
//
// KEY DESIGN DECISIONS (learned from crash analysis):
//
// (1) tmin_used = Mean("tmin") from DATA (tminRef[iq][iw]).
//     Do NOT override with the kinematic formula at Wc.
//     Reason: with a single wide W bin [2,10], Wc=6 GeV gives
//     tmin_kin ≈ -0.003 (nearly 0), which is wrong. The data mean
//     tminRef correctly reflects the actual event kinematics.
//
// (2) Q2 lower edge clamped to q2LoMin (default 0.5 GeV^2).
//     Reason: Q2_lo=0 causes SIGFPE in mdiffrad.f at line 178
//     (log(Q2) or 1/Q2 singularity). The analysis bin [0,2.27]
//     would pass Q2_lo=0 → Fortran crash → only 30/96 ALLmc rows.
//
// (3) W window: Wc ± wHalfWidth, centred on data Mean("W").
//     Reason: sending the full [2,10] GeV bin (W2 in [4,100])
//     makes DIFFRAD average over enormous phase space.
// ============================================================
inline static void
BuildInmdiForBins(const BinManager& bins,
                  const std::vector<std::vector<double>>& tminRef,
                  const std::vector<std::vector<double>>& WmeanRef,
                  const std::string& inmdiPath,
                  const DiffradRunParams& runParams,
                  double maxMtprime,
                  std::vector<SentBin>& sentBins)
{
  const auto& q2E  = bins.GetQ2Bins();
  const auto& wE   = bins.GetWBins();
  const auto& mtpE = bins.GetTprimeBins();

  const bool   hasW = !wE.empty();
  const size_t nQ   = (q2E.size() > 1) ? q2E.size() - 1 : 0;
  const size_t nW   = hasW ? wE.size() - 1 : 1;
  const size_t nT   = (mtpE.size() > 1) ? mtpE.size() - 1 : 0;

  if (nQ == 0 || nT == 0)
    throw std::runtime_error("BuildInmdiForBins: Q2 or tprime bins are empty");

  std::vector<double> W2min_v, W2max_v, Q2min_v, Q2max_v, tmin_v, tmax_v;
  sentBins.clear();

  std::cout << "\n[BuildInmdiForBins] wHalfWidth=" << runParams.wHalfWidth
            << " GeV  q2LoMin=" << runParams.q2LoMin
            << " GeV^2  maxMtprime=" << maxMtprime << " GeV^2\n\n";

  for (size_t iq = 0; iq < nQ; ++iq) {
    const double qLo_bin = q2E[iq], qHi_bin = q2E[iq + 1];
    const double Q2c     = 0.5*(qLo_bin + qHi_bin);

    // Clamp Q2 lower edge — prevents SIGFPE in Fortran when qLo=0
    const double qLo_send = std::max(qLo_bin, runParams.q2LoMin);
    const double qHi_send = qHi_bin;

    if (qLo_send != qLo_bin)
      std::cout << "[BuildInmdiForBins] Q2 lower edge clamped: "
                << qLo_bin << " → " << qLo_send
                << " (q2LoMin=" << runParams.q2LoMin << ")\n";

    for (size_t iw = 0; iw < nW; ++iw) {
      const double wLo_bin = hasW ? wE[iw]     : 2.0;
      const double wHi_bin = hasW ? wE[iw + 1] : 4.5;

      // Use data Mean("W") as the W window centre — not the arithmetic bin centre.
      // For a wide bin [2,10], arithmetic centre = 6 GeV gives wrong tmin.
      // Data mean is where the events actually are.
      const double Wc_data = (iq < WmeanRef.size() && iw < WmeanRef[iq].size())
                              ? WmeanRef[iq][iw]
                              : 0.5*(wLo_bin + wHi_bin);
      const double hw      = runParams.wHalfWidth;
      const double wLo_send = std::max(wLo_bin, Wc_data - hw);
      const double wHi_send = std::min(wHi_bin, Wc_data + hw);
      const double W2lo_send = wLo_send * wLo_send;
      const double W2hi_send = wHi_send * wHi_send;

      // tmin_used = Mean("tmin") from data — used directly, no kinematic override.
      // The kinematic override at Wc=6 GeV gives nearly zero tmin which is wrong.
      const double tmin_used = (iq < tminRef.size() && iw < tminRef[iq].size())
                                ? tminRef[iq][iw] : 0.0;

      if (!std::isfinite(tmin_used) || !(tmin_used < 0.0)) {
        std::cerr << "[BuildInmdiForBins] WARNING: non-physical tmin_used="
                  << tmin_used << " (iq=" << iq << ",iw=" << iw
                  << "). Skipping this (Q2,W) cell.\n";
        continue;
      }

      std::cout << std::fixed << std::setprecision(5)
                << "[BuildInmdiForBins] (iq=" << iq << ",iw=" << iw << ")"
                << "  Q2_bin=[" << qLo_bin << "," << qHi_bin
                << "]  Q2_send=[" << qLo_send << "," << qHi_send << "]"
                << "  Wc_data=" << Wc_data
                << "  W_send=[" << wLo_send << "," << wHi_send << "]"
                << "  tmin_used=" << tmin_used << "\n"
                << std::defaultfloat;

      for (size_t it = 0; it < nT; ++it) {
        const double mtp_lo = mtpE[it];
        const double mtp_hi = mtpE[it + 1];

        if (mtp_hi > maxMtprime) {
          std::cout << "  [skip] t' bin [" << mtp_lo << "," << mtp_hi
                    << "] > maxMtprime=" << maxMtprime << "\n";
          continue;
        }

        // Convert t-prime bin centre to a narrow physical t window.
        // Use bin MIDPOINT, not full edge-to-edge range.
        // DIFFRAD sig_rad_ diverges for wide t ranges -> RC>>1 and SIGFPE.
        // Matches backup style: t_centre +/- tHalfWidth (default 0.05 GeV^2).
        const double mtp_mid  = 0.5 * (mtp_lo + mtp_hi);
        const double t_centre = tmin_used - mtp_mid;  // physical t bin centre
        const double t_phys_lo = t_centre - runParams.tHalfWidth;  // more negative
        double       t_phys_hi = t_centre + runParams.tHalfWidth;  // less negative

        if (!std::isfinite(t_centre) || !(t_centre < 0.0)) {
          std::cerr << "[BuildInmdiForBins] WARNING: non-physical t_centre="
                    << t_centre << " tmin=" << tmin_used
                    << " mtp_mid=" << mtp_mid
                    << " (iq=" << iq << ",iw=" << iw << ",it=" << it << "). Skip.\n";
          continue;
        }
        // Clamp t_hi to stay below tmin (physical boundary)
        if (t_phys_hi >= tmin_used) t_phys_hi = tmin_used - 1e-6;

        // DIFFRAD convention: store negated W2 and Q2
        W2min_v.push_back(-W2hi_send);   // -W2_hi  (more negative)
        W2max_v.push_back(-W2lo_send);   // -W2_lo  (less negative)
        Q2min_v.push_back(-qHi_send);    // -Q2_hi
        Q2max_v.push_back(-qLo_send);    // -Q2_lo  (never 0 thanks to clamp)
        tmin_v.push_back(t_phys_lo);
        tmax_v.push_back(t_phys_hi);

        sentBins.push_back({iq, iw, it, tmin_used, mtp_lo, mtp_hi, Q2c, Wc_data});

        const double tc = 0.5*(t_phys_lo + t_phys_hi);
        std::cout << std::fixed << std::setprecision(5)
                  << "  [SEND #" << sentBins.size()-1 << "]"
                  << "  t'=[" << mtp_lo << "," << mtp_hi << "] mid=" << mtp_mid
                  << "  t_centre=" << t_centre
                  << "  t_send=[" << t_phys_lo << "," << t_phys_hi << "]"
                  << "  Q2c=" << Q2c
                  << "  Wc=" << Wc_data
                  << "  tmin=" << tmin_used
                  << "\n" << std::defaultfloat;
      }
    }
  }

  if (sentBins.empty())
    throw std::runtime_error(
        "BuildInmdiForBins: no valid bins (check maxMtprime / tmin sign)");

  std::cout << "\n[BuildInmdiForBins] Total sent to DIFFRAD: "
            << sentBins.size() << " points\n\n";

  WriteInmdiDat(inmdiPath, runParams,
                W2min_v, W2max_v, Q2min_v, Q2max_v, tmin_v, tmax_v);
}

// ============================================================
// RunDiffrad
// ============================================================
inline static void RunDiffrad(const DiffradCacheConfig& cfg)
{
  const std::string workDir = !cfg.diffradWorkDir.empty()
                              ? cfg.diffradWorkDir
                              : Dirname(cfg.runMacroPath);
  if (!fs::exists(workDir))
    throw std::runtime_error("RunDiffrad: workDir not found: " + workDir);

  const std::string inmdiSrc    = AbsPath(JoinPath(cfg.outDir, cfg.inmdiName));
  const std::string inmdiInWork = JoinPath(workDir, "inmdi.dat");
  if (!fs::exists(inmdiSrc))
    throw std::runtime_error("RunDiffrad: inmdi.dat not found: " + inmdiSrc);

  std::cout << "[RunDiffrad] " << inmdiSrc << " → " << inmdiInWork << "\n";
  fs::copy_file(inmdiSrc, inmdiInWork, fs::copy_options::overwrite_existing);

  TString old = gSystem->WorkingDirectory();
  gSystem->ChangeDirectory(workDir.c_str());

  if (!fs::exists(cfg.runMacroPath)) {
    gSystem->ChangeDirectory(old.Data());
    throw std::runtime_error("RunDiffrad: macro not found: " + cfg.runMacroPath);
  }

  Int_t err = 0;
  gROOT->ProcessLine(Form(".L %s+", cfg.runMacroPath.c_str()), &err);
  if (err) {
    gSystem->ChangeDirectory(old.Data());
    throw std::runtime_error("RunDiffrad: compile failed err=" + std::to_string(err));
  }

  std::cout << "[RunDiffrad] Calling RunMDiffrad(\""
            << inmdiInWork << "\", false, false)\n";
  gROOT->ProcessLine(
      Form("RunMDiffrad(\"%s\", false, false)", inmdiInWork.c_str()), &err);
  gSystem->ChangeDirectory(old.Data());

  const std::string produced = JoinPath(workDir, cfg.diffradOutRoot);
  const std::string target   = JoinPath(cfg.outDir, cfg.diffradOutRoot);
  if (!fs::exists(produced))
    throw std::runtime_error("RunDiffrad: output not found: " + produced);

  fs::copy_file(produced, target, fs::copy_options::overwrite_existing);
  std::cout << "[RunDiffrad] Output → " << target << "\n";
}

// ============================================================
// ConvertDiffradToCache
//
// Cache columns:
//   Q2          bin-centre Q^2 [GeV^2]
//   W           data Mean("W") for this bin [GeV]
//   t_centre    physical t centre from DIFFRAD [GeV^2, NEGATIVE]
//   tmin_used   Mean("tmin") from data [GeV^2, NEGATIVE]
//   mtprime     t' = tmin_used - t_phys [GeV^2, POSITIVE]
//   rad_corr    radiative correction factor
//   rad_corr_err sigma over nev samples
//   iq,iw,it    bin indices
// ============================================================
inline static void
ConvertDiffradToCache(const DiffradCacheConfig& cfg,
                      const std::vector<SentBin>& sentBins)
{
  const std::string diffradPath = JoinPath(cfg.outDir, cfg.diffradOutRoot);
  if (!fs::exists(diffradPath))
    throw std::runtime_error("ConvertDiffradToCache: not found: " + diffradPath);

  TFile fin(diffradPath.c_str(), "READ");
  if (fin.IsZombie())
    throw std::runtime_error("ConvertDiffradToCache: cannot open: " + diffradPath);

  TTree* tin = dynamic_cast<TTree*>(fin.Get(cfg.diffradOutTree.c_str()));
  if (!tin)
    throw std::runtime_error(
        "ConvertDiffradToCache: tree '" + cfg.diffradOutTree + "' not found");

  Double_t Q_in = 0, W_in = 0, t_in = 0, rc_in = 0, rcs_in = 0;
  tin->SetBranchAddress("Q",              &Q_in);
  tin->SetBranchAddress("W",              &W_in);
  tin->SetBranchAddress("t",              &t_in);
  tin->SetBranchAddress("rad_corr",       &rc_in);
  tin->SetBranchAddress("rad_corr_sigma", &rcs_in);

  const Long64_t N = tin->GetEntries();
  if (N != (Long64_t)sentBins.size())
    throw std::runtime_error(
        "ConvertDiffradToCache: rcsum entries=" + std::to_string(N) +
        " but sentBins=" + std::to_string(sentBins.size()) +
        " — mapping broken (DIFFRAD likely crashed for some points)");

  const std::string cachePath = JoinPath(cfg.outDir, cfg.cacheRoot);
  TFile fout(cachePath.c_str(), "RECREATE");
  if (fout.IsZombie())
    throw std::runtime_error("ConvertDiffradToCache: cannot create: " + cachePath);

  TTree tout(cfg.cacheTree.c_str(), "Phi rad-corr cache (DISANA_PhiDiffrad)");

  Double_t Q2_out, W_out, t_centre_out, tmin_out, mtprime_out, rc_out, rce_out;
  Int_t    iq_out, iw_out, it_out;

  tout.Branch("Q2",           &Q2_out,       "Q2/D");
  tout.Branch("W",            &W_out,        "W/D");
  tout.Branch("t_centre",     &t_centre_out, "t_centre/D");
  tout.Branch("tmin_used",    &tmin_out,     "tmin_used/D");
  tout.Branch("mtprime",      &mtprime_out,  "mtprime/D");
  tout.Branch("rad_corr",     &rc_out,       "rad_corr/D");
  tout.Branch("rad_corr_err", &rce_out,      "rad_corr_err/D");
  if (cfg.storeBinIndices) {
    tout.Branch("iq", &iq_out, "iq/I");
    tout.Branch("iw", &iw_out, "iw/I");
    tout.Branch("it", &it_out, "it/I");
  }

  std::cout << "\n[ConvertDiffradToCache] " << N << " entries\n";
  std::cout << std::fixed << std::setprecision(5);
  std::cout
    << "  entry | iq iw it | Q2c    Wc_data | tmin_used | t_phys   | mtprime  | rad_corr\n"
    << "  ------|----------|----------------|-----------|----------|----------|---------\n";

  for (Long64_t i = 0; i < N; ++i) {
    tin->GetEntry(i);

    const SentBin& sb = sentBins[i];
    iq_out = (Int_t)sb.iq;
    iw_out = (Int_t)sb.iw;
    it_out = (Int_t)sb.it;

    Q2_out       = sb.Q2c;
    W_out        = sb.Wc;
    tmin_out     = sb.tmin_used;
    t_centre_out = t_in;
    mtprime_out  = sb.tmin_used - t_in;   // t' = tmin_used - t_phys ≥ 0
    rc_out       = rc_in;
    rce_out      = rcs_in;

    std::cout
      << "  " << std::setw(5) << i
      << " | " << std::setw(2) << iq_out
      << " " << std::setw(2) << iw_out
      << " " << std::setw(2) << it_out
      << " | " << std::setw(6) << Q2_out
      << "  " << std::setw(6) << W_out
      << " | " << std::setw(9) << tmin_out
      << " | " << std::setw(8) << t_centre_out
      << " | " << std::setw(8) << mtprime_out
      << " | " << std::setw(8) << rc_out
      << "\n";

    tout.Fill();
  }

  tout.Write();
  fout.Close();
  fin.Close();

  std::cout << std::defaultfloat
            << "\n[PhiDiffrad] Cache written: " << cachePath
            << "  (" << N << " entries)\n"
            << "[PhiDiffrad] Columns: Q2, W, t_centre, tmin_used, mtprime, "
            << "rad_corr, rad_corr_err"
            << (cfg.storeBinIndices ? ", iq, iw, it" : "") << "\n";
}

}  // namespace PhiDiffrad
}  // namespace DISANA

#endif  // DISANA_PHI_DIFFRAD_H