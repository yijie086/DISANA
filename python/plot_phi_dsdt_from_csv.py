#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_phi_dsdt_from_csv.py
=========================

Complete  d/dt analysis  plots and dipole fits.

Fixes in this version:
  - Per-model beam energy auto-detection (Sp1810.6, Fall1810.594, Sp1910.2)
  - Q2_center / W_center / tprime_center derived from lo/hi if missing in CSV
  - NEW: If CSV omits W and/or Q2 (because binning is implicit), force:
      * W_center = hardcoded W_MEAN_GEV (requested: 2.8)
      * Q2_center = hardcoded per-iq mapping Q2_MEAN_BY_IQ_GEV2
    so Gamma_v does not become NaN just because W/Q2 columns are absent.
  - RadCorr matched by bin-containment using actual DIFFRAD column names
    (Q2_lo/hi/c, tprime_lo/hi/c, rad_corr) with comment='#' handling;
    falls back to nearest centroid with a printed warning if no bin matches
  - --beam-energy kept as optional fallback for unknown model names
  - --print-table / --table-outdir : prints + saves per-bin tables with
      raw_xs, Gamma_v, reduced_xs, Acceptance, Efficiency,
      xs_acc_corrected, xs_eff_corrected, RadCorr, xs_rad_corrected
  - Corrections panel ratio is a TRUE correction-factor plot:
      ratio = Final / Variant   (log-scale)
  - Acceptance / Efficiency group y-range forced for correction-factor plots
"""
#python3 ./../../../../../source/python/plot_phi_dsdt_from_csv.py --csv-root ./../ --models "Sp18_inb" "Sp18_outb" "Fall18_inb" "Fall18_outb" "Sp19_inb" --outdir plots_phi_dsdt --plot-reduced --use-external-radcorr --print-table
#python3 ./../../../../../source/python/plot_phi_dsdt_from_csv.py --csv-root ./../ --models "Sp18_inb" "Sp18_outb" "Fall18_inb" "Fall18_outb" "Sp19_inb" --outdir plots_phi_dsdt --plot-reduced --use-external-radcorr --print-table

import os
import re
import argparse
from glob import glob
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.ticker import AutoMinorLocator

# -----------------------------------------------------------------------------
#  External radiative correction tables (hard-coded paths)
# -----------------------------------------------------------------------------
K_RADCORR_FILE_10P2 = "/w/hallb-scshelf2102/clas12/singh/Softwares/Generators/PhiEventGen/DIFFRAD/diffradstn/outputs/10p2GeV/diffrad_rc_results.csv"
K_RADCORR_FILE_10P6 = "/w/hallb-scshelf2102/clas12/singh/Softwares/Generators/PhiEventGen/DIFFRAD/diffradstn/outputs/10p6GeV/diffrad_rc_results.csv"
# -----------------------------------------------------------------------------
W_MEAN_GEV = (1.8 + 2.5) / 2  # hardcoded mean W (requested)
# -----------------------------------------------------------------------------
#  Per-model beam energy lookup  (case-insensitive substring match on model name)
# -----------------------------------------------------------------------------
K_BEAM_ENERGY = {
    "sp18":    10.6,
    "sp2018":  10.6,
    "fall18":  10.594,
    "fall2018":10.594,
    "sp19":    10.2,
    "sp2019":  10.2,
}

def beam_energy_for_model(model_name):
    """Return beam energy [GeV] inferred from model name, else fall back to CLI global."""
    name_lower = model_name.lower()
    for key, energy in K_BEAM_ENERGY.items():
        if key in name_lower:
            return energy
    return G_BEAM_ENERGY_GEV   # CLI fallback (may be None)


# -----------------------------------------------------------------------------
#  Hatta–Strikman (HS) Ds(0) extraction (shape-template)
# -----------------------------------------------------------------------------
# HS defaults (Hatta–Strikman model; dipole A_s, tripole D_s)
HS_AS0_DEFAULT = 0.04
HS_mA_DEFAULT  = 1.13  # GeV
HS_mD_DEFAULT  = 0.76  # GeV

# -----------------------------------------------------------------------------
#  Lattice QCD reference values — Hackett, Pefkou, Shanahan
#  Phys. Rev. Lett. 132, 251904 (2024), Table I
#  MS-bar scheme at mu = 2 GeV, m_pi ~ 170 MeV ensemble
# -----------------------------------------------------------------------------
# Forward limits (dipole fit)
LAT_As0_DIPOLE  =  0.0257;  LAT_dAs0_DIPOLE  = 0.0095
LAT_Ds0_DIPOLE  = -0.18;    LAT_dDs0_DIPOLE  = 0.17
# Forward limits (z-expansion fit)
LAT_As0_ZEXP    =  0.032;   LAT_dAs0_ZEXP    = 0.012
LAT_Ds0_ZEXP    = -0.08;    LAT_dDs0_ZEXP    = 0.17
# Dipole mass^2 estimates for t-dependence of strange quark GFFs
# (exact values are in the Supplemental Material; these are estimated
#  from Fig. 1 of the paper — use with caution for shape comparison only)
LAT_mA2_LO  = 0.40   # GeV^2  lower bound on A_s dipole mass^2
LAT_mA2_HI  = 1.10   # GeV^2  upper bound
LAT_mD2_LO  = 0.35   # GeV^2  lower bound on D_s dipole mass^2
LAT_mD2_HI  = 0.85   # GeV^2  upper bound

def _lat_dipole(t_abs, G0, m2):
    """Dipole: G(t) = G0 / (1 + |t|/m^2)^2  (note: lattice uses dipole for D_s too)."""
    return G0 / (1.0 + t_abs / m2)**2

def _lat_band(t_arr, G0, dG0, m2_lo, m2_hi, n_sigma=1):
    """Return (central, lo_envelope, hi_envelope) for a lattice GFF band,
    varying G0 by +/-n_sigma*dG0 and m^2 across [m2_lo, m2_hi]."""
    curves = []
    for G0_v in [G0 - n_sigma*dG0, G0, G0 + n_sigma*dG0]:
        for m2 in np.linspace(m2_lo, m2_hi, 8):
            curves.append(_lat_dipole(t_arr, G0_v, m2))
    curves = np.array(curves)
    central = _lat_dipole(t_arr, G0, 0.5*(m2_lo + m2_hi))
    return central, np.min(curves, axis=0), np.max(curves, axis=0)

def hs_As(t, As0=HS_AS0_DEFAULT, mA=HS_mA_DEFAULT):
    """HS dipole form: A_s(t)=A_s(0)/(1 - t/mA^2)^2.  t is negative (spacelike)."""
    t = np.asarray(t, dtype=float)
    den = (1.0 - t/(mA*mA))
    return As0 / (den*den)

def hs_Ds(t, Ds0, mD=HS_mD_DEFAULT):
    """HS tripole form: D_s(t)=D_s(0)/(1 - t/mD^2)^3.  t is negative (spacelike)."""
    t = np.asarray(t, dtype=float)
    den = (1.0 - t/(mD*mD))
    return Ds0 / (den*den*den)

def hs_dsdt_template(t_abs, N, Ds0,
                     As0=HS_AS0_DEFAULT, mA=HS_mA_DEFAULT, mD=HS_mD_DEFAULT,
                     cAA=1.0, cAD=1.0, cDD=1.0):
    """HS-inspired shape template for dσ/dt in the same units as your data (nb/GeV^2).

    We evaluate form factors at spacelike t = -|t|:

        dσ/dt(|t|) = N * [ cAA*A_s(t)^2 + cAD*|t|*A_s(t)*D_s(t) + cDD*|t|^2*D_s(t)^2 ].

    N absorbs the overall normalization uncertainty; Ds0 is the extracted strange D-term at t=0.
    """
    t_abs = np.asarray(t_abs, dtype=float)
    t = -t_abs
    A = hs_As(t, As0=As0, mA=mA)
    D = hs_Ds(t, Ds0=Ds0, mD=mD)
    poly = (cAA*(A*A) + cAD*(t_abs)*(A*D) + cDD*(t_abs*t_abs)*(D*D))
    return N * np.maximum(poly, 0.0)

# Global run configuration set in main()
G_BEAM_ENERGY_GEV      = None   # float   fallback only; per-model takes priority
G_USE_EXTERNAL_RADCORR = False  # bool
G_RADCORR_FILE_OVERRIDE= None   # optional str
G_RADCORR_MATCH_MODE   = "nearest"


# -----------------------------------------------------------------------------
#  HS fit plotting helpers
# -----------------------------------------------------------------------------
def _make_hs_fit_plot(df, res, iq, iw, model_label, outdir, logy=True,
                      As0=HS_AS0_DEFAULT, mA=HS_mA_DEFAULT, mD=HS_mD_DEFAULT):
    """Make a 2-panel HS fit plot (data + fit, pull) for one (Q,W) bin.

    Plots are shown vs t' (experimental variable). The HS model is evaluated in
    |t| = t' + |t_min| (same convention used in the fit).
    """
    if df is None or res is None or "ReducedCrossSection" not in df.columns:
        return

    os.makedirs(outdir, exist_ok=True)

    tprime = df["tprime_center"].to_numpy(float)
    xs     = df["ReducedCrossSection"].to_numpy(float)
    xse    = (df["ReducedCrossSection_Err"].to_numpy(float)
              if "ReducedCrossSection_Err" in df.columns else np.zeros_like(xs))
    # 10% relative floor if missing/non-positive (avoids divide-by-zero in pulls)
    with np.errstate(invalid="ignore"):
        xse = np.where(np.isfinite(xse) & (xse > 0),
                       xse,
                       0.10 * np.where(np.isfinite(xs) & (xs > 0), xs, 1.0))

    # Use per-row tmin_abs if present, else fall back to 0
    tmin_abs = float(df["tmin_abs"].iloc[0]) if "tmin_abs" in df.columns else float(res.get("tmin_abs", 0.0))
    if not np.isfinite(tmin_abs):
        tmin_abs = 0.0

    # Apply same fit-range cuts (in t') to decide which points to show as 'fit points'
    tprime_min = float(res.get("tprime_min", np.nan))
    tprime_max = float(res.get("tprime_max", np.nan))

    mask = np.isfinite(tprime) & np.isfinite(xs) & (xs > 0) & np.isfinite(xse)
    if not mask.any():
        return

    # error floor (avoid 0 error -> inf pull)
    with np.errstate(invalid="ignore"):
        xse = np.where(np.isfinite(xse) & (xse > 0), xse, 0.10 * np.where(xs > 0, xs, 1.0))
    xse = _safe_sigma(xse)

    # x-errors from bin widths if available
    if "tprime_lo" in df.columns and "tprime_hi" in df.columns:
        xlo = tprime - df["tprime_lo"].to_numpy(float)
        xhi = df["tprime_hi"].to_numpy(float) - tprime
    else:
        xlo = np.zeros_like(tprime)
        xhi = np.zeros_like(tprime)

    # Dense curve in t' for display
    tp_min = max(0.0, float(np.nanmin(tprime[mask])) * 0.8)
    tp_max = float(np.nanmax(tprime[mask])) * 1.15
    tprime_dense = np.linspace(tp_min, tp_max, 600)
    t_abs_dense  = tprime_dense + tmin_abs

    y_dense = hs_dsdt_template(
        t_abs_dense,
        res["N"], res["Ds0"],
        As0=As0, mA=mA, mD=mD
    )

    # Predictions at data points for pulls
    t_abs = tprime + tmin_abs
    y_pred = hs_dsdt_template(
        t_abs,
        res["N"], res["Ds0"],
        As0=As0, mA=mA, mD=mD
    )
    pull = (xs - y_pred) / xse

    # Figure
    with mpl.rc_context(PANEL_RC):
        fig = plt.figure(figsize=(10, 12))
        gs  = mpl.gridspec.GridSpec(
            2, 1, figure=fig, height_ratios=[3, 1], hspace=0.0,
            left=0.12, right=0.97, top=0.92, bottom=0.09
        )
        ax_m = fig.add_subplot(gs[0])
        ax_r = fig.add_subplot(gs[1], sharex=ax_m)
        _style_ax(ax_m); _style_ax(ax_r)

        # Data
        ax_m.errorbar(
            tprime[mask], xs[mask],
            yerr=xse[mask],
            xerr=(xlo[mask], xhi[mask]),
            fmt="o", color="#111111", ecolor="#111111",
            elinewidth=1.8, capsize=4,
            markersize=8, mfc="white", mew=2.2, lw=1.8, zorder=6,
            label="data",
        )

        # Fit curve
        ax_m.plot(tprime_dense, y_dense, color="#CC1A33", lw=2.5, ls="--", zorder=5, label="HS fit")

        # Fit range shading (if available)
        if np.isfinite(tprime_min) and np.isfinite(tprime_max):
            ax_m.axvspan(tprime_min, tprime_max, alpha=0.06, color="#CC1A33", zorder=0)

        if logy:
            ax_m.set_yscale("log")
            ax_m.set_ylim(YLIM_XS)

        ax_m.set_ylabel(r"$\mathrm{d}\sigma/\mathrm{d}t'\ [\mathrm{nb/GeV}^2]$")
        ax_m.tick_params(labelbottom=False)
        ax_m.set_ylim(YLIM_XS)

        # Annotation
        chi2_str = (f"chi2/ndf = {res['chi2']:.1f}/{res['ndf']}" if res.get("ndf", 0) > 0 else "")
        txt = (
            f"Ds(0) = {_fmt_pm(res['Ds0'], res['dDs0'])}\n"
            f"N = {_fmt_pm(res['N'], res['dN'])}\n"
            f"|t_min| = {tmin_abs:.4f} GeV^2\n"
            + chi2_str
        ).strip()
        ax_m.text(
            0.97, 0.97, txt,
            transform=ax_m.transAxes,
            fontsize=14, va="top", ha="right", color="#111111",
            bbox=dict(boxstyle="round,pad=0.35", fc="white",
                      ec="#111111", alpha=0.90, lw=1.2),
        )

        # Pulls
        ax_r.axhline(0.0, color="gray", lw=1.2, ls="--", zorder=1)
        ax_r.axhline( 1.0, color="gray", lw=0.8, ls=":", zorder=1)
        ax_r.axhline(-1.0, color="gray", lw=0.8, ls=":", zorder=1)
        ax_r.errorbar(
            tprime[mask], pull[mask],
            yerr=np.ones(np.sum(mask)),
            xerr=(xlo[mask], xhi[mask]),
            fmt="o", color="#111111", ecolor="#111111",
            elinewidth=1.2, capsize=3,
            markersize=6, mfc="white", mew=2.0, zorder=5,
        )
        ax_r.set_ylim(-4.5, 4.5)
        ax_r.set_ylabel("Pull", fontsize=17)
        ax_r.set_xlabel(r"$-t'\ [\mathrm{GeV}^2]$")

        # Title
        title = f"HS Ds(0) fit  {model_label.replace('_',' ')}  Q{iq}" + (f" W{iw}" if iw is not None else "")
        fig.suptitle(title, fontsize=LEG_TITLE_FS + 2, y=0.965)

        # Save
        stem = f"hs_fit_{_sanitise(model_label)}_Q{iq}" + (f"_W{iw}" if iw is not None else "")
        _save(fig, outdir, stem)

ALPHA_EM  = 1/137.035999084
M_P_GEV   = 0.9382720813
M_PHI_GEV = 1.019461        # (1020) mass [GeV], PDG 2022

# -----------------------------------------------------------------------------
#  Hardcoded kinematics fallbacks (TOP, as requested)


# If Q2_center is missing, use an IQ->Q2 mean lookup (EDIT if needed)
Q2_MEAN_BY_IQ_GEV2 = {
    0: 1.0,
    1: 2.0,
    2: 3.0,
    3: 4.0,
}

# -----------------------------------------------------------------------------
#  Physics helpers
# -----------------------------------------------------------------------------

def _make_hs_gff_plot(df, res, iq, iw, model_label, outdir,
                      As0=HS_AS0_DEFAULT, mA=HS_mA_DEFAULT, mD=HS_mD_DEFAULT):
    """Plot the gravitational form factors A_s(t) (dipole) and D_s(t) (tripole).

    Uses the HS parameterizations:
        A_s(t) = A_s(0) / (1 - t/mA^2)^2
        D_s(t) = D_s(0) / (1 - t/mD^2)^3

    We plot vs |t| (GeV^2) with the physical spacelike substitution t = -|t|.
    The plotted |t| range is taken from the data points in the bin, extended slightly.
    """
    if df is None or res is None:
        return
    os.makedirs(outdir, exist_ok=True)

    # Determine |t| range from the bin
    if "tprime_center" in df.columns:
        tprime = df["tprime_center"].to_numpy(float)
    else:
        return
    tmin_abs = float(df["tmin_abs"].iloc[0]) if "tmin_abs" in df.columns else float(res.get("tmin_abs", 0.0))
    if not np.isfinite(tmin_abs):
        tmin_abs = 0.0
    t_abs_pts = tprime + tmin_abs
    t_abs_pts = t_abs_pts[np.isfinite(t_abs_pts)]
    if t_abs_pts.size == 0:
        return

    t_lo = max(0.0, float(np.nanmin(t_abs_pts)) * 0.8)
    t_hi = float(np.nanmax(t_abs_pts)) * 1.2
    if not np.isfinite(t_hi) or t_hi <= 0:
        return
    t_abs = np.linspace(t_lo, t_hi, 600)
    t = -t_abs  # spacelike Mandelstam t

    As = hs_As(t, As0=As0, mA=mA)
    Ds = hs_Ds(t, Ds0=float(res.get("Ds0", np.nan)), mD=mD)

    # Lattice bands — restrict to |t| <= 2 GeV^2 (lattice kinematic range)
    t_lat = np.linspace(0.0, min(t_hi, 2.0), 400)
    Ac, Alo, Ahi = _lat_band(t_lat, LAT_As0_DIPOLE, LAT_dAs0_DIPOLE,
                             LAT_mA2_LO, LAT_mA2_HI)
    Dc, Dlo, Dhi = _lat_band(t_lat, LAT_Ds0_DIPOLE, LAT_dDs0_DIPOLE,
                             LAT_mD2_LO, LAT_mD2_HI)

    stem_base = f"hs_gff_{_sanitise(model_label)}_Q{iq}" + (f"_W{iw}" if iw is not None else "")
    title = f"HS form factors  {model_label}  Q{iq}" + (f"_W{iw}" if iw is not None else "")

    for show_lattice in (False, True):
        with mpl.rc_context(PANEL_RC):
            fig, ax = plt.subplots(figsize=(9, 6))
            _style_ax(ax)

            if show_lattice:
                # Lattice bands drawn first, behind CLAS12 curves
                ax.fill_between(t_lat, Alo, Ahi,
                                color="steelblue", alpha=0.20, zorder=1,
                                label=r"$A_s^{\rm lat}\ 1\sigma$ (Hackett et al. 2024)")
                ax.plot(t_lat, Ac, color="steelblue", lw=1.5, ls="-", alpha=0.8, zorder=2)

                ax.fill_between(t_lat, Dlo, Dhi,
                                color="tomato", alpha=0.20, zorder=1,
                                label=r"$D_s^{\rm lat}\ 1\sigma$ (Hackett et al. 2024)")
                ax.plot(t_lat, Dc, color="tomato", lw=1.5, ls="--", alpha=0.8, zorder=2)

            # CLAS12 HS-extracted curves
            ax.plot(t_abs, As, lw=2.2, color="C0", zorder=4,
                    label=r"$A_s(t)$ CLAS12 (dipole)")
            ax.plot(t_abs, Ds, lw=2.2, ls="--", color="C1", zorder=4,
                    label=r"$D_s(t)$ CLAS12 (tripole)")

            ax.axhline(0.0, color="gray", lw=0.7, ls=":", zorder=0)

            ax.set_xlabel(r"$|t|\ [\mathrm{GeV}^2]$")
            ax.set_ylabel(r"Gravitational form factors")
            ax.set_title(title, fontsize=LEG_TITLE_FS)
            ax.legend(frameon=False, fontsize=LEG_FS - 2, loc="lower right")

            #if show_lattice:
                #ax.text(0.98, 0.2,
                 #       "Lattice: Hackett, Pefkou, Shanahan\nPRL 132, 251904 (2024), Table I",
                  #      transform=ax.transAxes, ha="right", va="bottom",
                   #     fontsize=max(LEG_FS - 4, 8), style="italic", color="gray")

            stem = stem_base + ("_with_lattice" if show_lattice else "")
            _save(fig, outdir, stem)

def xb_from_Q2_W(Q2, W, M=M_P_GEV):
    W2 = W * W
    return Q2 / (W2 - M*M + Q2)


def t_min_phi(Q2, W, M=M_P_GEV, m_phi=M_PHI_GEV):
    """Compute t_min for *p  p, i.e. the minimum kinematically allowed |t|.

    In the *p centre-of-mass frame the forward-scattering limit gives:

        t_min = ( (W + Q - m_)/(2W) - (W - M)/(2W) )    q_cm

    More explicitly, using 4-momentum conservation:

        t_min = [ (m_ + Q) / (2W)  (W  M)/(2W) ]    [q_*  q_]

    The standard compact formula (see e.g. Diehl 2003, or CLAS  papers) is:

        t_min = [ (m_ + Q) / (2W) ]  + [q_*  q__cm]   ... (signs)

    We use the fully explicit 4-vector approach to be unambiguous:

    In the *p CM frame:
        E_*  = (W  M  Q) / (2W)     [virtual photon energy, can be negative for large Q]
        |p_*| = sqrt(E_* + Q)           [3-momentum magnitude of virtual photon]
        E_   = (W + m_  M) / (2W)   [ energy at forward scattering]
        |p_|  = sqrt(E_  m_)

        t_min  = (p_*  p_)  [4-vector, evaluated at _cm = 0]
               = m_  Q  2(E_*  E_  |p_*|  |p_|)

    Returns
    -------
    float  : |t_min|  0  (we return the absolute value, consistent with t' = |t|  |t_min|)
             Returns np.nan if kinematics are forbidden (W < M + m_phi).
    """
    Q2  = float(Q2);  W = float(W)
    if W <= M + m_phi or Q2 < 0:
        return np.nan

    W2   = W * W
    mf2  = m_phi * m_phi
    M2   = M * M

    # virtual photon in CM frame
    E_gstar = (W2 - M2 - Q2) / (2.0 * W)
    p_gstar = np.sqrt(max(E_gstar**2 + Q2, 0.0))   # |p| = sqrt(E+Q) for spacelike

    #  meson in CM frame (forward, =0)
    E_phi_cm = (W2 + mf2 - M2) / (2.0 * W)
    p_phi_cm2 = E_phi_cm**2 - mf2
    if p_phi_cm2 < 0:
        return np.nan
    p_phi_cm = np.sqrt(p_phi_cm2)

    # t_min = (q_*  q_) at _cm = 0
    # = (E_*  E_)  (p_*  p_)   [q = Ep with metric +]
    # note: this is negative (spacelike), we return |t_min|
    t_min_val = (E_gstar - E_phi_cm)**2 - (p_gstar - p_phi_cm)**2
    return abs(float(t_min_val))


def compute_tmin_for_df(df):
    """Attach |t_min| column to df, computed row-by-row from Q2_center, W_center.

    Falls back to 0.0 if kinematics are not available (safe: then t  t').
    """
    df = df.copy()
    if "Q2_center" not in df.columns or "W_center" not in df.columns:
        df["tmin_abs"] = 0.0
        return df

    Q2 = df["Q2_center"].to_numpy(float)
    W  = df["W_center"].to_numpy(float)
    tmin = np.array([t_min_phi(q, w) for q, w in zip(Q2, W)], dtype=float)
    # fill NaN with 0 (safe fallback)
    tmin = np.where(np.isfinite(tmin), tmin, 0.0)
    df["tmin_abs"] = tmin
    return df

def hand_virtual_photon_flux(E, Q2, W, M=M_P_GEV, alpha=ALPHA_EM):
    """Hand convention virtual photon flux _v."""
    Q2 = float(Q2); W = float(W); E = float(E)
    if Q2 <= 0 or E <= 0 or W <= 0:
        return (float("nan"),) * 6

    W2  = W * W
    nu  = (W2 + Q2 - M*M) / (2*M)
    y   = nu / E
    Ep  = E - nu
    if Ep <= 0:
        return (float("nan"),) * 6

    K = (W2 - M*M) / (2*M)

    eps_num = 1 - y - Q2 / (4*E*E)
    eps_den = 1 - y + (y*y)/2 + Q2 / (4*E*E)
    eps = eps_num / eps_den if eps_den != 0 else float("nan")

    if (1 - eps) == 0 or Q2 == 0:
        Gamma_v = float("nan")
    else:
        Gamma_v = (alpha / (2*(np.pi**2))) * (Ep/E) * (K/Q2) * (1/(1 - eps))

    return Gamma_v, eps, y, nu, Ep, K

# -----------------------------------------------------------------------------
#  External RadCorr table loader + attacher (DIFFRAD)
# -----------------------------------------------------------------------------
def load_external_radcorr_table(path):
    """Load DIFFRAD CSV  standardised DataFrame with bin edges + crad."""
    tab = pd.read_csv(path, comment="#")
    tab.columns = [c.strip() for c in tab.columns]

    norm = {c: re.sub(r"[^a-z0-9]+", "", c.lower()) for c in tab.columns}

    def pick_col(*candidates):
        for cand in candidates:
            for col, nrm in norm.items():
                if nrm == cand:
                    return col
        for cand in candidates:
            for col, nrm in norm.items():
                if cand in nrm:
                    return col
        return None

    c_q2lo = pick_col("q2lo",  "q2_lo")
    c_q2hi = pick_col("q2hi",  "q2_hi")
    c_q2c  = pick_col("q2c",   "q2_c",  "q2center", "q2")
    c_tlo  = pick_col("tprimelo",  "tprime_lo",  "tlo")
    c_thi  = pick_col("tprimehi",  "tprime_hi",  "thi")
    c_tc   = pick_col("tprimec",   "tprime_c",   "tprimecenter", "tc")
    c_rc   = pick_col("radcorr",   "rad_corr",   "crad", "rc",
                      "radiativecorrection", "c_rad")

    if c_rc is None:
        raise ValueError(
            f"RadCorr table: cannot find rad_corr column. "
            f"Available columns: {list(tab.columns)}"
        )

    def to_f(col):
        return tab[col].astype(float) if col and col in tab.columns else None

    q2lo = to_f(c_q2lo)
    q2hi = to_f(c_q2hi)
    q2c  = to_f(c_q2c)
    tlo  = to_f(c_tlo)
    thi  = to_f(c_thi)
    tc   = to_f(c_tc)
    rc   = tab[c_rc].astype(float)

    out = pd.DataFrame({"crad": rc})
    if q2lo is not None: out["q2_lo"] = q2lo
    if q2hi is not None: out["q2_hi"] = q2hi
    if q2c  is not None: out["q2_c"]  = q2c
    elif q2lo is not None and q2hi is not None:
        out["q2_c"] = 0.5 * (q2lo + q2hi)
    if tlo  is not None: out["t_lo"]  = tlo
    if thi  is not None: out["t_hi"]  = thi
    if tc   is not None: out["t_c"]   = tc
    elif tlo is not None and thi is not None:
        out["t_c"] = 0.5 * (tlo + thi)

    out = out.replace([np.inf, -np.inf], np.nan).dropna(subset=["crad"])
    return out


def choose_radcorr_file_for_energy(E):
    if E is None:
        return None
    if abs(E - 10.6) <= abs(E - 10.2):
        return K_RADCORR_FILE_10P6
    return K_RADCORR_FILE_10P2


def _safe_col(df, col):
    return df[col].to_numpy(float) if col in df.columns else np.full(len(df), np.nan)


def attach_external_radcorr(df, E, path_override=None, mode="nearest"):
    path = path_override or choose_radcorr_file_for_energy(E)
    if path is None:
        return df

    try:
        tab = load_external_radcorr_table(path)
    except Exception as e:
        print(f"  [WARN] Could not load external RadCorr table '{path}': {e}")
        return df

    if "Q2_center" in df.columns:
        q2 = df["Q2_center"].to_numpy(float)
    elif "Q2_lo" in df.columns and "Q2_hi" in df.columns:
        q2 = 0.5 * (df["Q2_lo"].to_numpy(float) + df["Q2_hi"].to_numpy(float))
    else:
        q2 = np.full(len(df), np.nan)

    if "tprime_center" in df.columns:
        t = df["tprime_center"].to_numpy(float)
    elif "tprime_lo" in df.columns and "tprime_hi" in df.columns:
        t = 0.5 * (df["tprime_lo"].to_numpy(float) + df["tprime_hi"].to_numpy(float))
    else:
        t = np.full(len(df), np.nan)

    has_edges = ("q2_lo" in tab.columns and "q2_hi" in tab.columns
                 and "t_lo"  in tab.columns and "t_hi"  in tab.columns)
    has_center = "q2_c" in tab.columns and "t_c" in tab.columns

    tab_rc = tab["crad"].to_numpy(float)
    rc_out = np.full(len(df), np.nan)

    if has_edges:
        tab_q2lo = tab["q2_lo"].to_numpy(float)
        tab_q2hi = tab["q2_hi"].to_numpy(float)
        tab_tlo  = tab["t_lo"].to_numpy(float)
        tab_thi  = tab["t_hi"].to_numpy(float)

        for i in range(len(df)):
            if not (np.isfinite(q2[i]) and np.isfinite(t[i])):
                continue
            in_q2 = (q2[i] >= tab_q2lo) & (q2[i] < tab_q2hi)
            in_t  = (t[i]  >= tab_tlo)  & (t[i]  < tab_thi)
            match = np.where(in_q2 & in_t)[0]
            if len(match) > 0:
                rc_out[i] = tab_rc[match[0]]
            elif has_center:
                tab_q2c = tab["q2_c"].to_numpy(float)
                tab_tc  = tab["t_c"].to_numpy(float)
                dq2 = (tab_q2c - q2[i]) / np.maximum(tab_q2hi - tab_q2lo, 1e-9)
                dt  = (tab_tc  - t[i])  / np.maximum(tab_thi  - tab_tlo,  1e-9)
                j   = int(np.nanargmin(dq2**2 + dt**2))
                rc_out[i] = tab_rc[j]
                print(f"  [WARN] No exact RadCorr bin for Q2={q2[i]:.3f}, t={t[i]:.3f}"
                      f"  nearest centroid Q2={tab_q2c[j]:.3f}, t={tab_tc[j]:.3f}"
                      f", C_rad={tab_rc[j]:.4f}")

    elif has_center:
        tab_q2c = tab["q2_c"].to_numpy(float)
        tab_tc  = tab["t_c"].to_numpy(float)
        for i in range(len(df)):
            if not (np.isfinite(q2[i]) and np.isfinite(t[i])):
                continue
            dist2 = (tab_q2c - q2[i])**2 + (tab_tc - t[i])**2
            j = int(np.nanargmin(dist2))
            rc_out[i] = tab_rc[j]
    else:
        print("  [WARN] RadCorr table has neither bin edges nor centroid columns  skipping.")
        return df

    if "RadCorr" in df.columns:
        rc_existing = df["RadCorr"].to_numpy(float)
        rc_out = np.where(
            np.isfinite(rc_out) & (rc_out > 0), rc_out,
            np.where(np.isfinite(rc_existing) & (rc_existing > 0), rc_existing, np.nan)
        )

    df = df.copy()
    df["RadCorr"] = rc_out
    return df

def _attach_gamma_from_csv_columns(df):
    """If the CSV already contains a virtual photon flux column, standardize it.

    Accepts several common column names (depending on which macro produced the
    CSV). If found, the value is copied into a canonical column named
    `Gamma_v`.
    """
    gamma_cols = [
        "Gamma_v",
        "GammaV",
        "Gamma",
        "GammaV_mean",
        "Gamma_v_mean",
        "GammaV_center",
        "Gamma_v_center",
        "mean_GammaV",
        "GammaV_c",
        "Gammav",
    ]

    for c in gamma_cols:
        if c in df.columns and not df[c].isna().all():
            out = df.copy()
            out["Gamma_v"] = out[c].astype(float)
            return out
    return df


def attach_virtual_photon_flux_and_reduced_xs(df, E, prefer_csv_gamma=True):
    """Attach virtual photon flux (Gamma_v) and reduced cross-section.

    - If `prefer_csv_gamma=True` and the CSV already has a Gamma-like column,
      it will be used.
    - Otherwise (or if missing), Gamma_v is computed from (E, Q2_center, W_center).
    """
    df_in = df.copy()

    # 1) Prefer Gamma from CSV if provided
    if prefer_csv_gamma:
        df_in = _attach_gamma_from_csv_columns(df_in)

    # 2) Compute Gamma if still missing/invalid
    need_compute = (
        ("Gamma_v" not in df_in.columns)
        or df_in["Gamma_v"].isna().all()
        or (not np.isfinite(df_in["Gamma_v"].to_numpy(float)).any())
    )

    if need_compute:
        if E is None:
            return df_in
        if "Q2_center" not in df_in.columns:
            return df_in
        if "W_center" not in df_in.columns:
            return df_in

        Q2 = df_in["Q2_center"].to_numpy(float)
        W  = df_in["W_center"].to_numpy(float)

        Gamma = np.full(len(df_in), np.nan)
        eps   = np.full(len(df_in), np.nan)
        y_arr = np.full(len(df_in), np.nan)

        for i in range(len(df_in)):
            if not (np.isfinite(Q2[i]) and np.isfinite(W[i])):
                continue
            g, e, yy, nu, Ep, K = hand_virtual_photon_flux(E, Q2[i], W[i])
            Gamma[i] = g
            eps[i]   = e
            y_arr[i] = yy

        df_in["Gamma_v"] = Gamma
        df_in["epsilon"] = eps
        df_in["y"]       = y_arr

    # 3) Reduced cross-section
    if "CrossSection" in df_in.columns and "Gamma_v" in df_in.columns:
        xs   = df_in["CrossSection"].to_numpy(float)
        xerr = (
            df_in["CrossSection_Err"].to_numpy(float)
            if "CrossSection_Err" in df_in.columns else np.zeros_like(xs)
        )
        Gamma = df_in["Gamma_v"].to_numpy(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            df_in["ReducedCrossSection"]     = np.where(Gamma > 0, xs   / Gamma, np.nan)
            df_in["ReducedCrossSection_Err"] = np.where(Gamma > 0, xerr / Gamma, np.nan)

    return df_in

# -----------------------------------------------------------------------------
#  Derive center columns helper  (call before any physics attachment)
# -----------------------------------------------------------------------------
def _derive_centers(df):
    """Add *_center columns from lo/hi pairs if not already present."""
    df = df.copy()
    for base in ("Q2", "W", "tprime"):
        clo  = f"{base}_lo"
        chi  = f"{base}_hi"
        ccen = f"{base}_center"
        if ccen not in df.columns and clo in df.columns and chi in df.columns:
            df[ccen] = 0.5 * (df[clo].astype(float) + df[chi].astype(float))
    return df

def ensure_q2_w_centers(df, iq=None, w_mean=W_MEAN_GEV):
    """
    Ensure df has Q2_center and W_center so Gamma_v can be computed even when
    CSVs omit W/Q2 because they are implicit by bin/file.

    - W_center: if missing, fill with hardcoded W_MEAN_GEV
    - Q2_center: if missing, fill from Q2_MEAN_BY_IQ_GEV2 using iq
    """
    df = df.copy()

    # ---- W_center ----
    if "W_center" not in df.columns or df["W_center"].isna().all():
        for cand in ("W", "Wavg", "W_mean", "Wc", "W_c", "Wcenter"):
            if cand in df.columns and not df[cand].isna().all():
                df["W_center"] = df[cand].astype(float)
                break
        else:
            df["W_center"] = float(w_mean)

    # ---- Q2_center ----
    if "Q2_center" not in df.columns or df["Q2_center"].isna().all():
        for cand in ("Q2", "Q2avg", "Q2_mean", "Q2c", "Q2_c", "Q2center"):
            if cand in df.columns and not df[cand].isna().all():
                df["Q2_center"] = df[cand].astype(float)
                break
        else:
            if iq is not None and int(iq) in Q2_MEAN_BY_IQ_GEV2:
                df["Q2_center"] = float(Q2_MEAN_BY_IQ_GEV2[int(iq)])
            else:
                df["Q2_center"] = float(list(Q2_MEAN_BY_IQ_GEV2.values())[0])

    return df

# -----------------------------------------------------------------------------
#  Table dump
# -----------------------------------------------------------------------------
def dump_xs_table(df, variants, label, outdir):
    raw_xs, raw_err = variants["raw"]
    acc_xs, acc_err = variants["acc"]
    eff_xs, eff_err = variants["eff"]
    fin_xs, fin_err = variants["final"]

    n = len(df)

    def col_or_nan(name):
        return df[name].to_numpy(float) if name in df.columns else np.full(n, np.nan)

    tcent  = col_or_nan("tprime_center")
    A      = col_or_nan("Acceptance")
    Eff    = col_or_nan("Efficiency")
    Crad   = col_or_nan("RadCorr")
    Gamma  = col_or_nan("Gamma_v")

    with np.errstate(invalid="ignore", divide="ignore"):
        reduced = np.where(np.isfinite(Gamma) & (Gamma > 0), fin_xs / Gamma, np.nan)

    tab = pd.DataFrame({
        "tprime_center":         tcent,
        "raw_xs":                raw_xs,
        "Gamma_v":               Gamma,
        "reduced_xs":            reduced,
        "Acceptance":            A,
        "Efficiency":            Eff,
        "xs_acc_corrected":      acc_xs,
        "xs_eff_corrected":      eff_xs,
        "RadCorr":               Crad,
        "xs_rad_corrected":      fin_xs,
        "raw_xs_err":            raw_err,
        "xs_acc_corrected_err":  acc_err,
        "xs_eff_corrected_err":  eff_err,
        "xs_rad_corrected_err":  fin_err,
    })

    print(f"\n=== XS TABLE: {label} ===")
    print(tab.to_string(index=False))

    os.makedirs(outdir, exist_ok=True)
    tab.to_csv(os.path.join(outdir, f"table_{label}.csv"), index=False)

# -----------------------------------------------------------------------------
#  Style constants
# -----------------------------------------------------------------------------
MODEL_COLORS = ["#1F4FD8", "#E07010", "#009999", "#008800", "#7733CC", "#CC1A33", "#666666"]

CORR_COLORS = {"raw": "#000000", "acc": "#1F4FD8", "eff": "#7A3DB8", "final": "#009900"}
CORR_LABELS = {
    "raw":   "Raw  (no corr.)",
    "acc":   "Acceptance corr. only",
    "eff":   "Acc. + Eff. corr.",
    "final": "Acc. + Eff. + Rad. (final)",
}

PANEL_RC = {
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 22,
    "axes.labelsize": 22,
    "axes.titlesize": 22,
    "legend.fontsize": 20,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "xtick.major.size": 12,
    "xtick.major.width": 2,
    "ytick.major.size": 12,
    "ytick.major.width": 2,
    "xtick.minor.size": 6,
    "xtick.minor.width": 2,
    "ytick.minor.size": 6,
    "ytick.minor.width": 2,
    "axes.linewidth": 2,
    "lines.linewidth": 2,
}

LEG_FS       = 18
LEG_TITLE_FS = 18
YLIM_XS      = (1e-5, 80.0)
YLIM_XS_RAW      = (1e-8, 1e-4)

XLIM_T       = (0.0, 6.0)
YLIM_RED     = (1e-5, 50.0)  # reduced cross-section comparison y-limits
NCOLS_GRP    = 2

_CANVAS_CFG = {
    "xs": dict(
        col="CrossSection", err_col="CrossSection_Err",
        #ylabel=r"$\mathrm{d}\sigma/\mathrm{d}t'\ [\mathrm{nb/GeV}^2]$",
        ylabel=r"Yields (arb. units)",
        logy=True, tag="dsdt", title="Cross-section (final)"),
    "red": dict(
        col="ReducedCrossSection", err_col="ReducedCrossSection_Err",
        #ylabel=r"$\sigma_\mathrm{red} = (\mathrm{d}\sigma/\mathrm{d}t')/\Gamma_v\ [\mathrm{nb}]$",
        ylabel=r"Yields (arb. units)",
        logy=True, tag="dsdt_reduced", title="Reduced cross-section"),
    "acc": dict(
        col="Acceptance", err_col=None,
        ylabel=r"Acceptance $A(\varepsilon)$",
        logy=False, tag="acceptance", title="Acceptance"),
    "eff": dict(
        col="Efficiency", err_col=None,
        ylabel=r"Efficiency $\varepsilon$",
        logy=False, tag="efficiency", title="Efficiency"),
    "rad": dict(
        col="RadCorr", err_col=None,
        ylabel=r"$C_\mathrm{rad}$",
        logy=False, tag="radcorr", title="Radiative Correction"),
}

def _sanitise(s):
    return re.sub(r"[ /\\]", "_", s)

def _style_ax(ax):
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))

def _bin_label(df):
    q2lo = float(df["Q2_lo"].iloc[0]) if "Q2_lo" in df.columns else float("nan")
    q2hi = float(df["Q2_hi"].iloc[0]) if "Q2_hi" in df.columns else float("nan")
    lbl  = rf"$Q^2 \in [{q2lo:.2f},\,{q2hi:.2f}]\ \mathrm{{GeV}}^2$"
    if "W_lo" in df.columns and "W_hi" in df.columns and pd.notna(df["W_lo"].iloc[0]):
        wlo = float(df["W_lo"].iloc[0])
        whi = float(df["W_hi"].iloc[0])
        #lbl += rf"$\quad W \in [{wlo:.2f},\,{whi:.2f}]\ \mathrm{{GeV}}$"
    return lbl

def _save(fig, outdir, stem):
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"{stem}.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"    [OK] {stem}")

# -----------------------------------------------------------------------------
#  CSV loader  (per-model beam energy + center-column derivation + W/Q2 forcing)
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
#  Formatting helpers (avoid MathText parse issues in PDF backend)
# -----------------------------------------------------------------------------
def _fmt_sci(x, sig=3):
    """Format number with optional scientific notation using '×10^' (plain text)."""
    if x is None or not np.isfinite(x):
        return "nan"
    s = f"{float(x):.{sig}g}"
    if "e" in s or "E" in s:
        mant, exp = re.split(r"[eE]", s)
        exp = int(exp)
        return f"{mant}×10^{exp}"
    return s

def _fmt_pm(val, err, unit="", vsig=3, esig=2):
    """Format 'val ± err unit' as plain text (no $...$) to keep PDF mathtext happy."""
    v = _fmt_sci(val, sig=vsig)
    e = _fmt_sci(err, sig=esig)
    if unit:
        return f"{v} ± {e} {unit}"
    return f"{v} ± {e}"
def load_csv(csv_root, model, iq, iw=None):
    safe  = _sanitise(model)
    fname = os.path.join(
        csv_root, "CSVs", safe,
        f"dsdt_Q{iq}.csv" if iw is None else f"dsdt_Q{iq}_W{iw}.csv"
    )
    if not os.path.isfile(fname):
        return None

    df = pd.read_csv(fname)

    # -- derive center columns from lo/hi if possible ----------------------
    df = _derive_centers(df)

    # -- NEW: force W_center and Q2_center even if CSV omitted them --------
    df = ensure_q2_w_centers(df, iq=iq, w_mean=W_MEAN_GEV)

    # -- per-model beam energy --------------------------------------------
    E = beam_energy_for_model(model)

    if G_USE_EXTERNAL_RADCORR:
        df = attach_external_radcorr(
            df, E,
            path_override=G_RADCORR_FILE_OVERRIDE,
            mode=G_RADCORR_MATCH_MODE,
        )

    if E is not None:
        df = attach_virtual_photon_flux_and_reduced_xs(df, E)

    # -- attach |t_min| per row from *pp kinematics --------------------
    df = compute_tmin_for_df(df)

    return df

def discover_bin_keys(csv_root, models):
    keys = set()
    for m in models:
        for f in glob(os.path.join(csv_root, "CSVs", _sanitise(m), "dsdt_Q*.csv")):
            base  = os.path.basename(f).replace("dsdt_", "").replace(".csv", "")
            parts = base.split("_")
            iq    = int(parts[0][1:])
            iw    = int(parts[1][1:]) if len(parts) > 1 else None
            keys.add((iq, iw))
    return sorted(keys)

def _iw_unique(bin_keys):
    return sorted({iw for (_, iw) in bin_keys}, key=lambda x: (x is None, x))

# -----------------------------------------------------------------------------
#  Build correction variants
# -----------------------------------------------------------------------------
def _build_variants(df, luminosity=None, branching=None):
    xs   = _safe_col(df, "CrossSection")
    xerr = _safe_col(df, "CrossSection_Err")
    A    = _safe_col(df, "Acceptance")
    Ceff = _safe_col(df, "Efficiency")
    Crad = _safe_col(df, "RadCorr")
    Nsig = _safe_col(df, "RawCounts")
    Nerr = _safe_col(df, "RawCounts_Err")
    dT   = df["tprime_hi"].to_numpy(float) - df["tprime_lo"].to_numpy(float)

    A_s    = np.where(np.isfinite(A)    & (A    > 0), A,    1.0)
    # Clamp efficiency to (0, 1]: values >1 or <=0 are unphysical (low MC stats)
    # and are replaced with 1.0 so they don't corrupt xs_eff_corrected.
    Ceff_valid = np.isfinite(Ceff) & (Ceff > 0) & (Ceff <= 1.0)
    Ceff_s = np.where(Ceff_valid, Ceff, 1.0)
    Crad_s = np.where(np.isfinite(Crad) & (Crad > 0), Crad, 1.0)
    dT_s   = np.where(dT > 0, dT, 1.0)

    have_nsig = np.any(np.isfinite(Nsig) & (Nsig > 0))

    if have_nsig and luminosity is not None and branching is not None:
        L  = luminosity
        BR = branching
        denom_raw   = L * BR * dT_s
        denom_acc   = L * BR * dT_s * A_s
        denom_eff   = L * BR * dT_s * A_s * Ceff_s
        denom_final = L * BR * dT_s * A_s * Ceff_s * Crad_s

        raw_v   = np.where(denom_raw   > 0, Nsig / denom_raw,   np.nan)
        acc_v   = np.where(denom_acc   > 0, Nsig / denom_acc,   np.nan)
        eff_v   = np.where(denom_eff   > 0, Nsig / denom_eff,   np.nan)
        final_v = np.where(denom_final > 0, Nsig / denom_final, np.nan)
        raw_e   = np.where(denom_raw   > 0, Nerr / denom_raw,   np.nan)
        acc_e   = np.where(denom_acc   > 0, Nerr / denom_acc,   np.nan)
        eff_e   = np.where(denom_eff   > 0, Nerr / denom_eff,   np.nan)
        final_e = np.where(denom_final > 0, Nerr / denom_final, np.nan)
    else:
        # CrossSection already fully corrected — reconstruct intermediate variants
        final_v = xs.copy()
        final_e = xerr.copy()
        eff_v   = xs   * Crad_s                 # undo RadCorr only
        eff_e   = xerr * Crad_s
        acc_v   = xs   * Ceff_s * Crad_s       # undo Eff and RadCorr
        acc_e   = xerr * Ceff_s * Crad_s
        raw_v   = xs   * A_s * Ceff_s * Crad_s # undo Acceptance, Eff, RadCorr
        raw_e   = xerr * A_s * Ceff_s * Crad_s

    return {"raw": (raw_v, raw_e), "acc": (acc_v, acc_e), "eff": (eff_v, eff_e), "final": (final_v, final_e)}

# -----------------------------------------------------------------------------
#  Corrections canvas (before/after)
# -----------------------------------------------------------------------------
def make_corrections_canvas(csv_root, models, iq, iw, outdir, logy=True,
                             luminosity=None, branching=None,
                             print_table=False, table_outdir=None):
    dfs = {m: df for m in models if (df := load_csv(csv_root, m, iq, iw)) is not None}
    if not dfs:
        return

    n = len(dfs)
    with mpl.rc_context(PANEL_RC):
        fig   = plt.figure(figsize=(10 * n, 12))
        outer = mpl.gridspec.GridSpec(
            1, n, figure=fig,
            left=0.07, right=0.99, top=0.91, bottom=0.09, wspace=0.06,
        )
        main_axes, ratio_axes = [], []
        for i in range(n):
            inner = mpl.gridspec.GridSpecFromSubplotSpec(
                2, 1, subplot_spec=outer[i], height_ratios=[3, 1], hspace=0.0,
            )
            main_axes.append(fig.add_subplot(inner[0]))
            ratio_axes.append(fig.add_subplot(inner[1], sharex=main_axes[-1]))

        bin_header  = None
        glob_handles = glob_labels = None

        for im, (m, df) in enumerate(dfs.items()):
            ax_m = main_axes[im]
            ax_r = ratio_axes[im]
            _style_ax(ax_m); _style_ax(ax_r)

            t   = df["tprime_center"].to_numpy(float)
            xlo = t - df["tprime_lo"].to_numpy(float)
            xhi = df["tprime_hi"].to_numpy(float) - t

            variants = _build_variants(df, luminosity, branching)

            if print_table:
                tout  = table_outdir or os.path.join(os.path.dirname(outdir), "0_Tables")
                label = f"{m}_Q{iq}" + (f"_W{iw}" if iw is not None else "")
                dump_xs_table(df, variants, label, tout)

            loc_handles, loc_labels = [], []
            for vkey, (val, err) in variants.items():
                mask = np.isfinite(val) & (val > 0)
                if not mask.any():
                    continue
                colr = CORR_COLORS[vkey]
                h = ax_m.errorbar(
                    t[mask], val[mask], yerr=err[mask],
                    xerr=(xlo[mask], xhi[mask]),
                    fmt="o", color=colr, ecolor=colr, elinewidth=1.5, capsize=3,
                    markersize=7, mfc="white", mew=2.0, lw=1.8, zorder=5,
                    label=CORR_LABELS[vkey],
                )
                loc_handles.append(h)
                loc_labels.append(CORR_LABELS[vkey])

            final_v, final_e = variants["final"]
            ref_safe = np.where(np.isfinite(final_v) & (final_v > 0), final_v, np.nan)

            for vkey in ("raw", "acc", "eff"):
                val, err = variants[vkey]
                denom    = np.where(np.isfinite(val) & (val != 0), val, np.nan)
                ratio    = ref_safe / denom
                rmask    = np.isfinite(ratio) & (ratio > 0)
                if not rmask.any():
                    continue
                with np.errstate(invalid="ignore"):
                    ratio_err = ratio[rmask] * np.sqrt(
                        (final_e[rmask] / ref_safe[rmask])**2 +
                        (err[rmask]     / denom[rmask]   )**2
                    )
                ax_r.errorbar(
                    t[rmask], ratio[rmask], yerr=ratio_err,
                    xerr=(xlo[rmask], xhi[rmask]),
                    fmt="o", color=CORR_COLORS[vkey], ecolor=CORR_COLORS[vkey],
                    elinewidth=1.5, capsize=3,
                    markersize=6, mfc="white", mew=1.8, lw=1.5, zorder=5,
                )

            ax_r.axhline(1.0, color="gray", lw=1.4, ls="--", zorder=1)
            ax_r.set_yscale("log")
            ax_r.set_ylim(0.3, 500.0)
            ax_r.set_ylabel("Final / Variant", fontsize=16)
            ax_r.set_xlabel(r"$-t'\ [\mathrm{GeV}^2]$")
            ax_r.tick_params(axis="x", labelsize=18)
            ax_r.tick_params(axis="y", labelsize=16)

            ax_m.set_xlim(XLIM_T); ax_r.set_xlim(XLIM_T)
            ax_m.tick_params(labelbottom=False)
            if logy:
                ax_m.set_yscale("log")
                ax_m.set_ylim(YLIM_XS)
            ax_m.set_ylabel(_CANVAS_CFG["xs"]["ylabel"] if im == 0 else "")
            if im != 0:
                ax_m.tick_params(labelleft=False)
                ax_r.tick_params(labelleft=False)

            E_val = beam_energy_for_model(m)
            E_str = f"  (E = {E_val:.3g} GeV)" if E_val is not None else ""
            ax_m.annotate(
                m.replace("_", " ") + E_str,
                xy=(0.97, 0.97), xycoords="axes fraction",
                fontsize=LEG_FS, ha="right", va="top", color="black",
            )

            if bin_header is None:
                bin_header = _bin_label(df)
            if glob_handles is None:
                glob_handles = loc_handles
                glob_labels  = loc_labels

        if bin_header:
            fig.suptitle(bin_header, fontsize=LEG_TITLE_FS + 2, y=0.965)
        if glob_handles:
            main_axes[0].legend(
                glob_handles, glob_labels,
                loc="upper right", frameon=False, fontsize=LEG_FS - 2, ncol=1,
            )

        stem = f"corrections_Q{iq}" if iw is None else f"corrections_Q{iq}_W{iw}"
        _save(fig, outdir, stem)

# -----------------------------------------------------------------------------
#  Group / comparison plots
# -----------------------------------------------------------------------------
def _draw_panel(ax, models, csv_root, iq, iw, canvas_key, logy):
    cfg = _CANVAS_CFG[canvas_key]
    col = cfg["col"]
    ec  = cfg["err_col"]
    handles, labels = [], []
    header = None
    any_d  = False

    for im, m in enumerate(models):
        df = load_csv(csv_root, m, iq, iw)
        if df is None:
            continue
        if header is None:
            header = _bin_label(df)
        x    = df["tprime_center"].to_numpy(float)
        y    = _safe_col(df, col)
        # For efficiency, only plot physically valid values (0 < eff <= 1)
        if canvas_key == "eff":
            mask = np.isfinite(y) & (y > 0) & (y <= 1.0)
        else:
            mask = np.isfinite(y) & (y > 0)
        if not mask.any():
            continue
        color = MODEL_COLORS[im % len(MODEL_COLORS)]
        yerr  = _safe_col(df, ec) if ec else None
        h = ax.errorbar(
            x[mask], y[mask],
            yerr=yerr[mask] if yerr is not None else None,
            fmt="o", color=color, ecolor=color, elinewidth=1.5, capsize=3,
            markersize=7, mfc="white", mew=2.0, lw=1.8, zorder=5 + im,
            label=m.replace("_", " "),
        )
        handles.append(h); labels.append(m.replace("_", " "))
        any_d = True

    if cfg["logy"] and logy and any_d:
        ax.set_yscale("log")
    return handles, labels, header, any_d


def make_group_plot(csv_root, models, bin_keys, canvas_key, outdir, logy=True):
    cfg    = _CANVAS_CFG[canvas_key]
    tag    = cfg["tag"]
    ylabel = cfg["ylabel"]

    for iw_key in _iw_unique(bin_keys):
        iq_list = sorted(iq for (iq, iw) in bin_keys if iw == iw_key)
        if not iq_list:
            continue
        nQ    = len(iq_list)
        ncols = NCOLS_GRP + 1
        nrows = (nQ + NCOLS_GRP - 1) // NCOLS_GRP

        with mpl.rc_context(PANEL_RC):
            fig, axes_arr = plt.subplots(
                nrows, ncols, figsize=(10*ncols, 8*nrows), squeeze=False,
            )
            fig.subplots_adjust(left=0.07, right=0.99, top=0.93, bottom=0.09,
                                wspace=0.10, hspace=0.22)
            axes_flat = axes_arr.ravel()
            for ax in axes_flat:
                ax.axis("off")

            common_h = common_l = None
            w_label  = ""
            last_idx = -1

            for di, iq in enumerate(iq_list):
                slot = (di // NCOLS_GRP) * ncols + (di % NCOLS_GRP)
                ax   = axes_flat[slot]
                ax.axis("on"); _style_ax(ax)

                handles, labels_l, hdr, has_d = _draw_panel(
                    ax, models, csv_root, iq, iw_key, canvas_key, logy,
                )
                if not has_d:
                    ax.axis("off"); continue

                last_idx = di
                if common_h is None:
                    common_h = handles; common_l = labels_l

                ax.set_xlim(XLIM_T)
                if canvas_key == "xs":
                    ax.set_ylim(YLIM_XS)
                elif canvas_key == "red":
                    ax.set_ylim(YLIM_RED)
                elif canvas_key == "acc":
                    ax.set_ylim(0.0, 0.035)   # <-- requested
                elif canvas_key == "eff":
                    ax.set_ylim(0.6, 1.05)
                    ax.axhline(1.0, color="gray", lw=1.2, ls="--", zorder=1)
                elif canvas_key == "rad":
                    ax.set_ylim(0.6, 1.0)
                    ax.axhline(1.0, color="gray", lw=1.2, ls="--", zorder=1)
                # --- FORCE reduced comparison ymax to 30 (works even if --no-logy) ---

                df_ref = next(
                    (load_csv(csv_root, m, iq, iw_key) for m in models
                     if load_csv(csv_root, m, iq, iw_key) is not None), None
                )
                if df_ref is not None and "Q2_lo" in df_ref.columns and "Q2_hi" in df_ref.columns:
                    q2lo = float(df_ref["Q2_lo"].iloc[0])
                    q2hi = float(df_ref["Q2_hi"].iloc[0])
                    ax.annotate(
                        rf"$Q^2 \in [{q2lo:.2f},\,{q2hi:.2f}]\;\mathrm{{GeV}}^2$",
                        xy=(0.05, 0.97), xycoords="axes fraction",
                        fontsize=LEG_FS, ha="left", va="top",
                    )
                    if "W_lo" in df_ref.columns and pd.notna(df_ref["W_lo"].iloc[0]):
                        wlo = float(df_ref["W_lo"].iloc[0])
                        whi = float(df_ref["W_hi"].iloc[0])
                        w_label = rf"$W \in [{wlo:.2f},\,{whi:.2f}]\;\mathrm{{GeV}}$"

                cin = di % NCOLS_GRP
                ax.set_xlabel(r"$-t'\,[\mathrm{GeV}^2]$")
                ax.set_ylabel(ylabel if cin == 0 else "")
                if cin != 0:
                    ax.tick_params(labelleft=False)

            rows_used = (last_idx // NCOLS_GRP) + 1 if last_idx >= 0 else 0
            for row in range(rows_used):
                leg_slot = row * ncols + NCOLS_GRP
                if leg_slot < len(axes_flat) and common_h:
                    la = axes_flat[leg_slot]; la.axis("off")
                    la.legend(
                        common_h, common_l,
                        loc="center", frameon=False, fontsize=LEG_FS,
                        title=w_label or None, title_fontsize=LEG_TITLE_FS, ncol=1,
                    )

            title = cfg["title"] + r"  all $Q^2$ bins"
            if w_label:
                title += "     " + w_label
            fig.suptitle(title, fontsize=LEG_TITLE_FS + 2, y=0.97)
            stem = (f"{tag}_group_allQ2" if iw_key is None
                    else f"{tag}_group_allQ2_W{iw_key}")
            _save(fig, outdir, stem)


def make_radcorr_vs_Q2(csv_root, models, bin_keys, outdir):
    data_store   = defaultdict(dict)
    tprime_edges = None
    w_info       = {}

    for m in models:
        for (iq, iw) in bin_keys:
            df = load_csv(csv_root, m, iq, iw)
            if df is None:
                continue
            if iw not in w_info and "W_lo" in df.columns and pd.notna(df["W_lo"].iloc[0]):
                w_info[iw] = (float(df["W_lo"].iloc[0]), float(df["W_hi"].iloc[0]))
            if tprime_edges is None and "tprime_lo" in df.columns and "tprime_hi" in df.columns:
                tprime_edges = list(zip(
                    df["tprime_lo"].to_numpy(float),
                    df["tprime_hi"].to_numpy(float),
                ))
            rc  = _safe_col(df, "RadCorr")
            q2c = float(df["Q2_center"].iloc[0]) if "Q2_center" in df.columns else np.nan
            for it, rc_val in enumerate(rc):
                key = (iw, it)
                if m not in data_store[key]:
                    data_store[key][m] = ([], [])
                data_store[key][m][0].append(q2c)
                data_store[key][m][1].append(rc_val)

    if tprime_edges is None:
        return

    n_t = len(tprime_edges)
    for iw_key in _iw_unique(bin_keys):
        ncols = min(n_t, 4)
        nrows = (n_t + ncols - 1) // ncols

        with mpl.rc_context(PANEL_RC):
            fig, axes_arr = plt.subplots(
                nrows, ncols, figsize=(7*ncols, 6*nrows),
                constrained_layout=True, squeeze=False,
            )
            axes_flat = axes_arr.ravel()
            first_h, first_l = [], []

            for it in range(n_t):
                ax = axes_flat[it]; _style_ax(ax)
                tlo, thi = tprime_edges[it]
                ax.set_title(
                    rf"$-t' \in [{tlo:.3f},\,{thi:.3f}]\ \mathrm{{GeV}}^2$",
                    fontsize=LEG_FS, pad=4,
                )
                ax.set_xlabel(r"$Q^2\ [\mathrm{GeV}^2]$")
                ax.set_ylabel(r"$C_\mathrm{rad}$")
                ax.axhline(1.0, color="gray", lw=1.2, ls="--", zorder=1)

                for im, m in enumerate(models):
                    key = (iw_key, it)
                    if m not in data_store[key]:
                        continue
                    q2s, rcs = data_store[key][m]
                    if not q2s:
                        continue
                    order  = np.argsort(q2s)
                    q2s_s  = np.array(q2s)[order]
                    rcs_s  = np.array(rcs)[order]
                    mask   = np.isfinite(rcs_s) & (rcs_s != 0)
                    if not mask.any():
                        continue
                    color = MODEL_COLORS[im % len(MODEL_COLORS)]
                    h, = ax.plot(
                        q2s_s[mask], rcs_s[mask], "o-",
                        color=color, markersize=7, mfc="white", mew=2.0, lw=1.8,
                        zorder=5 + im, label=m.replace("_", " "),
                    )
                    if it == 0:
                        first_h.append(h); first_l.append(m.replace("_", " "))

            for it in range(n_t, len(axes_flat)):
                axes_flat[it].axis("off")

            if first_h:
                if n_t < len(axes_flat):
                    la = axes_flat[n_t]; la.axis("off")
                    la.legend(first_h, first_l, loc="center", frameon=False,
                              fontsize=LEG_FS, ncol=1)
                else:
                    fig.legend(first_h, first_l, loc="upper right",
                               frameon=False, fontsize=LEG_FS, ncol=1)

            w_title = ""
            if iw_key in w_info:
                wlo, whi = w_info[iw_key]
                w_title = rf"  $W \in [{wlo:.2f},\,{whi:.2f}]\ \mathrm{{GeV}}$"
            fig.suptitle(r"$C_\mathrm{rad}$ vs $Q^2$" + w_title,
                         fontsize=LEG_TITLE_FS + 2)
            stem = ("radcorr_vs_Q2" if iw_key is None
                    else f"radcorr_vs_Q2_W{iw_key}")
            _save(fig, outdir, stem)


def make_comparison_with_ratio(csv_root, models, iq, iw, canvas_key, outdir, logy=True, show_ratio=True):
    """
    Comparison plot across models for a single (Q,W) bin.

    If show_ratio=True:
        2-panel plot with lower panel showing Model / Reference ratio.
    If show_ratio=False:
        single-panel plot only (no ratio panel). This is used for Acceptance
        and Reduced-cross-section comparisons where the ratio panel is not needed.
    """
    cfg    = _CANVAS_CFG[canvas_key]
    col    = cfg["col"]
    ec     = cfg["err_col"]
    ylabel = cfg["ylabel"]

    dfs = {m: df for m in models if (df := load_csv(csv_root, m, iq, iw)) is not None}
    if not dfs:
        return

    ref_label = list(dfs.keys())[0]
    ref_df    = dfs[ref_label]
    t_ref     = ref_df["tprime_center"].to_numpy(float)
    y_ref     = _safe_col(ref_df, col)
    y_ref_s   = np.where(np.isfinite(y_ref) & (y_ref > 0), y_ref, np.nan)

    with mpl.rc_context(PANEL_RC):
        if show_ratio:
            fig = plt.figure(figsize=(10, 12))
            gs  = mpl.gridspec.GridSpec(
                2, 1, figure=fig, height_ratios=[3, 1], hspace=0.0,
                left=0.13, right=0.97, top=0.91, bottom=0.09,
            )
            ax_m = fig.add_subplot(gs[0])
            ax_r = fig.add_subplot(gs[1], sharex=ax_m)
            _style_ax(ax_m); _style_ax(ax_r)
        else:
            fig, ax_m = plt.subplots(figsize=(10, 8))
            ax_r = None
            _style_ax(ax_m)

        handles, labels_l = [], []
        for im, (m, df) in enumerate(dfs.items()):
            color = MODEL_COLORS[im % len(MODEL_COLORS)]
            t     = df["tprime_center"].to_numpy(float)
            y     = _safe_col(df, col)
            yerr  = _safe_col(df, ec) if ec else np.zeros_like(y)
            xlo   = t - df["tprime_lo"].to_numpy(float)
            xhi   = df["tprime_hi"].to_numpy(float) - t
            # For efficiency, mask out unphysical values (>1 or <=0)
            if canvas_key == "eff":
                mask = np.isfinite(y) & (y > 0) & (y <= 1.0)
            else:
                mask  = np.isfinite(y) & (y > 0)
            if not mask.any():
                continue

            h = ax_m.errorbar(
                t[mask], y[mask], yerr=yerr[mask] if yerr is not None else None,
                xerr=(xlo[mask], xhi[mask]),
                fmt="o", color=color, ecolor=color, elinewidth=1.5, capsize=3,
                markersize=7, mfc="white", mew=2.0, lw=1.8, zorder=5 + im,
                label=m.replace("_", " "),
            )
            handles.append(h); labels_l.append(m.replace("_", " "))

            if show_ratio and ax_r is not None:
                if im == 0:
                    ax_r.axhline(1.0, color=color, lw=1.5, ls="--", zorder=1)
                else:
                    # interpolate current model onto reference t' bins
                    y_int = np.interp(t_ref, t[mask], y[mask], left=np.nan, right=np.nan)
                    ratio = y_int / y_ref_s
                    rmask = np.isfinite(ratio)
                    if rmask.any():
                        ax_r.errorbar(
                            t_ref[rmask], ratio[rmask],
                            fmt="o", color=color, ecolor=color,
                            elinewidth=1.5, capsize=3,
                            markersize=6, mfc="white", mew=1.8, lw=1.5, zorder=5 + im,
                        )

        # --- Weighted average (only for reduced cross-section) ---
        if canvas_key == "red" and len(dfs) >= 2:
            avg_df = _weighted_average_reduced_xs(list(dfs.values()))
            if avg_df is not None and len(avg_df) > 0:
                t_a   = avg_df["tprime_center"].to_numpy(float)
                y_a   = avg_df["ReducedCrossSection"].to_numpy(float)
                ye_a  = avg_df["ReducedCrossSection_Err"].to_numpy(float)
                xlo_a = (t_a - avg_df["tprime_lo"].to_numpy(float)
                         if "tprime_lo" in avg_df.columns
                         else np.zeros_like(t_a))
                xhi_a = (avg_df["tprime_hi"].to_numpy(float) - t_a
                         if "tprime_hi" in avg_df.columns
                         else np.zeros_like(t_a))
                mask_a = np.isfinite(t_a) & np.isfinite(y_a) & (y_a > 0)
                if mask_a.any():
                    h_avg = ax_m.errorbar(
                        t_a[mask_a], y_a[mask_a],
                        yerr=ye_a[mask_a],
                        xerr=(np.where(np.isfinite(xlo_a[mask_a]), xlo_a[mask_a], 0.0),
                              np.where(np.isfinite(xhi_a[mask_a]), xhi_a[mask_a], 0.0)),
                        fmt="o", color="#111111", ecolor="#111111",
                        elinewidth=2.0, capsize=4,
                        markersize=8, mfc="white", mew=2.2, lw=2.0,
                        zorder=20,
                    )
                    handles.append(h_avg)
                    labels_l.append("Weighted avg.")

        # axes formatting
        ax_m.set_xlim(XLIM_T)

        if cfg["logy"] and logy and canvas_key not in ("acc", "rad", "eff"):
            ax_m.set_yscale("log")
            ax_m.set_ylim(YLIM_RED if canvas_key == "red" else YLIM_XS)
        elif canvas_key == "acc":
            ax_m.set_ylim(0.0, 0.035)
        elif canvas_key == "eff":
            ax_m.axhline(1.0, color="gray", lw=1.2, ls="--", zorder=0)
            ax_m.set_ylim(0.0, 1.05)
        elif canvas_key == "rad":
            ax_m.axhline(1.0, color="gray", lw=1.2, ls="--", zorder=0)
            ax_m.set_ylim(0.6, 1.0)

        ax_m.set_ylabel(ylabel)
        ax_m.set_xlabel(r"$-t'\ [\mathrm{GeV}^2]$")

        if show_ratio and ax_r is not None:
            ax_m.tick_params(labelbottom=False)
            ax_r.set_xlim(XLIM_T)
            ax_r.set_ylabel(rf"Model / {ref_label.replace('_', ' ')}", fontsize=16)
            ax_r.set_xlabel(r"$-t'\ [\mathrm{GeV}^2]$")
            ax_r.set_ylim(0.0, 2.5)

        fig.suptitle(_bin_label(ref_df), fontsize=LEG_TITLE_FS + 2, y=0.965)

        ax_m.legend(handles, labels_l, loc="upper right", frameon=False,
                    fontsize=LEG_FS, ncol=1)

        stem = (f"{cfg['tag']}_comparison_Q{iq}" if iw is None
                else f"{cfg['tag']}_comparison_Q{iq}_W{iw}")
        _save(fig, outdir, stem)

def make_final_xs_group(csv_root, models, bin_keys, outdir, logy=True):
    make_group_plot(csv_root, models, bin_keys, "xs", outdir, logy=logy)


def make_raw_yields_per_bin(csv_root, models, bin_keys, outdir,
                            logy=True, luminosity=None, branching=None):
    """Per-(Q²,W) bin plot of the raw (uncorrected) dσ/dt' vs t'.

    Plots the same data points as the black "Raw (no corr.)" series in the
    corrections canvas — i.e. the cross-section with both acceptance and
    radiative corrections undone (xs × A × C_rad), and without any virtual
    photon flux division.  One file per kinematic bin, all models overlaid.
    """
    os.makedirs(outdir, exist_ok=True)

    for (iq, iw) in bin_keys:
        dfs = {m: df for m in models
               if (df := load_csv(csv_root, m, iq, iw)) is not None}
        if not dfs:
            continue

        ref_df = next(iter(dfs.values()))

        with mpl.rc_context(PANEL_RC):
            fig, ax_m = plt.subplots(figsize=(10, 8))
            _style_ax(ax_m)

            handles, labels_l = [], []
            for im, (m, df) in enumerate(dfs.items()):
                color    = MODEL_COLORS[im % len(MODEL_COLORS)]
                t        = df["tprime_center"].to_numpy(float)
                xlo      = t - df["tprime_lo"].to_numpy(float)
                xhi      = df["tprime_hi"].to_numpy(float) - t

                variants = _build_variants(df, luminosity, branching)
                raw_v, raw_e = variants["raw"]

                mask = np.isfinite(raw_v) & (raw_v > 0)
                if not mask.any():
                    continue

                h = ax_m.errorbar(
                    t[mask], raw_v[mask],
                    yerr=raw_e[mask],
                    xerr=(xlo[mask], xhi[mask]),
                    fmt="o", color=color, ecolor=color,
                    elinewidth=1.5, capsize=3,
                    markersize=7, mfc="white", mew=2.0, lw=1.8,
                    zorder=5 + im,
                    label=m.replace("_", " "),
                )
                handles.append(h)
                labels_l.append(m.replace("_", " "))

            if not handles:
                plt.close(fig)
                continue

            ax_m.set_xlim(XLIM_T)
            if logy:
                ax_m.set_yscale("log")
                ax_m.set_ylim(YLIM_XS_RAW)
            ax_m.set_xlabel(r"$-t'\ [\mathrm{GeV}^2]$")
            ax_m.set_ylabel(
                #r"$\mathrm{d}\sigma/\mathrm{d}t'\ [\mathrm{nb/GeV}^2]$"
                r"Yields (arb. units)"
                "\n(raw: no acc. / rad. corr.)"
            )
            ax_m.legend(handles, labels_l, loc="upper right",
                        frameon=False, fontsize=LEG_FS, ncol=1)
            fig.suptitle(_bin_label(ref_df), fontsize=LEG_TITLE_FS + 2, y=0.97)

            stem = (f"raw_yields_Q{iq}" if iw is None
                    else f"raw_yields_Q{iq}_W{iw}")
            _save(fig, outdir, stem)


def make_reduced_xs_group(csv_root, models, bin_keys, outdir, logy=True):
    make_group_plot(csv_root, models, bin_keys, "red", outdir, logy=logy)


# -----------------------------------------------------------------------------
#  Dipole fit on REDUCED cross-section
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#  Dipole fit on REDUCED cross-section   model, fitter, per-bin plots,
#  residual panels, parameter-vs-Q/W summary plots
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#  Dipole fit (PAC39 / proposal-style)
# -----------------------------------------------------------------------------
# Fit is performed in absolute |t|, where:
#   |t| = t' + |t_min(Q2, W_mean)|
#
# Dipole model used (proposal / PAC39):
#   sigma_red(|t|) = sigma0 / (1 + |t| / m_g^2)^4
#
# The gluonic transverse radius (mean-square) is:
#   <b^2>_g [GeV^-2] = 8 / m_g^2
# and we convert with (0.197326 fm)^2 per GeV^-2.

GEV2_TO_FM2 = (0.197326 ** 2)

def _dipole_model_abs(t_abs, sigma0, mg2):
    """Proposal/PAC39 dipole: sigma0 / (1 + |t|/m_g^2)^4."""
    return sigma0 / (1.0 + (t_abs / mg2)) ** 4


def _safe_sigma(y_err):
    """Ensure strictly-positive sigmas for curve_fit."""
    y_err = np.asarray(y_err, dtype=float)
    finite = np.isfinite(y_err)
    if not finite.any():
        return np.ones_like(y_err, dtype=float)
    floor = np.nanmedian(y_err[finite]) * 1e-3
    if (not np.isfinite(floor)) or floor <= 0:
        floor = 1e-12
    return np.where((~finite) | (y_err <= 0), floor, y_err)



def _fit_hs_Ds0_from_dsdt(df, model_name,
                          tprime_min=0.0, tprime_max=1.0,
                          As0=HS_AS0_DEFAULT, mA=HS_mA_DEFAULT, mD=HS_mD_DEFAULT,
                          Ds0_bounds=(-5.0, 2.0),
                          cAA=1.0, cAD=1.0, cDD=1.0):
    """Fit Ds(0) using the HS-inspired template to the reduced cross-section σ_red = (dσ/dt')/Γ_v (nb).

    We convert each point's t' to |t| = t' + |t_min| (using df['tmin_abs'] already attached),
    then fit y(|t|) with free parameters (N, Ds0).

    Returns a dict with Ds0 and uncertainty (1σ from covariance) or None.
    """
    if "ReducedCrossSection" not in df.columns:
        return None

    if "tprime_center" not in df.columns:
        return None

    tprime = df["tprime_center"].to_numpy(float)

    # Robustly compute |t_min|. Prefer per-row column; fall back to
    # computing from Q2_center + W_MEAN_GEV so it is never silently 0.
    tmin_abs = 0.0
    if "tmin_abs" in df.columns:
        val = float(df["tmin_abs"].iloc[0])
        if np.isfinite(val) and val > 1e-6:
            tmin_abs = val
    if tmin_abs < 1e-6 and "Q2_center" in df.columns:
        q2c = float(df["Q2_center"].iloc[0])
        wc  = float(df["W_center"].iloc[0]) if "W_center" in df.columns else W_MEAN_GEV
        w_hs = min(wc, W_MEAN_GEV) if np.isfinite(wc) else W_MEAN_GEV
        tmin_computed = t_min_phi(q2c, w_hs)
        if np.isfinite(tmin_computed) and tmin_computed > 1e-6:
            tmin_abs = tmin_computed

    t_abs = tprime + tmin_abs

    # HS model validity: restrict to |t| < 1 GeV^2 (paper recommendation).
    hs_t_abs_max = 2.0
    tprime_max_eff = min(float(tprime_max), max(0.05, hs_t_abs_max - tmin_abs))

    y = df["ReducedCrossSection"].to_numpy(float)
    yerr = df["ReducedCrossSection_Err"].to_numpy(float) if "ReducedCrossSection_Err" in df.columns else np.zeros_like(y)

    mask = (
        np.isfinite(tprime) & np.isfinite(t_abs) &
        np.isfinite(y) & (y > 0) &
        (tprime >= float(tprime_min)) & (tprime <= tprime_max_eff)
    )
    if int(mask.sum()) < 3:
        print(f"  [WARN] {model_name}: only {int(mask.sum())} ds/dt points in t' "
              f"[{tprime_min},{tprime_max}] -> skipping Ds0 fit.")
        return None

    x_f = t_abs[mask]
    y_f = y[mask]

    with np.errstate(invalid="ignore"):
        yerr = np.where(np.isfinite(yerr) & (yerr > 0), yerr, 0.10 * y)
    s_f = _safe_sigma(yerr[mask])

    Ds0_0 = -0.5
    y0_model = hs_dsdt_template(
        np.array([x_f[0]]), 1.0, 0.0,
        As0=As0, mA=mA, mD=mD, cAA=cAA, cAD=cAD, cDD=cDD
    )[0]
    N0 = float(y_f[0] / y0_model) if np.isfinite(y0_model) and y0_model > 0 else float(np.nanmax(y_f))

    def fwrap(tabs, N, Ds0):
        return hs_dsdt_template(tabs, N, Ds0, As0=As0, mA=mA, mD=mD, cAA=cAA, cAD=cAD, cDD=cDD)

    try:
        popt, pcov = curve_fit(
            fwrap, x_f, y_f,
            sigma=s_f, absolute_sigma=True,
            p0=[N0, Ds0_0],
            bounds=([0.0, Ds0_bounds[0]], [np.inf, Ds0_bounds[1]]),
            maxfev=40000,
        )
    except Exception as exc:
        print(f"  [WARN] {model_name}: Ds0 fit failed: {exc}")
        return None

    perr = np.sqrt(np.diag(pcov)) if (pcov is not None and np.all(np.isfinite(pcov))) else [np.nan, np.nan]
    N_fit, Ds0_fit = [float(v) for v in popt]
    dN, dDs0 = [float(v) for v in perr]

    resid = y_f - fwrap(x_f, *popt)
    chi2 = float(np.sum((resid / s_f) ** 2))
    ndf = int(mask.sum()) - 2

    return dict(
        model=model_name,
        N=N_fit, dN=dN,
        Ds0=Ds0_fit, dDs0=dDs0,
        chi2=chi2, ndf=ndf,
        tmin_abs=tmin_abs,
        tprime_min=float(tprime_min),
        tprime_max=float(tprime_max_eff),
        t_abs_fit=x_f, y_fit=y_f, yerr_fit=s_f,
    )

def _fit_dipole_reduced(df, model_name, t_min=0.5, t_max=5.5):
    """Free fit (sigma0, m_g^2) to reduced cross-section, in |t|.

    We cut points by t' in [t_min, t_max] (same CLI meaning as before),
    but the fit x-variable is |t| = t' + |t_min|.
    """
    tprime = df["tprime_center"].to_numpy(float)
    y      = _safe_col(df, "ReducedCrossSection")
    yerr   = _safe_col(df, "ReducedCrossSection_Err")

    # Use mean-Q2 from the bin + hardcoded mean W for |t_min|
    Q2_mean = np.nan
    if "Q2_center" in df.columns:
        qvals = df["Q2_center"].dropna().to_numpy(float)
        if len(qvals) > 0:
            Q2_mean = float(np.nanmean(qvals))
    tmin_abs = float(t_min_phi(Q2_mean, W_MEAN_GEV)) if np.isfinite(Q2_mean) else 0.0
    if not np.isfinite(tmin_abs):
        tmin_abs = 0.0

    # Build |t| arrays
    t_abs = tprime + tmin_abs

    # Error handling: 10% relative floor if missing, then strictly-positive
    with np.errstate(invalid="ignore"):
        yerr = np.where(
            np.isfinite(yerr) & (yerr > 0),
            yerr,
            0.10 * np.where(np.isfinite(y) & (y > 0), y, 1.0),
        )
    yerr = _safe_sigma(yerr)

    mask_all = np.isfinite(tprime) & np.isfinite(t_abs) & np.isfinite(y) & (y > 0)
    mask     = mask_all & (tprime >= t_min) & (tprime <= t_max)

    if int(mask.sum()) < 3:
        print(f"  [WARN] {model_name}: only {int(mask.sum())} valid reduced-XS points "
              f"in t'  [{t_min},{t_max}] GeV^2  skipping dipole fit.")
        return None

    x_f = t_abs[mask]
    y_f = y[mask]
    s_f = yerr[mask]

    # Initial guess similar to the reference script
    sigma0_0 = float(y_f[0]) if np.isfinite(y_f[0]) else float(np.nanmedian(y_f))
    mg2_0    = 0.65

    try:
        popt, pcov = curve_fit(
            _dipole_model_abs,
            x_f, y_f,
            sigma=s_f,
            absolute_sigma=True,
            p0=[sigma0_0, mg2_0],
            bounds=([0.0, 1e-6], [np.inf, np.inf]),
            maxfev=40000,
        )
    except Exception as exc:
        print(f"  [WARN] {model_name}: dipole fit failed  {exc}")
        return None

    perr = (np.sqrt(np.diag(pcov))
            if (pcov is not None and np.all(np.isfinite(pcov)))
            else np.full(2, np.nan))

    sigma0, mg2 = [float(v) for v in popt]
    dsigma0, dmg2 = [float(v) for v in perr]

    # <b^2> in fm^2
    # The fit runs in |t| = t' + tmin_abs.  Re-expressing the dipole
    # σ₀/(1+|t|/m_g²)⁴ in t' space near t'=0 gives effective slope scale
    # m_eff² = m_g² + tmin_abs, so the physically correct formula is:
    #   <b²> = 8 / (m_g² + tmin_abs) × (ħc)²
    m_eff2 = mg2 + tmin_abs   # effective t'-slope scale [GeV²]
    if m_eff2 > 0 and np.isfinite(m_eff2):
        b2_gev  = 8.0 / m_eff2
        # uncertainty: only mg2 is a fit parameter; tmin_abs is fixed kinematics
        db2_gev = 8.0 * dmg2 / (m_eff2 ** 2) if np.isfinite(dmg2) else np.nan
        b2  = b2_gev  * GEV2_TO_FM2
        db2 = db2_gev * GEV2_TO_FM2 if np.isfinite(db2_gev) else np.nan
        b   = np.sqrt(b2) if (np.isfinite(b2) and b2 > 0) else np.nan
        db  = 0.5 * db2 / b if (np.isfinite(db2) and np.isfinite(b) and b > 0) else np.nan
    else:
        b2 = db2 = b = db = np.nan

    resid = y_f - _dipole_model_abs(x_f, *popt)
    chi2  = float(np.sum((resid / s_f) ** 2))
    ndf   = int(mask.sum()) - 2
    pull  = resid / s_f

    return dict(
        model=model_name,
        fit_mode="free",
        tmin_abs=tmin_abs,
        # Store the t' fit window used for this dipole fit
        tprime_min=float(t_min),
        tprime_max=float(t_max),
        Q2_mean=Q2_mean,
        # keep legacy key names used elsewhere in the script:
        A=sigma0, dA=dsigma0,
        t0=mg2,   dt0=dmg2,
        mg2=mg2,  dmg2=dmg2,
        b=b, db=db, b2=b2, db2=db2,
        chi2=chi2, ndf=ndf,
        t_fit=x_f, y_fit=y_f, yerr_fit=s_f,
        resid=resid, pull=pull,
        t_all=t_abs[mask_all],
        y_all=y[mask_all],
        yerr_all=yerr[mask_all],
        popt=np.array([sigma0, mg2], dtype=float),
    )


def _make_one_fit_canvas(ax_m, ax_r, df, model, color, t_min, t_max, logy, show_ylabel=True):
    """Draw one bin's reduced-XS dipole fit plotted vs t' (Mandelstam t-prime).

    Points are selected and displayed on the t' axis.  The dipole fit is
    performed internally in |t| = t' + |t_min(Q2, W_MEAN_GEV)|, and the
    resulting curve is back-transformed to t' for display.  This keeps the
    physics correct while giving the experimentally natural x-axis.
    """
    _style_ax(ax_m); _style_ax(ax_r)

    tprime = df["tprime_center"].to_numpy(float)
    y      = _safe_col(df, "ReducedCrossSection")
    yerr   = _safe_col(df, "ReducedCrossSection_Err")

    # Compute |t_min| for this kinematic bin (mean-Q2, hardcoded W_MEAN_GEV).
    # The fit uses |t| = t' + tmin_abs internally; we only display in t'.
    Q2_mean = np.nan
    if "Q2_center" in df.columns:
        qvals = df["Q2_center"].dropna().to_numpy(float)
        if len(qvals) > 0:
            Q2_mean = float(np.nanmean(qvals))
    tmin_abs = float(t_min_phi(Q2_mean, W_MEAN_GEV)) if np.isfinite(Q2_mean) else 0.0
    if not np.isfinite(tmin_abs):
        tmin_abs = 0.0

    # x-errors: bin half-widths in t' (tmin shift cancels)
    xlo = tprime - df["tprime_lo"].to_numpy(float)
    xhi = df["tprime_hi"].to_numpy(float) - tprime

    mask_plot = np.isfinite(tprime) & np.isfinite(y) & (y > 0)
    if not mask_plot.any():
        return None, None

    # yerr: 10% floor, strictly positive
    with np.errstate(invalid="ignore"):
        yerr_plot = np.where(
            np.isfinite(yerr) & (yerr > 0),
            yerr,
            0.10 * np.where(np.isfinite(y) & (y > 0), y, 1.0),
        )
    yerr_plot = _safe_sigma(yerr_plot)

    # data points plotted in t'
    h_data = ax_m.errorbar(
        tprime[mask_plot], y[mask_plot],
        yerr=yerr_plot[mask_plot],
        xerr=(xlo[mask_plot], xhi[mask_plot]),
        fmt="o", color=color, ecolor=color,
        elinewidth=1.5, capsize=3,
        markersize=7, mfc="white", mew=2.0, lw=1.8, zorder=5,
        label=model.replace("_", " "),
    )

    # dipole fit (done in |t| space inside _fit_dipole_reduced)
    res = _fit_dipole_reduced(df, model, t_min=t_min, t_max=t_max)

    if res is not None:
        # Evaluate fit curve in t' space: |t| = x_dense + tmin_abs
        x_dense_tp = np.linspace(
            max(0.0, np.nanmin(tprime[mask_plot]) * 0.7),
            np.nanmax(tprime[mask_plot]) * 1.10,
            600,
        )
        y_curve = _dipole_model_abs(x_dense_tp + tmin_abs, *res["popt"])

        chi2_str = (f"chi2/ndf = {res['chi2']:.1f}/{res['ndf']}"
                    if res["ndf"] > 0 else "")

        fit_txt = (
            f"|t_min| = {res['tmin_abs']:.4f} GeV^2\n"
            f"sigma0 = {_fmt_pm(res['A'], res['dA'])} nb\n"
            f"m_g^2 = {_fmt_pm(res['t0'], res['dt0'])} GeV^2\n"
            f"<b^2> = {_fmt_pm(res['b2'], res['db2'])} fm^2\n"
            f"{chi2_str}"
        ).strip()

        ax_m.plot(x_dense_tp, y_curve, color=color, lw=2.2, ls="--", zorder=4,
                  label="Dipole fit")

        # shade fit region in t' (no tmin shift needed on display axis)
        ax_m.axvspan(t_min, t_max, alpha=0.06, color=color, zorder=0)

        ax_m.text(
            0.97, 0.97, fit_txt,
            transform=ax_m.transAxes,
            fontsize=14, va="top", ha="right", color=color,
            bbox=dict(boxstyle="round,pad=0.35", fc="white",
                      ec=color, alpha=0.85, lw=1.2),
        )

        # pull panel: convert fit x-coords from |t| back to t'
        t_fit_tp = res["t_fit"] - tmin_abs
        ax_r.axhline(0.0, color="gray", lw=1.2, ls="--", zorder=1)
        ax_r.axhline( 1.0, color="gray", lw=0.8, ls=":", zorder=1)
        ax_r.axhline(-1.0, color="gray", lw=0.8, ls=":", zorder=1)
        ax_r.errorbar(
            t_fit_tp, res["pull"],
            yerr=np.ones_like(res["pull"]),
            fmt="o", color=color, ecolor=color,
            elinewidth=1.2, capsize=2,
            markersize=5, mfc="white", mew=1.8, zorder=5,
        )

    # axes formatting
    if logy:
        ax_m.set_yscale("log")
    if show_ylabel:
        ax_m.set_ylabel(r"$\sigma_\mathrm{red}\ [\mathrm{nb}]$", fontsize=20)
    ax_m.tick_params(labelbottom=False)
    ax_m.tick_params(axis="both", labelsize=18)
    ax_m.set_xlim(left=max(0.0, float(np.nanmin(tprime[mask_plot]) * 0.8)))

    ax_r.set_ylim(-4.5, 4.5)
    ax_r.set_xlabel(r"$-t'\ [\mathrm{GeV}^2]$", fontsize=20)
    if show_ylabel:
        ax_r.set_ylabel("Pull", fontsize=17)
    ax_r.tick_params(axis="both", labelsize=16)
    ax_r.set_xlim(ax_m.get_xlim())

    return res, h_data


# -----------------------------------------------------------------------------
#  Weighted-average reduced cross-section across datasets (inverse-variance)
# -----------------------------------------------------------------------------
def _weighted_average_reduced_xs(dfs_list):
    """
    Combine multiple DataFrames into a single inverse-variance weighted average
    reduced cross-section, bin-by-bin in t'.

    Weight per dataset point: w = 1 / sigma^2 where sigma is ReducedCrossSection_Err.
    If an error is missing/non-positive, we apply a conservative 10% relative floor.

    Returns
    -------
    pd.DataFrame or None
        Columns: tprime_center, tprime_lo, tprime_hi, ReducedCrossSection,
                 ReducedCrossSection_Err, n_datasets (+ kinematics from first df).
    """
    if not dfs_list:
        return None

    from collections import defaultdict

    num   = defaultdict(float)
    denom = defaultdict(float)
    count = defaultdict(int)
    t_lo_map, t_hi_map = {}, {}

    ref_df = dfs_list[0]

    for df in dfs_list:
        if df is None or len(df) == 0:
            continue

        t = df["tprime_center"].to_numpy(float) if "tprime_center" in df.columns else None
        if t is None:
            continue

        y = _safe_col(df, "ReducedCrossSection")
        ye = _safe_col(df, "ReducedCrossSection_Err")

        # 10% relative error floor where needed
        with np.errstate(invalid="ignore"):
            ye = np.where(
                np.isfinite(ye) & (ye > 0),
                ye,
                0.10 * np.where(np.isfinite(y) & (y > 0), y, np.nan),
            )

        tlo = _safe_col(df, "tprime_lo")
        thi = _safe_col(df, "tprime_hi")

        for i in range(len(df)):
            tv = float(t[i]); yv = float(y[i]); sv = float(ye[i])
            if not (np.isfinite(tv) and np.isfinite(yv) and yv > 0 and np.isfinite(sv) and sv > 0):
                continue
            tk = round(tv, 4)  # robust match across datasets
            w = 1.0 / (sv * sv)
            num[tk]   += yv * w
            denom[tk] += w
            count[tk] += 1
            if tk not in t_lo_map:
                t_lo_map[tk] = float(tlo[i]) if np.isfinite(tlo[i]) else np.nan
                t_hi_map[tk] = float(thi[i]) if np.isfinite(thi[i]) else np.nan

    if not denom:
        return None

    t_keys = sorted(denom.keys())
    t_cen  = np.array(t_keys, dtype=float)
    y_avg  = np.array([num[k] / denom[k] for k in t_keys], dtype=float)
    y_err  = np.array([1.0 / np.sqrt(denom[k]) for k in t_keys], dtype=float)
    n_ds   = np.array([count[k] for k in t_keys], dtype=int)
    t_lo   = np.array([t_lo_map.get(k, np.nan) for k in t_keys], dtype=float)
    t_hi   = np.array([t_hi_map.get(k, np.nan) for k in t_keys], dtype=float)

    out = pd.DataFrame({
        "tprime_center":           t_cen,
        "tprime_lo":               t_lo,
        "tprime_hi":               t_hi,
        "ReducedCrossSection":     y_avg,
        "ReducedCrossSection_Err": y_err,
        "n_datasets":              n_ds,
    })

    # Carry over kinematics from reference df (single value per bin)
    for col in ("Q2_lo", "Q2_hi", "Q2_center", "W_lo", "W_hi", "W_center"):
        if col in ref_df.columns and not ref_df[col].isna().all():
            out[col] = float(ref_df[col].iloc[0])

    return out


# -- combined fit canvas: individual datasets (faint) + weighted avg + fit -----


def _weighted_average_dsdt(dfs_list):
    """
    Combine multiple DataFrames into a single inverse-variance weighted average
    FINAL cross-section dσ/dt' (CrossSection), bin-by-bin in t'.

    Weight per dataset point: w = 1 / sigma^2 where sigma is CrossSection_Err.
    If an error is missing/non-positive, apply a conservative 10% relative floor.

    Returns
    -------
    pd.DataFrame or None
        Columns: tprime_center, tprime_lo, tprime_hi, CrossSection,
                 CrossSection_Err, n_datasets (+ kinematics from first df).
    """
    if not dfs_list:
        return None

    from collections import defaultdict

    num   = defaultdict(float)
    denom = defaultdict(float)
    count = defaultdict(int)
    t_lo_map, t_hi_map = {}, {}

    ref_df = dfs_list[0]

    for df in dfs_list:
        if df is None or len(df) == 0:
            continue
        if "tprime_center" not in df.columns:
            continue

        t  = df["tprime_center"].to_numpy(float)
        y  = _safe_col(df, "CrossSection")
        ye = _safe_col(df, "CrossSection_Err")

        # 10% relative error floor where needed
        with np.errstate(invalid="ignore"):
            ye = np.where(
                np.isfinite(ye) & (ye > 0),
                ye,
                0.10 * np.where(np.isfinite(y) & (y > 0), y, np.nan),
            )

        tlo = _safe_col(df, "tprime_lo")
        thi = _safe_col(df, "tprime_hi")

        for i in range(len(df)):
            tv = float(t[i]); yv = float(y[i]); sv = float(ye[i])
            if not (np.isfinite(tv) and np.isfinite(yv) and yv > 0 and np.isfinite(sv) and sv > 0):
                continue
            tk = round(tv, 4)
            w = 1.0 / (sv * sv)
            num[tk]   += yv * w
            denom[tk] += w
            count[tk] += 1
            if tk not in t_lo_map:
                t_lo_map[tk] = float(tlo[i]) if np.isfinite(tlo[i]) else np.nan
                t_hi_map[tk] = float(thi[i]) if np.isfinite(thi[i]) else np.nan

    if not denom:
        return None

    t_keys = sorted(denom.keys())
    t_cen  = np.array(t_keys, dtype=float)
    y_avg  = np.array([num[k] / denom[k] for k in t_keys], dtype=float)
    y_err  = np.array([1.0 / np.sqrt(denom[k]) for k in t_keys], dtype=float)
    n_ds   = np.array([count[k] for k in t_keys], dtype=int)
    t_lo   = np.array([t_lo_map.get(k, np.nan) for k in t_keys], dtype=float)
    t_hi   = np.array([t_hi_map.get(k, np.nan) for k in t_keys], dtype=float)

    out = pd.DataFrame({
        "tprime_center":     t_cen,
        "tprime_lo":         t_lo,
        "tprime_hi":         t_hi,
        "CrossSection":      y_avg,
        "CrossSection_Err":  y_err,
        "n_datasets":        n_ds,
    })

    # Carry over kinematics from reference df (single value per bin)
    for col in ("Q2_lo", "Q2_hi", "Q2_center", "W_lo", "W_hi", "W_center"):
        if col in ref_df.columns and not ref_df[col].isna().all():
            out[col] = float(ref_df[col].iloc[0])

    # Ensure tmin_abs exists for the fitter
    out = compute_tmin_for_df(out)

    return out
def _make_combined_fit_canvas(ax_m, ax_r, dfs_list, model_names,
                              avg_df, t_min, t_max, logy, show_ylabel=True):
    """
    Draw individual datasets (faint), the weighted-average spectrum (black),
    and a single dipole fit to the weighted average, using the ABS-|t| dipole:

        sigma(|t|) = sigma0 / (1 + |t|/m_g^2)^4

    where |t| = t' + |t_min|(Q2, W_MEAN_GEV).

    Returns the fit result dict (or None).
    """
    _style_ax(ax_m); _style_ax(ax_r)

    # Determine common |t_min| from Q2 mean of avg_df (fallback to first df)
    Q2_mean = np.nan
    if "Q2_center" in avg_df.columns:
        try:
            qvals = avg_df["Q2_center"].dropna().to_numpy(float)
            if len(qvals) > 0:
                Q2_mean = float(np.nanmean(qvals))
        except Exception:
            pass
    if not np.isfinite(Q2_mean) and dfs_list:
        df0 = dfs_list[0]
        if "Q2_center" in df0.columns:
            qvals = df0["Q2_center"].dropna().to_numpy(float)
            if len(qvals) > 0:
                Q2_mean = float(np.nanmean(qvals))
    tmin_abs = float(t_min_phi(Q2_mean, W_MEAN_GEV)) if np.isfinite(Q2_mean) else 0.0
    if not np.isfinite(tmin_abs):
        tmin_abs = 0.0

    # -- individual datasets in background (vs t') -------------------------
    for im, (m, df) in enumerate(zip(model_names, dfs_list)):
        tprime = df["tprime_center"].to_numpy(float)
        y      = _safe_col(df, "ReducedCrossSection")
        yerr   = _safe_col(df, "ReducedCrossSection_Err")
        xlo    = tprime - df["tprime_lo"].to_numpy(float)
        xhi    = df["tprime_hi"].to_numpy(float) - tprime
        mask   = np.isfinite(y) & (y > 0) & np.isfinite(tprime)
        if not mask.any():
            continue
        color = MODEL_COLORS[im % len(MODEL_COLORS)]
        ax_m.errorbar(
            tprime[mask], y[mask],
            yerr=np.where(np.isfinite(yerr[mask]) & (yerr[mask] > 0), yerr[mask], 0.0),
            xerr=(np.where(np.isfinite(xlo[mask]), xlo[mask], 0.0),
                  np.where(np.isfinite(xhi[mask]), xhi[mask], 0.0)),
            fmt="o", color=color, ecolor=color,
            alpha=0.30, elinewidth=1.0, capsize=2,
            markersize=5, mfc=color, mew=1.0, lw=0.8, zorder=3,
            label=m.replace("_", " "),
        )

    # -- weighted average (vs t') -----------------------------------------
    tprime_a = avg_df["tprime_center"].to_numpy(float)
    y_a      = avg_df["ReducedCrossSection"].to_numpy(float)
    ye_a     = avg_df["ReducedCrossSection_Err"].to_numpy(float)
    xlo_a    = tprime_a - avg_df["tprime_lo"].to_numpy(float)
    xhi_a    = avg_df["tprime_hi"].to_numpy(float) - tprime_a
    mask_a   = np.isfinite(tprime_a) & np.isfinite(y_a) & (y_a > 0)

    ax_m.errorbar(
        tprime_a[mask_a], y_a[mask_a],
        yerr=ye_a[mask_a],
        xerr=(np.where(np.isfinite(xlo_a[mask_a]), xlo_a[mask_a], 0.0),
              np.where(np.isfinite(xhi_a[mask_a]), xhi_a[mask_a], 0.0)),
        fmt="o", color="#111111", ecolor="#111111",
        elinewidth=2.0, capsize=4,
        markersize=8, mfc="white", mew=2.2, lw=1.8, zorder=8,
        label="Weighted average",
    )

    # -- fit to weighted average ------------------------------------------
    # Reuse the same fitter used for per-model plots (absolute |t| model).
    res = _fit_dipole_reduced(avg_df, "combined", t_min=t_min, t_max=t_max)
    if res is not None:
        # Curve displayed in t': evaluate dipole at (x_dense + tmin_abs)
        t_dense_tp = np.linspace(
            max(0.0, np.nanmin(tprime_a[mask_a]) * 0.95),
            np.nanmax(tprime_a[mask_a]) * 1.10,
            600,
        )
        y_curve = _dipole_model_abs(t_dense_tp + tmin_abs, *res["popt"])
        ax_m.plot(t_dense_tp, y_curve, color="#111111", lw=2.5, ls="--",
                  zorder=9, label="Dipole fit (combined)")

        # fit range shading in t' (no tmin shift on display axis)
        ax_m.axvspan(t_min, t_max, alpha=0.035,
                     color="gray", zorder=0)

        # pull panel: convert |t| fit coords back to t'
        t_fit_tp = res["t_fit"] - tmin_abs
        ax_r.axhline(0.0, color="gray", lw=1.2, ls="--", zorder=1)
        ax_r.axhline( 1.0, color="gray", lw=0.8, ls=":", zorder=1)
        ax_r.axhline(-1.0, color="gray", lw=0.8, ls=":", zorder=1)
        ax_r.errorbar(
            t_fit_tp, res["pull"],
            yerr=np.ones_like(res["pull"]),
            fmt="o", color="#111111", ecolor="#111111",
            elinewidth=1.2, capsize=2,
            markersize=5, mfc="white", mew=1.8, zorder=5,
        )

        chi2_str = (rf"$\chi^2/\nu = {res['chi2']:.1f}/{res['ndf']}$"
                    if res["ndf"] > 0 else "")
        # Avoid any LaTeX parse issues with scientific notation by using plain text
        fit_txt = (
            f"|t_min|={res['tmin_abs']:.4f}  "
            f"m_g^2={res['mg2']:.3f}±{res['dmg2']:.3f}  "
            f"<b^2>={res['b2']:.3f}±{res['db2']:.3f} fm^2  "
            + (f"  chi2/ndf={res['chi2']:.1f}/{res['ndf']}" if res["ndf"] > 0 else "")
        )
        ax_m.text(
            0.97, 0.97, fit_txt,
            transform=ax_m.transAxes,
            fontsize=13, va="top", ha="right", color="#111111",
            bbox=dict(boxstyle="round,pad=0.35", fc="white",
                      ec="#111111", alpha=0.88, lw=1.2),
        )

    # axes
    ax_m.set_xlim(0.0, 6.0)
    if logy:
        ax_m.set_yscale("log")
    if show_ylabel:
        ax_m.set_ylabel(r"$\sigma_\mathrm{red}\ [\mathrm{nb}]$", fontsize=20)
    ax_m.tick_params(labelbottom=False, labelsize=18)
    ax_m.set_xlabel("")

    ax_r.set_xlim(0.0, 6.0)
    ax_r.set_ylim(-4.5, 4.5)
    ax_r.set_xlabel(r"$-t'\ [\mathrm{GeV}^2]$", fontsize=20)
    if show_ylabel:
        ax_r.set_ylabel("Pull", fontsize=17)
    ax_r.tick_params(labelsize=16)

    # Top x label on main panel
    ax_m.set_xlabel(r"$-t'\ [\mathrm{GeV}^2]$", fontsize=20)

    return res


def make_dipole_fit_plots(csv_root, models, bin_keys, outdir,
                          t_min=0.1, t_max=5.5, logy=True):
    """
    Two sets of output per (Q, W) bin:

    A) Per-model canvases  (unchanged  individual fit to each dataset)
        dipole_reduced_Q{iq}[_W{iw}].pdf

    B) Combined canvas  (NEW  weighted average of all datasets + single fit)
        dipole_combined_Q{iq}[_W{iw}].pdf

    The combined fit results (one row per kinematic bin) are used for
    the b vs x_B plot and summary CSV.
    """
    os.makedirs(outdir, exist_ok=True)

    per_model_results  = []   # per-model fit rows  (for parameter-vs-Q plots)
    combined_results   = []   # combined fit rows   (for b vs x_B)

    for (iq, iw) in bin_keys:
        dfs_m = [(m, load_csv(csv_root, m, iq, iw)) for m in models]
        dfs_m = [(m, df) for m, df in dfs_m if df is not None]
        if not dfs_m:
            continue

        m_names = [m  for m, _  in dfs_m]
        dfs     = [df for _, df in dfs_m]
        n_m     = len(dfs_m)

        # -- kinematics from first available df ----------------------------
        ref_df = dfs[0]
        q2c = (float(ref_df["Q2_center"].iloc[0])
               if "Q2_center" in ref_df.columns else np.nan)
        wc  = (float(ref_df["W_center"].iloc[0])
               if "W_center"  in ref_df.columns else np.nan)
        bin_title = _bin_label(ref_df)

        # -- A) per-model fit canvases -------------------------------------
        with mpl.rc_context(PANEL_RC):
            fig = plt.figure(figsize=(10 * n_m, 12))
            outer = mpl.gridspec.GridSpec(
                1, n_m, figure=fig,
                left=0.09, right=0.98, top=0.90, bottom=0.09,
                wspace=0.06,
            )
            for im, (m, df) in enumerate(dfs_m):
                inner = mpl.gridspec.GridSpecFromSubplotSpec(
                    2, 1, subplot_spec=outer[im],
                    height_ratios=[3, 1], hspace=0.0,
                )
                ax_m = fig.add_subplot(inner[0])
                ax_r = fig.add_subplot(inner[1], sharex=ax_m)

                color     = MODEL_COLORS[im % len(MODEL_COLORS)]
                res, h_d  = _make_one_fit_canvas(
                    ax_m, ax_r, df, m, color,
                    t_min, t_max, logy, show_ylabel=(im == 0),
                )
                if h_d is not None:
                    ax_m.legend(loc="lower left", frameon=False, fontsize=LEG_FS - 2)
                if im != 0:
                    ax_m.tick_params(labelleft=False)
                    ax_r.tick_params(labelleft=False)
                ax_m.annotate(
                    m.replace("_", " "),
                    xy=(0.03, 0.03), xycoords="axes fraction",
                    fontsize=LEG_FS, ha="left", va="bottom", color=color,
                )
                if res is not None:
                    per_model_results.append(dict(
                        model=m, iq=iq, iw=iw,
                        Q2_center=q2c, W_center=wc,
                        **{k: res[k] for k in
                           ("A","dA","b","db","t0","dt0","chi2","ndf","tmin_abs")},
                    ))

            fig.suptitle(bin_title, fontsize=LEG_TITLE_FS + 2, y=0.96)
            stem = (f"dipole_reduced_Q{iq}" if iw is None
                    else f"dipole_reduced_Q{iq}_W{iw}")
            _save(fig, outdir, stem)

        # -- B) combined weighted-average canvas ---------------------------
        avg_df = _weighted_average_reduced_xs(dfs)
        if avg_df is None:
            print(f"  [WARN] Could not build weighted average for Q{iq} W{iw}.")
            continue

        # save weighted-average data as CSV
        avg_csv = os.path.join(
            outdir,
            f"combined_avg_Q{iq}.csv" if iw is None
            else f"combined_avg_Q{iq}_W{iw}.csv"
        )
        avg_df.to_csv(avg_csv, index=False)

        with mpl.rc_context(PANEL_RC):
            fig2 = plt.figure(figsize=(10, 12))
            inner2 = mpl.gridspec.GridSpec(
                2, 1, figure=fig2,
                height_ratios=[3, 1], hspace=0.0,
                left=0.12, right=0.97, top=0.90, bottom=0.09,
            )
            ax_m2 = fig2.add_subplot(inner2[0])
            ax_r2 = fig2.add_subplot(inner2[1], sharex=ax_m2)

            res_comb = _make_combined_fit_canvas(
                ax_m2, ax_r2,
                dfs, m_names, avg_df,
                t_min, t_max, logy, show_ylabel=True,
            )

            # legend: individual models + average + fit
            handles_l, labels_l = ax_m2.get_legend_handles_labels()
            ax_m2.legend(handles_l, labels_l,
                         loc="lower left", frameon=True, framealpha=0.85,
                         edgecolor="#cccccc", fontsize=LEG_FS - 3, ncol=1)

            # annotate number of datasets
            n_ds_min = int(avg_df["n_datasets"].min())
            n_ds_max = int(avg_df["n_datasets"].max())
            ds_str   = (f"{n_ds_min} dataset" + ("s" if n_ds_min != 1 else "")
                        if n_ds_min == n_ds_max
                        else f"{n_ds_min}{n_ds_max} datasets")
            ax_m2.annotate(
                f"Weighted average ({ds_str})",
                xy=(0.03, 0.03), xycoords="axes fraction",
                fontsize=LEG_FS - 2, ha="left", va="bottom", color="#111111",
            )

            fig2.suptitle(bin_title + "   [combined]",
                          fontsize=LEG_TITLE_FS + 2, y=0.96)
            stem2 = (f"dipole_combined_Q{iq}" if iw is None
                     else f"dipole_combined_Q{iq}_W{iw}")
            _save(fig2, outdir, stem2)

        # store combined result row
        if res_comb is not None:
            combined_results.append(dict(
                model="combined", iq=iq, iw=iw,
                Q2_center=q2c, W_center=wc,
                n_datasets=int(avg_df["n_datasets"].mean().round()),
                **{k: res_comb[k] for k in
                   ("A","dA","b","db","t0","dt0","chi2","ndf","tmin_abs")},
            ))

    # -- save CSVs ---------------------------------------------------------
    if not per_model_results and not combined_results:
        print("  [WARN] No dipole fit results  skipping summary plots.")
        return

    if per_model_results:
        sumtab = pd.DataFrame(per_model_results)
        sumtab.to_csv(os.path.join(outdir, "dipole_fit_summary_permodel.csv"),
                      index=False)
        print("\n=== Per-model Dipole Fit Summary ===")
        print(sumtab.to_string(index=False))

    if combined_results:
        comb_tab = pd.DataFrame(combined_results)
        comb_path = os.path.join(outdir, "dipole_fit_summary_combined.csv")
        comb_tab.to_csv(comb_path, index=False)
        print("\n=== Combined Dipole Fit Summary ===")
        print(comb_tab.to_string(index=False))
        print(f"   saved to {comb_path}")
    else:
        comb_tab = pd.DataFrame(per_model_results)   # fallback

    # -- parameter-vs-Q / W plots (per-model) -----------------------------
    if per_model_results:
        d_par = os.path.join(outdir, "fit_parameters")
        sumtab_all = pd.DataFrame(per_model_results)
        _make_param_vs_Q2_plots(sumtab_all, models, d_par)
        _make_param_vs_W_plots( sumtab_all, models, d_par)
        _make_param_summary_per_model(sumtab_all, models, d_par)
        _make_param_correlation_plots(sumtab_all, models, d_par)

    # -- b vs x_B  use combined fit results ----------------------------
    make_b2_vs_xB_plot(comb_tab, ["combined"],
                       os.path.join(outdir, "b2_vs_xB"))


# -- helper: group the summary table by W-bin and plot param vs Q ------------
def _make_param_vs_Q2_plots(sumtab, models, outdir):
    """b, A, t0 vs Q for each W-bin and each model, with error bars."""
    os.makedirs(outdir, exist_ok=True)
    params = [
        ("A",  "dA",  r"$A\ [\mathrm{nb}]$",           "norm_A"),
        ("t0", "dt0", r"$t_0\ [\mathrm{GeV}^2]$",      "t0"),
    ]
    iw_vals = sorted(sumtab["iw"].unique(), key=lambda x: (x is None, x))

    for par, dpar, ylabel, ptag in params:
        for iw_key in iw_vals:
            sub = sumtab[sumtab["iw"] == iw_key].copy()
            if sub.empty:
                continue

            with mpl.rc_context(PANEL_RC):
                fig, ax = plt.subplots(figsize=(9, 6))
                _style_ax(ax)

                for im, m in enumerate(models):
                    ms = sub[sub["model"] == m].sort_values("Q2_center")
                    if ms.empty:
                        continue
                    color = MODEL_COLORS[im % len(MODEL_COLORS)]
                    ax.errorbar(
                        ms["Q2_center"].to_numpy(float),
                        ms[par].to_numpy(float),
                        yerr=ms[dpar].to_numpy(float),
                        fmt="o-", color=color, ecolor=color,
                        elinewidth=1.8, capsize=4,
                        markersize=8, mfc="white", mew=2.2, lw=1.8,
                        label=m.replace("_", " "), zorder=5 + im,
                    )

                ax.set_xlabel(r"$Q^2\ [\mathrm{GeV}^2]$")
                ax.set_ylabel(ylabel)
                w_str = (f"  (W-bin {iw_key})" if iw_key is not None else "")
                ax.set_title(
                    rf"Dipole fit parameter {ylabel.split('[')[0].strip()}"
                    rf" vs $Q^2${w_str}",
                    fontsize=LEG_TITLE_FS,
                )
                ax.legend(frameon=False, fontsize=LEG_FS - 2)

                stem = (f"{ptag}_vs_Q2" if iw_key is None
                        else f"{ptag}_vs_Q2_W{iw_key}")
                _save(fig, outdir, stem)


def _make_param_vs_W_plots(sumtab, models, outdir):
    """b, A, t0 vs W for each Q-bin and each model, with error bars."""
    os.makedirs(outdir, exist_ok=True)
    params = [
        ("A",  "dA",  r"$A\ [\mathrm{nb}]$",           "norm_A"),
        ("t0", "dt0", r"$t_0\ [\mathrm{GeV}^2]$",      "t0"),
    ]
    iq_vals = sorted(sumtab["iq"].unique())

    for par, dpar, ylabel, ptag in params:
        for iq_key in iq_vals:
            sub = sumtab[sumtab["iq"] == iq_key].copy()
            if sub.empty or sub["W_center"].isna().all():
                continue

            with mpl.rc_context(PANEL_RC):
                fig, ax = plt.subplots(figsize=(9, 6))
                _style_ax(ax)

                for im, m in enumerate(models):
                    ms = sub[sub["model"] == m].sort_values("W_center")
                    if ms.empty:
                        continue
                    color = MODEL_COLORS[im % len(MODEL_COLORS)]
                    ax.errorbar(
                        ms["W_center"].to_numpy(float),
                        ms[par].to_numpy(float),
                        yerr=ms[dpar].to_numpy(float),
                        fmt="o-", color=color, ecolor=color,
                        elinewidth=1.8, capsize=4,
                        markersize=8, mfc="white", mew=2.2, lw=1.8,
                        label=m.replace("_", " "), zorder=5 + im,
                    )

                ax.set_xlabel(r"$W\ [\mathrm{GeV}]$")
                ax.set_ylabel(ylabel)
                ax.set_title(
                    rf"Dipole fit parameter {ylabel.split('[')[0].strip()}"
                    rf" vs $W$  (Q-bin {iq_key})",
                    fontsize=LEG_TITLE_FS,
                )
                ax.legend(frameon=False, fontsize=LEG_FS - 2)

                stem = f"{ptag}_vs_W_Q{iq_key}"
                _save(fig, outdir, stem)


def _make_param_summary_per_model(sumtab, models, outdir):
    """3-panel (A, b, t0) summary figure per model  all bins together."""
    os.makedirs(outdir, exist_ok=True)
    params = [
        ("A",  "dA",  r"$A\ [\mathrm{nb}]$"),
        ("t0", "dt0", r"$t_0\ [\mathrm{GeV}^2]$"),
    ]

    for m in models:
        ms = sumtab[sumtab["model"] == m].sort_values(["iw", "Q2_center"])
        if ms.empty:
            continue

        color = MODEL_COLORS[models.index(m) % len(MODEL_COLORS)]

        # Build a human-readable bin label: "Q[lo,hi]  W[lo,hi]"
        bin_labels = []
        for _, row in ms.iterrows():
            lbl = rf"$Q^2$={row['Q2_center']:.2f}"
            if pd.notna(row["W_center"]):
                lbl += rf"  $W$={row['W_center']:.2f}"
            bin_labels.append(lbl)

        x = np.arange(len(ms))

        with mpl.rc_context(PANEL_RC):
            fig, axes = plt.subplots(
                2, 1, figsize=(max(8, 1.2 * len(ms) + 3), 12),
                sharex=True,
            )
            fig.subplots_adjust(hspace=0.08, left=0.12, right=0.97,
                                top=0.93, bottom=0.18)

            for ax, (par, dpar, ylabel) in zip(axes, params):
                _style_ax(ax)
                ax.errorbar(
                    x, ms[par].to_numpy(float),
                    yerr=ms[dpar].to_numpy(float),
                    fmt="o", color=color, ecolor=color,
                    elinewidth=1.8, capsize=4,
                    markersize=9, mfc="white", mew=2.2, lw=0, zorder=5,
                )
                # horizontal reference lines
                mean_val = float(np.nanmean(ms[par]))
                ax.axhline(mean_val, color=color, lw=1.2, ls="--",
                           alpha=0.5, zorder=2,
                           label=rf"mean = {mean_val:.3g}")
                ax.set_ylabel(ylabel, fontsize=20)
                ax.legend(frameon=False, fontsize=LEG_FS - 4, loc="upper right")
                ax.tick_params(axis="both", labelsize=17)

            axes[-1].set_xticks(x)
            axes[-1].set_xticklabels(bin_labels, rotation=35,
                                     ha="right", fontsize=13)
            axes[-1].set_xlabel("Kinematic bin", fontsize=20)

            # /ndf on top panel
            chi2_str = "  ".join(
                rf"$\chi^2/\nu={row['chi2']:.1f}/{row['ndf']}$"
                for _, row in ms.iterrows()
                if row["ndf"] > 0
            )
            if chi2_str:
                axes[0].annotate(
                    chi2_str,
                    xy=(0.02, 0.02), xycoords="axes fraction",
                    fontsize=11, va="bottom", color="gray",
                )

            fig.suptitle(
                rf"Dipole fit parameters  {m.replace('_', ' ')}",
                fontsize=LEG_TITLE_FS + 2, y=0.96,
            )
            _save(fig, outdir, f"param_summary_{_sanitise(m)}")


# -----------------------------------------------------------------------------
#  Dipole fit parameter correlation plots
# -----------------------------------------------------------------------------
def _make_param_correlation_plots(sumtab, models, outdir):
    """
    2-D scatter / correlation plots for every pair of dipole fit parameters
    (A, t0) with 1- error ellipses.

    Pairs plotted: (A, t0)

    One figure is produced that overlays all models (colour-coded) and all
    kinematic bins (marker shape cycles over Q-bins).

    Error ellipses are drawn using the marginal uncertainties (dA, dt0) as
    semi-axes  they represent the 1 box projected as an ellipse; for the
    true joint confidence region you would need the full covariance matrix,
    but the marginal representation is informative for quick comparison.
    """
    import matplotlib.patches as mpatches
    from matplotlib.patches import Ellipse

    os.makedirs(outdir, exist_ok=True)

    # The pure-dipole model has only 2 free params: A and t0.
    # (b is NaN in the pure-dipole case  silently skip it.)
    param_pairs = [
        ("A",  "dA",  r"$A\ [\mathrm{nb}]$",
         "t0", "dt0", r"$t_0\ [\mathrm{GeV}^2]$",
         "corr_A_vs_t0"),
    ]

    iq_vals = sorted(sumtab["iq"].unique())
    markers = ["o", "s", "^", "D", "v", "P", "*", "h"]

    for (px, dpx, xlbl, py, dpy, ylbl, stem) in param_pairs:

        # Skip if both columns are all-NaN (e.g. b when using pure dipole)
        if sumtab[px].isna().all() or sumtab[py].isna().all():
            continue

        with mpl.rc_context(PANEL_RC):
            fig, ax = plt.subplots(figsize=(9, 7))
            _style_ax(ax)

            for im, m in enumerate(models):
                ms = sumtab[sumtab["model"] == m]
                if ms.empty:
                    continue
                color = MODEL_COLORS[im % len(MODEL_COLORS)]

                for iq_idx, iq_key in enumerate(iq_vals):
                    row = ms[ms["iq"] == iq_key]
                    if row.empty:
                        continue

                    mkr = markers[iq_idx % len(markers)]

                    for _, r in row.iterrows():
                        xv  = float(r[px]);   xe  = float(r[dpx])
                        yv  = float(r[py]);   ye  = float(r[dpy])
                        if not (np.isfinite(xv) and np.isfinite(yv)):
                            continue

                        # central marker
                        ax.plot(xv, yv, mkr,
                                color=color, mfc="white", mew=2.0,
                                markersize=9, zorder=6)
                        # error bars
                        ax.errorbar(xv, yv,
                                    xerr=xe if np.isfinite(xe) else None,
                                    yerr=ye if np.isfinite(ye) else None,
                                    fmt="none", color=color, ecolor=color,
                                    elinewidth=1.6, capsize=4, zorder=5)
                        # 1- error ellipse
                        if np.isfinite(xe) and np.isfinite(ye) and xe > 0 and ye > 0:
                            ell = Ellipse(
                                xy=(xv, yv),
                                width=2 * xe, height=2 * ye,
                                angle=0,
                                edgecolor=color, facecolor=color,
                                alpha=0.12, lw=1.2, zorder=4,
                            )
                            ax.add_patch(ell)

            # -- legend: models (colour) + Q-bins (shape) --------------
            legend_handles = []
            for im, m in enumerate(models):
                legend_handles.append(
                    mlines.Line2D([], [], color=MODEL_COLORS[im % len(MODEL_COLORS)],
                                  marker="o", linestyle="None",
                                  markersize=9, mfc="white", mew=2.0,
                                  label=m.replace("_", " "))
                )
            for iq_idx, iq_key in enumerate(iq_vals):
                mkr = markers[iq_idx % len(markers)]
                q2c_vals = sumtab[sumtab["iq"] == iq_key]["Q2_center"].dropna()
                q2_lbl   = (rf"$Q^2 \approx {q2c_vals.mean():.2f}\ \mathrm{{GeV}}^2$"
                            if not q2c_vals.empty else rf"$Q^2$-bin {iq_key}")
                legend_handles.append(
                    mlines.Line2D([], [], color="gray",
                                  marker=mkr, linestyle="None",
                                  markersize=9, mfc="white", mew=1.8,
                                  label=q2_lbl)
                )

            ax.legend(handles=legend_handles,
                      frameon=True, framealpha=0.90, edgecolor="#cccccc",
                      fontsize=LEG_FS - 3, loc="best")

            ax.set_xlabel(xlbl, fontsize=22)
            ax.set_ylabel(ylbl, fontsize=22)
            ax.set_title(
                rf"Dipole fit parameter correlation: {xlbl.split('[')[0].strip()} vs "
                rf"{ylbl.split('[')[0].strip()}"
                "\n(error bars = 1 marginal; ellipse = 1 marginal area)",
                fontsize=14, pad=8,
            )
            _save(fig, outdir, stem)


# -----------------------------------------------------------------------------
#  <b>_g  vs  x_B  plot  (from dipole fit t0 parameter)
# -----------------------------------------------------------------------------
# Conversion factor:  1 GeV = 0.389379 fm
_GEVMINUSSQ_TO_FM2 = (0.197326 ** 2)

# Literature points digitized from PR12-12-007 (proposal) Fig. (gluonic radius)
# Columns: xB, <b^2>_g [fm^2]. Used for overlay comparison only.
LITERATURE_B2_FROM_PROPOSAL = {
    # Digitized from PR12-12-007 proposal Fig. 6 (world data overlay).
    # Units: xB (dimensionless), <b^2>_g in fm^2.
    # NOTE: Published uncertainties are not all explicitly tabulated in the proposal figure.
    # The db2 values below are conservative placeholders; replace with the exact values
    # if/when you have them in a table.
    "FNAL 82": {
        "xB":     [0.04944],
        "b2_fm2": [0.22945],
        "db2_fm2":[0.1010],
        "dxB_lo": [0.03229],   # optional, if you decide to draw x-errors
        "dxB_hi": [0.03612],
        "color": "C1",
        "marker": "s",
    },
    "H1 05": {
        "xB":     [0.0001492937173643, 0.0002720869272003, 0.0018121176816118],
        "b2_fm2": [0.3908317494795909, 0.3705470572984184, 0.328262339418526],
        "db2_fm2":[0.030, 0.030, 0.030],  # placeholder
        "color":  "C0",
        "marker": "o",
    },
    "ZEUS 02 e+e-": {
        "xB":     [0.0001389289146607, 0.000229400884545, 0.0004013424777244, 0.0007171492515423, 0.0016662323472533, 0.0062681039747397],
        "b2_fm2": [0.2920949902407286, 0.3363805076628352, 0.3036484566705411, 0.317098972637575, 0.2802185960591133, 0.2499575496342738],
        "db2_fm2":[0.030, 0.030, 0.030, 0.030, 0.030, 0.030],  # placeholder
        "color":  "C2",
        "marker": "^",
    },
    "ZEUS 02 mu+mu-": {
        "xB":     [0.0003457313937778, 0.0004419840398351, 0.0006016795339795, 0.000892858329558, 0.003671795459816],
        "b2_fm2": [0.3439308363059849, 0.3338610568929371, 0.2981057471264367, 0.3008911274699729, 0.2792016154085119],
        "db2_fm2":[0.030, 0.030, 0.030, 0.030, 0.030],  # placeholder
        "color":  "C4",
        "marker": "v",
    },
}



def _b2_from_params(mg2_fit, dmg2_fit, tmin_abs=0.0):
    """Compute <b^2>_g (fm^2) from fitted m_g^2 and |t_min|.

    The fit is performed in |t| = t' + tmin_abs.  Re-expressing the dipole
    σ₀/(1 + |t|/m_g²)⁴ as a function of t' near t'=0 gives effective slope
    scale m_eff² = m_g² + tmin_abs.  The physically correct formula is therefore:

        <b²>_g [GeV⁻²] = 8 / (m_g² + tmin_abs)

    Uncertainty (tmin_abs is treated as exact kinematics, no uncertainty):

        δ<b²>_g = 8 × δm_g² / (m_g² + tmin_abs)²

    Conversion: 1 GeV⁻² = (0.197326 fm)² ≈ 0.038938 fm².
    """
    mg2_fit  = float(mg2_fit)
    tmin_abs = float(tmin_abs)
    m_eff2   = mg2_fit + tmin_abs
    if not np.isfinite(m_eff2) or m_eff2 <= 0:
        return np.nan, np.nan
    b2_gev  = 8.0 / m_eff2
    db2_gev = 8.0 * float(dmg2_fit) / (m_eff2 ** 2) if np.isfinite(dmg2_fit) else np.nan
    b2  = b2_gev  * _GEVMINUSSQ_TO_FM2
    db2 = db2_gev * _GEVMINUSSQ_TO_FM2 if np.isfinite(db2_gev) else np.nan
    return b2, db2

def make_b2_vs_xB_plot(sumtab, models, outdir):
    """
    Plot <b^2>_g (fm^2) vs x_B derived from the dipole fit t0 parameter (CLAS12),
    overlay literature points from the proposal, and draw ONE fit line to the
    world (literature) data only (NOT CLAS12). The fit is extrapolated to xB=1.

    Fit form (proposal-style straight-line in ln(1/xB)):
        <b^2> = A + B * ln(1/xB)
    """
    os.makedirs(outdir, exist_ok=True)

    # ---- CLAS12 points: compute xB and <b^2> from dipole fit parameters ----
    tab = sumtab.copy()

    # x_B (keep consistent with the |t_min| convention used elsewhere in this script)
    Q2 = tab["Q2_center"].to_numpy(float)
    W_used = float(W_MEAN_GEV)
    with np.errstate(invalid="ignore"):
        tab["xB"] = Q2 / ((W_used * W_used) - M_P_GEV**2 + Q2)

    # <b^2>_g = 8/(m_g^2 + |t_min|) * (ħc)^2  [fm^2]
    b2_vals  = np.full(len(tab), np.nan)
    db2_vals = np.full(len(tab), np.nan)
    for i, row in tab.iterrows():
        mg2_fit  = float(row.get("t0", np.nan))
        dmg2_fit = float(row.get("dt0", np.nan))

        if "tmin_abs" in tab.columns and np.isfinite(row.get("tmin_abs", np.nan)):
            tmin_row = float(row["tmin_abs"])
        else:
            q2_row = float(row.get("Q2_center", np.nan))
            tmin_row = float(t_min_phi(q2_row, W_MEAN_GEV)) if np.isfinite(q2_row) else 0.0
            if not np.isfinite(tmin_row):
                tmin_row = 0.0

        if np.isfinite(mg2_fit) and mg2_fit > 0:
            b2, db2 = _b2_from_params(mg2_fit, dmg2_fit, tmin_abs=tmin_row)
            b2_vals[i]  = b2
            db2_vals[i] = db2

    tab["b2"]  = b2_vals
    tab["db2"] = db2_vals
    tab = tab[np.isfinite(tab["xB"]) & (tab["xB"] > 0) & np.isfinite(tab["b2"])].copy()

    # ---- marker cycle for multiple W-bins (CLAS12) ----
    iw_vals = sorted(tab["iw"].unique(), key=lambda x: (x is None, x)) if (not tab.empty and "iw" in tab.columns) else [None]
    markers = ["o", "s", "^", "D", "v", "P", "*"]

    with mpl.rc_context(PANEL_RC):
        fig, ax = plt.subplots(figsize=(9, 6.5))
        _style_ax(ax)

        # -----------------------------
        # 1) Literature points (colored + y error bars)
        # -----------------------------
        lit_handles, lit_labels = [], []
        fit_x, fit_y = [], []

        for exp, d in LITERATURE_B2_FROM_PROPOSAL.items():
            xb_lit  = np.asarray(d.get("xB", []), dtype=float)
            b2_lit  = np.asarray(d.get("b2_fm2", []), dtype=float)
            db2_lit = d.get("db2_fm2", None)
            db2_lit = np.asarray(db2_lit, dtype=float) if db2_lit is not None else None

            msk = np.isfinite(xb_lit) & (xb_lit > 0) & np.isfinite(b2_lit)
            if db2_lit is not None:
                msk = msk & np.isfinite(db2_lit)

            if not np.any(msk):
                continue

            # collect for world-data fit (NOT CLAS12)
            fit_x.extend(xb_lit[msk].tolist())
            fit_y.extend(b2_lit[msk].tolist())

            col = d.get("color", "black")
            mkr = d.get("marker", "s")

            h = ax.errorbar(
                xb_lit[msk], b2_lit[msk],
                yerr=(db2_lit[msk] if db2_lit is not None else None),
                fmt=mkr, linestyle="None",
                color=col, ecolor=col,
                elinewidth=1.6, capsize=4,
                markersize=8, mfc="white", mew=2.0,
                alpha=0.95, zorder=20,
                label=exp,
            )
            lit_handles.append(h)
            lit_labels.append(exp)

        # -----------------------------
        # 2) One straight-line fit to WORLD data only:  A + B ln(1/xB)
        # -----------------------------
        h_fit = None
        fit_label = None
        fit_x = np.asarray(fit_x, dtype=float)
        fit_y = np.asarray(fit_y, dtype=float)

        fit_mask = np.isfinite(fit_x) & (fit_x > 0) & np.isfinite(fit_y)
        fit_x = fit_x[fit_mask]
        fit_y = fit_y[fit_mask]

        if fit_x.size >= 2:
            X = np.log(1.0 / fit_x)
            # y = A + B*X  => polyfit returns [B, A]
            B, A = np.polyfit(X, fit_y, 1)

            print(f"[b2_vs_xB] world fit: <b^2> = {A:.6f} + {B:.6f} * ln(1/xB)")

            x_curve = np.logspace(np.log10(max(np.min(fit_x) * 0.8, 1e-6)), 0.0, 500)  # up to xB=1
            y_curve = A + B * np.log(1.0 / x_curve)

            h_fit, = ax.plot(
                x_curve, y_curve,
                ls="--", lw=2.4, color="black", zorder=18,
                label=r"$\langle b^2\rangle = A + B\,\ln(1/x_B)$ (world fit)",
            )
            fit_label = r"$\langle b^2\rangle = A + B\,\ln(1/x_B)$ fit"

        # -----------------------------
        # 3) CLAS12 points (NOT used in the fit)
        # -----------------------------
        clas_handles, clas_labels = [], []
        if not tab.empty:
            for im, m in enumerate(models):
                color = MODEL_COLORS[im % len(MODEL_COLORS)]
                ms = tab[tab["model"] == m]
                if ms.empty:
                    continue

                for iw_idx, iw_key in enumerate(iw_vals):
                    mw = ms if iw_key is None else ms[ms["iw"] == iw_key]
                    mw = mw.sort_values("xB")
                    if mw.empty:
                        continue

                    mkr = markers[iw_idx % len(markers)]
                    lbl = m.replace("_", " ")
                    if len(iw_vals) > 1:
                        lbl += f"  W-bin {iw_key}" if iw_key is not None else "  (all W)"

                    h = ax.errorbar(
                        mw["xB"].to_numpy(float),
                        mw["b2"].to_numpy(float),
                        yerr=mw["db2"].to_numpy(float),
                        fmt=mkr, color=color, ecolor=color,
                        elinewidth=2.0, capsize=5,
                        markersize=9, mfc="white", mew=2.2, lw=1.5,
                        zorder=5 + im * 10 + iw_idx,
                        label=lbl,
                    )
                    clas_handles.append(h)
                    clas_labels.append(lbl)

                    ax.plot(
                        mw["xB"].to_numpy(float),
                        mw["b2"].to_numpy(float),
                        color=color, lw=1.2, ls="-", alpha=0.5, zorder=4,
                    )

        # -----------------------------
        # 4) Axis + decorations
        # -----------------------------
        ax.set_xlabel(r"$x_B$", fontsize=22)
        ax.set_ylabel(r"$\langle b^2 \rangle_g \ [\mathrm{fm}^2]$", fontsize=22)
        ax.set_title(
            r"Mean-Square Gluonic Transverse Radius $\langle b^2\rangle_g$ vs $x_B$"
            "\n" r"(from dipole fit: $\langle b^2\rangle_g = 8/(t_0^\mathrm{fit} + |t_\mathrm{min}|) \times 0.389\ \mathrm{fm}^2$)",
            fontsize=15, pad=8,
        )

        # Use log-x if wide range OR if literature included
        all_x = []
        if not tab.empty:
            all_x.append(tab["xB"].to_numpy(float))
        if fit_x.size > 0:
            all_x.append(fit_x)
        all_x = np.concatenate(all_x) if len(all_x) else np.array([0.1, 1.0])
        all_x = all_x[np.isfinite(all_x) & (all_x > 0)]
        xb_range = (np.nanmax(all_x) / max(np.nanmin(all_x), 1e-6)) if all_x.size else 1.0
        if xb_range > 5:
            ax.set_xscale("log")

        # proton charge radius reference line
        r_Ep_fm = 0.8414
        h_ref = ax.axhline(
            r_Ep_fm**2,
            color="#008800",
            lw=1.5,
            ls="--",
            zorder=1,
            label=rf"$r_{{Ep}}^2 = {r_Ep_fm**2:.3f}\ \mathrm{{fm}}^2$",
        )

        ax.set_ylim(bottom=0.0, top= 0.85)
        ax.set_xlim(left=min(1e-5, float(np.nanmin(all_x))*0.8) if all_x.size else 1e-5, right=1.0)

        # -----------------------------
        # 5) Legend: upper-left (includes ONE fit line)
        # -----------------------------
        leg_h, leg_l = [], []

        # CLAS12 series first (if present)
        leg_h += clas_handles
        leg_l += clas_labels

        # literature series
        leg_h += lit_handles
        leg_l += lit_labels

        # fit line (one only)
        if h_fit is not None and fit_label is not None:
            leg_h += [h_fit]
            leg_l += [fit_label]

        # reference line
        leg_h += [h_ref]
        leg_l += [rf"$r_{{Ep}}^2 = {r_Ep_fm**2:.3f}\ \mathrm{{fm}}^2$"]

        ax.legend(
            leg_h, leg_l,
            loc="upper left",
            frameon=True,
            framealpha=0.92,
            edgecolor="#cccccc",
            fontsize=LEG_FS - 2,
        )

        _save(fig, outdir, "b2_vs_xB")

    # Save CLAS12 derived table (literature is separate)
    if not tab.empty:
        out_cols = ["model", "iq", "iw", "Q2_center", "W_center",
                    "xB", "b2", "db2", "t0", "dt0", "tmin_abs", "chi2", "ndf"]
        out_cols = [c for c in out_cols if c in tab.columns]
        tab[out_cols].to_csv(os.path.join(outdir, "b2_vs_xB_table.csv"), index=False)
        print("\n=== <b^2> vs x_B table (CLAS12 only) ===")
        print(tab[out_cols].to_string(index=False))


# -----------------------------------------------------------------------------
#  HS summary plots: Ds(0) vs Q2 and vs xB
# -----------------------------------------------------------------------------
def _compute_xB_from_Q2(Q2, W_used=W_MEAN_GEV, M=M_P_GEV):
    """Compute Bjorken x_B from Q2 and a chosen W.

    This script's other summary plot (b2 vs xB) uses W_MEAN_GEV to stay
    consistent with the |t_min| convention. We do the same here.
    """
    Q2 = np.asarray(Q2, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return Q2 / ((W_used * W_used) - (M * M) + Q2)


def make_hs_ds0_summary_plots(ds0_table, outdir, prefer_combined=True):
    """Make final summary plots for the extracted strange D-term Ds(0).

    Produces:
      - hs_Ds0_vs_Q2.(pdf|png)
      - hs_Ds0_vs_xB.(pdf|png)
      - hs_Ds0_with_xB.csv

    Parameters
    ----------
    ds0_table : pd.DataFrame
        Table produced by the HS extraction step.
        Must contain columns: Ds0, dDs0, Q2_center.
    outdir : str
        Output directory (typically <outdir>/HS_fits)
    prefer_combined : bool
        If True, plot the rows with model == 'combined' when available.
    """
    if ds0_table is None or len(ds0_table) == 0:
        return

    os.makedirs(outdir, exist_ok=True)

    tab = ds0_table.copy()
    if "Q2_center" not in tab.columns or "Ds0" not in tab.columns:
        return

    # Choose which rows represent the "final" result.
    if prefer_combined and ("model" in tab.columns) and (tab["model"] == "combined").any():
        tab_plot = tab[tab["model"] == "combined"].copy()
    else:
        tab_plot = tab.copy()

    # Compute xB (consistent with other script summary: uses W_MEAN_GEV)
    tab_plot["xB"] = _compute_xB_from_Q2(tab_plot["Q2_center"].to_numpy(float))

    # Save a convenience table with xB
    out_csv = os.path.join(outdir, "hs_Ds0_with_xB.csv")
    tab_plot.to_csv(out_csv, index=False)

    # Marker cycle for multiple W-bins
    iw_vals = sorted(tab_plot["iw"].unique(), key=lambda x: (x is None, x)) if "iw" in tab_plot.columns else [None]
    markers = ["o", "s", "^", "D", "v", "P", "*", "h"]

    def _plot_one(xcol, xlabel, stem, xscale=None):
        # Pre-compute lattice x-range for the band
        xcol_vals = tab_plot[xcol].to_numpy(float)
        xcol_finite = xcol_vals[np.isfinite(xcol_vals) & (xcol_vals > 0)]
        if xcol_finite.size > 0:
            xlim_lo = float(np.nanmin(xcol_finite)) * 0.7
            xlim_hi = float(np.nanmax(xcol_finite)) * 1.3
        else:
            xlim_lo, xlim_hi = 0.0, 1.0
        xrng = np.array([xlim_lo, xlim_hi])

        for show_lattice in (False, True):
            with mpl.rc_context(PANEL_RC):
                fig, ax = plt.subplots(figsize=(9, 6))
                _style_ax(ax)

                for iw_idx, iw_key in enumerate(iw_vals):
                    sub = tab_plot if iw_key is None else tab_plot[tab_plot["iw"] == iw_key]
                    if sub.empty:
                        continue
                    sub = sub.sort_values(xcol)
                    x = sub[xcol].to_numpy(float)
                    y = sub["Ds0"].to_numpy(float)
                    ye = sub["dDs0"].to_numpy(float) if "dDs0" in sub.columns else None
                    mkr = markers[iw_idx % len(markers)]
                    lbl = ("combined" if prefer_combined else "all")
                    if len(iw_vals) > 1:
                        lbl += f"  W-bin {iw_key}" if iw_key is not None else "  (all W)"

                    mask = np.isfinite(x) & np.isfinite(y)
                    if ye is not None:
                        mask = mask & np.isfinite(ye)
                    if not np.any(mask):
                        continue

                    ax.errorbar(
                        x[mask], y[mask],
                        yerr=(ye[mask] if ye is not None else None),
                        fmt=mkr,
                        color="#111111", ecolor="#111111",
                        elinewidth=2.0, capsize=5,
                        markersize=9, mfc="white", mew=2.2,
                        lw=1.2, zorder=5 + iw_idx,
                        label=lbl,
                    )
                    ax.plot(x[mask], y[mask], color="#111111", lw=1.0, alpha=0.45, zorder=3)

                if show_lattice:
                    # Lattice QCD reference band (Q^2-independent: it's a ground-state property)
                    ax.fill_between(xrng,
                                    [LAT_Ds0_DIPOLE - LAT_dDs0_DIPOLE] * 2,
                                    [LAT_Ds0_DIPOLE + LAT_dDs0_DIPOLE] * 2,
                                    color="steelblue", alpha=0.18, zorder=0,
                                    label=rf"Lat. dipole $D_s(0)={LAT_Ds0_DIPOLE:.2f}\pm{LAT_dDs0_DIPOLE:.2f}$")
                    ax.axhline(LAT_Ds0_DIPOLE, color="steelblue", lw=1.6, ls="-", alpha=0.85, zorder=1)

                    ax.fill_between(xrng,
                                    [LAT_Ds0_ZEXP - LAT_dDs0_ZEXP] * 2,
                                    [LAT_Ds0_ZEXP + LAT_dDs0_ZEXP] * 2,
                                    color="tomato", alpha=0.12, zorder=0,
                                    label=rf"Lat. z-exp $D_s(0)={LAT_Ds0_ZEXP:.2f}\pm{LAT_dDs0_ZEXP:.2f}$")
                    ax.axhline(LAT_Ds0_ZEXP, color="tomato", lw=1.6, ls="--", alpha=0.85, zorder=1)

                ax.axhline(0.0, color="gray", lw=1.2, ls="--", zorder=1)
                ax.set_xlabel(xlabel)
                ax.set_ylabel(r"$D_s(0)$")
                ax.set_title(r"Strange D-term extraction (Hatta–Strikman Fit)")
                if xscale:
                    ax.set_xscale(xscale)
                ax.legend(frameon=False, fontsize=LEG_FS - 2, loc="center left")

                #if show_lattice:
                 #   ax.text(0.98, 0.97,
                  #          "Lattice: Hackett, Pefkou, Shanahan\nPRL 132, 251904 (2024), Table I",
                   #         transform=ax.transAxes, ha="right", va="top",
                    #        fontsize=max(LEG_FS - 4, 8), style="italic", color="gray")

                out_stem = stem + ("_with_lattice" if show_lattice else "")
                _save(fig, outdir, out_stem)

    _plot_one("Q2_center", r"$Q^2\ [\mathrm{GeV}^2]$", "hs_Ds0_vs_Q2")

    # xB range can span orders of magnitude in some datasets; use log if wide.
    xb_arr = tab_plot["xB"].to_numpy(float)
    xb_mask = np.isfinite(xb_arr) & (xb_arr > 0)
    xb_scale = None
    if np.any(xb_mask):
        xb_min = float(np.nanmin(xb_arr[xb_mask]))
        xb_max = float(np.nanmax(xb_arr[xb_mask]))
        if xb_min > 0 and (xb_max / xb_min) > 5:
            xb_scale = "log"
    _plot_one("xB", r"$x_B$", "hs_Ds0_vs_xB", xscale=xb_scale)


# -----------------------------------------------------------------------------
#  main
# -----------------------------------------------------------------------------

def _run_hs_ds0_workflow(*, csv_root, models, bin_keys, outdir,
                         tprime_min, tprime_max, As0, mA, mD, Ds0_bounds, logy=True):
    """Run HS Ds(0) extraction on *reduced* cross-section and produce fit plots + summary plots.

    This is called both from the full pipeline and from --fits-only mode.
    """
    print("\n=== HS Ds(0) extraction (template fit on reduced XS) ===")

    hs_outdir = os.path.join(outdir, "HS_fits")
    os.makedirs(hs_outdir, exist_ok=True)

    ds_rows = []

    # Per-model Ds0 fits
    for (iq, iw) in bin_keys:
        for m in models:
            df = load_csv(csv_root, m, iq, iw)
            if df is None:
                continue

            res = _fit_hs_Ds0_from_dsdt(
                df, m,
                tprime_min=tprime_min,
                tprime_max=tprime_max,
                As0=As0, mA=mA, mD=mD,
                Ds0_bounds=tuple(Ds0_bounds),
            )
            if res is None:
                continue

            _make_hs_fit_plot(df, res, iq, iw, m, hs_outdir, logy=logy, As0=As0, mA=mA, mD=mD)
            _make_hs_gff_plot(df, res, iq, iw, m, hs_outdir,
                              As0=As0, mA=mA, mD=mD)

            q2c = float(df["Q2_center"].iloc[0]) if "Q2_center" in df.columns else np.nan
            wc  = float(df["W_center"].iloc[0])  if "W_center"  in df.columns else np.nan

            ds_rows.append(dict(
                model=m, iq=iq, iw=iw,
                Q2_center=q2c, W_center=wc,
                Ds0=res["Ds0"], dDs0=res["dDs0"],
                N=res["N"], dN=res["dN"],
                chi2=res["chi2"], ndf=res["ndf"],
                tmin_abs=res["tmin_abs"],
            ))

    # Combined (weighted-average) Ds0 fits (uses inverse-variance weighted reduced XS)
    for (iq, iw) in bin_keys:
        dfs = []
        for m in models:
            df = load_csv(csv_root, m, iq, iw)
            if df is not None:
                dfs.append(df)

        if len(dfs) < 2:
            continue

        avg_df = _weighted_average_reduced_xs(dfs)
        if avg_df is None or len(avg_df) == 0:
            continue

        comb_csv = os.path.join(
            hs_outdir,
            f"combined_red_avg_Q{iq}.csv" if iw is None else f"combined_red_avg_Q{iq}_W{iw}.csv"
        )
        avg_df.to_csv(comb_csv, index=False)

        res_c = _fit_hs_Ds0_from_dsdt(
            avg_df, "combined",
            tprime_min=tprime_min,
            tprime_max=tprime_max,
            As0=As0, mA=mA, mD=mD,
            Ds0_bounds=tuple(Ds0_bounds),
        )
        if res_c is None:
            continue

        _make_hs_fit_plot(avg_df, res_c, iq, iw, "combined", hs_outdir, logy=logy, As0=As0, mA=mA, mD=mD)
        _make_hs_gff_plot(avg_df, res_c, iq, iw, 'combined', hs_outdir,
                          As0=As0, mA=mA, mD=mD)

        q2c = float(avg_df["Q2_center"].iloc[0]) if "Q2_center" in avg_df.columns else np.nan
        wc  = float(avg_df["W_center"].iloc[0])  if "W_center"  in avg_df.columns else np.nan

        ds_rows.append(dict(
            model="combined", iq=iq, iw=iw,
            Q2_center=q2c, W_center=wc,
            Ds0=res_c["Ds0"], dDs0=res_c["dDs0"],
            N=res_c["N"], dN=res_c["dN"],
            chi2=res_c["chi2"], ndf=res_c["ndf"],
            tmin_abs=res_c["tmin_abs"],
            n_datasets=int(np.nanmin(avg_df["n_datasets"])) if "n_datasets" in avg_df.columns else np.nan,
            combined_spectrum=comb_csv,
        ))

    if ds_rows:
        out = pd.DataFrame(ds_rows)
        out = out.sort_values(["model", "iw", "iq"], key=lambda s: s.map(lambda x: (x is None, x)))
        out_csv = os.path.join(hs_outdir, "hs_Ds0_summary.csv")
        out.to_csv(out_csv, index=False)
        print("\n=== HS Ds(0) Summary (including combined) ===")
        print(out.to_string(index=False))
        print(f"\n[OK] wrote {out_csv}")
        make_hs_ds0_summary_plots(out, hs_outdir, prefer_combined=True)
    else:
        print("  [WARN] no successful Ds(0) fits.")

def main():
    ap = argparse.ArgumentParser(
        description=" d/dt analysis: corrections, acceptance, efficiency, RadCorr, "
                    "final cross-section, dipole fits."
    )
    ap.add_argument("--csv-root",  default=".",
                    help="Root dir containing CSVs/ sub-tree.")
    ap.add_argument("--models",    nargs="+", required=True,
                    help='Model names, e.g. "Sp18_inb" "Sp19_inb"')
    ap.add_argument("--outdir",    default="phi_results",
                    help="Top-level output directory.")
    ap.add_argument("--luminosity", type=float, default=None,
                    help="Integrated luminosity [nb] for variant reconstruction "
                         "from RawCounts. If omitted, scaling from CrossSection used.")
    ap.add_argument("--branching", type=float, default=0.492,
                    help="KK branching ratio.")
    ap.add_argument("--beam-energy", type=float, default=None,
                    help="Fallback beam energy E [GeV] for model names not matched "
                         "by auto-detection. Auto: Sp1810.6, Fall1810.594, Sp1910.2.")
    ap.add_argument("--plot-reduced", action="store_true",
                    help="Also produce the group-overview reduced cross-section plots "
                         "(dipole fits on reduced XS are always produced).")
    ap.add_argument("--plot-raw-yields", action="store_true",
                    help="Produce a group-overview plot of raw (uncorrected) yields "
                         "N_raw vs t' for all Q² bins. Saved to <outdir>/0a_RawYields_group. "
                         "Use this for presentations: show raw data first, then corrections, "
                         "then the final reduced cross-sections.")
    ap.add_argument("--fits-only", action="store_true",
                    help="Skip all correction/overview plots and run only the dipole and HS fit stages using the existing CSVs.")
    ap.add_argument("--use-external-radcorr", action="store_true",
                    help="Overwrite RadCorr using DIFFRAD CSV tables.")
    ap.add_argument("--radcorr-file", default=None,
                    help="Override external RadCorr CSV path (skips auto-selection).")
    ap.add_argument("--print-table", action="store_true",
                    help="Print + save cross-check tables per (Q,W) bin.")
    ap.add_argument("--table-outdir", default=None,
                    help="Where to save tables (default: <outdir>/0_Tables).")
    ap.add_argument("--no-logy", action="store_true",
                    help="Disable log-y on cross-section canvases.")
    ap.add_argument("--t-min", type=float, default=0.2,
                    help="Min tprime for fitting [GeV^2].")
    ap.add_argument("--t-max", type=float, default=3.0,
                    help="Max |t| for fitting [GeV].")

    ap.add_argument("--extract-ds0", action="store_true",
                    help="Extract strange D-term Ds(0) using HS-inspired template fit to final dσ/dt'.")
    ap.add_argument("--ds0-tprime-min", type=float, default=0.20,
                    help="Min t' for Ds0 fit [GeV^2].")
    ap.add_argument("--ds0-tprime-max", type=float, default=3.0,
                    help="Max t' for Ds0 fit [GeV^2].")
    ap.add_argument("--hs-As0", type=float, default=HS_AS0_DEFAULT,
                    help="HS As(0) (default 0.04).")
    ap.add_argument("--hs-mA", type=float, default=HS_mA_DEFAULT,
                    help="HS mA [GeV] (default 1.13).")
    ap.add_argument("--hs-mD", type=float, default=HS_mD_DEFAULT,
                    help="HS mD [GeV] (default 0.76).")
    ap.add_argument("--ds0-bounds", nargs=2, type=float, default=[-5.0, 2.0],
                    metavar=("DS0_MIN", "DS0_MAX"),
                    help="Bounds for Ds(0) fit (default: -5 to 2).")
    args = ap.parse_args()

    global G_BEAM_ENERGY_GEV, G_USE_EXTERNAL_RADCORR, G_RADCORR_FILE_OVERRIDE
    G_BEAM_ENERGY_GEV       = args.beam_energy          # fallback only
    G_USE_EXTERNAL_RADCORR  = bool(args.use_external_radcorr)
    G_RADCORR_FILE_OVERRIDE = args.radcorr_file

    csv_root     = args.csv_root
    models       = args.models
    outdir       = args.outdir
    logy         = not args.no_logy
    lumi         = args.luminosity
    br           = args.branching
    plot_red     = args.plot_reduced
    plot_raw     = args.plot_raw_yields
    table_outdir = args.table_outdir or os.path.join(outdir, "0_Tables")

    print("Beam energy assignment:")
    for m in models:
        E = beam_energy_for_model(m)
        print(f"  {m:20s}  {E} GeV  "
              f"[RadCorr file: {os.path.basename(choose_radcorr_file_for_energy(E)) if E else 'N/A'}]")

    present = [
        m for m in models
        if os.path.isdir(os.path.join(csv_root, "CSVs", _sanitise(m)))
        and glob(os.path.join(csv_root, "CSVs", _sanitise(m), "dsdt_Q*.csv"))
    ]
    if not present:
        print("[ERROR] No models found.")
        return

    bin_keys = discover_bin_keys(csv_root, present)

    # -----------------------------------------------------------------
    #  Fast path: fits only (dipole + optional HS Ds(0))
    # -----------------------------------------------------------------
    if args.fits_only:
        make_dipole_fit_plots(
            csv_root, present, bin_keys,
            os.path.join(outdir, "7_Dipolefit_ReducedXS"),
            t_min=args.t_min, t_max=args.t_max, logy=logy,
        )

        if args.extract_ds0:
            _run_hs_ds0_workflow(
                csv_root=csv_root, models=present, bin_keys=bin_keys, outdir=outdir,
                tprime_min=args.ds0_tprime_min, tprime_max=args.ds0_tprime_max,
                As0=args.hs_As0, mA=args.hs_mA, mD=args.hs_mD,
                Ds0_bounds=tuple(args.ds0_bounds), logy=logy,
            )

        print(f"\n[DONE] Outputs in {outdir}/")
        return

    d1 = os.path.join(outdir, "1_Corrections_beforeafter")
    # ------------------------------------------------------------------
    #  0a) Raw yields overview  (show this first in presentations)
    # ------------------------------------------------------------------
    if plot_raw:
        make_raw_yields_per_bin(csv_root, present, bin_keys,
                                os.path.join(outdir, "0a_RawYields"),
                                logy=logy, luminosity=lumi, branching=br)

    for (iq, iw) in bin_keys:
        make_corrections_canvas(
            csv_root, present, iq, iw, d1, logy,
            luminosity=lumi, branching=br,
            print_table=args.print_table, table_outdir=table_outdir,
        )

    make_group_plot(csv_root, present, bin_keys, "acc",
                    os.path.join(outdir, "2_Acceptance_group"), logy=False)

    make_group_plot(csv_root, present, bin_keys, "eff",
                    os.path.join(outdir, "2b_Efficiency_group"), logy=False)

    make_group_plot(csv_root, present, bin_keys, "rad",
                    os.path.join(outdir, "3_RadCorr_group"), logy=False)

    make_radcorr_vs_Q2(csv_root, present, bin_keys,
                       os.path.join(outdir, "4_RadCorr_vs_Q2"))

    d5xs  = os.path.join(outdir, "5_Comparison_xs")
    d5acc = os.path.join(outdir, "5_Comparison_acc")
    d5eff = os.path.join(outdir, "5_Comparison_eff")
    d5rad = os.path.join(outdir, "5_Comparison_rad")
    d5red = os.path.join(outdir, "5_Comparison_reduced")
    for (iq, iw) in bin_keys:
        make_comparison_with_ratio(csv_root, present, iq, iw, "xs",  d5xs,  logy)
        make_comparison_with_ratio(csv_root, present, iq, iw, "acc", d5acc, False, show_ratio=False)
        make_comparison_with_ratio(csv_root, present, iq, iw, "eff", d5eff, False, show_ratio=False)
        make_comparison_with_ratio(csv_root, present, iq, iw, "rad", d5rad, False)
        if plot_red:
            make_comparison_with_ratio(csv_root, present, iq, iw, "red", d5red, logy, show_ratio=False)

    make_final_xs_group(csv_root, present, bin_keys,
                        os.path.join(outdir, "6_FinalCrossSection"), logy)

    if plot_red:
        make_reduced_xs_group(csv_root, present, bin_keys,
                              os.path.join(outdir, "6b_ReducedCrossSection_group"),
                              logy=logy)

    # -- Dipole fits on REDUCED cross-section (always run) ------------------
    make_dipole_fit_plots(
        csv_root, present, bin_keys,
        os.path.join(outdir, "7_Dipolefit_ReducedXS"),
        t_min=args.t_min, t_max=args.t_max, logy=logy,
    )

    # -- HS Ds(0) extraction (optional) --------------------------------------
    if args.extract_ds0:
        _run_hs_ds0_workflow(
            csv_root=csv_root, models=present, bin_keys=bin_keys, outdir=outdir,
            tprime_min=args.ds0_tprime_min, tprime_max=args.ds0_tprime_max,
            As0=args.hs_As0, mA=args.hs_mA, mD=args.hs_mD,
            Ds0_bounds=tuple(args.ds0_bounds), logy=logy,
        )

    print(f"\n[DONE] Outputs in {outdir}/")

if __name__ == "__main__":
    main()