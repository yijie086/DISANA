#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hs_only_fit_phi.py
==================

Fast HS-only fitter for strangeness D-term Ds(0) using Hatta–Strikman-inspired
shape template, *always* fitting the REDUCED cross section:

    sigma_red(t') = (dσ/dt') / Γ_v

Inputs: the same per-bin CSVs your main workflow writes, e.g.
  <csv-root>/CSVs/<model>/dsdt_Q<i>.csv
  <csv-root>/CSVs/<model>/dsdt_Q<i>_W<j>.csv

Outputs (in <outdir>/HS_fits):
  - hs_fit_<model>_Q<i>[_W<j>].pdf/.png
  - hs_fit_combined_Q<i>[_W<j>].pdf/.png
  - combined_red_avg_Q<i>[_W<j>].csv   (weighted-average reduced-XSec spectrum)
  - hs_Ds0_summary.csv                 (per-model + combined fit results)
  - hs_Ds0_vs_Q2.pdf/.png
  - hs_Ds0_vs_xB.pdf/.png
  - hs_Ds0_with_xB.csv

This script intentionally does NOT do any acceptance/radiative corrections or
overview plots. It reads whatever is already stored in the CSVs.
"""

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
from matplotlib.ticker import AutoMinorLocator

# -------------------------
# Constants / kinematics
# -------------------------
ALPHA_EM  = 1/137.035999084
M_P_GEV   = 0.9382720813
M_PHI_GEV = 1.019461

# Consistent with your main script default
W_MEAN_GEV = 2.8

# Optional fallback if Q2_center missing and you know iq->Q2 mapping
Q2_MEAN_BY_IQ_GEV2 = {0: 1.0, 1: 2.0, 2: 3.0, 3: 4.0}

# -------------------------
# HS defaults (Eq. 20)
# -------------------------
HS_AS0_DEFAULT = 0.04
HS_mA_DEFAULT  = 1.13  # GeV
HS_mD_DEFAULT  = 0.76  # GeV

# -------------------------
# Plot style
# -------------------------
PANEL_RC = {
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 18,
    "axes.labelsize": 18,
    "axes.titlesize": 18,
    "legend.fontsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "xtick.major.size": 10,
    "xtick.major.width": 2,
    "ytick.major.size": 10,
    "ytick.major.width": 2,
    "xtick.minor.size": 5,
    "xtick.minor.width": 1.8,
    "ytick.minor.size": 5,
    "ytick.minor.width": 1.8,
    "axes.linewidth": 2,
    "lines.linewidth": 2,
}
MODEL_COLORS = ["#1F4FD8", "#E07010", "#009999", "#008800", "#7733CC", "#CC1A33", "#666666"]

def _style_ax(ax):
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))

def _sanitise(s: str) -> str:
    return re.sub(r"[ /\\]", "_", s)

# -------------------------
# Physics helpers
# -------------------------
def t_min_phi(Q2, W, M=M_P_GEV, m_phi=M_PHI_GEV):
    """Return |t_min| for gamma* p -> phi p. Returns NaN if forbidden."""
    Q2  = float(Q2);  W = float(W)
    if W <= M + m_phi or Q2 < 0:
        return np.nan
    W2   = W * W
    mf2  = m_phi * m_phi
    M2   = M * M
    E_gstar = (W2 - M2 - Q2) / (2.0 * W)
    p_gstar = np.sqrt(max(E_gstar**2 + Q2, 0.0))
    E_phi_cm = (W2 + mf2 - M2) / (2.0 * W)
    p_phi_cm2 = E_phi_cm**2 - mf2
    if p_phi_cm2 < 0:
        return np.nan
    p_phi_cm = np.sqrt(p_phi_cm2)
    t_min_val = (E_gstar - E_phi_cm)**2 - (p_gstar - p_phi_cm)**2
    return abs(float(t_min_val))

def xb_from_Q2_W(Q2, W, M=M_P_GEV):
    W2 = float(W) * float(W)
    Q2 = float(Q2)
    return Q2 / (W2 - M*M + Q2)

def hand_virtual_photon_flux(E, Q2, W, M=M_P_GEV, alpha=ALPHA_EM):
    """Hand convention Γ_v. Returns (Gamma_v, eps, y, nu, Ep, K)."""
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

# -------------------------
# HS model functions
# -------------------------
def hs_As(t, As0=HS_AS0_DEFAULT, mA=HS_mA_DEFAULT):
    t = np.asarray(t, dtype=float)
    den = (1.0 - t/(mA*mA))
    return As0 / (den*den)

def hs_Ds(t, Ds0, mD=HS_mD_DEFAULT):
    t = np.asarray(t, dtype=float)
    den = (1.0 - t/(mD*mD))
    return Ds0 / (den*den*den)

def hs_template_sigma_red(t_abs, N, Ds0,
                          As0=HS_AS0_DEFAULT, mA=HS_mA_DEFAULT, mD=HS_mD_DEFAULT,
                          cAA=1.0, cAD=1.0, cDD=1.0):
    """
    HS-inspired shape template for *reduced* cross section (units: nb):
      sigma_red(|t|) = N * [ cAA*As^2 + cAD*|t|*As*Ds + cDD*|t|^2*Ds^2 ]
    with spacelike t = -|t|.
    """
    t_abs = np.asarray(t_abs, dtype=float)
    t = -t_abs
    A = hs_As(t, As0=As0, mA=mA)
    D = hs_Ds(t, Ds0=Ds0, mD=mD)
    poly = (cAA*(A*A) + cAD*(t_abs)*(A*D) + cDD*(t_abs*t_abs)*(D*D))
    return N * np.maximum(poly, 0.0)

# -------------------------
# Data helpers
# -------------------------
def _derive_centers(df):
    df = df.copy()
    for base in ("Q2", "W", "tprime"):
        clo, chi, ccen = f"{base}_lo", f"{base}_hi", f"{base}_center"
        if ccen not in df.columns and clo in df.columns and chi in df.columns:
            df[ccen] = 0.5 * (df[clo].astype(float) + df[chi].astype(float))
    return df

def ensure_q2_w_centers(df, iq=None):
    df = df.copy()
    if "W_center" not in df.columns or df["W_center"].isna().all():
        df["W_center"] = float(W_MEAN_GEV)
    if "Q2_center" not in df.columns or df["Q2_center"].isna().all():
        if iq is not None and int(iq) in Q2_MEAN_BY_IQ_GEV2:
            df["Q2_center"] = float(Q2_MEAN_BY_IQ_GEV2[int(iq)])
        else:
            df["Q2_center"] = float(list(Q2_MEAN_BY_IQ_GEV2.values())[0])
    return df

def compute_tmin_for_df(df):
    df = df.copy()
    if "Q2_center" not in df.columns or "W_center" not in df.columns:
        df["tmin_abs"] = 0.0
        return df
    Q2 = df["Q2_center"].to_numpy(float)
    W  = df["W_center"].to_numpy(float)
    tmin = np.array([t_min_phi(q, w) for q, w in zip(Q2, W)], dtype=float)
    tmin = np.where(np.isfinite(tmin), tmin, 0.0)
    df["tmin_abs"] = tmin
    return df

def _attach_gamma_from_csv_columns(df):
    gamma_cols = ["Gamma_v","GammaV","Gamma","Gamma_v_center","GammaV_center","GammaV_mean","Gamma_v_mean"]
    for c in gamma_cols:
        if c in df.columns and not df[c].isna().all():
            out = df.copy()
            out["Gamma_v"] = out[c].astype(float)
            return out
    return df

def ensure_reduced_xs(df, beam_energy=None, prefer_csv_gamma=True):
    """
    Ensure ReducedCrossSection exists.
    If already present -> keep.
    Else if CrossSection + Gamma_v exist -> compute.
    Else if CrossSection + Q2_center + W_center + beam_energy -> compute Gamma_v then reduced.
    """
    df = df.copy()
    if "ReducedCrossSection" in df.columns and not df["ReducedCrossSection"].isna().all():
        return df

    if prefer_csv_gamma:
        df = _attach_gamma_from_csv_columns(df)

    if "Gamma_v" not in df.columns or df["Gamma_v"].isna().all():
        if beam_energy is None:
            return df
        if "Q2_center" not in df.columns or "W_center" not in df.columns:
            return df
        Q2 = df["Q2_center"].to_numpy(float)
        W  = df["W_center"].to_numpy(float)
        Gamma = np.full(len(df), np.nan)
        for i in range(len(df)):
            if np.isfinite(Q2[i]) and np.isfinite(W[i]):
                Gamma[i] = hand_virtual_photon_flux(beam_energy, Q2[i], W[i])[0]
        df["Gamma_v"] = Gamma

    if "CrossSection" in df.columns and "Gamma_v" in df.columns:
        xs = df["CrossSection"].to_numpy(float)
        xe = df["CrossSection_Err"].to_numpy(float) if "CrossSection_Err" in df.columns else np.zeros_like(xs)
        G  = df["Gamma_v"].to_numpy(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            df["ReducedCrossSection"] = np.where(G > 0, xs/G, np.nan)
            df["ReducedCrossSection_Err"] = np.where(G > 0, xe/G, np.nan)
    return df

def load_csv(csv_root, model, iq, iw=None, beam_energy=None):
    safe = _sanitise(model)
    fname = os.path.join(csv_root, "CSVs", safe,
                         f"dsdt_Q{iq}.csv" if iw is None else f"dsdt_Q{iq}_W{iw}.csv")
    if not os.path.isfile(fname):
        return None
    df = pd.read_csv(fname)
    df = _derive_centers(df)
    df = ensure_q2_w_centers(df, iq=iq)
    df = compute_tmin_for_df(df)
    df = ensure_reduced_xs(df, beam_energy=beam_energy)
    return df

def discover_bin_keys(csv_root, models):
    keys = set()
    for m in models:
        safe = _sanitise(m)
        for f in glob(os.path.join(csv_root, "CSVs", safe, "dsdt_Q*.csv")):
            base = os.path.basename(f).replace("dsdt_", "").replace(".csv", "")
            parts = base.split("_")
            iq = int(parts[0][1:])
            iw = int(parts[1][1:]) if len(parts) > 1 else None
            keys.add((iq, iw))
    return sorted(keys)

# -------------------------
# Weighted average reduced xs (same logic as main script)
# -------------------------
def _safe_sigma(yerr):
    yerr = np.asarray(yerr, dtype=float)
    finite = np.isfinite(yerr)
    if not finite.any():
        return np.ones_like(yerr)
    floor = np.nanmedian(yerr[finite]) * 1e-3
    if (not np.isfinite(floor)) or floor <= 0:
        floor = 1e-12
    return np.where((~finite) | (yerr <= 0), floor, yerr)

def weighted_average_reduced_xs(dfs_list):
    if not dfs_list:
        return None
    num = defaultdict(float)
    den = defaultdict(float)
    count = defaultdict(int)
    t_lo_map, t_hi_map = {}, {}
    ref_df = dfs_list[0]

    for df in dfs_list:
        if df is None or len(df) == 0:
            continue
        if "tprime_center" not in df.columns or "ReducedCrossSection" not in df.columns:
            continue
        t = df["tprime_center"].to_numpy(float)
        y = df["ReducedCrossSection"].to_numpy(float)
        ye = df["ReducedCrossSection_Err"].to_numpy(float) if "ReducedCrossSection_Err" in df.columns else np.full_like(y, np.nan)

        with np.errstate(invalid="ignore"):
            ye = np.where(np.isfinite(ye) & (ye > 0), ye, 0.10 * np.where(np.isfinite(y) & (y > 0), y, np.nan))
        ye = _safe_sigma(ye)

        tlo = df["tprime_lo"].to_numpy(float) if "tprime_lo" in df.columns else np.full_like(t, np.nan)
        thi = df["tprime_hi"].to_numpy(float) if "tprime_hi" in df.columns else np.full_like(t, np.nan)

        for i in range(len(df)):
            tv, yv, sv = float(t[i]), float(y[i]), float(ye[i])
            if not (np.isfinite(tv) and np.isfinite(yv) and yv > 0 and np.isfinite(sv) and sv > 0):
                continue
            tk = round(tv, 4)
            w = 1.0/(sv*sv)
            num[tk] += yv*w
            den[tk] += w
            count[tk] += 1
            if tk not in t_lo_map:
                t_lo_map[tk] = float(tlo[i]) if np.isfinite(tlo[i]) else np.nan
                t_hi_map[tk] = float(thi[i]) if np.isfinite(thi[i]) else np.nan

    if not den:
        return None

    t_keys = sorted(den.keys())
    out = pd.DataFrame({
        "tprime_center": np.array(t_keys, dtype=float),
        "tprime_lo": np.array([t_lo_map.get(k, np.nan) for k in t_keys], dtype=float),
        "tprime_hi": np.array([t_hi_map.get(k, np.nan) for k in t_keys], dtype=float),
        "ReducedCrossSection": np.array([num[k]/den[k] for k in t_keys], dtype=float),
        "ReducedCrossSection_Err": np.array([1.0/np.sqrt(den[k]) for k in t_keys], dtype=float),
        "n_datasets": np.array([count[k] for k in t_keys], dtype=int),
    })
    # carry kinematics (single-valued) from ref
    for col in ("Q2_center","W_center","Q2_lo","Q2_hi","W_lo","W_hi"):
        if col in ref_df.columns and not ref_df[col].isna().all():
            out[col] = float(ref_df[col].iloc[0])
    out = compute_tmin_for_df(out)
    return out

# -------------------------
# HS fit routine (reduced xs)
# -------------------------
def fit_hs_Ds0(df, label,
               tprime_min=0.0, tprime_max=1.0,
               As0=HS_AS0_DEFAULT, mA=HS_mA_DEFAULT, mD=HS_mD_DEFAULT,
               Ds0_bounds=(-5.0, 2.0),
               cAA=1.0, cAD=1.0, cDD=1.0):
    if "ReducedCrossSection" not in df.columns:
        return None

    tprime = df["tprime_center"].to_numpy(float)
    tmin_abs = float(df["tmin_abs"].iloc[0]) if "tmin_abs" in df.columns else 0.0
    if not np.isfinite(tmin_abs):
        tmin_abs = 0.0
    t_abs = tprime + tmin_abs

    y = df["ReducedCrossSection"].to_numpy(float)
    yerr = df["ReducedCrossSection_Err"].to_numpy(float) if "ReducedCrossSection_Err" in df.columns else np.zeros_like(y)

    mask = (np.isfinite(tprime) & np.isfinite(t_abs) & np.isfinite(y) & (y > 0) &
            (tprime >= float(tprime_min)) & (tprime <= float(tprime_max)))
    if int(mask.sum()) < 3:
        print(f"  [WARN] {label}: only {int(mask.sum())} points in t' [{tprime_min},{tprime_max}] -> skip")
        return None

    x_f = t_abs[mask]
    y_f = y[mask]
    with np.errstate(invalid="ignore"):
        yerr = np.where(np.isfinite(yerr) & (yerr > 0), yerr, 0.10*y)
    s_f = _safe_sigma(yerr[mask])

    Ds0_0 = -0.5
    y0_model = hs_template_sigma_red(np.array([x_f[0]]), 1.0, 0.0, As0=As0, mA=mA, mD=mD, cAA=cAA, cAD=cAD, cDD=cDD)[0]
    N0 = float(y_f[0]/y0_model) if np.isfinite(y0_model) and y0_model > 0 else float(np.nanmax(y_f))

    def fwrap(tabs, N, Ds0):
        return hs_template_sigma_red(tabs, N, Ds0, As0=As0, mA=mA, mD=mD, cAA=cAA, cAD=cAD, cDD=cDD)

    try:
        popt, pcov = curve_fit(
            fwrap, x_f, y_f,
            sigma=s_f, absolute_sigma=True,
            p0=[N0, Ds0_0],
            bounds=([0.0, Ds0_bounds[0]], [np.inf, Ds0_bounds[1]]),
            maxfev=40000,
        )
    except Exception as exc:
        print(f"  [WARN] {label}: HS fit failed: {exc}")
        return None

    perr = np.sqrt(np.diag(pcov)) if (pcov is not None and np.all(np.isfinite(pcov))) else [np.nan, np.nan]
    N_fit, Ds0_fit = [float(v) for v in popt]
    dN, dDs0 = [float(v) for v in perr]

    resid = y_f - fwrap(x_f, *popt)
    chi2 = float(np.sum((resid/s_f)**2))
    ndf = int(mask.sum()) - 2

    return dict(
        label=label,
        N=N_fit, dN=dN,
        Ds0=Ds0_fit, dDs0=dDs0,
        chi2=chi2, ndf=ndf,
        tmin_abs=tmin_abs,
        t_abs_fit=x_f, y_fit=y_f, yerr_fit=s_f,
    )

# -------------------------
# Plot HS fit
# -------------------------
def plot_hs_fit(df, res, iq, iw, label, outdir):
    os.makedirs(outdir, exist_ok=True)

    tprime = df["tprime_center"].to_numpy(float)
    y = df["ReducedCrossSection"].to_numpy(float)
    yerr = df["ReducedCrossSection_Err"].to_numpy(float) if "ReducedCrossSection_Err" in df.columns else np.zeros_like(y)
    with np.errstate(invalid="ignore"):
        yerr = np.where(np.isfinite(yerr) & (yerr > 0), yerr, 0.10*y)
    yerr = _safe_sigma(yerr)

    tmin = float(res["tmin_abs"])
    t_abs = tprime + tmin

    x_dense = np.linspace(max(0.0, np.nanmin(tprime)*0.7) + tmin, np.nanmax(tprime)*1.1 + tmin, 500)
    y_dense = hs_template_sigma_red(x_dense, res["N"], res["Ds0"])

    with mpl.rc_context(PANEL_RC):
        fig = plt.figure(figsize=(9, 11))
        gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.0,
                                   left=0.13, right=0.97, top=0.92, bottom=0.09)
        ax = fig.add_subplot(gs[0])
        axr = fig.add_subplot(gs[1], sharex=ax)
        _style_ax(ax); _style_ax(axr)

        # plot in t' on x-axis, but evaluate model at |t| = t'+tmin
        ax.errorbar(tprime, y, yerr=yerr, fmt="o", color="black",
                    ecolor="black", elinewidth=1.5, capsize=3,
                    markersize=7, mfc="white", mew=2.0, zorder=5)
        ax.plot(x_dense - tmin, y_dense, color="#CC1A33", lw=2.2, ls="--", zorder=4)

        ax.set_yscale("log")
        ax.set_ylabel(r"$\sigma_{\rm red} = (d\sigma/dt')/\Gamma_v\ [\mathrm{nb}]$")
        ax.set_xlabel("")  # bottom panel

        txt = (f"{label}\n"
               f"Ds(0) = {res['Ds0']:.3f} ± {res['dDs0']:.3f}\n"
               f"N = {res['N']:.3g} ± {res['dN']:.3g}\n"
               f"|t_min| = {res['tmin_abs']:.4f} GeV$^2$\n"
               f"chi2/ndf = {res['chi2']:.1f}/{res['ndf']}")
        ax.text(0.97, 0.97, txt, transform=ax.transAxes,
                ha="right", va="top", fontsize=13,
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#888888", alpha=0.9))

        y_fit_pts = hs_template_sigma_red(t_abs, res["N"], res["Ds0"])
        pull = (y - y_fit_pts) / yerr
        axr.axhline(0.0, color="gray", lw=1.2, ls="--")
        axr.axhline( 1.0, color="gray", lw=0.8, ls=":")
        axr.axhline(-1.0, color="gray", lw=0.8, ls=":")
        axr.errorbar(tprime, pull, yerr=np.ones_like(pull),
                     fmt="o", color="black", ecolor="black",
                     elinewidth=1.2, capsize=2, markersize=5,
                     mfc="white", mew=1.6)
        axr.set_ylim(-4.5, 4.5)
        axr.set_ylabel("Pull")
        axr.set_xlabel(r"$-t'\ [\mathrm{GeV}^2]$")

        stem = f"hs_fit_{_sanitise(label)}_Q{iq}" + (f"_W{iw}" if iw is not None else "")
        for ext in ("pdf", "png"):
            fig.savefig(os.path.join(outdir, f"{stem}.{ext}"), dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"    [OK] {stem}")

# -------------------------
# Summary plots Ds0 vs Q2 / xB
# -------------------------
def plot_ds0_summaries(tab, outdir):
    os.makedirs(outdir, exist_ok=True)
    if tab.empty:
        return

    # compute xB using W_MEAN_GEV for consistency
    t = tab.copy()
    t["xB"] = [xb_from_Q2_W(row["Q2_center"], W_MEAN_GEV) if np.isfinite(row.get("Q2_center", np.nan)) else np.nan
               for _, row in t.iterrows()]

    # save augmented table
    t.to_csv(os.path.join(outdir, "hs_Ds0_with_xB.csv"), index=False)

    # prefer combined if present
    if "model" in t.columns and (t["model"] == "combined").any():
        plot_tab = t[t["model"] == "combined"].copy()
        title_suffix = " (combined)"
    else:
        plot_tab = t.copy()
        title_suffix = ""

    plot_tab = plot_tab[np.isfinite(plot_tab["Ds0"]) & np.isfinite(plot_tab["dDs0"])].copy()
    if plot_tab.empty:
        return

    with mpl.rc_context(PANEL_RC):
        # Ds0 vs Q2
        fig, ax = plt.subplots(figsize=(9, 6.5))
        _style_ax(ax)
        plot_tab = plot_tab.sort_values("Q2_center")
        ax.errorbar(plot_tab["Q2_center"], plot_tab["Ds0"], yerr=plot_tab["dDs0"],
                    fmt="o", color="black", ecolor="black",
                    elinewidth=1.8, capsize=4, markersize=8, mfc="white", mew=2.0)
        ax.set_xlabel(r"$Q^2\ [\mathrm{GeV}^2]$")
        ax.set_ylabel(r"$D_s(0)$")
        ax.set_title(r"$D_s(0)$ vs $Q^2$" + title_suffix)
        for ext in ("pdf","png"):
            fig.savefig(os.path.join(outdir, f"hs_Ds0_vs_Q2.{ext}"), dpi=200, bbox_inches="tight")
        plt.close(fig)

        # Ds0 vs xB
        fig, ax = plt.subplots(figsize=(9, 6.5))
        _style_ax(ax)
        plot_tab = plot_tab.sort_values("xB")
        ax.errorbar(plot_tab["xB"], plot_tab["Ds0"], yerr=plot_tab["dDs0"],
                    fmt="o", color="black", ecolor="black",
                    elinewidth=1.8, capsize=4, markersize=8, mfc="white", mew=2.0)
        ax.set_xlabel(r"$x_B$ (computed with $W_\mathrm{mean}$)")
        ax.set_ylabel(r"$D_s(0)$")
        ax.set_title(r"$D_s(0)$ vs $x_B$" + title_suffix)
        # log-x if wide range
        xb = plot_tab["xB"].to_numpy(float)
        if np.nanmax(xb) / max(np.nanmin(xb), 1e-6) > 5:
            ax.set_xscale("log")
        for ext in ("pdf","png"):
            fig.savefig(os.path.join(outdir, f"hs_Ds0_vs_xB.{ext}"), dpi=200, bbox_inches="tight")
        plt.close(fig)

# -------------------------
# Main
# -------------------------
def main():
    ap = argparse.ArgumentParser(description="Fast HS-only Ds(0) fitter on reduced cross sections.")
    ap.add_argument("--csv-root", default=".", help="Root dir containing CSVs/ sub-tree.")
    ap.add_argument("--models", nargs="+", required=True, help="Model names (subdirs under CSVs/).")
    ap.add_argument("--outdir", default="phi_results", help="Output top directory (HS_fits will be created inside).")

    ap.add_argument("--beam-energy", type=float, default=None,
                    help="Optional: beam energy [GeV] used ONLY if ReducedCrossSection is missing and we need to compute Γv.")
    ap.add_argument("--tprime-min", type=float, default=0.0, help="Min t' for HS fit [GeV^2].")
    ap.add_argument("--tprime-max", type=float, default=1.0, help="Max t' for HS fit [GeV^2].")

    ap.add_argument("--hs-As0", type=float, default=HS_AS0_DEFAULT, help="HS As(0) (default 0.04).")
    ap.add_argument("--hs-mA", type=float, default=HS_mA_DEFAULT, help="HS mA [GeV] (default 1.13).")
    ap.add_argument("--hs-mD", type=float, default=HS_mD_DEFAULT, help="HS mD [GeV] (default 0.76).")
    ap.add_argument("--ds0-bounds", nargs=2, type=float, default=[-5.0, 2.0],
                    metavar=("DS0_MIN","DS0_MAX"), help="Bounds for Ds(0) fit.")
    ap.add_argument("--no-combined", action="store_true",
                    help="Skip weighted-average combined HS fit (per-model only).")
    args = ap.parse_args()

    models = args.models
    csv_root = args.csv_root
    outdir = args.outdir
    hs_dir = os.path.join(outdir, "HS_fits")
    os.makedirs(hs_dir, exist_ok=True)

    # Filter models that exist
    present = []
    for m in models:
        d = os.path.join(csv_root, "CSVs", _sanitise(m))
        if os.path.isdir(d) and glob(os.path.join(d, "dsdt_Q*.csv")):
            present.append(m)
    if not present:
        print("[ERROR] No model CSV directories found under CSVs/.")
        return

    bin_keys = discover_bin_keys(csv_root, present)
    if not bin_keys:
        print("[ERROR] No dsdt_Q*.csv files found.")
        return

    print(f"Found {len(bin_keys)} (Q,W) bins. Running HS-only fits...")

    rows = []

    for (iq, iw) in bin_keys:
        dfs = []
        mnames = []
        for m in present:
            df = load_csv(csv_root, m, iq, iw, beam_energy=args.beam_energy)
            if df is None:
                continue
            if "ReducedCrossSection" not in df.columns or df["ReducedCrossSection"].isna().all():
                print(f"  [WARN] Missing ReducedCrossSection for {m} Q{iq} W{iw}; skipping this dataset.")
                continue

            res = fit_hs_Ds0(
                df, m,
                tprime_min=args.tprime_min,
                tprime_max=args.tprime_max,
                As0=args.hs_As0, mA=args.hs_mA, mD=args.hs_mD,
                Ds0_bounds=tuple(args.ds0_bounds),
            )
            if res is None:
                continue

            # plot
            plot_hs_fit(df, res, iq, iw, m, hs_dir)

            q2c = float(df["Q2_center"].iloc[0]) if "Q2_center" in df.columns else np.nan
            wc  = float(df["W_center"].iloc[0]) if "W_center" in df.columns else np.nan
            rows.append(dict(model=m, iq=iq, iw=iw, Q2_center=q2c, W_center=wc,
                             Ds0=res["Ds0"], dDs0=res["dDs0"], N=res["N"], dN=res["dN"],
                             chi2=res["chi2"], ndf=res["ndf"], tmin_abs=res["tmin_abs"]))

            dfs.append(df)
            mnames.append(m)

        # Combined
        if (not args.no_combined) and dfs:
            avg_df = weighted_average_reduced_xs(dfs)
            if avg_df is not None:
                # save spectrum used
                avg_csv = os.path.join(hs_dir, f"combined_red_avg_Q{iq}.csv" if iw is None else f"combined_red_avg_Q{iq}_W{iw}.csv")
                avg_df.to_csv(avg_csv, index=False)

                res_c = fit_hs_Ds0(
                    avg_df, "combined",
                    tprime_min=args.tprime_min,
                    tprime_max=args.tprime_max,
                    As0=args.hs_As0, mA=args.hs_mA, mD=args.hs_mD,
                    Ds0_bounds=tuple(args.ds0_bounds),
                )
                if res_c is not None:
                    plot_hs_fit(avg_df, res_c, iq, iw, "combined", hs_dir)

                    q2c = float(avg_df["Q2_center"].iloc[0]) if "Q2_center" in avg_df.columns else np.nan
                    wc  = float(avg_df["W_center"].iloc[0]) if "W_center" in avg_df.columns else np.nan
                    rows.append(dict(model="combined", iq=iq, iw=iw, Q2_center=q2c, W_center=wc,
                                     Ds0=res_c["Ds0"], dDs0=res_c["dDs0"], N=res_c["N"], dN=res_c["dN"],
                                     chi2=res_c["chi2"], ndf=res_c["ndf"], tmin_abs=res_c["tmin_abs"]))
            else:
                print(f"  [WARN] Could not build weighted average for Q{iq} W{iw}.")

    if not rows:
        print("[WARN] No successful HS fits.")
        return

    tab = pd.DataFrame(rows)
    tab_path = os.path.join(hs_dir, "hs_Ds0_summary.csv")
    tab.to_csv(tab_path, index=False)
    print(f"\n[OK] wrote {tab_path}")
    print(tab.to_string(index=False))

    # Summary plots
    plot_ds0_summaries(tab, hs_dir)
    print(f"\n[DONE] HS outputs in {hs_dir}/")

if __name__ == "__main__":
    main()
