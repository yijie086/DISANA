#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_phi_dsdt_from_csv.py
=========================

Complete φ dσ/dt analysis — plots and dipole fits.

Fixes in this version:
  - Per-model beam energy auto-detection (Sp18→10.6, Fall18→10.594, Sp19→10.2)
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
      raw_xs, Gamma_v, reduced_xs, Acceptance, xs_acc_corrected, RadCorr, xs_rad_corrected
  - Corrections panel ratio is a TRUE correction-factor plot:
      ratio = Final / Variant   (log-scale)
  - Acceptance group y-range forced to 0.0–0.05
"""
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

# ─────────────────────────────────────────────────────────────────────────────
#  External radiative correction tables (hard-coded paths)
# ─────────────────────────────────────────────────────────────────────────────
K_RADCORR_FILE_10P2 = "/w/hallb-scshelf2102/clas12/singh/Softwares/Generators/PhiEventGen/DIFFRAD/diffradstn/outputs/10p2GeV/diffrad_rc_results.csv"
K_RADCORR_FILE_10P6 = "/w/hallb-scshelf2102/clas12/singh/Softwares/Generators/PhiEventGen/DIFFRAD/diffradstn/outputs/10p6GeV/diffrad_rc_results.csv"

# ─────────────────────────────────────────────────────────────────────────────
#  Per-model beam energy lookup  (case-insensitive substring match on model name)
# ─────────────────────────────────────────────────────────────────────────────
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

# Global run configuration set in main()
G_BEAM_ENERGY_GEV      = None   # float  – fallback only; per-model takes priority
G_USE_EXTERNAL_RADCORR = False  # bool
G_RADCORR_FILE_OVERRIDE= None   # optional str
G_RADCORR_MATCH_MODE   = "nearest"

ALPHA_EM  = 1/137.035999084
M_P_GEV   = 0.9382720813

# ─────────────────────────────────────────────────────────────────────────────
#  Hardcoded kinematics fallbacks (TOP, as requested)
# ─────────────────────────────────────────────────────────────────────────────
W_MEAN_GEV = 2.8  # hardcoded mean W (requested)

# If Q2_center is missing, use an IQ->Q2 mean lookup (EDIT if needed)
Q2_MEAN_BY_IQ_GEV2 = {
    0: 1.0,
    1: 2.0,
    2: 3.0,
    3: 4.0,
}

# ─────────────────────────────────────────────────────────────────────────────
#  Physics helpers
# ─────────────────────────────────────────────────────────────────────────────
def xb_from_Q2_W(Q2, W, M=M_P_GEV):
    W2 = W * W
    return Q2 / (W2 - M*M + Q2)

def hand_virtual_photon_flux(E, Q2, W, M=M_P_GEV, alpha=ALPHA_EM):
    """Hand convention virtual photon flux Γ_v."""
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

# ─────────────────────────────────────────────────────────────────────────────
#  External RadCorr table loader + attacher (DIFFRAD)
# ─────────────────────────────────────────────────────────────────────────────
def load_external_radcorr_table(path):
    """Load DIFFRAD CSV → standardised DataFrame with bin edges + crad."""
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
                      f" → nearest centroid Q2={tab_q2c[j]:.3f}, t={tab_tc[j]:.3f}"
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
        print("  [WARN] RadCorr table has neither bin edges nor centroid columns — skipping.")
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

# ─────────────────────────────────────────────────────────────────────────────
#  Derive center columns helper  (call before any physics attachment)
# ─────────────────────────────────────────────────────────────────────────────
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

# ─────────────────────────────────────────────────────────────────────────────
#  Table dump
# ─────────────────────────────────────────────────────────────────────────────
def dump_xs_table(df, variants, label, outdir):
    raw_xs, raw_err = variants["raw"]
    acc_xs, acc_err = variants["acc"]
    fin_xs, fin_err = variants["final"]

    n = len(df)

    def col_or_nan(name):
        return df[name].to_numpy(float) if name in df.columns else np.full(n, np.nan)

    tcent  = col_or_nan("tprime_center")
    A      = col_or_nan("Acceptance")
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
        "xs_acc_corrected":      acc_xs,
        "RadCorr":               Crad,
        "xs_rad_corrected":      fin_xs,
        "raw_xs_err":            raw_err,
        "xs_acc_corrected_err":  acc_err,
        "xs_rad_corrected_err":  fin_err,
    })

    print(f"\n=== XS TABLE: {label} ===")
    print(tab.to_string(index=False))

    os.makedirs(outdir, exist_ok=True)
    tab.to_csv(os.path.join(outdir, f"table_{label}.csv"), index=False)

# ─────────────────────────────────────────────────────────────────────────────
#  Style constants
# ─────────────────────────────────────────────────────────────────────────────
MODEL_COLORS = ["#1F4FD8", "#E07010", "#009999", "#008800", "#7733CC", "#CC1A33", "#666666"]

CORR_COLORS = {"raw": "#000000", "acc": "#1F4FD8", "final": "#009900"}
CORR_LABELS = {
    "raw":   "Raw  (no corr.)",
    "acc":   "Acceptance corr. only",
    "final": "Acc. + Rad. (final)",
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
YLIM_XS      = (1e-5, 20.0)
XLIM_T       = (0.0, 6.0)
NCOLS_GRP    = 2

_CANVAS_CFG = {
    "xs": dict(
        col="CrossSection", err_col="CrossSection_Err",
        ylabel=r"$\mathrm{d}\sigma/\mathrm{d}t'\ [\mathrm{nb/GeV}^2]$",
        logy=True, tag="dsdt", title="Cross-section (final)"),
    "red": dict(
        col="ReducedCrossSection", err_col="ReducedCrossSection_Err",
        ylabel=r"$\sigma_\mathrm{red} = (\mathrm{d}\sigma/\mathrm{d}t')/\Gamma_v\ [\mathrm{nb}]$",
        logy=True, tag="dsdt_reduced", title="Reduced cross-section"),
    "acc": dict(
        col="Acceptance", err_col=None,
        ylabel=r"Acceptance $A(\varepsilon)$",
        logy=False, tag="acceptance", title="Acceptance"),
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
        lbl += rf"$\quad W \in [{wlo:.2f},\,{whi:.2f}]\ \mathrm{{GeV}}$"
    return lbl

def _save(fig, outdir, stem):
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"{stem}.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"    [OK] {stem}")

# ─────────────────────────────────────────────────────────────────────────────
#  CSV loader  (per-model beam energy + center-column derivation + W/Q2 forcing)
# ─────────────────────────────────────────────────────────────────────────────
def load_csv(csv_root, model, iq, iw=None):
    safe  = _sanitise(model)
    fname = os.path.join(
        csv_root, "CSVs", safe,
        f"dsdt_Q{iq}.csv" if iw is None else f"dsdt_Q{iq}_W{iw}.csv"
    )
    if not os.path.isfile(fname):
        return None

    df = pd.read_csv(fname)

    # ── derive center columns from lo/hi if possible ──────────────────────
    df = _derive_centers(df)

    # ── NEW: force W_center and Q2_center even if CSV omitted them ────────
    df = ensure_q2_w_centers(df, iq=iq, w_mean=W_MEAN_GEV)

    # ── per-model beam energy ────────────────────────────────────────────
    E = beam_energy_for_model(model)

    if G_USE_EXTERNAL_RADCORR:
        df = attach_external_radcorr(
            df, E,
            path_override=G_RADCORR_FILE_OVERRIDE,
            mode=G_RADCORR_MATCH_MODE,
        )

    if E is not None:
        df = attach_virtual_photon_flux_and_reduced_xs(df, E)

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

# ─────────────────────────────────────────────────────────────────────────────
#  Build correction variants
# ─────────────────────────────────────────────────────────────────────────────
def _build_variants(df, luminosity=None, branching=None):
    xs   = _safe_col(df, "CrossSection")
    xerr = _safe_col(df, "CrossSection_Err")
    A    = _safe_col(df, "Acceptance")
    Crad = _safe_col(df, "RadCorr")
    Nsig = _safe_col(df, "RawCounts")
    Nerr = _safe_col(df, "RawCounts_Err")
    dT   = df["tprime_hi"].to_numpy(float) - df["tprime_lo"].to_numpy(float)

    A_s    = np.where(np.isfinite(A)    & (A    > 0), A,    1.0)
    Crad_s = np.where(np.isfinite(Crad) & (Crad > 0), Crad, 1.0)
    dT_s   = np.where(dT > 0, dT, 1.0)

    have_nsig = np.any(np.isfinite(Nsig) & (Nsig > 0))

    if have_nsig and luminosity is not None and branching is not None:
        L  = luminosity
        BR = branching
        denom_raw   = L * BR * dT_s
        denom_acc   = L * BR * dT_s * A_s
        denom_final = L * BR * dT_s * A_s * Crad_s

        raw_v   = np.where(denom_raw   > 0, Nsig / denom_raw,   np.nan)
        acc_v   = np.where(denom_acc   > 0, Nsig / denom_acc,   np.nan)
        final_v = np.where(denom_final > 0, Nsig / denom_final, np.nan)
        raw_e   = np.where(denom_raw   > 0, Nerr / denom_raw,   np.nan)
        acc_e   = np.where(denom_acc   > 0, Nerr / denom_acc,   np.nan)
        final_e = np.where(denom_final > 0, Nerr / denom_final, np.nan)
    else:
        # CrossSection already fully corrected → reconstruct variants
        final_v = xs.copy()
        final_e = xerr.copy()
        acc_v   = xs   * Crad_s          # undo RadCorr
        acc_e   = xerr * Crad_s
        raw_v   = xs   * A_s * Crad_s    # undo both
        raw_e   = xerr * A_s * Crad_s

    return {"raw": (raw_v, raw_e), "acc": (acc_v, acc_e), "final": (final_v, final_e)}

# ─────────────────────────────────────────────────────────────────────────────
#  Corrections canvas (before/after)
# ─────────────────────────────────────────────────────────────────────────────
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

            for vkey in ("raw", "acc"):
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

# ─────────────────────────────────────────────────────────────────────────────
#  Group / comparison plots
# ─────────────────────────────────────────────────────────────────────────────
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
                elif canvas_key == "acc":
                    ax.set_ylim(0.0, 0.05)   # <-- requested
                elif canvas_key == "rad":
                    ax.set_ylim(0.6, 1.0)
                    ax.axhline(1.0, color="gray", lw=1.2, ls="--", zorder=1)

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

            title = cfg["title"] + r" — all $Q^2$ bins"
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
                w_title = rf" — $W \in [{wlo:.2f},\,{whi:.2f}]\ \mathrm{{GeV}}$"
            fig.suptitle(r"$C_\mathrm{rad}$ vs $Q^2$" + w_title,
                         fontsize=LEG_TITLE_FS + 2)
            stem = ("radcorr_vs_Q2" if iw_key is None
                    else f"radcorr_vs_Q2_W{iw_key}")
            _save(fig, outdir, stem)


def make_comparison_with_ratio(csv_root, models, iq, iw, canvas_key, outdir, logy=True):
    cfg    = _CANVAS_CFG[canvas_key]
    col    = cfg["col"]; ec = cfg["err_col"]; ylabel = cfg["ylabel"]

    dfs = {m: df for m in models if (df := load_csv(csv_root, m, iq, iw)) is not None}
    if not dfs:
        return

    ref_label = list(dfs.keys())[0]
    ref_df    = dfs[ref_label]
    t_ref     = ref_df["tprime_center"].to_numpy(float)
    y_ref     = _safe_col(ref_df, col)
    y_ref_s   = np.where(np.isfinite(y_ref) & (y_ref > 0), y_ref, np.nan)

    with mpl.rc_context(PANEL_RC):
        fig = plt.figure(figsize=(10, 12))
        gs  = mpl.gridspec.GridSpec(
            2, 1, figure=fig, height_ratios=[3, 1], hspace=0.0,
            left=0.13, right=0.97, top=0.91, bottom=0.09,
        )
        ax_m = fig.add_subplot(gs[0])
        ax_r = fig.add_subplot(gs[1], sharex=ax_m)
        _style_ax(ax_m); _style_ax(ax_r)

        handles, labels_l = [], []
        for im, (m, df) in enumerate(dfs.items()):
            color = MODEL_COLORS[im % len(MODEL_COLORS)]
            t     = df["tprime_center"].to_numpy(float)
            y     = _safe_col(df, col)
            yerr  = _safe_col(df, ec) if ec else np.zeros_like(y)
            xlo   = t - df["tprime_lo"].to_numpy(float)
            xhi   = df["tprime_hi"].to_numpy(float) - t
            mask  = np.isfinite(y) & (y > 0)
            if not mask.any():
                continue
            h = ax_m.errorbar(
                t[mask], y[mask], yerr=yerr[mask],
                xerr=(xlo[mask], xhi[mask]),
                fmt="o", color=color, ecolor=color, elinewidth=1.5, capsize=3,
                markersize=7, mfc="white", mew=2.0, lw=1.8, zorder=5 + im,
                label=m.replace("_", " "),
            )
            handles.append(h); labels_l.append(m.replace("_", " "))
            if im == 0:
                ax_r.axhline(1.0, color=color, lw=1.5, ls="--", zorder=1)
            else:
                y_int = np.interp(t_ref, t[mask], y[mask],
                                  left=np.nan, right=np.nan)
                ratio = y_int / y_ref_s
                rmask = np.isfinite(ratio)
                if rmask.any():
                    ax_r.errorbar(
                        t_ref[rmask], ratio[rmask],
                        fmt="o", color=color, ecolor=color,
                        elinewidth=1.5, capsize=3,
                        markersize=6, mfc="white", mew=1.8, lw=1.5, zorder=5 + im,
                    )

        ax_m.set_xlim(XLIM_T); ax_r.set_xlim(XLIM_T)
        ax_m.tick_params(labelbottom=False)
        if cfg["logy"] and logy:
            ax_m.set_yscale("log"); ax_m.set_ylim(YLIM_XS)
        elif canvas_key == "acc":
            ax_m.set_ylim(0.0, 0.05)  # <-- requested
        elif canvas_key == "rad":
            ax_m.axhline(1.0, color="gray", lw=1.2, ls="--", zorder=0)
            ax_m.set_ylim(0.6, 1.0)
        ax_m.set_ylabel(ylabel)
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


def make_reduced_xs_group(csv_root, models, bin_keys, outdir, logy=True):
    make_group_plot(csv_root, models, bin_keys, "red", outdir, logy=logy)


# ─────────────────────────────────────────────────────────────────────────────
#  Dipole fit on REDUCED cross-section
# ─────────────────────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────────────────────
#  Dipole fit on REDUCED cross-section  — model, fitter, per-bin plots,
#  residual panels, parameter-vs-Q²/W summary plots
# ─────────────────────────────────────────────────────────────────────────────

def _dipole_model(t, A, b, t0):
    """Dipole parameterisation: A * exp(-b*t) / (1 + t/t0)^4.

    Parameters
    ----------
    t  : |t'| values (GeV²)
    A  : overall normalisation (nb)
    b  : exponential slope (GeV⁻²)
    t0 : dipole mass-squared scale (GeV²), must be > 0
    """
    return A * np.exp(-b * t) / (1.0 + t / t0) ** 4


def _fit_dipole_reduced(df, model_name, t_min=0.1, t_max=5.5):
    """Fit the reduced cross-section in *df* with the dipole model.

    Returns a result dict, or None if the fit cannot be performed.
    All arrays in the dict cover only the fitted range (mask applied).
    """
    t    = df["tprime_center"].to_numpy(float)
    y    = _safe_col(df, "ReducedCrossSection")
    yerr = _safe_col(df, "ReducedCrossSection_Err")

    # Replace zero / nan errors with a 10 % relative floor so curve_fit
    # does not blow up on missing uncertainties.
    with np.errstate(invalid="ignore"):
        yerr = np.where(
            np.isfinite(yerr) & (yerr > 0),
            yerr,
            0.10 * np.where(np.isfinite(y) & (y > 0), y, 1.0),
        )

    # full-range mask (for plotting outside the fit window)
    mask_all = np.isfinite(t) & np.isfinite(y) & (y > 0)
    # fit-range mask
    mask = mask_all & np.isfinite(yerr) & (yerr > 0) & (t >= t_min) & (t <= t_max)

    if mask.sum() < 3:
        print(f"  [WARN] {model_name}: only {mask.sum()} valid reduced-XS points "
              f"in [{t_min},{t_max}] GeV² — skipping dipole fit.")
        return None

    t_f = t[mask];  y_f = y[mask];  s_f = yerr[mask]

    # ── initial guesses ───────────────────────────────────────────────────
    A0  = float(np.nanmedian(y_f)) if np.isfinite(y_f).any() else 1.0
    b0  = 2.0
    t00 = 0.71   # standard dipole ≃ m_ρ²

    try:
        popt, pcov = curve_fit(
            _dipole_model, t_f, y_f, sigma=s_f, absolute_sigma=True,
            p0=[A0, b0, t00],
            bounds=([0.0, 0.0, 1e-3], [1e9, 100.0, 100.0]),
            maxfev=40000,
        )
    except (RuntimeError, ValueError) as exc:
        print(f"  [WARN] {model_name}: dipole fit failed — {exc}")
        return None

    perr = (np.sqrt(np.diag(pcov))
            if (pcov is not None and np.all(np.isfinite(pcov)))
            else np.full(3, np.nan))
    A, b, t0 = popt
    dA, db, dt0 = perr

    # ── goodness-of-fit ───────────────────────────────────────────────────
    resid = y_f - _dipole_model(t_f, *popt)
    chi2  = float(np.sum((resid / s_f) ** 2))
    ndf   = int(mask.sum()) - 3

    # pull = normalised residual per point
    pull  = resid / s_f

    return dict(
        model=model_name,
        # parameters
        A=A, dA=dA, b=b, db=db, t0=t0, dt0=dt0,
        # goodness-of-fit
        chi2=chi2, ndf=ndf,
        # fit-range arrays
        t_fit=t_f, y_fit=y_f, yerr_fit=s_f,
        resid=resid, pull=pull,
        # full data arrays (for plotting outside fit window)
        t_all=t[mask_all], y_all=y[mask_all], yerr_all=yerr[mask_all],
        popt=popt,
    )


# ── per-bin fit plot  (data + fit curve + residuals panel) ────────────────────
def _make_one_fit_canvas(ax_m, ax_r, df, model, color, t_min, t_max, logy,
                         show_ylabel=True):
    """Draw one bin's reduced-XS fit onto a pre-created (ax_m, ax_r) pair.

    Returns the fit result dict (or None) and a legend handle for the data.
    """
    _style_ax(ax_m); _style_ax(ax_r)

    t    = df["tprime_center"].to_numpy(float)
    y    = _safe_col(df, "ReducedCrossSection")
    yerr = _safe_col(df, "ReducedCrossSection_Err")
    xlo  = t - df["tprime_lo"].to_numpy(float)
    xhi  = df["tprime_hi"].to_numpy(float) - t

    mask_plot = np.isfinite(y) & (y > 0)
    if not mask_plot.any():
        return None, None

    # ── data points (all in range) ────────────────────────────────────────
    h_data = ax_m.errorbar(
        t[mask_plot], y[mask_plot],
        yerr=np.where(np.isfinite(yerr[mask_plot]) & (yerr[mask_plot] > 0),
                      yerr[mask_plot], 0.0),
        xerr=(xlo[mask_plot], xhi[mask_plot]),
        fmt="o", color=color, ecolor=color,
        elinewidth=1.5, capsize=3,
        markersize=7, mfc="white", mew=2.0, lw=1.8, zorder=5,
        label=model.replace("_", " "),
    )

    # ── dipole fit ────────────────────────────────────────────────────────
    res = _fit_dipole_reduced(df, model, t_min=t_min, t_max=t_max)
    if res is not None:
        t_dense  = np.linspace(max(0.0, t_min * 0.5), t_max * 1.05, 600)
        y_curve  = _dipole_model(t_dense, *res["popt"])
        chi2_str = (rf"$\chi^2/\nu = {res['chi2']:.1f}/{res['ndf']}$"
                    if res["ndf"] > 0 else "")
        fit_txt  = (
            rf"$A = {res['A']:.3g} \pm {res['dA']:.1g}$ nb" "\n"
            rf"$b = {res['b']:.3f} \pm {res['db']:.3f}$ GeV$^{{-2}}$" "\n"
            rf"$t_0 = {res['t0']:.3f} \pm {res['dt0']:.3f}$ GeV$^2$" "\n"
            + chi2_str
        )

        ax_m.plot(t_dense, y_curve, color=color, lw=2.2, ls="--", zorder=4,
                  label="Dipole fit")

        # shade the fit region faintly
        ax_m.axvspan(t_min, t_max, alpha=0.06, color=color, zorder=0)

        # fit-parameter box inside the plot
        ax_m.text(
            0.97, 0.97, fit_txt,
            transform=ax_m.transAxes,
            fontsize=14, va="top", ha="right", color=color,
            bbox=dict(boxstyle="round,pad=0.35", fc="white",
                      ec=color, alpha=0.85, lw=1.2),
        )

        # ── residuals (pull = (data−fit)/σ) ──────────────────────────────
        ax_r.axhline(0.0, color="gray", lw=1.2, ls="--", zorder=1)
        ax_r.axhline( 1.0, color="gray", lw=0.8, ls=":", zorder=1)
        ax_r.axhline(-1.0, color="gray", lw=0.8, ls=":", zorder=1)
        ax_r.errorbar(
            res["t_fit"], res["pull"],
            yerr=np.ones_like(res["pull"]),
            fmt="o", color=color, ecolor=color,
            elinewidth=1.2, capsize=2,
            markersize=5, mfc="white", mew=1.8, zorder=5,
        )

    ax_m.set_xlim(XLIM_T)
    if logy:
        ax_m.set_yscale("log")
    if show_ylabel:
        ax_m.set_ylabel(
            r"$\sigma_\mathrm{red}\ [\mathrm{nb}]$", fontsize=20,
        )
    ax_m.tick_params(labelbottom=False)
    ax_m.tick_params(axis="both", labelsize=18)

    ax_r.set_xlim(XLIM_T)
    ax_r.set_ylim(-4.5, 4.5)
    ax_r.set_xlabel(r"$-t'\ [\mathrm{GeV}^2]$", fontsize=20)
    if show_ylabel:
        ax_r.set_ylabel("Pull", fontsize=17)
    ax_r.tick_params(axis="both", labelsize=16)

    return res, h_data



# ─────────────────────────────────────────────────────────────────────────────
#  Weighted-average reduced cross-section across datasets
# ─────────────────────────────────────────────────────────────────────────────
def _weighted_average_reduced_xs(dfs_list):
    """
    Combine multiple DataFrames into a single inverse-variance weighted
    average reduced cross-section.

    For every unique t-bin (tprime_center rounded to 4 dp):
        y_avg  = sum(y_i / sig_i^2) / sum(1 / sig_i^2)
        sig_avg = 1 / sqrt(sum(1 / sig_i^2))

    Returns a single DataFrame with the same column layout as the input
    DataFrames (ReducedCrossSection, ReducedCrossSection_Err, tprime_*,
    Q2/W kinematics from the first dataset), plus n_datasets per bin.
    """
    if not dfs_list:
        return None

    from collections import defaultdict
    num    = defaultdict(float)
    denom  = defaultdict(float)
    count  = defaultdict(int)
    t_lo_map, t_hi_map = {}, {}

    # reference df for kinematics
    ref_df = dfs_list[0]

    for df in dfs_list:
        t_col    = df["tprime_center"].to_numpy(float)
        y_col    = _safe_col(df, "ReducedCrossSection")
        yerr_col = _safe_col(df, "ReducedCrossSection_Err")

        # 10 % floor if stat error missing
        with np.errstate(invalid="ignore"):
            yerr_col = np.where(
                np.isfinite(yerr_col) & (yerr_col > 0),
                yerr_col,
                0.10 * np.where(np.isfinite(y_col) & (y_col > 0), y_col, np.nan),
            )

        for i in range(len(df)):
            t_v = float(t_col[i]);  y_v = float(y_col[i]);  s_v = float(yerr_col[i])
            if not (np.isfinite(t_v) and np.isfinite(y_v)
                    and np.isfinite(s_v) and y_v > 0 and s_v > 0):
                continue
            tk = round(t_v, 4)
            w  = 1.0 / s_v ** 2
            num[tk]   += y_v * w
            denom[tk] += w
            count[tk] += 1
            if tk not in t_lo_map and "tprime_lo" in df.columns:
                t_lo_map[tk] = float(df["tprime_lo"].iloc[i])
                t_hi_map[tk] = float(df["tprime_hi"].iloc[i])

    if not denom:
        return None

    t_keys = sorted(denom.keys())
    t_cen  = np.array(t_keys, dtype=float)
    y_avg  = np.array([num[tk] / denom[tk]      for tk in t_keys])
    y_err  = np.array([1.0 / np.sqrt(denom[tk]) for tk in t_keys])
    n_ds   = np.array([count[tk]                for tk in t_keys], dtype=int)
    t_lo   = np.array([t_lo_map.get(tk, np.nan) for tk in t_keys])
    t_hi   = np.array([t_hi_map.get(tk, np.nan) for tk in t_keys])

    out = pd.DataFrame({
        "tprime_center":           t_cen,
        "tprime_lo":               t_lo,
        "tprime_hi":               t_hi,
        "ReducedCrossSection":     y_avg,
        "ReducedCrossSection_Err": y_err,
        "n_datasets":              n_ds,
    })

    # carry over Q2 / W kinematics from reference dataset
    for col in ("Q2_lo", "Q2_hi", "Q2_center", "W_lo", "W_hi", "W_center"):
        if col in ref_df.columns:
            out[col] = float(ref_df[col].iloc[0])

    return out


# ── combined fit canvas: individual datasets (faint) + weighted avg + fit ─────
def _make_combined_fit_canvas(ax_m, ax_r, dfs_list, model_names,
                               avg_df, t_min, t_max, logy, show_ylabel=True):
    """
    Draw individual datasets as faint coloured markers, the weighted average
    as filled black circles, and a single dipole fit to the average.
    Returns the fit result dict (or None).
    """
    _style_ax(ax_m);  _style_ax(ax_r)
    COLOR_AVG = "#111111"

    # ── individual datasets in the background ─────────────────────────────
    for im, (m, df) in enumerate(zip(model_names, dfs_list)):
        t    = df["tprime_center"].to_numpy(float)
        y    = _safe_col(df, "ReducedCrossSection")
        yerr = _safe_col(df, "ReducedCrossSection_Err")
        xlo  = t - _safe_col(df, "tprime_lo")
        xhi  = _safe_col(df, "tprime_hi") - t
        mask = np.isfinite(y) & (y > 0)
        if not mask.any():
            continue
        color = MODEL_COLORS[im % len(MODEL_COLORS)]
        ax_m.errorbar(
            t[mask], y[mask],
            yerr=np.where(np.isfinite(yerr[mask]) & (yerr[mask] > 0), yerr[mask], 0.0),
            xerr=(np.where(np.isfinite(xlo[mask]), xlo[mask], 0.0),
                  np.where(np.isfinite(xhi[mask]), xhi[mask], 0.0)),
            fmt="o", color=color, ecolor=color,
            alpha=0.35, elinewidth=1.0, capsize=2,
            markersize=5, mfc=color, mew=1.0, lw=0.8, zorder=3,
            label=m.replace("_", " "),
        )

    # ── weighted average ──────────────────────────────────────────────────
    t_a   = avg_df["tprime_center"].to_numpy(float)
    y_a   = avg_df["ReducedCrossSection"].to_numpy(float)
    ye_a  = avg_df["ReducedCrossSection_Err"].to_numpy(float)
    xlo_a = t_a - avg_df["tprime_lo"].to_numpy(float)
    xhi_a = avg_df["tprime_hi"].to_numpy(float) - t_a
    mask_a = np.isfinite(y_a) & (y_a > 0)

    ax_m.errorbar(
        t_a[mask_a], y_a[mask_a],
        yerr=ye_a[mask_a],
        xerr=(np.where(np.isfinite(xlo_a[mask_a]), xlo_a[mask_a], 0.0),
              np.where(np.isfinite(xhi_a[mask_a]), xhi_a[mask_a], 0.0)),
        fmt="o", color=COLOR_AVG, ecolor=COLOR_AVG,
        elinewidth=2.0, capsize=4,
        markersize=8, mfc="white", mew=2.2, lw=1.8, zorder=8,
        label="Weighted average",
    )

    # ── dipole fit ────────────────────────────────────────────────────────
    res = _fit_dipole_reduced(avg_df, "combined", t_min=t_min, t_max=t_max)
    if res is not None:
        t_dense = np.linspace(max(0.0, t_min * 0.5), t_max * 1.05, 600)
        y_curve = _dipole_model(t_dense, *res["popt"])
        chi2_str = (rf"$\chi^2/\nu = {res['chi2']:.1f}/{res['ndf']}$"
                    if res["ndf"] > 0 else "")
        fit_txt = (
            rf"$A = {res['A']:.3g} \pm {res['dA']:.1g}$ nb"      "\n"
            rf"$b = {res['b']:.3f} \pm {res['db']:.3f}$ GeV$^{{-2}}$" "\n"
            rf"$t_0 = {res['t0']:.3f} \pm {res['dt0']:.3f}$ GeV$^2$"  "\n"
            + chi2_str
        )
        ax_m.plot(t_dense, y_curve,
                  color=COLOR_AVG, lw=2.5, ls="--", zorder=9,
                  label="Dipole fit (combined)")
        ax_m.axvspan(t_min, t_max, alpha=0.05, color="gray", zorder=0)
        ax_m.text(
            0.97, 0.97, fit_txt,
            transform=ax_m.transAxes,
            fontsize=14, va="top", ha="right", color=COLOR_AVG,
            bbox=dict(boxstyle="round,pad=0.35", fc="white",
                      ec=COLOR_AVG, alpha=0.88, lw=1.2),
        )
        # pull panel
        ax_r.axhline( 0.0, color="gray", lw=1.2, ls="--", zorder=1)
        ax_r.axhline( 1.0, color="gray", lw=0.8, ls=":",  zorder=1)
        ax_r.axhline(-1.0, color="gray", lw=0.8, ls=":",  zorder=1)
        ax_r.errorbar(
            res["t_fit"], res["pull"],
            yerr=np.ones_like(res["pull"]),
            fmt="o", color=COLOR_AVG, ecolor=COLOR_AVG,
            elinewidth=1.2, capsize=2,
            markersize=5, mfc="white", mew=1.8, zorder=5,
        )

    ax_m.set_xlim(XLIM_T)
    if logy:
        ax_m.set_yscale("log")
    if show_ylabel:
        ax_m.set_ylabel(r"$\sigma_\mathrm{red}\ [\mathrm{nb}]$", fontsize=20)
    ax_m.tick_params(labelbottom=False, labelsize=18)

    ax_r.set_xlim(XLIM_T)
    ax_r.set_ylim(-4.5, 4.5)
    ax_r.set_xlabel(r"$-t'\ [\mathrm{GeV}^2]$", fontsize=20)
    if show_ylabel:
        ax_r.set_ylabel("Pull", fontsize=17)
    ax_r.tick_params(labelsize=16)

    return res

def make_dipole_fit_plots(csv_root, models, bin_keys, outdir,
                          t_min=0.1, t_max=5.5, logy=True):
    """
    Two sets of output per (Q², W) bin:

    A) Per-model canvases  (unchanged — individual fit to each dataset)
       → dipole_reduced_Q{iq}[_W{iw}].pdf

    B) Combined canvas  (NEW — weighted average of all datasets + single fit)
       → dipole_combined_Q{iq}[_W{iw}].pdf

    The combined fit results (one row per kinematic bin) are used for
    the ⟨b²⟩ vs x_B plot and summary CSV.
    """
    os.makedirs(outdir, exist_ok=True)

    per_model_results  = []   # per-model fit rows  (for parameter-vs-Q² plots)
    combined_results   = []   # combined fit rows   (for b² vs x_B)

    for (iq, iw) in bin_keys:
        dfs_m = [(m, load_csv(csv_root, m, iq, iw)) for m in models]
        dfs_m = [(m, df) for m, df in dfs_m if df is not None]
        if not dfs_m:
            continue

        m_names = [m  for m, _  in dfs_m]
        dfs     = [df for _, df in dfs_m]
        n_m     = len(dfs_m)

        # ── kinematics from first available df ────────────────────────────
        ref_df = dfs[0]
        q2c = (float(ref_df["Q2_center"].iloc[0])
               if "Q2_center" in ref_df.columns else np.nan)
        wc  = (float(ref_df["W_center"].iloc[0])
               if "W_center"  in ref_df.columns else np.nan)
        bin_title = _bin_label(ref_df)

        # ── A) per-model fit canvases ─────────────────────────────────────
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
                           ("A","dA","b","db","t0","dt0","chi2","ndf")},
                    ))

            fig.suptitle(bin_title, fontsize=LEG_TITLE_FS + 2, y=0.96)
            stem = (f"dipole_reduced_Q{iq}" if iw is None
                    else f"dipole_reduced_Q{iq}_W{iw}")
            _save(fig, outdir, stem)

        # ── B) combined weighted-average canvas ───────────────────────────
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
                        else f"{n_ds_min}–{n_ds_max} datasets")
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
                   ("A","dA","b","db","t0","dt0","chi2","ndf")},
            ))

    # ── save CSVs ─────────────────────────────────────────────────────────
    if not per_model_results and not combined_results:
        print("  [WARN] No dipole fit results — skipping summary plots.")
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
        print(f"  → saved to {comb_path}")
    else:
        comb_tab = pd.DataFrame(per_model_results)   # fallback

    # ── parameter-vs-Q² / W plots (per-model) ─────────────────────────────
    if per_model_results:
        d_par = os.path.join(outdir, "fit_parameters")
        sumtab_all = pd.DataFrame(per_model_results)
        _make_param_vs_Q2_plots(sumtab_all, models, d_par)
        _make_param_vs_W_plots( sumtab_all, models, d_par)
        _make_param_summary_per_model(sumtab_all, models, d_par)

    # ── ⟨b²⟩ vs x_B — use combined fit results ────────────────────────────
    make_b2_vs_xB_plot(comb_tab, ["combined"],
                       os.path.join(outdir, "b2_vs_xB"))


# ── helper: group the summary table by W-bin and plot param vs Q² ────────────
def _make_param_vs_Q2_plots(sumtab, models, outdir):
    """b, A, t0 vs Q² for each W-bin and each model, with error bars."""
    os.makedirs(outdir, exist_ok=True)
    params = [
        ("b",  "db",  r"$b\ [\mathrm{GeV}^{-2}]$",   "slope_b"),
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
                        ms["Q2_center"], ms[par], yerr=ms[dpar],
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
    """b, A, t0 vs W for each Q²-bin and each model, with error bars."""
    os.makedirs(outdir, exist_ok=True)
    params = [
        ("b",  "db",  r"$b\ [\mathrm{GeV}^{-2}]$",   "slope_b"),
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
                        ms["W_center"], ms[par], yerr=ms[dpar],
                        fmt="o-", color=color, ecolor=color,
                        elinewidth=1.8, capsize=4,
                        markersize=8, mfc="white", mew=2.2, lw=1.8,
                        label=m.replace("_", " "), zorder=5 + im,
                    )

                ax.set_xlabel(r"$W\ [\mathrm{GeV}]$")
                ax.set_ylabel(ylabel)
                ax.set_title(
                    rf"Dipole fit parameter {ylabel.split('[')[0].strip()}"
                    rf" vs $W$  (Q²-bin {iq_key})",
                    fontsize=LEG_TITLE_FS,
                )
                ax.legend(frameon=False, fontsize=LEG_FS - 2)

                stem = f"{ptag}_vs_W_Q{iq_key}"
                _save(fig, outdir, stem)


def _make_param_summary_per_model(sumtab, models, outdir):
    """3-panel (A, b, t0) summary figure per model — all bins together."""
    os.makedirs(outdir, exist_ok=True)
    params = [
        ("A",  "dA",  r"$A\ [\mathrm{nb}]$"),
        ("b",  "db",  r"$b\ [\mathrm{GeV}^{-2}]$"),
        ("t0", "dt0", r"$t_0\ [\mathrm{GeV}^2]$"),
    ]

    for m in models:
        ms = sumtab[sumtab["model"] == m].sort_values(["iw", "Q2_center"])
        if ms.empty:
            continue

        color = MODEL_COLORS[models.index(m) % len(MODEL_COLORS)]

        # Build a human-readable bin label: "Q²∈[lo,hi]  W∈[lo,hi]"
        bin_labels = []
        for _, row in ms.iterrows():
            lbl = rf"$Q^2$={row['Q2_center']:.2f}"
            if pd.notna(row["W_center"]):
                lbl += rf"  $W$={row['W_center']:.2f}"
            bin_labels.append(lbl)

        x = np.arange(len(ms))

        with mpl.rc_context(PANEL_RC):
            fig, axes = plt.subplots(
                3, 1, figsize=(max(8, 1.2 * len(ms) + 3), 16),
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

            # χ²/ndf on top panel
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
                rf"Dipole fit parameters — {m.replace('_', ' ')}",
                fontsize=LEG_TITLE_FS + 2, y=0.96,
            )
            _save(fig, outdir, f"param_summary_{_sanitise(m)}")


# ─────────────────────────────────────────────────────────────────────────────
#  <b²>_g  vs  x_B  plot  (from dipole fit t0 parameter)
# ─────────────────────────────────────────────────────────────────────────────
# Conversion factor:  1 GeV⁻² = 0.389379 fm²
_GEVMINUSSQ_TO_FM2 = 0.389379

def _b2_from_t0(t0, dt0):
    """
    For a dipole form factor  F(t) = (1 - t/t0)^{-2},
    the mean transverse radius is:
        <b²> = 4 / t0   [in units of t0's inverse, i.e. GeV⁻²]
    Convert to fm²:
        <b²> [fm²] = (4 / t0 [GeV²]) * 0.389379
    Propagated error (dt0 → db2):
        db2 = (4 / t0²) * dt0 * 0.389379
    """
    b2  = (4.0 / t0)  * _GEVMINUSSQ_TO_FM2
    db2 = (4.0 / t0**2) * dt0 * _GEVMINUSSQ_TO_FM2
    return b2, db2


def make_b2_vs_xB_plot(sumtab, models, outdir):
    """
    Plot <b²>_g (fm²) vs x_B derived from the dipole fit t0 parameter.

    x_B is computed per bin from Q2_center and W_center (already in sumtab).
    One series per model, all W-bins overlaid (distinguished by marker shape
    if multiple W-bins exist).

    Parameters
    ----------
    sumtab : pd.DataFrame
        Output of make_dipole_fit_plots — must have columns
        model, Q2_center, W_center, t0, dt0.
    models : list of str
    outdir : str
    """
    os.makedirs(outdir, exist_ok=True)

    # ── compute x_B and <b²> for every row ──────────────────────────────────
    tab = sumtab.copy()

    # x_B = Q2 / (W² - Mp² + Q2)
    Q2 = tab["Q2_center"].to_numpy(float)
    W  = tab["W_center"].to_numpy(float)
    with np.errstate(invalid="ignore"):
        W2  = W * W
        xb  = Q2 / (W2 - M_P_GEV**2 + Q2)
    tab["xB"] = xb

    # <b²> and its error from t0
    b2_vals  = np.full(len(tab), np.nan)
    db2_vals = np.full(len(tab), np.nan)
    for i, row in tab.iterrows():
        t0  = float(row["t0"])
        dt0 = float(row["dt0"])
        if np.isfinite(t0) and t0 > 0:
            b2, db2 = _b2_from_t0(t0, dt0)
            b2_vals[i]  = b2
            db2_vals[i] = db2
    tab["b2"]  = b2_vals
    tab["db2"] = db2_vals

    # drop rows with no valid b² or x_B
    tab = tab[np.isfinite(tab["xB"]) & np.isfinite(tab["b2"])].copy()
    if tab.empty:
        print("  [WARN] make_b2_vs_xB_plot: no valid (xB, b²) pairs — skipping.")
        return

    # ── marker cycle for multiple W-bins ────────────────────────────────────
    iw_vals   = sorted(tab["iw"].unique(), key=lambda x: (x is None, x))
    markers   = ["o", "s", "^", "D", "v", "P", "*"]

    with mpl.rc_context(PANEL_RC):
        fig, ax = plt.subplots(figsize=(9, 6.5))
        _style_ax(ax)

        handles_legend, labels_legend = [], []

        for im, m in enumerate(models):
            color = MODEL_COLORS[im % len(MODEL_COLORS)]
            ms = tab[tab["model"] == m]
            if ms.empty:
                continue

            for iw_idx, iw_key in enumerate(iw_vals):
                mw = ms[ms["iw"] == iw_key].sort_values("xB")
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
                handles_legend.append(h)
                labels_legend.append(lbl)

                # connect points with a thin line to guide the eye
                ax.plot(
                    mw["xB"].to_numpy(float),
                    mw["b2"].to_numpy(float),
                    color=color, lw=1.2, ls="-", alpha=0.5, zorder=4,
                )

        # ── axis decoration ──────────────────────────────────────────────────
        ax.set_xlabel(r"$x_B$", fontsize=22)
        ax.set_ylabel(
            r"$\langle b^2 \rangle_g \ [\mathrm{fm}^2]$", fontsize=22,
        )
        ax.set_title(
            r"Mean-Square Gluonic Transverse Radius $\langle b^2\rangle_g$ vs $x_B$"
            "\n" r"(from dipole fit: $\langle b^2\rangle_g = 4/t_0 \times 0.389\ \mathrm{fm}^2$)",
            fontsize=15, pad=8,
        )

        # x-axis: log scale since x_B spans ~0.1–0.4 but may vary; keep linear if range small
        xb_arr = tab["xB"].to_numpy(float)
        xb_range = np.nanmax(xb_arr) / max(np.nanmin(xb_arr), 1e-6)
        if xb_range > 5:
            ax.set_xscale("log")

        # proton charge radius² reference line
        r_Ep_fm = 0.8414       # CODATA 2018 value, fm
        ax.axhline(r_Ep_fm**2, color="#008800", lw=1.5, ls="--", zorder=1,
                   label=rf"$r_{{Ep}}^2 = {r_Ep_fm**2:.3f}\ \mathrm{{fm}}^2$")

        ax.set_ylim(bottom=0.0)
        ax.legend(handles_legend + [ax.lines[-1]],
                  labels_legend  + [rf"$r_{{Ep}}^2 = {r_Ep_fm**2:.3f}\ \mathrm{{fm}}^2$"],
                  loc="upper right", frameon=True, framealpha=0.92,
                  edgecolor="#cccccc", fontsize=LEG_FS - 2)

        _save(fig, outdir, "b2_vs_xB")

    # ── also save the derived table ──────────────────────────────────────────
    out_cols = ["model", "iq", "iw", "Q2_center", "W_center",
                "xB", "b2", "db2", "t0", "dt0", "chi2", "ndf"]
    out_cols = [c for c in out_cols if c in tab.columns]
    tab[out_cols].to_csv(os.path.join(outdir, "b2_vs_xB_table.csv"), index=False)
    print("\n=== <b²> vs x_B table ===")
    print(tab[out_cols].to_string(index=False))


# ─────────────────────────────────────────────────────────────────────────────
#  main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description="φ dσ/dt analysis: corrections, acceptance, RadCorr, "
                    "final cross-section, dipole fits."
    )
    ap.add_argument("--csv-root",  default=".",
                    help="Root dir containing CSVs/ sub-tree.")
    ap.add_argument("--models",    nargs="+", required=True,
                    help='Model names, e.g. "Sp18_inb" "Sp19_inb"')
    ap.add_argument("--outdir",    default="phi_results",
                    help="Top-level output directory.")
    ap.add_argument("--luminosity", type=float, default=None,
                    help="Integrated luminosity [nb⁻¹] for variant reconstruction "
                         "from RawCounts. If omitted, scaling from CrossSection used.")
    ap.add_argument("--branching", type=float, default=0.492,
                    help="φ→K⁺K⁻ branching ratio.")
    ap.add_argument("--beam-energy", type=float, default=None,
                    help="Fallback beam energy E [GeV] for model names not matched "
                         "by auto-detection. Auto: Sp18→10.6, Fall18→10.594, Sp19→10.2.")
    ap.add_argument("--plot-reduced", action="store_true",
                    help="Also produce the group-overview reduced cross-section plots "
                         "(dipole fits on reduced XS are always produced).")
    ap.add_argument("--use-external-radcorr", action="store_true",
                    help="Overwrite RadCorr using DIFFRAD CSV tables.")
    ap.add_argument("--radcorr-file", default=None,
                    help="Override external RadCorr CSV path (skips auto-selection).")
    ap.add_argument("--print-table", action="store_true",
                    help="Print + save cross-check tables per (Q²,W) bin.")
    ap.add_argument("--table-outdir", default=None,
                    help="Where to save tables (default: <outdir>/0_Tables).")
    ap.add_argument("--no-logy", action="store_true",
                    help="Disable log-y on cross-section canvases.")
    ap.add_argument("--t-min", type=float, default=0.1,
                    help="Min |t′| for fitting [GeV²].")
    ap.add_argument("--t-max", type=float, default=5.5,
                    help="Max |t′| for fitting [GeV²].")
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
    table_outdir = args.table_outdir or os.path.join(outdir, "0_Tables")

    print("Beam energy assignment:")
    for m in models:
        E = beam_energy_for_model(m)
        print(f"  {m:20s} → {E} GeV  "
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

    d1 = os.path.join(outdir, "1_Corrections_beforeafter")
    for (iq, iw) in bin_keys:
        make_corrections_canvas(
            csv_root, present, iq, iw, d1, logy,
            luminosity=lumi, branching=br,
            print_table=args.print_table, table_outdir=table_outdir,
        )

    make_group_plot(csv_root, present, bin_keys, "acc",
                    os.path.join(outdir, "2_Acceptance_group"), logy=False)

    make_group_plot(csv_root, present, bin_keys, "rad",
                    os.path.join(outdir, "3_RadCorr_group"), logy=False)

    make_radcorr_vs_Q2(csv_root, present, bin_keys,
                       os.path.join(outdir, "4_RadCorr_vs_Q2"))

    d5xs  = os.path.join(outdir, "5_Comparison_xs")
    d5acc = os.path.join(outdir, "5_Comparison_acc")
    d5rad = os.path.join(outdir, "5_Comparison_rad")
    for (iq, iw) in bin_keys:
        make_comparison_with_ratio(csv_root, present, iq, iw, "xs",  d5xs,  logy)
        make_comparison_with_ratio(csv_root, present, iq, iw, "acc", d5acc, False)
        make_comparison_with_ratio(csv_root, present, iq, iw, "rad", d5rad, False)

    make_final_xs_group(csv_root, present, bin_keys,
                        os.path.join(outdir, "6_FinalCrossSection"), logy)

    if plot_red:
        make_reduced_xs_group(csv_root, present, bin_keys,
                              os.path.join(outdir, "6b_ReducedCrossSection_group"),
                              logy=logy)

    # ── Dipole fits on REDUCED cross-section (always run) ──────────────────
    make_dipole_fit_plots(
        csv_root, present, bin_keys,
        os.path.join(outdir, "7_Dipolefit_ReducedXS"),
        t_min=args.t_min, t_max=args.t_max, logy=logy,
    )

    print(f"\n[DONE] Outputs in {outdir}/")


if __name__ == "__main__":
    main()