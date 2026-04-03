#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_phi_rlt_from_csv.py
========================
Plot R = sigma_L / sigma_T and r^04_00 (SDME) vs -t' from the
RLT_summary.csv produced by DISANAcomparer::PlotPhiRLT_FromCache().

Also computes an inverse-variance weighted average across all datasets
per (Q2, W, t') bin, with a chi^2 consistency check.

CSV schema (one row per model / Q2 / W / t' bin):
  model, Q2_lo, Q2_hi, W_lo, W_hi, tprime_center,
  r04_00, r04_00_err, R, R_err

Usage examples
--------------
# Single CSV, all models, with CLAS6 reference band:
python3 plot_phi_rlt_from_csv.py \
    --csv PhiRLT/RLT_summary.csv --outdir RLT_plots --show-clas6

# Merge several CSVs, keep specific models:
python3 plot_phi_rlt_from_csv.py \
    --csv Sp18/PhiRLT/RLT_summary.csv Fall18/PhiRLT/RLT_summary.csv \
    --models Sp18_inb Fall18_inb \
    --outdir RLT_combined --show-clas6

# Skip weighted average (individual datasets only):
python3 plot_phi_rlt_from_csv.py --csv ... --no-wavg
"""

import os
import re
import argparse
from itertools import cycle

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# ---------------------------------------------------------------------------
#  Plot style  (mirrors plot_phi_dsdt_from_csv.py)
# ---------------------------------------------------------------------------
PANEL_RC = {
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 22,
    "axes.labelsize": 22,
    "axes.titlesize": 20,
    "legend.fontsize": 18,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "xtick.major.size": 12,
    "xtick.major.width": 2,
    "ytick.major.size": 12,
    "ytick.major.width": 2,
    "xtick.minor.size":  6,
    "xtick.minor.width": 2,
    "ytick.minor.size":  6,
    "ytick.minor.width": 2,
    "axes.linewidth": 2,
    "lines.linewidth": 2,
}

MODEL_COLORS  = ["#1F4FD8", "#E07010", "#009999", "#008800",
                 "#7733CC", "#CC1A33", "#666666", "#BB6600"]
MODEL_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "h"]

WAVG_COLOR  = "#000000"
WAVG_MARKER = "*"
WAVG_SIZE   = 14

LEG_FS       = 16
LEG_TITLE_FS = 18

CLAS6_R_VALUE = 0.7
CLAS6_R_ERR   = 0.2

# ---------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------

def _sanitise(s):
    return re.sub(r"[ /\\]", "_", s)

def _style_ax(ax):
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))

def _save(fig, outdir, stem):
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"{stem}.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"    [OK] {stem}")

def _bin_label(q2lo, q2hi, wlo, whi):
    lbl = rf"$Q^2 \in [{q2lo:.2f},\,{q2hi:.2f}]\ \mathrm{{GeV}}^2$"
    try:
        if np.isfinite(float(wlo)) and np.isfinite(float(whi)) and float(wlo) > 0:
            lbl += rf"$\quad W \in [{wlo:.2f},\,{whi:.2f}]\ \mathrm{{GeV}}$"
    except (TypeError, ValueError):
        pass
    return lbl

def _isnan_safe(x):
    try:
        return bool(np.isnan(x))
    except (TypeError, ValueError):
        return x is None

def _select_bin(df, q2lo, q2hi, wlo, whi):
    """Filter dataframe to one (Q2, W) bin, NaN-safe."""
    mask = np.isclose(df["Q2_lo"], q2lo) & np.isclose(df["Q2_hi"], q2hi)
    if not _isnan_safe(wlo) and not _isnan_safe(whi) and float(wlo) > 0:
        mask &= (np.isclose(df["W_lo"], wlo, equal_nan=True) &
                 np.isclose(df["W_hi"], whi, equal_nan=True))
    return df[mask]

# ---------------------------------------------------------------------------
#  CSV loading
# ---------------------------------------------------------------------------

def load_summary(csv_paths, model_filter):
    frames = []
    for path in csv_paths:
        if not os.path.isfile(path):
            print(f"[WARN] CSV not found: {path}")
            continue
        df = pd.read_csv(path, comment="#")
        print(f"  Loaded {len(df)} rows from {path}")
        frames.append(df)
    if not frames:
        raise FileNotFoundError("No valid CSV files found.")

    df = pd.concat(frames, ignore_index=True)

    required = {"model", "Q2_lo", "Q2_hi", "tprime_center",
                "r04_00", "r04_00_err", "R", "R_err"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing columns: {missing}")

    for col in ("W_lo", "W_hi"):
        if col not in df.columns:
            df[col] = np.nan

    if model_filter:
        before = len(df)
        df = df[df["model"].isin(model_filter)]
        print(f"  After model filter: {len(df)}/{before} rows "
              f"({sorted(df['model'].unique())})")

    if df.empty:
        raise ValueError("No data remaining after model filter.")
    return df

def discover_bins(df):
    cols = ["Q2_lo", "Q2_hi", "W_lo", "W_hi"]
    return list(df[cols].drop_duplicates().sort_values(cols).to_records(index=False))

# ---------------------------------------------------------------------------
#  Weighted average
# ---------------------------------------------------------------------------

def compute_weighted_average(df, obs_col, err_col):
    """
    Inverse-variance weighted average across all models per
    (Q2_lo, Q2_hi, W_lo, W_hi, tprime_center) bin.

    Weighted mean:
        y_avg = sum(y_i / sigma_i^2) / sum(1 / sigma_i^2)
        sigma_avg = 1 / sqrt( sum(1 / sigma_i^2) )

    Chi2 consistency (Birge ratio test):
        chi2 = sum( (y_i - y_avg)^2 / sigma_i^2 )
        ndf  = N - 1

    If chi2/ndf > 1 the datasets are in tension and the uncertainty
    can be inflated by sqrt(chi2/ndf) (Particle Data Group prescription).
    The inflated error is stored as <obs>_wavg_err_pdg.
    """
    group_cols = ["Q2_lo", "Q2_hi", "W_lo", "W_hi", "tprime_center"]
    rows = []

    for keys, grp in df.groupby(group_cols, dropna=False):
        y  = grp[obs_col].to_numpy(float)
        ye = grp[err_col].to_numpy(float)

        mask = np.isfinite(y) & np.isfinite(ye) & (ye > 0) & (y > 0)
        y, ye = y[mask], ye[mask]
        n = len(y)
        if n == 0:
            continue

        w     = 1.0 / ye**2
        wsum  = w.sum()
        y_avg = (w * y).sum() / wsum
        y_err = 1.0 / np.sqrt(wsum)          # stat. uncertainty

        if n > 1:
            chi2 = float(np.sum(w * (y - y_avg)**2))
            ndf  = n - 1
            # PDG scale factor: inflate if chi2/ndf > 1
            scale = max(1.0, np.sqrt(chi2 / ndf))
        else:
            chi2  = np.nan
            ndf   = 0
            scale = 1.0

        if isinstance(keys, tuple):
            row = dict(zip(group_cols, keys))
        else:
            row = {group_cols[0]: keys}

        row[f"{obs_col}_wavg"]         = y_avg
        row[f"{obs_col}_wavg_err"]     = y_err           # stat only
        row[f"{obs_col}_wavg_err_pdg"] = y_err * scale   # PDG inflated
        row[f"{obs_col}_chi2"]         = chi2
        row[f"{obs_col}_ndf"]          = ndf
        row["n_datasets"]              = n
        rows.append(row)

    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(group_cols).reset_index(drop=True)


def save_weighted_average_csv(df_R, df_r04, outdir):
    group_cols = ["Q2_lo", "Q2_hi", "W_lo", "W_hi", "tprime_center"]
    os.makedirs(outdir, exist_ok=True)

    if df_R.empty and df_r04.empty:
        return
    if df_R.empty:
        merged = df_r04
    elif df_r04.empty:
        merged = df_R
    else:
        merged = pd.merge(df_R, df_r04,
                          on=group_cols + ["n_datasets"], how="outer")

    path = os.path.join(outdir, "RLT_weighted_average.csv")
    merged.to_csv(path, index=False, float_format="%.6e")
    print(f"\n[CSV] Weighted average → {path}")
    print(merged.to_string(index=False))


# ---------------------------------------------------------------------------
#  Per-bin plot  (models + weighted average overlay)
# ---------------------------------------------------------------------------

def plot_bin(df_bin, models, obs_col, err_col, ylabel, yrange,
             bin_label, outdir, stem,
             df_wavg_bin=None, show_clas6=False, draw_zero_line=False):

    with mpl.rc_context(PANEL_RC):
        fig, ax = plt.subplots(figsize=(10, 8),
                               gridspec_kw=dict(left=0.13, right=0.97,
                                                top=0.91, bottom=0.11))
        _style_ax(ax)
        cc = cycle(MODEL_COLORS)
        mc = cycle(MODEL_MARKERS)

        plotted = 0
        for model in models:
            dfm = df_bin[df_bin["model"] == model].sort_values("tprime_center")
            if dfm.empty:
                continue
            tp = dfm["tprime_center"].to_numpy(float)
            y  = dfm[obs_col].to_numpy(float)
            ye = dfm[err_col].to_numpy(float)
            xlo = xhi = np.zeros_like(tp)
            if "tprime_lo" in dfm.columns and "tprime_hi" in dfm.columns:
                xlo = tp - dfm["tprime_lo"].to_numpy(float)
                xhi = dfm["tprime_hi"].to_numpy(float) - tp
            mask = np.isfinite(tp) & np.isfinite(y) & (y > 0)
            if not mask.any():
                continue
            col = next(cc); mrk = next(mc)
            ax.errorbar(tp[mask], y[mask], yerr=ye[mask],
                        xerr=(xlo[mask], xhi[mask]),
                        fmt=mrk, color=col, ecolor=col,
                        elinewidth=1.8, capsize=4,
                        markersize=8, mfc="white", mew=2.2,
                        label=model.replace("_", " "), zorder=5)
            plotted += 1

        if plotted == 0:
            plt.close(fig)
            return

        # Weighted average overlay
        if df_wavg_bin is not None and not df_wavg_bin.empty:
            wav_col  = f"{obs_col}_wavg"
            werr_col = f"{obs_col}_wavg_err_pdg"   # PDG-inflated
            chi2_col = f"{obs_col}_chi2"
            ndf_col  = f"{obs_col}_ndf"

            if wav_col in df_wavg_bin.columns:
                dfw = df_wavg_bin.sort_values("tprime_center")
                tp_w  = dfw["tprime_center"].to_numpy(float)
                y_w   = dfw[wav_col].to_numpy(float)
                ye_w  = dfw[werr_col].to_numpy(float)
                mw = np.isfinite(tp_w) & np.isfinite(y_w) & (y_w > 0)
                if mw.any():
                    ax.errorbar(tp_w[mw], y_w[mw], yerr=ye_w[mw],
                                fmt=WAVG_MARKER, color=WAVG_COLOR,
                                ecolor=WAVG_COLOR,
                                elinewidth=2.5, capsize=5,
                                markersize=WAVG_SIZE, mfc=WAVG_COLOR,
                                mew=1.5, label="Weighted avg. (PDG)", zorder=10)

                    # chi2/ndf annotation
                    if chi2_col in dfw.columns and ndf_col in dfw.columns:
                        chi2v = dfw.loc[mw, chi2_col].to_numpy(float)
                        ndfv  = dfw.loc[mw, ndf_col].to_numpy(float)
                        valid = np.isfinite(chi2v) & (ndfv > 0)
                        if valid.any():
                            mean_chi2ndf = np.mean(chi2v[valid] / ndfv[valid])
                            ax.text(0.03, 0.96,
                                    f"$\\langle\\chi^2/\\mathrm{{ndf}}\\rangle = {mean_chi2ndf:.2f}$",
                                    transform=ax.transAxes, fontsize=15,
                                    va="top", ha="left", color="#444444",
                                    bbox=dict(boxstyle="round,pad=0.3",
                                              fc="white", ec="#AAAAAA", alpha=0.85))

        if show_clas6 and obs_col == "R":
            ax.axhspan(CLAS6_R_VALUE - CLAS6_R_ERR,
                       CLAS6_R_VALUE + CLAS6_R_ERR,
                       color="#AAAAAA", alpha=0.25, zorder=0,
                       label=f"CLAS6 R={CLAS6_R_VALUE}±{CLAS6_R_ERR}")

        if draw_zero_line:
            ax.axhline(0.0, color="gray", lw=1.2, ls="--", zorder=0)

        ax.set_xlabel(r"$-t'\ [\mathrm{GeV}^2]$")
        ax.set_ylabel(ylabel)
        ax.set_ylim(*yrange)
        ax.legend(loc="upper right", framealpha=0.85,
                  edgecolor="#999999", fontsize=LEG_FS)
        fig.suptitle(bin_label, fontsize=LEG_TITLE_FS, y=0.97)
        _save(fig, outdir, stem)


# ---------------------------------------------------------------------------
#  Weighted-average-only plot with chi2/ndf colorbar
# ---------------------------------------------------------------------------

def plot_wavg_only(df_wavg, bins, obs_col, ylabel, yrange, outdir, stem_prefix,
                   show_clas6=False):
    wav_col  = f"{obs_col}_wavg"
    werr_col = f"{obs_col}_wavg_err_pdg"
    chi2_col = f"{obs_col}_chi2"
    ndf_col  = f"{obs_col}_ndf"

    if wav_col not in df_wavg.columns:
        return

    cmap = mpl.cm.RdYlGn_r
    norm = mpl.colors.Normalize(vmin=0.0, vmax=3.0)

    for ib, (q2lo, q2hi, wlo, whi) in enumerate(bins):
        dfw = _select_bin(df_wavg, q2lo, q2hi, wlo, whi).sort_values("tprime_center")
        if dfw.empty:
            continue

        tp   = dfw["tprime_center"].to_numpy(float)
        y    = dfw[wav_col].to_numpy(float)
        ye   = dfw[werr_col].to_numpy(float)
        mask = np.isfinite(tp) & np.isfinite(y) & (y > 0)
        if not mask.any():
            continue

        has_w = not _isnan_safe(wlo) and float(wlo) > 0
        stem  = f"{stem_prefix}_Q{ib}" + (f"_W{ib}" if has_w else "")

        # Compute chi2/ndf per point for colour coding
        if chi2_col in dfw.columns and ndf_col in dfw.columns:
            chi2v = dfw.loc[mask, chi2_col].to_numpy(float)
            ndfv  = dfw.loc[mask, ndf_col].to_numpy(float)
            with np.errstate(invalid="ignore", divide="ignore"):
                chi2ndf = np.where(ndfv > 0, chi2v / ndfv, np.nan)
            colors = [cmap(norm(c)) if np.isfinite(c) else "#888888"
                      for c in chi2ndf]
            have_chi2 = not np.all(np.isnan(chi2ndf))
        else:
            colors = [WAVG_COLOR] * int(mask.sum())
            have_chi2 = False

        with mpl.rc_context(PANEL_RC):
            fig, ax = plt.subplots(
                figsize=(10, 8),
                gridspec_kw=dict(left=0.13, right=0.93 if have_chi2 else 0.97,
                                 top=0.91, bottom=0.11))
            _style_ax(ax)

            for x_, y_, ye_, c_ in zip(tp[mask], y[mask], ye[mask], colors):
                ax.errorbar(x_, y_, yerr=ye_,
                            fmt=WAVG_MARKER, color=c_, ecolor=c_,
                            elinewidth=2.5, capsize=5,
                            markersize=WAVG_SIZE, mfc=c_, mew=1.5, zorder=5)

            if have_chi2:
                sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
                sm.set_array([])
                cbar = fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.035)
                cbar.set_label(r"$\chi^2/\mathrm{ndf}$  (dataset consistency)",
                               fontsize=16)
                cbar.ax.tick_params(labelsize=14)

            if show_clas6 and obs_col == "R":
                ax.axhspan(CLAS6_R_VALUE - CLAS6_R_ERR,
                           CLAS6_R_VALUE + CLAS6_R_ERR,
                           color="#AAAAAA", alpha=0.25, zorder=0,
                           label=f"CLAS6 R={CLAS6_R_VALUE}±{CLAS6_R_ERR}")
                ax.legend(fontsize=LEG_FS)

            ax.set_xlabel(r"$-t'\ [\mathrm{GeV}^2]$")
            ax.set_ylabel(ylabel)
            ax.set_ylim(*yrange)
            fig.suptitle(_bin_label(q2lo, q2hi, wlo, whi),
                         fontsize=LEG_TITLE_FS, y=0.97)
            _save(fig, outdir, stem)


# ---------------------------------------------------------------------------
#  Multi-panel summary
# ---------------------------------------------------------------------------

def plot_summary_panel(df, models, bins, obs_col, err_col,
                       ylabel, yrange, outdir, stem,
                       df_wavg=None, show_clas6=False):
    n = len(bins)
    if n == 0:
        return

    ncols = min(3, n)
    nrows = (n + ncols - 1) // ncols
    wav_col  = f"{obs_col}_wavg"
    werr_col = f"{obs_col}_wavg_err_pdg"

    color_map  = {m: MODEL_COLORS[i % len(MODEL_COLORS)] for i, m in enumerate(models)}
    marker_map = {m: MODEL_MARKERS[i % len(MODEL_MARKERS)] for i, m in enumerate(models)}

    with mpl.rc_context(PANEL_RC):
        fig, axes = plt.subplots(
            nrows, ncols, figsize=(7 * ncols, 6 * nrows),
            squeeze=False,
            gridspec_kw=dict(hspace=0.38, wspace=0.30,
                             left=0.10, right=0.98,
                             top=0.93, bottom=0.09))
        axes_flat = axes.flatten()

        for idx, (q2lo, q2hi, wlo, whi) in enumerate(bins):
            ax = axes_flat[idx]
            _style_ax(ax)
            df_bin  = _select_bin(df, q2lo, q2hi, wlo, whi)
            plotted = 0

            for model in models:
                dfm = df_bin[df_bin["model"] == model].sort_values("tprime_center")
                if dfm.empty:
                    continue
                tp   = dfm["tprime_center"].to_numpy(float)
                y    = dfm[obs_col].to_numpy(float)
                ye   = dfm[err_col].to_numpy(float)
                mask = np.isfinite(tp) & np.isfinite(y) & (y > 0)
                if not mask.any():
                    continue
                ax.errorbar(tp[mask], y[mask], yerr=ye[mask],
                            fmt=marker_map[model], color=color_map[model],
                            ecolor=color_map[model],
                            elinewidth=1.5, capsize=3,
                            markersize=6, mfc="white", mew=2.0,
                            label=model.replace("_", " "), zorder=5)
                plotted += 1

            if df_wavg is not None and not df_wavg.empty and wav_col in df_wavg.columns:
                dfw = _select_bin(df_wavg, q2lo, q2hi, wlo, whi).sort_values("tprime_center")
                tp_w  = dfw["tprime_center"].to_numpy(float)
                y_w   = dfw[wav_col].to_numpy(float)
                ye_w  = dfw[werr_col].to_numpy(float)
                mw = np.isfinite(tp_w) & np.isfinite(y_w) & (y_w > 0)
                if mw.any():
                    ax.errorbar(tp_w[mw], y_w[mw], yerr=ye_w[mw],
                                fmt=WAVG_MARKER, color=WAVG_COLOR,
                                ecolor=WAVG_COLOR,
                                elinewidth=2.0, capsize=4,
                                markersize=WAVG_SIZE, mfc=WAVG_COLOR,
                                mew=1.5, label="Weighted avg.", zorder=10)

            if show_clas6 and obs_col == "R" and plotted > 0:
                ax.axhspan(CLAS6_R_VALUE - CLAS6_R_ERR,
                           CLAS6_R_VALUE + CLAS6_R_ERR,
                           color="#AAAAAA", alpha=0.20, zorder=0)

            ax.set_ylim(*yrange)
            ax.set_xlabel(r"$-t'\ [\mathrm{GeV}^2]$", fontsize=17)
            ax.set_ylabel(ylabel, fontsize=17)
            ax.set_title(_bin_label(q2lo, q2hi, wlo, whi), fontsize=14)
            if plotted > 0 and idx == 0:
                ax.legend(fontsize=12, framealpha=0.85)

        for idx in range(len(bins), len(axes_flat)):
            axes_flat[idx].set_visible(False)

        fig.suptitle(ylabel, fontsize=LEG_TITLE_FS + 2, y=0.98)
        _save(fig, outdir, stem)


# ---------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Plot R=sigma_L/sigma_T and r04_00 vs -t', "
                    "with inverse-variance weighted average across datasets.")
    ap.add_argument("--csv", nargs="+", required=True,
                    help="Path(s) to RLT_summary.csv file(s).")
    ap.add_argument("--models", nargs="+", default=None,
                    help="Subset of model names to include (default: all).")
    ap.add_argument("--outdir", default="RLT_plots",
                    help="Output directory for plots and CSV.")
    ap.add_argument("--show-clas6", action="store_true",
                    help="Overlay CLAS6 reference band on R plots.")
    ap.add_argument("--R-ymax", type=float, default=5.0,
                    help="Upper y-limit for R plots (default 5.0).")
    ap.add_argument("--no-wavg", action="store_true",
                    help="Skip weighted average computation and plots.")
    ap.add_argument("--no-panel", action="store_true",
                    help="Skip the multi-panel summary figures.")
    args = ap.parse_args()

    print("=" * 60)
    print("  plot_phi_rlt_from_csv.py")
    print("=" * 60)

    df     = load_summary(args.csv, args.models)
    models = args.models if args.models else sorted(df["model"].unique())
    bins   = discover_bins(df)

    print(f"\nModels  : {models}")
    print(f"Bins    : {len(bins)} (Q2, W) combination(s)")
    print(f"Outdir  : {args.outdir}\n")

    # Weighted averages
    df_wavg_R = df_wavg_r04 = pd.DataFrame()
    if not args.no_wavg:
        print("--- Computing inverse-variance weighted averages ---")
        df_wavg_R   = compute_weighted_average(df, "R",      "R_err")
        df_wavg_r04 = compute_weighted_average(df, "r04_00", "r04_00_err")
        if not df_wavg_R.empty or not df_wavg_r04.empty:
            save_weighted_average_csv(
                df_wavg_R, df_wavg_r04,
                os.path.join(args.outdir, "3_WeightedAverage"))

    obs_cfg = [
        dict(obs_col="R",      err_col="R_err",
             ylabel=r"$R = \sigma_L / \sigma_T$",
             yrange=(0.0, args.R_ymax), tag="R",
             df_wavg=df_wavg_R, show_clas6=args.show_clas6),
        dict(obs_col="r04_00", err_col="r04_00_err",
             ylabel=r"$r^{04}_{00}$",
             yrange=(0.0, 1.0), tag="r04_00",
             df_wavg=df_wavg_r04, show_clas6=False),
    ]

    # Per-bin plots
    for cfg in obs_cfg:
        print(f"\n--- Plotting {cfg['tag']} ---")
        per_bin_dir = os.path.join(args.outdir, f"1_PerBin_{cfg['tag']}")
        for ib, (q2lo, q2hi, wlo, whi) in enumerate(bins):
            has_w  = not _isnan_safe(wlo) and float(wlo) > 0
            stem   = f"{cfg['tag']}_Q{ib}" + (f"_W{ib}" if has_w else "")
            db     = _select_bin(df, q2lo, q2hi, wlo, whi)
            dfw_b  = (_select_bin(cfg["df_wavg"], q2lo, q2hi, wlo, whi)
                      if not cfg["df_wavg"].empty else None)
            plot_bin(df_bin=db, models=models,
                     obs_col=cfg["obs_col"], err_col=cfg["err_col"],
                     ylabel=cfg["ylabel"], yrange=cfg["yrange"],
                     bin_label=_bin_label(q2lo, q2hi, wlo, whi),
                     outdir=per_bin_dir, stem=stem,
                     df_wavg_bin=dfw_b,
                     show_clas6=cfg["show_clas6"])

    # Weighted-average-only plots
    if not args.no_wavg:
        for cfg in obs_cfg:
            if cfg["df_wavg"].empty:
                continue
            print(f"\n--- Weighted-average plot: {cfg['tag']} ---")
            plot_wavg_only(
                df_wavg=cfg["df_wavg"], bins=bins,
                obs_col=cfg["obs_col"], ylabel=cfg["ylabel"],
                yrange=cfg["yrange"],
                outdir=os.path.join(args.outdir, "3_WeightedAverage"),
                stem_prefix=f"wavg_{cfg['tag']}",
                show_clas6=cfg["show_clas6"])

    # Multi-panel summary
    if not args.no_panel:
        for cfg in obs_cfg:
            print(f"\n--- Summary panel: {cfg['tag']} ---")
            plot_summary_panel(
                df=df, models=models, bins=bins,
                obs_col=cfg["obs_col"], err_col=cfg["err_col"],
                ylabel=cfg["ylabel"], yrange=cfg["yrange"],
                outdir=os.path.join(args.outdir, "2_Summary"),
                stem=f"summary_{cfg['tag']}",
                df_wavg=cfg["df_wavg"] if not args.no_wavg else None,
                show_clas6=cfg["show_clas6"])

    print(f"\n[DONE] All outputs saved to {args.outdir}/")


if __name__ == "__main__":
    main()