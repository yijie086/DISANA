#!/usr/bin/env python3
"""
plot_radcorr.py
---------------
Reads the phi_radcorr_*.root cache files produced by DISANA_PhiDiffrad
and plots rad_corr vs mtprime for every (Q2, W) bin.

Usage:
    python3 plot_radcorr.py --files /path/to/phi_radcorr_rgasp19_inb.root \
                                    /path/to/phi_radcorr_rgafall18_inb.root \
                            --tree  radcorr \
                            --q2bins 0.9 1.2 1.733 2.533 8.0 \
                            --wbins  2.0 10.0 \
                            --outdir RadCorrPlots

All arguments have sensible defaults — edit the DEFAULTS section below
if you prefer not to pass flags every time.
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ── optional: uproot for pure-Python ROOT reading ────────────────────────────
try:
    import uproot
    HAS_UPROOT = True
except ImportError:
    HAS_UPROOT = False

# ── optional: PyROOT fallback ─────────────────────────────────────────────────
try:
    import ROOT
    ROOT.gROOT.SetBatch(True)
    HAS_ROOT = True
except ImportError:
    HAS_ROOT = False

if not HAS_UPROOT and not HAS_ROOT:
    sys.exit("ERROR: neither uproot nor PyROOT is available. "
             "Install uproot with:  pip install uproot awkward")

# =============================================================================
# DEFAULTS — edit here instead of passing CLI flags every time
# =============================================================================
DEFAULT_FILES   = []          # e.g. ["phi_radcorr_rgasp19_inb.root"]
DEFAULT_TREE    = "radcorr"
DEFAULT_Q2BINS  = [0.9, 1.2, 1.733, 2.533, 8.0]   # GeV^2 edges
DEFAULT_WBINS   = [2.0, 10.0]                        # GeV edges (1 bin = no W split)
DEFAULT_OUTDIR  = "RadCorrPlots"
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# I/O helpers
# ─────────────────────────────────────────────────────────────────────────────

def load_tree_uproot(path, tree):
    with uproot.open(path) as f:
        t = f[tree]
        return {k: t[k].array(library="np") for k in t.keys()}

def load_tree_pyroot(path, tree):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {path}")
    t = f.Get(tree)
    if not t:
        raise RuntimeError(f"Tree '{tree}' not found in {path}")
    arrays = {}
    for br in t.GetListOfBranches():
        name = br.GetName()
        arrays[name] = []
    for entry in t:
        for name in arrays:
            arrays[name].append(getattr(entry, name))
    f.Close()
    return {k: np.array(v) for k, v in arrays.items()}

def load_tree(path, tree):
    if HAS_UPROOT:
        return load_tree_uproot(path, tree)
    return load_tree_pyroot(path, tree)


# ─────────────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────────────

# Distinct colours / markers for multiple files
COLORS  = ["#e6194b","#3cb44b","#4363d8","#f58231","#911eb4",
           "#42d4f4","#f032e6","#bfef45","#fabed4","#469990"]
MARKERS = ["o","s","^","D","v","P","X","*","h","8"]

def label_from_path(path):
    """Extract a short label from the file path."""
    base = os.path.basename(path)
    # strip phi_radcorr_ prefix and .root suffix if present
    label = base.replace("phi_radcorr_", "").replace(".root", "")
    return label


def plot_radcorr(all_data, q2bins, wbins, outdir):
    """
    all_data : list of (label, arrays_dict)
    q2bins   : Q^2 bin edges [GeV^2]
    wbins    : W   bin edges [GeV]
    outdir   : output directory for PDFs/PNGs
    """
    os.makedirs(outdir, exist_ok=True)

    nQ = len(q2bins) - 1
    nW = len(wbins)  - 1
    hasW = nW > 1

    # ── one figure per (Q2, W) bin ──────────────────────────────────────────
    for iq in range(nQ):
        q2lo, q2hi = q2bins[iq], q2bins[iq + 1]

        for iw in range(nW):
            wlo = wbins[iw]
            whi = wbins[iw + 1]

            fig, ax = plt.subplots(figsize=(7, 5))

            any_data = False

            for idx, (label, arr) in enumerate(all_data):
                Q2  = arr.get("Q2",  None)
                W   = arr.get("W",   None)
                mtp = arr.get("mtprime",  None)
                rc  = arr.get("rad_corr", None)
                rce = arr.get("rad_corr_err", None)

                if Q2 is None or mtp is None or rc is None:
                    print(f"  [WARN] {label}: missing required columns, skipping.")
                    continue

                # ── select this (Q2,W) bin ──────────────────────────────────
                mask = (Q2 > q2lo) & (Q2 <= q2hi)
                if W is not None and hasW:
                    mask &= (W > wlo) & (W <= whi)

                mtp_sel = mtp[mask]
                rc_sel  = rc[mask]
                rce_sel = rce[mask] if rce is not None else None

                if len(mtp_sel) == 0:
                    continue

                # sort by mtprime for a clean line
                order = np.argsort(mtp_sel)
                mtp_sel = mtp_sel[order]
                rc_sel  = rc_sel[order]
                if rce_sel is not None:
                    rce_sel = rce_sel[order]

                color  = COLORS[idx % len(COLORS)]
                marker = MARKERS[idx % len(MARKERS)]

                if rce_sel is not None:
                    ax.errorbar(mtp_sel, rc_sel, yerr=rce_sel,
                                fmt=marker+"-", color=color,
                                markersize=6, linewidth=1.4,
                                capsize=3, label=label)
                else:
                    ax.plot(mtp_sel, rc_sel,
                            marker+"-", color=color,
                            markersize=6, linewidth=1.4, label=label)
                any_data = True

            # ── reference line at 1 ─────────────────────────────────────────
            ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.9,
                       label="No correction (= 1)")

            # ── labels & cosmetics ──────────────────────────────────────────
            q2_title = f"$Q^2 \\in [{q2lo:.3g},\\ {q2hi:.3g}]$ GeV$^2$"
            w_title  = (f"$W \\in [{wlo:.2g},\\ {whi:.2g}]$ GeV"
                        if hasW else "")
            title = q2_title + (f",  {w_title}" if w_title else "")

            ax.set_title(title, fontsize=11)
            ax.set_xlabel("$-t'$ (= $|t - t_{\\min}|$)  [GeV$^2$]", fontsize=11)
            ax.set_ylabel("$C_{\\rm rad}$", fontsize=12)
            ax.legend(fontsize=8, framealpha=0.7)
            ax.grid(True, alpha=0.3)
            ax.set_ylim(bottom=0)

            if not any_data:
                ax.text(0.5, 0.5, "No data in this bin",
                        ha="center", va="center", transform=ax.transAxes,
                        fontsize=13, color="gray")

            fig.tight_layout()

            tag = f"Q{iq}_W{iw}"
            png_path = os.path.join(outdir, f"radcorr_{tag}.png")
            pdf_path = os.path.join(outdir, f"radcorr_{tag}.pdf")
            fig.savefig(png_path, dpi=150)
            fig.savefig(pdf_path)
            plt.close(fig)
            print(f"  Saved: {png_path}")

    # ── summary canvas: all (Q2,W) bins on one page ─────────────────────────
    fig_all, axes = plt.subplots(nQ, max(nW, 1),
                                 figsize=(5 * max(nW, 1), 4 * nQ),
                                 squeeze=False)

    for iq in range(nQ):
        q2lo, q2hi = q2bins[iq], q2bins[iq + 1]

        for iw in range(nW):
            wlo = wbins[iw]
            whi = wbins[iw + 1]
            ax  = axes[iq][iw]

            for idx, (label, arr) in enumerate(all_data):
                Q2  = arr.get("Q2",  None)
                W   = arr.get("W",   None)
                mtp = arr.get("mtprime",  None)
                rc  = arr.get("rad_corr", None)
                rce = arr.get("rad_corr_err", None)

                if Q2 is None or mtp is None or rc is None:
                    continue

                mask = (Q2 > q2lo) & (Q2 <= q2hi)
                if W is not None and hasW:
                    mask &= (W > wlo) & (W <= whi)

                mtp_sel = mtp[mask]
                rc_sel  = rc[mask]
                rce_sel = rce[mask] if rce is not None else None

                if len(mtp_sel) == 0:
                    continue

                order   = np.argsort(mtp_sel)
                mtp_sel = mtp_sel[order]
                rc_sel  = rc_sel[order]
                if rce_sel is not None:
                    rce_sel = rce_sel[order]

                color  = COLORS[idx % len(COLORS)]
                marker = MARKERS[idx % len(MARKERS)]

                if rce_sel is not None:
                    ax.errorbar(mtp_sel, rc_sel, yerr=rce_sel,
                                fmt=marker+"-", color=color,
                                markersize=4, linewidth=1.2,
                                capsize=2, label=label)
                else:
                    ax.plot(mtp_sel, rc_sel, marker+"-",
                            color=color, markersize=4,
                            linewidth=1.2, label=label)

            ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8)
            ax.set_title(
                f"$Q^2\\in[{q2lo:.2g},{q2hi:.2g}]$"
                + (f"\n$W\\in[{wlo:.2g},{whi:.2g}]$" if hasW else ""),
                fontsize=8)
            ax.set_xlabel("$-t'$ [GeV$^2$]", fontsize=7)
            ax.set_ylabel("$C_{\\rm rad}$", fontsize=8)
            ax.grid(True, alpha=0.25)
            ax.set_ylim(bottom=0)
            if iq == 0 and iw == 0:
                ax.legend(fontsize=6, framealpha=0.7)

    fig_all.suptitle("Radiative Correction $C_{\\rm rad}$ vs $-t'$",
                     fontsize=13, y=1.01)
    fig_all.tight_layout()
    summary_path = os.path.join(outdir, "radcorr_summary.pdf")
    fig_all.savefig(summary_path, bbox_inches="tight")
    plt.close(fig_all)
    print(f"\n  Summary saved: {summary_path}")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--files",   nargs="+", default=DEFAULT_FILES,
                   help="One or more phi_radcorr_*.root cache files")
    p.add_argument("--tree",    default=DEFAULT_TREE,
                   help="TTree name inside the ROOT files (default: radcorr)")
    p.add_argument("--q2bins",  nargs="+", type=float, default=DEFAULT_Q2BINS,
                   help="Q^2 bin edges in GeV^2")
    p.add_argument("--wbins",   nargs="+", type=float, default=DEFAULT_WBINS,
                   help="W bin edges in GeV")
    p.add_argument("--outdir",  default=DEFAULT_OUTDIR,
                   help="Output directory for plots")
    return p.parse_args()


def main():
    args = parse_args()

    if not args.files:
        sys.exit("ERROR: no ROOT files specified. "
                 "Use --files or edit DEFAULT_FILES in the script.")

    print(f"Loading {len(args.files)} file(s) from tree '{args.tree}' ...")
    all_data = []
    for path in args.files:
        if not os.path.exists(path):
            print(f"  [WARN] File not found, skipping: {path}")
            continue
        label = label_from_path(path)
        print(f"  Reading: {path}  →  label='{label}'")
        arr = load_tree(path, args.tree)
        print(f"    Columns : {list(arr.keys())}")
        print(f"    Entries : {len(next(iter(arr.values())))}")
        all_data.append((label, arr))

    if not all_data:
        sys.exit("ERROR: no files could be loaded.")

    print(f"\nPlotting  nQ={len(args.q2bins)-1}  nW={len(args.wbins)-1}  "
          f"bins → {args.outdir}/")
    plot_radcorr(all_data, args.q2bins, args.wbins, args.outdir)
    print("\nDone.")


if __name__ == "__main__":
    main()