from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

src = Path('datapoint_SDHEP_DIS_Cross_Section.txt')
out = Path('./')
out.mkdir(parents=True, exist_ok=True)

groups = defaultdict(list)
with src.open() as f:
    next(f)
    for line in f:
        fields = line.rstrip().split('\t')
        values = list(map(float, fields[:12]))
        key = tuple(values[6:12])
        groups[key].append(values[:6])


def draw(ax, key, rows):
    data = np.asarray(sorted(rows, key=lambda r: r[3]))
    phi, alu, err = data[:, 3], data[:, 4], data[:, 5]
    finite = np.isfinite(phi) & np.isfinite(alu) & np.isfinite(err) & (err > 0)
    ax.errorbar(phi[finite], alu[finite], yerr=err[finite], fmt='o', ms=3.2,
                lw=0.8, capsize=1.5, color='#1769aa', ecolor='#6f8797')
    ax.axhline(0, color='0.45', lw=0.7)
    a = b = a_err = b_err = np.nan
    if np.count_nonzero(finite) >= 2:
        cos_phi = np.cos(np.deg2rad(phi[finite]))
        design = np.column_stack((np.ones_like(cos_phi), cos_phi))
        weighted_design = design / err[finite, None]
        weighted_alu = alu[finite] / err[finite]
        if np.linalg.matrix_rank(weighted_design) == 2:
            covariance = np.linalg.inv(weighted_design.T @ weighted_design)
            a, b = covariance @ (weighted_design.T @ weighted_alu)
            a_err, b_err = np.sqrt(np.diag(covariance))
            grid = np.linspace(0, 360, 361)
            grid_design = np.column_stack(
                (np.ones_like(grid), np.cos(np.deg2rad(grid))))
            curve = grid_design @ np.array([a, b])
            curve_error = np.sqrt(np.einsum(
                'ij,jk,ik->i', grid_design, covariance, grid_design))
            ax.fill_between(grid, curve - curve_error, curve + curve_error,
                            color='#c23b22', alpha=0.22, linewidth=0)
            ax.plot(grid, curve, color='#c23b22', lw=1.2)
    xi0, xi1, c0, c1, t0, t1 = key
    fit_text = (rf'$a={a:.3f}\pm {a_err:.3f},\ b={b:.3f}\pm {b_err:.3f}$'
                if np.all(np.isfinite([a, b, a_err, b_err]))
                else 'fit unavailable')
    ax.set_title(rf'$\xi$ [{xi0:g},{xi1:g}), $\cos\theta$ [{c0:g},{c1:g})' + '\n' +
                 rf'$|t|$ [{t0:g},{t1:g}), ' + fit_text, fontsize=8)
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.grid(alpha=0.18, lw=0.5)
    return a, b, a_err, b_err


ordered = sorted(groups, key=lambda k: (k[4], k[0], k[2]))
global_xi_ranges = sorted({(k[0], k[1]) for k in ordered})
# Row zero is the page top, so descending order makes cos(theta) increase
# from the bottom of every page to the top.
global_cos_ranges = sorted({(k[2], k[3]) for k in ordered}, reverse=True)
pdf_path = out / 'sdhep_xsec_all_bins.pdf'
with PdfPages(pdf_path) as pdf:
    for t_range in sorted({(k[4], k[5]) for k in ordered}):
        keys = [k for k in ordered if (k[4], k[5]) == t_range]
        xi_ranges = global_xi_ranges
        cos_ranges = global_cos_ranges
        ncols = len(global_xi_ranges)
        nrows = len(global_cos_ranges)
        fig, axes = plt.subplots(nrows, ncols, figsize=(13, 3.2 * nrows), squeeze=False)
        key_lookup = {(k[0], k[1], k[2], k[3]): k for k in keys}
        for row, cos_range in enumerate(cos_ranges):
            for col, xi_range in enumerate(xi_ranges):
                ax = axes[row, col]
                key = key_lookup.get((*xi_range, *cos_range))
                if key is None:
                    ax.axhline(0, color='0.45', lw=0.7)
                    ax.set_xlim(0, 360)
                    ax.set_xticks([0, 90, 180, 270, 360])
                    ax.grid(alpha=0.18, lw=0.5)
                    ax.set_title(
                        rf'$\xi$ [{xi_range[0]:g},{xi_range[1]:g}), '
                        rf'$\cos\theta$ [{cos_range[0]:g},{cos_range[1]:g})' + '\n' +
                        rf'$|t|$ [{t_range[0]:g},{t_range[1]:g}), No data',
                        fontsize=8)
                    ax.text(0.5, 0.5, 'No data', transform=ax.transAxes,
                            ha='center', va='center', color='0.5', fontsize=9)
                else:
                    draw(ax, key, groups[key])
                if row == nrows - 1:
                    ax.set_xlabel(rf'$\xi\in[{xi_range[0]:g},{xi_range[1]:g})$')
                if col == 0:
                    ax.set_ylabel(rf'$\cos\theta\in[{cos_range[0]:g},{cos_range[1]:g})$' + '\n' + 'xsec')
        fig.supxlabel(r'$\phi_{\rm SDHEP}$ [deg]')
        fig.suptitle(rf'SDHEP BSA — $|t|\in[{t_range[0]:g},{t_range[1]:g})$ GeV$^2$', fontsize=14)
        fig.tight_layout(rect=(0.025, 0.025, 1, 0.97))
        pdf.savefig(fig)
        plt.close(fig)

# A single high-resolution overview image.
ncols = 5
nrows = int(np.ceil(len(ordered) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(16, 2.8 * nrows), squeeze=False)
for ax, key in zip(axes.flat, ordered):
    draw(ax, key, groups[key])
for ax in axes.flat[len(ordered):]:
    ax.axis('off')
fig.supxlabel(r'$\phi_{\rm SDHEP}$ [deg]')
fig.supylabel('xsec')
fig.suptitle('SDHEP beam-spin asymmetry — RGK 6.5 GeV', fontsize=16)
fig.tight_layout(rect=(0.02, 0.02, 1, 0.985))
fig.savefig(out / 'sdhep_bsa_overview.png', dpi=180)
plt.close(fig)

print(pdf_path)
