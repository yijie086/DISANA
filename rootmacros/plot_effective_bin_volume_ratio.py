#!/usr/bin/env python3
"""Plot Veff / (Delta xB Delta Q2 Delta |t|) in xB-Q2 bins.

By default, all effective_bin_volume_Eb*.csv files in the current directory
are collected into one multipage PDF, with one page for every |t| bin.
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path
from typing import Iterable

try:
    import matplotlib

    matplotlib.use("Agg")

    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import Normalize
    from matplotlib.patches import Rectangle
    import numpy as np
    import pandas as pd
except ModuleNotFoundError as exc:
    raise SystemExit(
        f"Missing Python package '{exc.name}'. Install dependencies with: "
        "python -m pip install matplotlib pandas numpy"
    ) from exc


REQUIRED_COLUMNS = {
    "xB_min",
    "xB_max",
    "Q2_min",
    "Q2_max",
    "t_min",
    "t_max",
    "Veff_over_Vnominal",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Draw xB versus Q2 maps of effective-bin-volume ratio; "
            "each |t| bin is written as a separate PDF page."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        default=["effective_bin_volume_Eb*.csv"],
        help="input CSV file(s) or glob(s) (default: effective_bin_volume_Eb*.csv)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="effective_bin_volume_ratio_xB_Q2_by_t.pdf",
        help="output multipage PDF",
    )
    parser.add_argument(
        "--no-labels",
        action="store_true",
        help="do not print numerical ratios inside the bins",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=160,
        help="PDF figure DPI (default: 160)",
    )
    return parser.parse_args()


def expand_inputs(patterns: Iterable[str]) -> list[Path]:
    paths: list[Path] = []
    for pattern in patterns:
        matches = [Path(p) for p in glob.glob(pattern)]
        if matches:
            paths.extend(matches)
        elif Path(pattern).is_file():
            paths.append(Path(pattern))

    unique = sorted({p.resolve() for p in paths})
    if not unique:
        raise SystemExit(
            "No input CSV files found. Run from the directory containing "
            "effective_bin_volume_Eb*.csv, or pass CSV paths explicitly."
        )
    return unique


def load_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = sorted(REQUIRED_COLUMNS.difference(df.columns))
    if missing:
        raise SystemExit(f"{path}: missing required column(s): {', '.join(missing)}")

    numeric = list(REQUIRED_COLUMNS)
    if "beam_energy" in df.columns:
        numeric.append("beam_energy")
    for col in numeric:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    bad_edges = df[list(REQUIRED_COLUMNS - {"Veff_over_Vnominal"})].isna().any(axis=1)
    if bad_edges.any():
        print(f"[warning] {path.name}: skipped {int(bad_edges.sum())} row(s) with invalid bin edges")
        df = df.loc[~bad_edges].copy()

    df["_source"] = path.name
    return df


def cut_summary(df: pd.DataFrame) -> str:
    fields = [
        ("gen_xB_min", "xBmin"),
        ("gen_xB_max", "xBmax"),
        ("gen_Q2_min", "Q2min"),
        ("gen_Q2_max", "Q2max"),
        ("gen_y_min", "ymin"),
        ("gen_y_max", "ymax"),
        ("gen_W2_min", "W2min"),
        ("gen_t_min", "|t|min"),
        ("gen_t_max", "|t|max"),
    ]
    pieces = []
    for col, label in fields:
        if col in df.columns:
            values = pd.to_numeric(df[col], errors="coerce").dropna().unique()
            if len(values) == 1:
                pieces.append(f"{label}={values[0]:g}")
    return "generator: " + ", ".join(pieces) if pieces else ""


def draw_page(
    pdf: PdfPages,
    group: pd.DataFrame,
    source: str,
    t_min: float,
    t_max: float,
    beam_energy: float | None,
    labels: bool,
    dpi: int,
) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 6.5), constrained_layout=True)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("lightgray")
    norm = Normalize(vmin=0.0, vmax=1.0, clip=True)

    for row in group.itertuples(index=False):
        ratio = float(row.Veff_over_Vnominal)
        valid = np.isfinite(ratio)
        color = cmap(norm(ratio)) if valid else "lightgray"
        width = float(row.xB_max - row.xB_min)
        height = float(row.Q2_max - row.Q2_min)
        ax.add_patch(
            Rectangle(
                (row.xB_min, row.Q2_min),
                width,
                height,
                facecolor=color,
                edgecolor="black",
                linewidth=0.7,
            )
        )
        if labels:
            text = f"{ratio:.3f}" if valid else "nan"
            text_color = "white" if valid and ratio < 0.55 else "black"
            ax.text(
                row.xB_min + width / 2.0,
                row.Q2_min + height / 2.0,
                text,
                ha="center",
                va="center",
                fontsize=7,
                color=text_color,
            )

    x_min, x_max = group["xB_min"].min(), group["xB_max"].max()
    q_min, q_max = group["Q2_min"].min(), group["Q2_max"].max()
    x_pad = max(0.01 * (x_max - x_min), 0.002)
    q_pad = max(0.01 * (q_max - q_min), 0.02)
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(q_min - q_pad, q_max + q_pad)
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"$Q^2\;[\mathrm{GeV}^2]$")

    energy_text = "" if beam_energy is None else rf", $E_b={beam_energy:g}$ GeV"
    ax.set_title(
        rf"$V_{{\rm eff}}/(\Delta x_B\,\Delta Q^2\,\Delta |t|)$"
        + "\n"
        + rf"${t_min:g} \leq |t| < {t_max:g}$ GeV$^2${energy_text}"
    )
    subtitle = cut_summary(group)
    if subtitle:
        ax.text(0.5, 1.005, subtitle, transform=ax.transAxes, ha="center", va="bottom", fontsize=7)
    ax.text(0.01, 0.01, source, transform=ax.transAxes, fontsize=7, color="0.35")
    ax.grid(False)

    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar.set_array([])
    cbar = fig.colorbar(scalar, ax=ax, pad=0.025)
    cbar.set_label(r"$V_{\rm eff}/V_{\rm nominal}$")
    pdf.savefig(fig, dpi=dpi)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    input_paths = expand_inputs(args.inputs)
    output = Path(args.output).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    page_count = 0
    with PdfPages(output) as pdf:
        for path in input_paths:
            df = load_csv(path)
            if df.empty:
                print(f"[warning] {path.name}: no valid rows")
                continue

            group_columns = ["t_min", "t_max"]
            has_energy = "beam_energy" in df.columns
            if has_energy:
                group_columns.insert(0, "beam_energy")

            for keys, group in df.sort_values(group_columns).groupby(
                group_columns, sort=True, dropna=False
            ):
                keys = keys if isinstance(keys, tuple) else (keys,)
                if has_energy:
                    beam_energy, t_min, t_max = keys
                    beam_energy = None if pd.isna(beam_energy) else float(beam_energy)
                else:
                    t_min, t_max = keys
                    beam_energy = None
                draw_page(
                    pdf,
                    group,
                    path.name,
                    float(t_min),
                    float(t_max),
                    beam_energy,
                    not args.no_labels,
                    args.dpi,
                )
                page_count += 1

    if page_count == 0:
        output.unlink(missing_ok=True)
        raise SystemExit("No valid t-bin pages were produced.")
    print(f"[done] wrote {page_count} page(s) to {output}")


if __name__ == "__main__":
    main()
