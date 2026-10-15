#!/usr/bin/env python3
"""Auto-detect and compare DIS corrections, one PDF for each kinematic bin."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# ================= 用户配置：直接修改这里，无需命令行参数 =================
SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIRS = [
    SCRIPT_DIR,
    SCRIPT_DIR / "XSEC6535_0",
]
OUTPUT_DIR = SCRIPT_DIR / "correction_comparison_pdfs"  # 可改为 Path("/你的/输出文件夹")
LABELS = ["Spring24", "Fall18"]  # 可改为 ["version1", "version2"]，顺序和数量与 INPUT_DIRS 一致
CORRECTION = "all"  # 可选："all", "acc", "pi0", "eff", "rad", "p1"
# ======================================================================

try:
    from plot_current_datapoints import dataset_labels, draw_series, plot_sources, resolve_inputs
except ModuleNotFoundError as error:
    if error.name != "plot_current_datapoints":
        raise
    # Some remote mounts permit file reads but return an empty directory listing.
    import importlib.util

    helper_path = Path(__file__).resolve().with_name("plot_current_datapoints.py")
    spec = importlib.util.spec_from_file_location("plot_current_datapoints", helper_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load plotting helpers from {helper_path}") from error
    helpers = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helpers)
    dataset_labels = helpers.dataset_labels
    draw_series = helpers.draw_series
    plot_sources = helpers.plot_sources
    resolve_inputs = helpers.resolve_inputs


CORRECTIONS = {
    "acc": {
        "filename": "datapoint_DIS_accCorr.txt",
        "observable": "accCorr",
        "title": "Acceptance correction",
        "ylabel": r"$A_{\mathrm{acc}}$",
        "ylim": (0.01, 1.0),
        "logy": True,
    },
    "pi0": {
        "filename": "datapoint_DIS_pi0Corr.txt",
        "observable": "pi0Corr",
        "title": r"$\pi^0$ correction",
        "ylabel": r"$\eta^{\pi^0}$",
        "ylim": (0.0, 1.0),
        "logy": False,
    },
    "eff": {
        "filename": "datapoint_DIS_effCorr.txt",
        "observable": "effCorr",
        "title": "Efficiency correction",
        "ylabel": r"$A_{\mathrm{eff}}$",
        "ylim": (0.1, 1.1),
        "logy": False,
    },
    "rad": {
        "filename": "datapoint_DIS_radCorr.txt",
        "observable": "radCorr",
        "title": "Radiative correction",
        "ylabel": r"$C_{\mathrm{rad}}$",
        "ylim": (0.5, 1.5),
        "logy": False,
    },
    "p1": {
        "filename": "datapoint_DIS_P1Cut.txt",
        "observable": "P1Cut",
        "title": "P1-cut correction",
        "ylabel": r"$C_{P1}$",
        "ylim": (0.0, 1.2),
        "logy": False,
    },
}


def output_name(observable: str, points: np.ndarray) -> str:
    x_b, q2, t_abs = np.mean(points[:, :3], axis=0)
    return (
        f"correction_{observable}_xB{x_b:.3f}_Q2{q2:.3f}_t{t_abs:.3f}.pdf"
    )


def draw_bin(
    config: dict, points: np.ndarray, destination: Path,
    series: list[tuple[str, np.ndarray]] | None = None, compare: bool = False,
) -> None:
    x_b, q2, t_abs = np.mean(points[:, :3], axis=0)
    fig, ax = plt.subplots()
    try:
        draw_series(
            ax, series if series is not None else [("Correction", points)], compare,
            "Correction, stat. unc.",
        )
        ax.set_xlabel("φ (deg)")
        ax.set_ylabel(config["ylabel"])
        if config["logy"]:
            ax.set_yscale("log")
        ax.set_ylim(*config["ylim"])
        ax.set_title(f"{config['title']}  xB={x_b:.3f}, Q2={q2:.3f}, t={t_abs:.3f}")
        ax.legend()
        fig.tight_layout()
        fig.savefig(destination)
    finally:
        plt.close(fig)


def plot_correction(
    correction: str, source: Path | list[Path], output_dir: Path,
    labels: list[str] | None = None,
) -> tuple[int, int]:
    config = CORRECTIONS[correction]
    sources = [source] if isinstance(source, Path) else list(source)
    return plot_sources(
        sources, dataset_labels(sources, labels), output_dir, config["logy"],
        lambda points: output_name(config["observable"], points),
        lambda points, destination, series, compare: draw_bin(
            config, points, destination, series, compare
        ),
    )


def main() -> None:
    script_dir = SCRIPT_DIR
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folders", nargs="*", type=Path, help="input folders to compare")
    parser.add_argument("--correction", choices=("all", *CORRECTIONS), default=CORRECTION,
                        help="plot available corrections, or select one type (default: all)")
    parser.add_argument("--input", "--input-dirs", dest="inputs", nargs="+", type=Path,
                        help="input folders, or TXT files for one explicit correction type")
    parser.add_argument("--labels", nargs="+", default=LABELS, help="legend labels in input order")
    parser.add_argument("--output", type=Path, default=OUTPUT_DIR, help="override the output directory")
    args = parser.parse_args()
    if args.output is not None:
        args.output = Path(args.output).expanduser()
    if args.folders and args.inputs:
        parser.error("use positional folders or --input, not both")
    paths = args.inputs or args.folders or [Path(folder).expanduser() for folder in INPUT_DIRS] or [script_dir]
    for path in paths:
        if not path.exists():
            parser.error(f"input path not found: {path}")
        if (args.correction == "all" or args.folders) and not path.is_dir():
            parser.error("TXT files require --input and one explicit --correction type")
    selected = list(CORRECTIONS) if args.correction == "all" else [args.correction]
    plotted = []
    try:
        input_labels = dataset_labels(paths, args.labels)
        for correction in selected:
            config = CORRECTIONS[correction]
            sources, labels = resolve_inputs(
                paths, input_labels, config["filename"], allow_missing=args.correction == "all"
            )
            if not sources:
                continue
            if args.output is None:
                output_dir = script_dir / f"datapoint_{correction}_correction_pdfs"
            elif args.correction == "all":
                output_dir = args.output / correction
            else:
                output_dir = args.output
            written, skipped = plot_correction(correction, sources, output_dir, labels)
            plotted.append(correction)
            print(f"{correction}: {written} PDFs in {output_dir} ({skipped} empty bins skipped)")
        if not plotted:
            expected = ", ".join(config["filename"] for config in CORRECTIONS.values())
            parser.error(f"no correction datapoint TXT found; expected one of: {expected}")
    except (ValueError, OSError) as error:
        parser.error(str(error))
    if args.correction == "all":
        print(f"auto-detected correction types: {', '.join(plotted)}")
        missing = [name for name in CORRECTIONS if name not in plotted]
        if missing:
            print(f"skipped missing correction types: {', '.join(missing)}")


if __name__ == "__main__":
    main()
