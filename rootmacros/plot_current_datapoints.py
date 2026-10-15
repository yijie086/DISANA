#!/usr/bin/env python3
"""Plot and compare DIS xsec and BSA data, one PDF for each kinematic bin."""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# ================= 用户配置：直接修改这里，无需命令行参数 =================
SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIRS = [
    SCRIPT_DIR,
    SCRIPT_DIR / "XSEC6535_0",
]
OUTPUT_DIR = SCRIPT_DIR / "datapoint_comparison_pdfs"  # 可改为 Path("/你的/输出文件夹")
LABELS = ["Fall18", "Spring24"]  # 可改为 ["version1", "version2"]，顺序和数量与 INPUT_DIRS 一致
CORRECTION = "all"  # 可选："all", "acc", "pi0", "eff", "rad", "p1"
# ======================================================================


BinKey = tuple[float, float, float, float, float, float]
Point = tuple[float, float, float, float, float, float]


def read_datapoints(path: Path) -> dict[BinKey, list[Point]]:
    groups: dict[BinKey, list[Point]] = defaultdict(list)
    with path.open(encoding="utf-8") as stream:
        for line_number, raw_line in enumerate(stream, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                raise ValueError(
                    f"{path}:{line_number}: expected at least 12 tab-separated fields"
                )
            try:
                values = tuple(float(field) for field in fields[:12])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid numeric field") from error
            groups[values[6:12]].append(values[:6])
    if not groups:
        raise ValueError(f"No datapoints found in {path}")
    return groups


def finite_points(rows: list[Point], require_positive: bool) -> np.ndarray:
    points = [
        row
        for row in rows
        if all(math.isfinite(value) for value in row)
        and (not require_positive or row[4] > 0.0)
    ]
    return np.asarray(sorted(points, key=lambda row: row[3]), dtype=float)


def output_name(kind: str, points: np.ndarray) -> str:
    x_b, q2, t_abs = np.mean(points[:, :3], axis=0)
    observable = "XS" if kind == "xsec" else "BSA"
    return f"{kind}_{observable}_xB{x_b:.3f}_Q2{q2:.3f}_t{t_abs:.3f}.pdf"


def dataset_labels(paths: list[Path], labels: list[str] | None = None) -> list[str]:
    if labels is not None:
        if len(labels) != len(paths):
            raise ValueError("--labels must contain one label per input")
        return labels
    folders = [path if path.is_dir() else path.parent for path in paths]
    names = [folder.resolve().name for folder in folders]
    return [
        name if names.count(name) == 1 else str(path.resolve())
        for name, path in zip(names, paths)
    ]


def resolve_inputs(
    paths: list[Path], labels: list[str], filename: str, allow_missing: bool
) -> tuple[list[Path], list[str]]:
    sources = []
    available_labels = []
    for path, label in zip(paths, labels):
        source = path / filename if path.is_dir() else path
        if not source.is_file():
            if allow_missing:
                print(f"Skipping missing input: {source}")
                continue
            raise ValueError(f"Input file not found: {source}")
        sources.append(source)
        available_labels.append(label)
    return sources, available_labels


def draw_series(ax, series: list[tuple[str, np.ndarray]], compare: bool, default_label: str) -> None:
    markers = ("s", "o", "^", "D", "v", "P", "X", "<", ">", "*")
    colors = plt.get_cmap("tab20").colors
    colors = colors[::2] + colors[1::2]
    for index, (label, points) in enumerate(series):
        if not points.size:
            continue
        color = colors[index % len(colors)] if compare else "black"
        ax.errorbar(
            points[:, 3], points[:, 4], yerr=np.abs(points[:, 5]),
            fmt=markers[index % len(markers)] if compare else "s",
            color=color, ecolor=color, zorder=3,
            label=f"{label}, stat. unc." if compare else f"{label}: {default_label}",
        )


def plot_sources(
    sources: list[Path], labels: list[str], output_dir: Path,
    require_positive: bool, name_for_points, draw,
) -> tuple[int, int]:
    if not sources or len(sources) != len(labels):
        raise ValueError("Provide one label for each input file")
    groups = [read_datapoints(source) for source in sources]
    keys = set().union(*(group.keys() for group in groups))
    output_dir.mkdir(parents=True, exist_ok=True)
    written = skipped = 0
    for key in sorted(keys, key=lambda item: (item[4], item[0], item[2], item)):
        series = []
        for label, group in zip(labels, groups):
            points = finite_points(group.get(key, []), require_positive)
            series.append((label, points))
        if not any(points.size for _, points in series):
            skipped += 1
            continue
        points = np.concatenate([points for _, points in series if points.size])
        destination = output_dir / name_for_points(points)
        draw(points, destination, series, len(sources) > 1)
        written += 1
    return written, skipped


def draw_bin(
    kind: str, points: np.ndarray, destination: Path,
    series: list[tuple[str, np.ndarray]] | None = None, compare: bool = False,
) -> None:
    x_b, q2, t_abs = np.mean(points[:, :3], axis=0)
    fig, ax = plt.subplots()
    try:
        draw_series(
            ax, series if series is not None else [("Data", points)], compare,
            "Data (acceptance correction only; not normalized), stat. unc.",
        )
        ax.set_xlabel("φ (deg)")
        if kind == "xsec":
            ax.set_ylabel("Cross section (nb/GeV^4)")
            ax.set_yscale("log")
            observable = "XS"
        else:
            ax.set_ylabel(r"$A_{LU}$")
            ax.axhline(0.0, color="black", linestyle="--", linewidth=1)
            ax.set_ylim(-0.6, 0.6)
            observable = "BSA"
        ax.set_title(
            f"{observable}  xB={x_b:.3f}, Q2={q2:.3f}, t={t_abs:.3f}\n"
            "Acceptance correction only; not normalized"
        )
        ax.legend()
        fig.tight_layout()
        fig.savefig(destination)
    finally:
        plt.close(fig)


def plot_file(
    source: Path | list[Path], output_dir: Path, kind: str,
    labels: list[str] | None = None,
) -> tuple[int, int]:
    sources = [source] if isinstance(source, Path) else list(source)
    return plot_sources(
        sources, dataset_labels(sources, labels), output_dir, kind == "xsec",
        lambda points: output_name(kind, points),
        lambda points, destination, series, compare: draw_bin(
            kind, points, destination, series, compare
        ),
    )


def main() -> None:
    script_dir = SCRIPT_DIR
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folders", nargs="*", type=Path, help="input folders to compare")
    parser.add_argument("--input", "--input-dirs", dest="input_dirs", nargs="+", type=Path,
                        help="input folders to compare (alternative to positional folders)")
    parser.add_argument("--labels", nargs="+", default=LABELS, help="legend labels in input order")
    parser.add_argument("--xsec-input", nargs="+", type=Path,
                        help="cross-section TXT files or folders")
    parser.add_argument("--bsa-input", nargs="+", type=Path,
                        help="BSA TXT files or folders")
    parser.add_argument("--output", type=Path, default=OUTPUT_DIR, help="output root, with xsec/ and bsa/ subfolders")
    parser.add_argument("--xsec-output", type=Path, help="cross-section output folder")
    parser.add_argument("--bsa-output", type=Path, help="BSA output folder")
    args = parser.parse_args()
    if args.output is not None:
        args.output = Path(args.output).expanduser()
    if args.folders and args.input_dirs:
        parser.error("use positional folders or --input-dirs, not both")
    folders = args.input_dirs or args.folders or [Path(folder).expanduser() for folder in INPUT_DIRS]
    for folder in folders:
        if not folder.is_dir():
            parser.error(f"input folder not found: {folder}")
    total_written = 0
    try:
        for kind, filename in (("xsec", "datapoint_DIS_Cross_Section.txt"),
                               ("bsa", "datapoint_DIS_BSA.txt")):
            explicit = getattr(args, f"{kind}_input")
            paths = explicit or folders or [script_dir / filename]
            labels = dataset_labels(paths, args.labels)
            sources, labels = resolve_inputs(paths, labels, filename,
                                             allow_missing=bool(folders) and explicit is None)
            if not sources:
                print(f"{kind}: no available input files, skipped")
                continue
            output_dir = getattr(args, f"{kind}_output") or (
                args.output / kind if args.output else script_dir / f"datapoint_{kind}_pdfs"
            )
            written, skipped = plot_file(sources, output_dir, kind, labels)
            total_written += written
            print(f"{kind}: {written} PDFs in {output_dir} ({skipped} empty bins skipped)")
        if not total_written:
            parser.error("no plottable datapoints found in the supplied inputs")
    except (ValueError, OSError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
