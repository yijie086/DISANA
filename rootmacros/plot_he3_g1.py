#!/usr/bin/env python3
"""Plot free and bound proton/neutron x*g1 data in polarized He-3.

Expected input columns:
    x, g1p_free, g1p_He3, g1p_He3/free,
       g1n_free, g1n_He3, g1n_He3/free

Examples:
    python plot_he3_g1.py g1p_g1n_in_He3.dat --show
    python plot_he3_g1.py g1p_g1n_in_He3.dat -o he3_g1.pdf
    python plot_he3_g1.py g1p_g1n_in_He3.dat --xmin 0.2 --xmax 0.8 \
        --ratio-scale linear --ratio-ymin 0.85 --ratio-ymax 1.15 --show
"""

from __future__ import annotations

import argparse
import ctypes
import os
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    import numpy as np
except ModuleNotFoundError as exc:
    raise SystemExit(
        "Missing plotting dependency. Install it with:\n"
        "  python3 -m pip install numpy matplotlib"
    ) from exc


# -----------------------------------------------------------------------------
# CONFIG: edit this block for the usual plotting changes
# -----------------------------------------------------------------------------
CONFIG = {
    "input": Path("g1p_g1n_in_He3.dat"),
    "output": Path("he3_xg1_all_data.png"),

    # Shared x axis for all three panels.
    # scale choices: "linear" or "log". Log requires x > 0.
    "x_axis": {
        "min": None,       # e.g. 0.02
        "max": None,       # e.g. 1.00
        "scale": "linear",
    },

    # Each panel has an independent y range and scale.
    # scale choices: "linear", "log", or "symlog".
    # Use symlog when data contain both positive and negative values.
    "proton_y_axis": {
        "min": None,
        "max": None,
        "scale": "linear",
        "linthresh": 0.01,
    },
    "neutron_y_axis": {
        "min": None,
        "max": None,
        "scale": "linear",
        "linthresh": 0.005,
    },
    "ratio_y_axis": {
        "min": 0.8,       # e.g. 0.85 for a zoomed linear plot
        "max": 1.2,       # e.g. 1.15
        "scale": "linear",
        "linthresh": 0.2,
    },
    "difference_y_axis": {
        "min": None,
        "max": None,
        "scale": "linear",
        "linthresh": 0.005,
    },

    # Wide canvas leaves a dedicated right margin for legends. Keeping all
    # legends outside the axes prevents them from hiding physics curves.
    "figsize": (14.5, 8.5),
    "font_size": 12,
    "line_width": 1.8,
    "marker_size": 2.5,
    "markers": False,
    "mark_zero_crossings": True,
    # The input table contains x but no Q2 column, so the JAM25 comparison
    # must be evaluated at one explicitly chosen fixed scale.
    "jam25": True,
    "jam25_q2": 5.0,     # GeV^2
    "dpi": 220,
    "title": r"Free and bound nucleon $xg_1$ in $^3$He",
    "q2_label": "",       # e.g. r"$Q^2=5\,\mathrm{GeV}^2$"
}

COLORS = {
    "p_free": "#1f77b4",
    "p_he3": "#ff7f0e",
    "p_jam25": "#111111",
    "n_free": "#2ca02c",
    "n_he3": "#d62728",
    "n_jam25": "#111111",
    "ratio_p": "#9467bd",
    "ratio_n": "#8c564b",
    "difference_p": "#ff7f0e",
    "difference_n": "#d62728",
}

FIGSIZE = CONFIG["figsize"]
LINEWIDTH = CONFIG["line_width"]
MARKER_SIZE = CONFIG["marker_size"]
FONT_SIZE = CONFIG["font_size"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot x*g1p and x*g1n from a He-3 model table."
    )
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=CONFIG["input"],
        help=f"Input .dat file (default: {CONFIG['input']})",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=CONFIG["output"],
        help=f"Output image/PDF path (default: {CONFIG['output']})",
    )
    parser.add_argument(
        "--xmin", type=float, default=CONFIG["x_axis"]["min"],
        help="Minimum x to display",
    )
    parser.add_argument(
        "--xmax", type=float, default=CONFIG["x_axis"]["max"],
        help="Maximum x to display",
    )
    parser.add_argument(
        "--xscale", choices=("linear", "log"),
        default=CONFIG["x_axis"]["scale"],
        help="Shared x-axis scale",
    )
    parser.add_argument(
        "--proton-yscale", choices=("linear", "log", "symlog"),
        default=CONFIG["proton_y_axis"]["scale"],
    )
    parser.add_argument("--proton-ymin", type=float, default=CONFIG["proton_y_axis"]["min"])
    parser.add_argument("--proton-ymax", type=float, default=CONFIG["proton_y_axis"]["max"])
    parser.add_argument(
        "--neutron-yscale", choices=("linear", "log", "symlog"),
        default=CONFIG["neutron_y_axis"]["scale"],
    )
    parser.add_argument("--neutron-ymin", type=float, default=CONFIG["neutron_y_axis"]["min"])
    parser.add_argument("--neutron-ymax", type=float, default=CONFIG["neutron_y_axis"]["max"])
    parser.add_argument(
        "--ratio-scale",
        choices=("linear", "log", "symlog"),
        default=CONFIG["ratio_y_axis"]["scale"],
        help="Vertical scale for the ratio panel",
    )
    parser.add_argument(
        "--linthresh",
        type=float,
        default=CONFIG["ratio_y_axis"]["linthresh"],
        help="Linear threshold for the symlog ratio scale",
    )
    parser.add_argument("--ratio-ymin", type=float, default=CONFIG["ratio_y_axis"]["min"])
    parser.add_argument("--ratio-ymax", type=float, default=CONFIG["ratio_y_axis"]["max"])
    parser.add_argument(
        "--difference-yscale", choices=("linear", "log", "symlog"),
        default=CONFIG["difference_y_axis"]["scale"],
    )
    parser.add_argument(
        "--difference-ymin", type=float,
        default=CONFIG["difference_y_axis"]["min"],
    )
    parser.add_argument(
        "--difference-ymax", type=float,
        default=CONFIG["difference_y_axis"]["max"],
    )
    parser.add_argument(
        "--markers",
        action=argparse.BooleanOptionalAction,
        default=CONFIG["markers"],
        help="Draw a marker at every input row in addition to lines",
    )
    parser.add_argument(
        "--zero-markers",
        action=argparse.BooleanOptionalAction,
        default=CONFIG["mark_zero_crossings"],
        help="Mark neutron zero-crossing positions",
    )
    parser.add_argument(
        "--jam25",
        action=argparse.BooleanOptionalAction,
        default=CONFIG["jam25"],
        help=("Overlay the generator's JAMpol25 g1p and g1n curves "
              "(default: enabled)"),
    )
    parser.add_argument(
        "--jam25-q2",
        type=float,
        default=CONFIG["jam25_q2"],
        help="Fixed Q2 in GeV^2 used for the JAM25 curves (default: 5)",
    )
    parser.add_argument(
        "--title",
        default=CONFIG["title"],
        help="Figure title",
    )
    parser.add_argument(
        "--q2-label",
        default=CONFIG["q2_label"],
        help=r"Optional text such as '$Q^2=5\,\mathrm{GeV}^2$'",
    )
    parser.add_argument("--dpi", type=int, default=CONFIG["dpi"], help="Raster output DPI")
    parser.add_argument("--show", action="store_true", help="Open an interactive window")
    return parser.parse_args()


def read_table(path: Path) -> dict[str, np.ndarray]:
    table = np.loadtxt(path, comments="#", ndmin=2)
    if table.shape[1] < 7:
        raise ValueError(
            f"Expected at least 7 columns, but {path} has {table.shape[1]}."
        )

    return {
        "x": table[:, 0],
        "g1p_free": table[:, 1],
        "g1p_he3": table[:, 2],
        "ratio_p": table[:, 3],
        "g1n_free": table[:, 4],
        "g1n_he3": table[:, 5],
        "ratio_n": table[:, 6],
    }


def load_jam25_nucleon_function():
    """Load the exact JAM25 interface used by the He3DIS generator.

    The shared generator library currently returns the JAM25 PID-908 value
    divided by x.  JAM25's own plotting utility identifies xfxQ2(908,x,Q2)
    itself as g1, so evaluate_jam25 multiplies the library result by x to undo
    that extra division locally.  The generator source/library is not changed.
    """
    candidates: list[Path] = []
    configured = os.environ.get("HE3DIS_SPIN_MODEL_LIB")
    if configured:
        candidates.append(Path(configured).expanduser())

    generator_root = (
        Path(__file__).resolve().parent
        / "../../../../Simulation/He3DIS/He3DIS"
    ).resolve()
    candidates.extend(
        [
            generator_root / "build-ifarm/libhe3dis_spin_model.so",
            generator_root / "build/libhe3dis_spin_model.dylib",
            generator_root / "build/libhe3dis_spin_model.so",
        ]
    )

    errors: list[str] = []
    for path in candidates:
        if not path.is_file():
            continue
        try:
            library = ctypes.CDLL(str(path))
            function = library.he3dis_nucleon_observables
            function.argtypes = [
                ctypes.c_int,
                ctypes.c_double,
                ctypes.c_double,
                ctypes.c_double,
            ] + [ctypes.POINTER(ctypes.c_double)] * 7
            function.restype = ctypes.c_int
            # Retain the library for as long as the function is in use.
            function._he3dis_library = library
            print(f"Loaded JAM25 generator model: {path}")
            return function
        except (OSError, AttributeError) as exc:
            errors.append(f"  {path}: {exc}")

    details = "\n".join(errors) if errors else "  no candidate library exists"
    raise RuntimeError(
        "Could not load the He3DIS JAM25 model library. Build it with the "
        "He3DIS setup/build command or set HE3DIS_SPIN_MODEL_LIB.\n" + details
    )


def evaluate_jam25(x: np.ndarray, q2: float) -> dict[str, np.ndarray]:
    """Evaluate generator JAM25 g1p and g1n at alpha=1 and fixed Q2."""
    if not np.isfinite(q2) or q2 <= 0.0:
        raise ValueError("--jam25-q2 must be positive and finite.")

    function = load_jam25_nucleon_function()

    def evaluate_nucleon(struck_neutron: bool) -> np.ndarray:
        result = np.full(x.shape, np.nan, dtype=float)
        failed_statuses: dict[int, int] = {}
        for index, x_value in enumerate(x):
            if not np.isfinite(x_value) or not 0.0 < x_value < 1.0:
                continue
            outputs = [ctypes.c_double(np.nan) for _ in range(7)]
            status = function(
                int(struck_neutron),
                float(x_value),
                1.0,  # alpha_active=1 makes x_active=x_b=x for a free nucleon
                float(q2),
                *(ctypes.byref(value) for value in outputs),
            )
            if status == 0:
                result[index] = outputs[2].value
            else:
                failed_statuses[status] = failed_statuses.get(status, 0) + 1
        if failed_statuses:
            particle = "neutron" if struck_neutron else "proton"
            print(f"Warning: JAM25 {particle} failed statuses: {failed_statuses}")
        return result

    return {
        "g1p_jam25": evaluate_nucleon(False),
        "g1n_jam25": evaluate_nucleon(True),
    }


def restrict_x(data: dict[str, np.ndarray], xmin: float | None, xmax: float | None):
    x = data["x"]
    mask = np.ones_like(x, dtype=bool)
    if xmin is not None:
        mask &= x >= xmin
    if xmax is not None:
        mask &= x <= xmax
    if not np.any(mask):
        raise ValueError("The selected x range contains no input rows.")
    return {name: values[mask] for name, values in data.items()}


def zero_crossings(x: np.ndarray, y: np.ndarray) -> list[float]:
    """Find zero crossings using linear interpolation between adjacent rows."""
    crossings: list[float] = []
    for x0, x1, y0, y1 in zip(x[:-1], x[1:], y[:-1], y[1:]):
        if y0 == 0:
            crossings.append(float(x0))
        elif np.isfinite(y0) and np.isfinite(y1) and y0 * y1 < 0:
            crossings.append(float(x0 - y0 * (x1 - x0) / (y1 - y0)))
    if y[-1] == 0:
        crossings.append(float(x[-1]))
    return crossings


def plot_line(ax, x, y, *, label, color, linestyle="-", markers=False):
    ax.plot(
        x,
        y,
        color=color,
        linestyle=linestyle,
        linewidth=LINEWIDTH,
        marker="o" if markers else None,
        markersize=MARKER_SIZE if markers else None,
        markeredgewidth=0,
        label=label,
    )


def format_axis(ax, ylabel: str):
    ax.set_ylabel(ylabel)
    # Publication-style axes: no background grid, inward ticks on all sides,
    # and minor ticks for reading values without competing with the data.
    ax.grid(False)
    ax.minorticks_on()
    ax.tick_params(
        which="major", direction="in", top=True, right=True,
        length=6, width=1.0,
    )
    ax.tick_params(
        which="minor", direction="in", top=True, right=True,
        length=3, width=0.8,
    )
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
    ax.margins(x=0)


def add_external_legends(fig, axes):
    """Place one compact legend group per panel in the right-hand margin."""
    specs = (
        (axes[0], r"(a) Proton $xg_1$", 0.83),
        (axes[1], r"(b) Neutron $xg_1$", 0.61),
        (axes[2], r"(c) Nuclear ratios", 0.38),
        (axes[3], r"(d) Nuclear shifts", 0.17),
    )
    for ax, title, y_anchor in specs:
        handles, labels = ax.get_legend_handles_labels()
        legend = fig.legend(
            handles,
            labels,
            title=title,
            loc="center left",
            bbox_to_anchor=(0.795, y_anchor),
            frameon=False,
            handlelength=2.8,
            handletextpad=0.7,
            labelspacing=0.55,
            borderaxespad=0.0,
        )
        legend._legend_box.align = "left"


def set_y_axis(ax, *, scale: str, ymin, ymax, linthresh: float):
    """Apply one panel's y scale and optional manual limits."""
    if scale == "symlog":
        ax.set_yscale("symlog", linthresh=linthresh)
    else:
        ax.set_yscale(scale)
    if ymin is not None or ymax is not None:
        ax.set_ylim(bottom=ymin, top=ymax)


def make_figure(data: dict[str, np.ndarray], args: argparse.Namespace):
    plt.rcParams.update(
        {
            "font.size": FONT_SIZE,
            "font.family": "serif",
            "mathtext.fontset": "stix",
            "axes.labelsize": FONT_SIZE + 1,
            "axes.titlesize": FONT_SIZE,
            "legend.fontsize": FONT_SIZE - 1,
            "legend.title_fontsize": FONT_SIZE,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(
        2, 2, figsize=FIGSIZE, sharex=True, constrained_layout=False
    )
    ax_p, ax_n, ax_r, ax_d = axes.ravel()
    x = data["x"]
    # The input table and generator API contain g1.  Multiply every absolute
    # free, bound, and JAM25 curve by the same Bjorken x before plotting.
    xg1p_free = x * data["g1p_free"]
    xg1p_he3 = x * data["g1p_he3"]
    xg1n_free = x * data["g1n_free"]
    xg1n_he3 = x * data["g1n_he3"]

    plot_line(
        ax_p,
        x,
        xg1p_free,
        label=r"$xg_1^p$ free",
        color=COLORS["p_free"],
        markers=args.markers,
    )
    if "g1p_jam25" in data:
        plot_line(
            ax_p,
            x,
            x * data["g1p_jam25"],
            label=rf"JAM25 $xg_1^p$ ($Q^2={args.jam25_q2:g}$ GeV$^2$)",
            color=COLORS["p_jam25"],
            linestyle=":",
            markers=False,
        )
    plot_line(
        ax_p,
        x,
        xg1p_he3,
        label=r"$xg_1^p$ in $^3$He",
        color=COLORS["p_he3"],
        linestyle="--",
        markers=args.markers,
    )
    format_axis(ax_p, r"$xg_1^p(x)$")
    set_y_axis(
        ax_p,
        scale=args.proton_yscale,
        ymin=args.proton_ymin,
        ymax=args.proton_ymax,
        linthresh=CONFIG["proton_y_axis"]["linthresh"],
    )

    plot_line(
        ax_n,
        x,
        xg1n_free,
        label=r"$xg_1^n$ free",
        color=COLORS["n_free"],
        markers=args.markers,
    )
    if "g1n_jam25" in data:
        plot_line(
            ax_n,
            x,
            x * data["g1n_jam25"],
            label=rf"JAM25 $xg_1^n$ ($Q^2={args.jam25_q2:g}$ GeV$^2$)",
            color=COLORS["n_jam25"],
            linestyle=":",
            markers=False,
        )
    plot_line(
        ax_n,
        x,
        xg1n_he3,
        label=r"$xg_1^n$ in $^3$He",
        color=COLORS["n_he3"],
        linestyle="--",
        markers=args.markers,
    )
    ax_n.axhline(0.0, color="0.4", linewidth=1.0)

    if args.zero_markers:
        crossing_specs = (
            ("g1n_free", COLORS["n_free"], "free"),
            ("g1n_he3", COLORS["n_he3"], r"$^3$He"),
        )
        for key, color, label in crossing_specs:
            for crossing in zero_crossings(x, data[key]):
                # Do not label the trivial x=1 endpoint.
                if crossing >= x.max() - 1e-12:
                    continue
                ax_n.axvline(crossing, color=color, linewidth=1.0, linestyle=":")
                ax_n.annotate(
                    rf"{label}: $x_0={crossing:.4f}$",
                    xy=(crossing, 0),
                    xytext=(5, 8 if key == "g1n_free" else -18),
                    textcoords="offset points",
                    color=color,
                    fontsize=FONT_SIZE - 1,
                )
    format_axis(ax_n, r"$xg_1^n(x)$")
    set_y_axis(
        ax_n,
        scale=args.neutron_yscale,
        ymin=args.neutron_ymin,
        ymax=args.neutron_ymax,
        linthresh=CONFIG["neutron_y_axis"]["linthresh"],
    )

    # Matplotlib cannot draw infinity. Keeping NaN at those locations makes the
    # line stop honestly at the last finite point; the endpoint is annotated.
    ratio_p = np.where(np.isfinite(data["ratio_p"]), data["ratio_p"], np.nan)
    ratio_n = np.where(np.isfinite(data["ratio_n"]), data["ratio_n"], np.nan)
    plot_line(
        ax_r,
        x,
        ratio_p,
        label=r"$g_{1p}^{^3\mathrm{He}}/g_{1p}^{\mathrm{free}}$",
        color=COLORS["ratio_p"],
        markers=args.markers,
    )
    plot_line(
        ax_r,
        x,
        ratio_n,
        label=r"$g_{1n}^{^3\mathrm{He}}/g_{1n}^{\mathrm{free}}$",
        color=COLORS["ratio_n"],
        linestyle="--",
        markers=args.markers,
    )
    ax_r.axhline(1.0, color="0.4", linewidth=1.0, linestyle=":")

    set_y_axis(
        ax_r,
        scale=args.ratio_scale,
        ymin=args.ratio_ymin,
        ymax=args.ratio_ymax,
        linthresh=args.linthresh,
    )

    if np.any(~np.isfinite(data["ratio_p"])) or np.any(~np.isfinite(data["ratio_n"])):
        ax_r.annotate(
            r"Ratios $\rightarrow\infty$ at $x=1$ because free $g_1\rightarrow0$",
            xy=(1.0, 1.0),
            xycoords=("data", "axes fraction"),
            xytext=(-5, -5),
            textcoords="offset points",
            ha="right",
            va="top",
            fontsize=FONT_SIZE - 1,
        )

    format_axis(ax_r, r"$g_1(^3\mathrm{He})/g_1(\mathrm{free})$")
    ax_r.set_xlabel(r"$x$")

    difference_p = x * (data["g1p_he3"] - data["g1p_free"])
    difference_n = x * (data["g1n_he3"] - data["g1n_free"])
    plot_line(
        ax_d,
        x,
        difference_p,
        label=r"$x(g_{1p}^{*}-g_{1p})$",
        color=COLORS["difference_p"],
        markers=args.markers,
    )
    plot_line(
        ax_d,
        x,
        difference_n,
        label=r"$x(g_{1n}^{*}-g_{1n})$",
        color=COLORS["difference_n"],
        linestyle="--",
        markers=args.markers,
    )
    ax_d.axhline(0.0, color="0.4", linewidth=1.0, linestyle=":")
    format_axis(ax_d, r"$x\,\Delta g_1=x(g_1^{*}-g_1)$")
    set_y_axis(
        ax_d,
        scale=args.difference_yscale,
        ymin=args.difference_ymin,
        ymax=args.difference_ymax,
        linthresh=CONFIG["difference_y_axis"]["linthresh"],
    )
    ax_d.set_xlabel(r"$x$")

    for ax in axes.ravel():
        ax.set_xscale(args.xscale)

    # Reserve the right 21% for legends and keep comfortable panel spacing.
    fig.subplots_adjust(
        left=0.085, right=0.77, bottom=0.105, top=0.88,
        wspace=0.22, hspace=0.12,
    )
    add_external_legends(fig, axes.ravel())

    # Keep bookkeeping such as the number of table rows out of the figure;
    # publication captions can carry that information if it is needed.
    heading = args.title
    if args.q2_label:
        heading += "\n" + args.q2_label
    if heading:
        fig.suptitle(heading, fontsize=FONT_SIZE + 2, y=0.965)
    return fig


def main() -> None:
    args = parse_args()
    data = restrict_x(read_table(args.input), args.xmin, args.xmax)
    if args.jam25:
        data.update(evaluate_jam25(data["x"], args.jam25_q2))
    fig = make_figure(data, args)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    print(f"Saved: {args.output.resolve()}")

    if args.show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
