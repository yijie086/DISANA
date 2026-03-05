import os
import glob
import re
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# --- Physics Constants from PAC39 / CLAS12 ---
M_P_GEV = 0.93827
M_PHI_GEV = 1.019461
W_MEAN = 2.8

# Map Q-bin index (from filenames) -> <Q^2> [GeV^2]
Q2_MAP = {0: 1.0, 1: 2.0, 2: 3.0, 3: 4.0}

# Conversion: 1 GeV^-2 = (0.197326 fm)^2
GEV2_TO_FM2 = 0.197326**2


def get_t_min(Q2: float, W: float) -> float:
    """Calculates the |t_min| boundary to shift t' to absolute |t|."""
    W2, M2, mf2 = W**2, M_P_GEV**2, M_PHI_GEV**2
    E_gstar = (W2 - M2 - Q2) / (2.0 * W)
    p_gstar = np.sqrt(max(E_gstar**2 + Q2, 0.0))
    E_phi_cm = (W2 + mf2 - M2) / (2.0 * W)
    p_phi_cm = np.sqrt(max(E_phi_cm**2 - mf2, 0.0))
    return abs((E_gstar - E_phi_cm)**2 - (p_gstar - p_phi_cm)**2)


def dipole_model(t_abs: np.ndarray, sigma0: float, mg2: float) -> np.ndarray:
    """Dipole form factor: sigma(|t|) = sigma0 / (1 + |t|/mg2)^4"""
    return sigma0 / (1.0 + t_abs / mg2) ** 4


def _safe_sigma(y_err: np.ndarray) -> np.ndarray:
    """curve_fit requires strictly positive sigmas; replace non-positive with a small floor."""
    y_err = np.asarray(y_err, dtype=float)
    finite = np.isfinite(y_err)
    if not finite.any():
        return np.ones_like(y_err)
    floor = np.nanmedian(y_err[finite]) * 1e-3
    if not np.isfinite(floor) or floor <= 0:
        floor = 1e-12
    y_err = np.where((~finite) | (y_err <= 0), floor, y_err)
    return y_err


def run_analysis(input_dir: str, output_dir: str) -> None:
    os.makedirs(output_dir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(input_dir, "table_*_Q*_W*.csv")))
    results = []
    combined = []  # store per-Q data+fit info for a combined overlay plot

    for fpath in files:
        base = os.path.basename(fpath)
        match = re.search(r"_Q(\d+)", base)
        if not match:
            continue

        q_idx = int(match.group(1))
        q2 = float(Q2_MAP.get(q_idx, 1.0))

        df = pd.read_csv(fpath)
        df = df.replace([np.inf, -np.inf], np.nan).dropna()

        # Only keep physical (positive) reduced cross section
        df = df[df["reduced_xs"] > 0].copy()
        if len(df) < 3:
            continue

        t_min = get_t_min(q2, W_MEAN)
        t_abs = df["tprime_center"].to_numpy(dtype=float) + t_min

        # x_B using mean W for the bin set (as used in your script)
        xb = q2 / (W_MEAN**2 + q2 - M_P_GEV**2)

        y = df["reduced_xs"].to_numpy(dtype=float)

        # y_err = (rad-corrected xs err) / Gamma_v
        y_err = (df["xs_rad_corrected_err"].to_numpy(dtype=float) /
                 df["Gamma_v"].to_numpy(dtype=float))
        y_err = _safe_sigma(y_err)

        # Initial guess: sigma0 ~ first y-point, mg2 ~ 0.65 GeV^2
        p0 = [float(y[0]), 0.65]

        try:
            # Constrain mg2 to positive, keep sigma0 non-negative.
            popt, pcov = curve_fit(
                dipole_model,
                t_abs,
                y,
                sigma=y_err,
                p0=p0,
                absolute_sigma=True,
                bounds=([0.0, 1e-6], [np.inf, np.inf]),
                maxfev=10000,
            )

            sigma0, mg2 = popt
            mg2_err = float(np.sqrt(pcov[1, 1])) if np.isfinite(pcov[1, 1]) else np.nan

            # PAC39 Relation: <b^2> = 8/m_g^2  [GeV^-2]
            b2_gev = 8.0 / mg2
            b2_err_gev = 8.0 * mg2_err / (mg2**2) if np.isfinite(mg2_err) else np.nan

            # Convert to fm^2
            b2_fm = b2_gev * GEV2_TO_FM2
            b2_err_fm = b2_err_gev * GEV2_TO_FM2 if np.isfinite(b2_err_gev) else np.nan

            results.append(
                {"x_B": xb, "Q2": q2, "b2_fm": b2_fm, "b2_err_fm": b2_err_fm,
                 "mg2": mg2, "mg2_err": mg2_err, "Q_idx": q_idx, "file": base}
            )

            # --- Individual fit plot ---
            plt.figure(figsize=(8, 6))
            plt.errorbar(t_abs, y, yerr=y_err, fmt="ok", label="Data")

            t_fine = np.linspace(min(t_abs), max(t_abs) * 1.1, 200)
            plt.plot(
                t_fine,
                dipole_model(t_fine, *popt),
                "-",
                label=rf"$m_g^2={mg2:.3f}$ GeV$^2$",
            )

            plt.yscale("log")
            plt.xlabel(r"$|t|$ [GeV$^2$]")
            plt.ylabel(r"Reduced $\sigma$ [nb/GeV$^2$]")
            plt.title(rf"Q-bin {q_idx} ($Q^2={q2:.1f}$ GeV$^2$)")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"fit_Q{q_idx}.png"), dpi=200)
            plt.close()

            combined.append(
                {
                    "q_idx": q_idx,
                    "Q2": q2,
                    "xB": xb,
                    "t_abs": t_abs,
                    "y": y,
                    "y_err": y_err,
                    "popt": popt,
                }
            )

        except Exception as e:
            print(f"Fit failed for {base} (Q{q_idx}): {e}")

    # --- Combined overlay plot (all Q bins on one plot) ---
    if combined:
        plt.figure(figsize=(8, 6))

        for item in sorted(combined, key=lambda d: d["q_idx"]):
            q_idx = item["q_idx"]
            q2 = item["Q2"]
            xb = item["xB"]
            t_abs = item["t_abs"]
            y = item["y"]
            y_err = item["y_err"]
            sigma0, mg2 = item["popt"]

            plt.errorbar(
                t_abs,
                y,
                yerr=y_err,
                fmt="o",
                ms=4,
                capsize=2,
                label=rf"Q{q_idx}: $Q^2={q2:.1f}$, $x_B={xb:.3f}$",
            )

            t_fine = np.linspace(min(t_abs), max(t_abs) * 1.1, 200)
            plt.plot(
                t_fine,
                dipole_model(t_fine, sigma0, mg2),
                "-",
                alpha=0.9,
            )

        plt.yscale("log")
        plt.xlabel(r"$|t|$ [GeV$^2$]")
        plt.ylabel(r"Reduced $\sigma$ [nb/GeV$^2$]")
        plt.title("All Q bins: data + dipole fits")
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=9)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "combined_fits_allQ.png"), dpi=200)
        plt.close()

    # --- Summary plot: gluonic radius squared in fm^2 vs x_B ---
    if results:
        res_df = pd.DataFrame(results).sort_values("x_B")
        res_df.to_csv(os.path.join(output_dir, "gluon_radius_fm.csv"), index=False)

        plt.figure(figsize=(8, 6))
        plt.errorbar(
            res_df["x_B"],
            res_df["b2_fm"],
            yerr=res_df["b2_err_fm"],
            fmt="o--",
            capsize=5,
            label=r"$\langle b^2 \rangle$",
        )
        plt.xlabel(r"$x_B$")
        plt.ylabel(r"$\langle b^2 \rangle$ [fm$^2$]")
        plt.title("Proton gluonic radius squared (from dipole fits)")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "summary_radius_fm.png"), dpi=200)
        plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, help="Directory with CSV tables")
    parser.add_argument("--output-dir", required=True, help="Directory for plots")
    args = parser.parse_args()
    run_analysis(args.input_dir, args.output_dir)