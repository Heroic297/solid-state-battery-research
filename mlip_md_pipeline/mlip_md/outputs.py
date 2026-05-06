"""
Output module: CSV summary table and publication-quality plots.

Outputs
-------
outputs/
  screening_results.csv        — one row per candidate, all key metrics
  {candidate_id}_msd.png       — MSD(t) at all temperatures
  {candidate_id}_arrhenius.png — Arrhenius plot ln(D) vs 1000/T
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")   # non-interactive backend for headless/server use
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

logger = logging.getLogger(__name__)

# Colour palette consistent across temperatures
_TEMP_COLORS = {
    600:  "#2196F3",   # blue
    800:  "#FF9800",   # orange
    1000: "#F44336",   # red
}
_DEFAULT_COLOR = "#607D8B"


def save_results_csv(
    arrhenius_results: list,   # list of ArrheniusResult
    output_dir: Path,
    filename: str = "screening_results.csv",
) -> Path:
    """
    Compile all ArrheniusResult objects into a flat CSV.

    Columns (matching downstream stability/report agent schema)
    ----------------------------------------------------------
    candidate_id, mobile_species,
    T600_D_cm2s, T800_D_cm2s, T1000_D_cm2s,
    T600_D_std_cm2s, T800_D_std_cm2s, T1000_D_std_cm2s,
    T600_R2_msd, T800_R2_msd, T1000_R2_msd,
    T600_sigma_S_cm, T800_sigma_S_cm, T1000_sigma_S_cm,
    Ea_eV, Ea_std_eV, D0_cm2s, D0_std_cm2s,
    arrhenius_R2,
    sigma_300K_extrap_S_cm,
    n_temperatures_fit, status
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for res in arrhenius_results:
        row: dict = {
            "candidate_id": res.candidate_id,
            "mobile_species": res.mobile_species,
            "Ea_eV": _safe_round(res.Ea_eV, 4),
            "Ea_std_eV": _safe_round(res.Ea_std_eV, 4),
            "D0_cm2s": f"{res.D0_cm2s:.4e}",
            "D0_std_cm2s": f"{res.D0_std_cm2s:.4e}",
            "arrhenius_R2": _safe_round(res.arrhenius_r_squared, 4),
            "sigma_300K_extrap_S_cm": f"{res.sigma_S_cm_300K_extrap:.4e}",
            "n_temperatures_fit": len(res.temperatures_K),
            "status": "ok",
        }

        for i, T in enumerate(res.temperatures_K):
            t_int = int(T)
            row[f"T{t_int}_D_cm2s"] = f"{res.diffusivities_cm2s[i]:.4e}"
            row[f"T{t_int}_D_std_cm2s"] = f"{res.diffusivities_std_cm2s[i]:.4e}"
            row[f"T{t_int}_R2_msd"] = _safe_round(res.r_squared_msd[i], 4)
            row[f"T{t_int}_sigma_S_cm"] = f"{res.sigma_S_cm[i]:.4e}"

        rows.append(row)

    df = pd.DataFrame(rows)
    out_path = output_dir / filename
    df.to_csv(out_path, index=False)
    logger.info("Saved CSV: %s  (%d rows)", out_path, len(df))
    return out_path


def plot_msd(
    msd_results: list,    # list of MSDResult for ONE candidate
    output_dir: Path,
    candidate_id: str,
) -> Path:
    """
    Plot MSD(τ) curves for all temperatures on a single figure.

    Returns path to the saved PNG.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 5))

    sorted_results = sorted(msd_results, key=lambda r: r.temperature_K)

    for res in sorted_results:
        T = int(res.temperature_K)
        color = _TEMP_COLORS.get(T, _DEFAULT_COLOR)
        label = f"{T} K"

        ax.plot(
            res.lag_times_ps,
            res.msd_angstrom2,
            color=color,
            linewidth=1.5,
            label=label,
            alpha=0.9,
        )

        # Highlight the fitting window
        fit_mask = (res.lag_times_ps >= res.fit_start_ps) & (res.lag_times_ps <= res.fit_end_ps)
        ax.plot(
            res.lag_times_ps[fit_mask],
            res.msd_angstrom2[fit_mask],
            color=color,
            linewidth=3.0,
            alpha=0.5,
        )

        # Linear fit line
        if res.diffusivity_cm2s > 0:
            D_Aps = res.diffusivity_cm2s / 1e-4   # convert cm²/s → Å²/ps
            t_range = np.linspace(res.fit_start_ps, res.fit_end_ps, 100)
            msd_fit = 6 * D_Aps * t_range
            ax.plot(t_range, msd_fit, "--", color=color, linewidth=1.2,
                    label=f"fit {T} K (D={res.diffusivity_cm2s:.2e} cm²/s)")

    ax.set_xlabel("Lag time (ps)", fontsize=12)
    ax.set_ylabel("MSD (Å²)", fontsize=12)
    ax.set_title(f"MSD — {candidate_id}", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    plt.tight_layout()

    out_path = output_dir / f"{candidate_id}_msd.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    logger.info("Saved MSD plot: %s", out_path)
    return out_path


def plot_arrhenius(
    arrhenius_result,   # ArrheniusResult
    output_dir: Path,
    extrap_T_range_K: Optional[tuple[float, float]] = (300.0, 1100.0),
) -> Path:
    """
    Arrhenius plot: ln(D) vs 1000/T with fit line.

    Returns path to the saved PNG.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    res = arrhenius_result
    T_arr = res.temperatures_K
    D_arr = res.diffusivities_cm2s
    D_std = res.diffusivities_std_cm2s

    # 1000/T for x-axis
    x_data = 1000.0 / T_arr

    fig, ax = plt.subplots(figsize=(6, 5))

    # Data points with error bars
    ax.errorbar(
        x_data,
        np.log(D_arr),
        yerr=D_std / D_arr,    # σ(ln D) ≈ σ(D)/D
        fmt="o",
        color="#1565C0",
        markersize=8,
        elinewidth=1.5,
        capsize=4,
        label="MD data",
        zorder=5,
    )

    # Fit line (extrapolated over requested range)
    from mlip_md.arrhenius import kB_eV_per_K
    if extrap_T_range_K:
        T_fit = np.linspace(extrap_T_range_K[0], extrap_T_range_K[1], 200)
        D_fit = res.D0_cm2s * np.exp(-res.Ea_eV / (kB_eV_per_K * T_fit))
        x_fit = 1000.0 / T_fit
        ax.plot(
            x_fit,
            np.log(D_fit),
            "--",
            color="#F44336",
            linewidth=1.8,
            label=(
                f"Arrhenius fit\n"
                f"$E_a$ = {res.Ea_eV:.3f} ± {res.Ea_std_eV:.3f} eV\n"
                f"$R^2$ = {res.arrhenius_r_squared:.3f}"
            ),
            zorder=4,
        )

    ax.set_xlabel("1000 / T  (K⁻¹)", fontsize=12)
    ax.set_ylabel("ln D  (D in cm²/s)", fontsize=12)
    ax.set_title(
        f"Arrhenius — {res.candidate_id}  [{res.mobile_species}⁺]",
        fontsize=13, fontweight="bold"
    )
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3)

    # Secondary x-axis: temperature in K
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    # Tick positions in 1000/T corresponding to nice T values
    T_ticks = np.array([400, 500, 600, 700, 800, 900, 1000])
    x_ticks = 1000.0 / T_ticks
    x_ticks_in_range = x_ticks[(x_ticks >= ax.get_xlim()[0]) &
                                (x_ticks <= ax.get_xlim()[1])]
    T_labels = (1000.0 / x_ticks_in_range).astype(int)
    ax2.set_xticks(x_ticks_in_range)
    ax2.set_xticklabels([f"{T}" for T in T_labels], fontsize=9)
    ax2.set_xlabel("T (K)", fontsize=11)

    plt.tight_layout()

    out_path = output_dir / f"{res.candidate_id}_arrhenius.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    logger.info("Saved Arrhenius plot: %s", out_path)
    return out_path


def _safe_round(value, decimals: int):
    try:
        return round(float(value), decimals)
    except (TypeError, ValueError):
        return None
