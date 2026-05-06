#!/usr/bin/env python
"""
Post-hoc MSD convergence diagnostic.

Plots MSD(t) in sliding windows to assess whether 200 ps trajectories
have converged to a stable diffusion coefficient.

Usage:
    python scripts/plot_convergence.py \
        --traj outputs/LGPS/LGPS_T600K.traj \
        --mobile Li \
        --timestep 2.0 \
        --write-interval 10 \
        --n-windows 5
"""

from __future__ import annotations

import sys
from pathlib import Path

# Make `mlip_md` importable when running this script directly from inside
# `mlip_md_pipeline/` without an editable install or PYTHONPATH set.
_PKG_ROOT = Path(__file__).resolve().parent.parent
if str(_PKG_ROOT) not in sys.path:
    sys.path.insert(0, str(_PKG_ROOT))

import click
import logging
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


@click.command()
@click.option("--traj", "traj_path", required=True, help="ASE .traj file path.")
@click.option("--mobile", default="Li", show_default=True,
              help="Mobile species symbol.")
@click.option("--timestep", "timestep_fs", default=2.0, show_default=True,
              help="MD timestep in fs.")
@click.option("--write-interval", "write_interval", default=10, show_default=True,
              help="Trajectory write interval (steps).")
@click.option("--n-windows", "n_windows", default=5, show_default=True,
              help="Number of time windows to compare.")
@click.option("--output", "output_path", default=None,
              help="Output PNG path (default: next to traj file).")
def main(traj_path, mobile, timestep_fs, write_interval, n_windows, output_path):
    logging.basicConfig(level=logging.INFO)

    traj_path = Path(traj_path)
    if output_path is None:
        output_path = traj_path.with_suffix("_convergence.png")

    from mlip_md.md_runner import read_trajectory_positions
    from mlip_md.msd import compute_msd, fit_diffusivity

    positions, _ = read_trajectory_positions(traj_path, [mobile])
    n_frames = positions.shape[0]

    window_size = n_frames // n_windows
    if window_size < 50:
        logger.warning("Trajectory too short for %d windows; using 2.", n_windows)
        n_windows = 2
        window_size = n_frames // n_windows

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    D_windows = []

    for i in range(n_windows):
        start = i * window_size
        end = start + window_size
        pos_win = positions[start:end]

        lag_times, msd = compute_msd(
            pos_win, timestep_fs=timestep_fs, write_interval=write_interval
        )
        D, _, R2 = fit_diffusivity(lag_times, msd)
        D_windows.append(D)

        label = f"Window {i+1}: D={D:.2e} cm²/s"
        ax1.plot(lag_times, msd, label=label, alpha=0.8)

    ax1.set_xlabel("Lag time (ps)")
    ax1.set_ylabel("MSD (Å²)")
    ax1.set_title(f"MSD windows — {traj_path.stem}")
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    ax2.bar(range(1, n_windows + 1), D_windows, color="#1565C0", alpha=0.8)
    ax2.axhline(np.mean(D_windows), color="#F44336", linewidth=2,
                linestyle="--", label=f"Mean D = {np.mean(D_windows):.2e}")
    ax2.set_xlabel("Window index")
    ax2.set_ylabel("D (cm²/s)")
    ax2.set_title("D convergence across windows")
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3, axis="y")

    # Coefficient of variation
    cv = np.std(D_windows) / np.mean(D_windows) if np.mean(D_windows) > 0 else float("nan")
    fig.suptitle(
        f"{traj_path.stem} | CV(D) = {cv:.1%} "
        f"({'converged' if cv < 0.3 else 'NOT converged — extend trajectory'})",
        fontsize=11,
    )

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=150)
    plt.close(fig)
    logger.info("Convergence plot saved: %s", output_path)
    print(f"D values per window (cm²/s): {[f'{d:.3e}' for d in D_windows]}")
    print(f"CV = {cv:.1%}")


if __name__ == "__main__":
    main()
