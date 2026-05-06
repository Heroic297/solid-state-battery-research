"""
Mean-squared displacement (MSD) extraction and diffusivity estimation.

MSD definition
--------------
  MSD(τ) = < |r(t + τ) − r(t)|² >

averaged over all mobile ions and all time origins t.

The Einstein relation gives the self-diffusion coefficient D:

  MSD(τ) = 6 D τ   (3-D isotropic diffusion)

We fit a linear regression to the diffusive regime (beyond the ballistic
onset, typically τ > ~1 ps) to extract D.

Window-averaged MSD
-------------------
We use the standard overlapping time-origin (OTO) algorithm, which is
equivalent to the windowed correlation approach but O(N log N) via FFT.
This is implemented in ase.utils.msd or manually here for robustness.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.stats import linregress

logger = logging.getLogger(__name__)


@dataclass
class MSDResult:
    """Container for MSD analysis at a single temperature."""
    temperature_K: float
    candidate_id: str
    lag_times_ps: np.ndarray           # (n_lags,) array of lag times in ps
    msd_angstrom2: np.ndarray          # (n_lags,) MSD in Å²
    diffusivity_cm2s: float            # D in cm²/s from linear fit
    diffusivity_std_cm2s: float        # 1-σ uncertainty on D (from fit residuals)
    r_squared: float                   # goodness-of-fit R² of MSD linear regression
    fit_start_ps: float                # start of the fitting window (ps)
    fit_end_ps: float                  # end of the fitting window (ps)
    n_mobile_ions: int


def compute_msd(
    positions: np.ndarray,
    timestep_fs: float,
    write_interval: int = 1,
    max_lag_fraction: float = 0.5,
    fit_start_fraction: float = 0.1,
    fit_end_fraction: float = 0.9,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute MSD using overlapping time origins (OTO) via FFT.

    Parameters
    ----------
    positions : np.ndarray, shape (n_frames, n_atoms, 3)
        Unwrapped Cartesian positions in Ångström.
    timestep_fs : float
        MD timestep in femtoseconds.
    write_interval : int
        Frames were written every `write_interval` MD steps.
    max_lag_fraction : float
        Fraction of total trajectory to use as maximum lag (default 0.5
        ensures reasonable statistics for all time origins).
    fit_start_fraction : float
        Fraction of max lag time at which to start linear fit.
    fit_end_fraction : float
        Fraction of max lag time at which to end linear fit.

    Returns
    -------
    lag_times_ps : np.ndarray, shape (n_lags,)
    msd_angstrom2 : np.ndarray, shape (n_lags,)
    """
    n_frames, n_atoms, _ = positions.shape
    frame_time_ps = timestep_fs * write_interval / 1000.0  # ps per frame

    max_lag_frames = max(1, int(n_frames * max_lag_fraction))

    # OTO MSD via FFT (autocorrelation trick)
    # For each atom: MSD(τ) = < |r(t+τ) - r(t)|² >_t
    msd_sum = np.zeros(max_lag_frames)
    counts = np.zeros(max_lag_frames)

    for atom_idx in range(n_atoms):
        r = positions[:, atom_idx, :]  # (n_frames, 3)
        atom_msd = _msd_fft_1d(r, max_lag_frames)
        msd_sum += atom_msd
        counts += np.arange(n_frames, n_frames - max_lag_frames, -1)

    # Actually counts is the same for all atoms; recompute correctly
    n_origins = np.arange(n_frames, n_frames - max_lag_frames, -1)
    msd = msd_sum / (n_atoms * n_origins)

    lag_times_ps = np.arange(max_lag_frames) * frame_time_ps

    return lag_times_ps, msd


def _msd_fft_1d(r: np.ndarray, max_lag: int) -> np.ndarray:
    """
    Compute the sum (over time origins) of |r(t+τ) - r(t)|² for one atom.

    Uses the FFT autocorrelation identity:
      MSD(τ) = <r²(t+τ)> + <r²(t)> − 2 <r(t+τ)·r(t)>

    Parameters
    ----------
    r : (n_frames, 3) array of positions
    max_lag : int

    Returns
    -------
    msd_sum : (max_lag,) array — sum over time origins (not normalised)
    """
    n = len(r)
    # Squared norm per frame
    r2 = np.sum(r ** 2, axis=1)  # (n_frames,)

    # Sum of squared norms for each lag τ:
    # S1(τ) = Σ_{t=0}^{N-τ-1} [r²(t+τ) + r²(t)]
    cum_r2 = np.cumsum(r2)
    cum_r2 = np.concatenate([[0], cum_r2])

    # Autocorrelation via FFT for dot product term
    # Pad to next power of 2 for efficiency
    nfft = 2 ** int(np.ceil(np.log2(2 * n - 1)))
    msd_sum = np.zeros(max_lag)

    for dim in range(r.shape[1]):
        f = np.fft.rfft(r[:, dim], n=nfft)
        ac = np.fft.irfft(f * np.conj(f))[:n]  # autocorrelation
        for tau in range(max_lag):
            n_origins = n - tau
            s1 = (cum_r2[n] - cum_r2[n - n_origins] +
                  cum_r2[n_origins] - cum_r2[0])
            msd_sum[tau] += s1 - 2 * ac[tau]

    return msd_sum


def fit_diffusivity(
    lag_times_ps: np.ndarray,
    msd_angstrom2: np.ndarray,
    fit_start_fraction: float = 0.1,
    fit_end_fraction: float = 0.9,
) -> tuple[float, float, float]:
    """
    Fit MSD(τ) = 6D·τ in the diffusive regime.

    Parameters
    ----------
    lag_times_ps : np.ndarray
    msd_angstrom2 : np.ndarray
    fit_start_fraction : float
        Fraction of max lag at which to start fit (skip ballistic regime).
    fit_end_fraction : float
        Fraction of max lag at which to end fit.

    Returns
    -------
    D_cm2s : float
        Self-diffusion coefficient in cm²/s.
    D_std_cm2s : float
        Standard error of slope converted to cm²/s.
    r_squared : float
        R² of the linear fit.

    Notes
    -----
    MSD [Å²] vs t [ps]  →  slope [Å²/ps]
    D [cm²/s] = slope / 6  [Å²/ps] × (1e-8)² cm²/Å² × 1e12 ps/s
             = slope / 6 × 1e-4  [cm²/s]
    """
    n = len(lag_times_ps)
    i_start = max(1, int(n * fit_start_fraction))
    i_end = max(i_start + 2, int(n * fit_end_fraction))

    t_fit = lag_times_ps[i_start:i_end]
    msd_fit = msd_angstrom2[i_start:i_end]

    result = linregress(t_fit, msd_fit)
    slope = result.slope        # Å²/ps
    slope_stderr = result.stderr

    # Conversion: Å²/ps → cm²/s
    # 1 Å = 1e-8 cm  →  1 Å² = 1e-16 cm²
    # 1 ps = 1e-12 s
    # D [cm²/s] = slope [Å²/ps] × (1e-16 cm²/Å²) / (1e-12 s/ps)
    #           = slope × 1e-4  [cm²/s]
    conversion = 1e-4 / 6.0

    D_cm2s = max(0.0, slope * conversion)
    D_std_cm2s = abs(slope_stderr * conversion)
    r_squared = result.rvalue ** 2

    logger.info(
        "MSD fit: slope=%.4e Å²/ps  D=%.4e cm²/s  R²=%.4f",
        slope, D_cm2s, r_squared,
    )
    return D_cm2s, D_std_cm2s, r_squared


def analyse_trajectory(
    traj_path,
    temperature_K: float,
    candidate_id: str,
    mobile_species: list[str],
    timestep_fs: float,
    write_interval: int = 10,
    fit_start_fraction: float = 0.10,
    fit_end_fraction: float = 0.90,
) -> MSDResult:
    """
    Full pipeline: trajectory → MSDResult.

    Calls read_trajectory_positions from md_runner, then compute_msd and
    fit_diffusivity.
    """
    from mlip_md.md_runner import read_trajectory_positions

    logger.info(
        "Analysing trajectory: %s (T=%.0f K)", traj_path, temperature_K
    )

    positions, _ = read_trajectory_positions(traj_path, mobile_species)
    n_frames, n_mobile, _ = positions.shape
    logger.info("Trajectory: %d frames, %d mobile ions", n_frames, n_mobile)

    lag_times_ps, msd_a2 = compute_msd(
        positions=positions,
        timestep_fs=timestep_fs,
        write_interval=write_interval,
    )

    D, D_std, R2 = fit_diffusivity(
        lag_times_ps, msd_a2,
        fit_start_fraction=fit_start_fraction,
        fit_end_fraction=fit_end_fraction,
    )

    max_lag = lag_times_ps[-1] if len(lag_times_ps) else 0.0
    fit_start = max_lag * fit_start_fraction
    fit_end = max_lag * fit_end_fraction

    return MSDResult(
        temperature_K=temperature_K,
        candidate_id=candidate_id,
        lag_times_ps=lag_times_ps,
        msd_angstrom2=msd_a2,
        diffusivity_cm2s=D,
        diffusivity_std_cm2s=D_std,
        r_squared=R2,
        fit_start_ps=fit_start,
        fit_end_ps=fit_end,
        n_mobile_ions=n_mobile,
    )
