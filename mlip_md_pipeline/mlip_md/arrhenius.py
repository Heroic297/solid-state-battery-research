"""
Arrhenius fitting and Nernst-Einstein ionic conductivity conversion.

Arrhenius model
---------------
  D(T) = D₀ exp(−Eₐ / kB T)

Linearised as:
  ln D = ln D₀ − Eₐ / (kB T)

A weighted linear regression of ln(D) vs 1/T gives:
  slope = −Eₐ / kB  →  Eₐ = −slope × kB
  intercept = ln D₀  →  D₀ = exp(intercept)

Nernst-Einstein conductivity
----------------------------
  σ = n q² D / (kB T)

where:
  n   = number density of mobile ions [m⁻³]
  q   = ionic charge [C]  (e = 1.602e-19 C for monovalent Li+/Na+)
  D   = diffusivity [m²/s]
  kB  = Boltzmann constant [J/K]
  T   = temperature [K]

Unit pipeline:
  D [cm²/s] → D [m²/s] by × 1e-4
  n [Å⁻³]   → n [m⁻³] by × 1e30
  σ [S/m]   → σ [S/cm] by × 1e-2

Important caveats
-----------------
- Nernst-Einstein assumes uncorrelated ion motion (Haven ratio H_R = 1).
  Real solid electrolytes exhibit H_R < 1 due to cooperative hopping.
  This pipeline does NOT apply Haven ratio corrections.
- Diffusivities from short (200 ps) MLIP-MD are likely under-converged
  for materials with Eₐ > ~0.4 eV. Treat values as screening estimates.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.stats import linregress

logger = logging.getLogger(__name__)

kB_eV_per_K = 8.617333262e-5   # eV K⁻¹
kB_J_per_K  = 1.380649e-23     # J K⁻¹
e_charge_C  = 1.602176634e-19  # C


@dataclass
class ArrheniusResult:
    """Arrhenius fit and Nernst-Einstein conductivity results for one candidate."""
    candidate_id: str
    mobile_species: str           # "Li" or "Na"
    temperatures_K: np.ndarray
    diffusivities_cm2s: np.ndarray
    diffusivities_std_cm2s: np.ndarray
    r_squared_msd: np.ndarray     # per-temperature MSD fit R²

    # Arrhenius fit
    Ea_eV: float
    Ea_std_eV: float              # propagated from linregress stderr
    D0_cm2s: float
    D0_std_cm2s: float
    arrhenius_r_squared: float    # R² of ln(D) vs 1/T fit

    # Nernst-Einstein at each temperature
    sigma_S_cm: np.ndarray        # σ [S/cm] at each T
    sigma_S_cm_300K_extrap: float # extrapolated to 300 K using Arrhenius fit


def fit_arrhenius(
    temperatures_K: np.ndarray,
    diffusivities_cm2s: np.ndarray,
    diffusivities_std_cm2s: Optional[np.ndarray] = None,
) -> tuple[float, float, float, float, float]:
    """
    Weighted Arrhenius regression on ln(D) vs 1/(kB T).

    Parameters
    ----------
    temperatures_K : array-like, shape (N,)
    diffusivities_cm2s : array-like, shape (N,)
    diffusivities_std_cm2s : array-like, shape (N,) or None
        If provided, used as weights (w = 1/σ²) in the regression.

    Returns
    -------
    Ea_eV : float
    Ea_std_eV : float
    D0_cm2s : float
    D0_std_cm2s : float
    r_squared : float
    """
    T = np.asarray(temperatures_K, dtype=float)
    D = np.asarray(diffusivities_cm2s, dtype=float)

    # Guard: remove non-positive D values (unphysical, MSD fit failed)
    valid = D > 0
    if valid.sum() < 2:
        raise ValueError(
            f"Need at least 2 positive diffusivities for Arrhenius fit. "
            f"Got D={D}"
        )

    T = T[valid]
    D = D[valid]

    x = 1.0 / (kB_eV_per_K * T)   # 1/(kB T) in eV⁻¹
    y = np.log(D)                  # ln(D [cm²/s])

    if diffusivities_std_cm2s is not None:
        std = np.asarray(diffusivities_std_cm2s, dtype=float)[valid]
        std = np.where(std > 0, std, 1e-30)
        # Weight by 1/σ² in log space: σ_lnD ≈ σ_D / D
        sigma_ln = std / D
        weights = 1.0 / sigma_ln ** 2
        # Weighted least-squares via numpy
        X = np.column_stack([np.ones_like(x), x])
        W = np.diag(weights)
        XtWX = X.T @ W @ X
        XtWy = X.T @ W @ y
        try:
            coeffs = np.linalg.solve(XtWX, XtWy)
        except np.linalg.LinAlgError:
            coeffs, _, _, _ = np.linalg.lstsq(XtWX, XtWy, rcond=None)
        intercept, slope = coeffs
        # Propagate errors
        cov = np.linalg.inv(XtWX)
        intercept_std = np.sqrt(cov[0, 0])
        slope_std = np.sqrt(cov[1, 1])
    else:
        result = linregress(x, y)
        slope = result.slope
        intercept = result.intercept
        slope_std = result.stderr
        intercept_std = result.intercept_stderr if hasattr(result, 'intercept_stderr') else abs(slope_std)

    # Arrhenius: ln D = ln D₀ − Eₐ/(kB T)  →  slope = −Eₐ [eV]
    Ea_eV = -slope
    Ea_std_eV = slope_std

    D0_cm2s = np.exp(intercept)
    D0_std_cm2s = D0_cm2s * intercept_std  # propagated from ln(D₀) uncertainty

    # R² of fit
    y_pred = intercept + slope * x
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    logger.info(
        "Arrhenius: Eₐ=%.4f ± %.4f eV  D₀=%.4e cm²/s  R²=%.4f",
        Ea_eV, Ea_std_eV, D0_cm2s, r_squared,
    )
    return Ea_eV, Ea_std_eV, D0_cm2s, D0_std_cm2s, r_squared


def nernst_einstein_conductivity(
    diffusivity_cm2s: float,
    temperature_K: float,
    n_mobile_ions: int,
    cell_volume_angstrom3: float,
    ion_charge: int = 1,
) -> float:
    """
    Compute ionic conductivity via Nernst-Einstein equation.

      σ = n q² D / (kB T)

    Parameters
    ----------
    diffusivity_cm2s : float
        Self-diffusion coefficient [cm²/s].
    temperature_K : float
    n_mobile_ions : int
        Number of mobile ions in the simulation cell.
    cell_volume_angstrom3 : float
        Simulation cell volume [Å³].
    ion_charge : int
        Formal charge of mobile ion (1 for Li⁺/Na⁺).

    Returns
    -------
    sigma_S_cm : float
        Ionic conductivity in S/cm.
    """
    # Number density: n [m⁻³] = N / V [Å⁻³] × 1e30 [Å³/m³]
    n_per_angstrom3 = n_mobile_ions / cell_volume_angstrom3
    n_m3 = n_per_angstrom3 * 1e30

    # Diffusivity: [m²/s]
    D_m2s = diffusivity_cm2s * 1e-4

    # Charge
    q_C = ion_charge * e_charge_C

    # σ [S/m]
    sigma_S_m = (n_m3 * q_C ** 2 * D_m2s) / (kB_J_per_K * temperature_K)

    # Convert to S/cm (1 S/m = 0.01 S/cm)
    sigma_S_cm = sigma_S_m * 1e-2

    logger.info(
        "Nernst-Einstein: T=%.0f K  D=%.4e cm²/s  n=%.4e Å⁻³  σ=%.4e S/cm",
        temperature_K, diffusivity_cm2s, n_per_angstrom3, sigma_S_cm,
    )
    return sigma_S_cm


def full_arrhenius_analysis(
    candidate_id: str,
    mobile_species: str,
    msd_results,           # list of MSDResult
    cell_volume_angstrom3: float,
    ion_charge: int = 1,
    extrap_temperature_K: float = 300.0,
) -> ArrheniusResult:
    """
    Run full Arrhenius + Nernst-Einstein analysis for one candidate.

    Parameters
    ----------
    candidate_id : str
    mobile_species : str
    msd_results : list of MSDResult
        One per temperature, from msd.analyse_trajectory().
    cell_volume_angstrom3 : float
    ion_charge : int
    extrap_temperature_K : float
        Temperature for extrapolated σ (typically 300 K for RT comparison).

    Returns
    -------
    ArrheniusResult
    """
    from mlip_md.msd import MSDResult

    results_sorted = sorted(msd_results, key=lambda r: r.temperature_K)
    T_arr = np.array([r.temperature_K for r in results_sorted])
    D_arr = np.array([r.diffusivity_cm2s for r in results_sorted])
    D_std_arr = np.array([r.diffusivity_std_cm2s for r in results_sorted])
    R2_arr = np.array([r.r_squared for r in results_sorted])
    n_mobile = results_sorted[0].n_mobile_ions

    Ea, Ea_std, D0, D0_std, arr_R2 = fit_arrhenius(T_arr, D_arr, D_std_arr)

    # Nernst-Einstein at each simulation temperature
    sigma_arr = np.array([
        nernst_einstein_conductivity(
            diffusivity_cm2s=D,
            temperature_K=T,
            n_mobile_ions=n_mobile,
            cell_volume_angstrom3=cell_volume_angstrom3,
            ion_charge=ion_charge,
        )
        for T, D in zip(T_arr, D_arr)
    ])

    # Extrapolated D and σ at extrap_temperature_K
    D_extrap = D0 * np.exp(-Ea / (kB_eV_per_K * extrap_temperature_K))
    sigma_extrap = nernst_einstein_conductivity(
        diffusivity_cm2s=D_extrap,
        temperature_K=extrap_temperature_K,
        n_mobile_ions=n_mobile,
        cell_volume_angstrom3=cell_volume_angstrom3,
        ion_charge=ion_charge,
    )

    return ArrheniusResult(
        candidate_id=candidate_id,
        mobile_species=mobile_species,
        temperatures_K=T_arr,
        diffusivities_cm2s=D_arr,
        diffusivities_std_cm2s=D_std_arr,
        r_squared_msd=R2_arr,
        Ea_eV=Ea,
        Ea_std_eV=Ea_std,
        D0_cm2s=D0,
        D0_std_cm2s=D0_std,
        arrhenius_r_squared=arr_R2,
        sigma_S_cm=sigma_arr,
        sigma_S_cm_300K_extrap=sigma_extrap,
    )
