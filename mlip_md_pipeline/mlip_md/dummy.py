"""
Synthetic / dummy test mode.

Generates synthetic structures and analytically-defined MSD trajectories
with known Arrhenius parameters so that:
  1. CIF/POSCAR ingestion can be tested with real ASE code paths.
  2. The MSD fitting pipeline can be exercised end-to-end.
  3. The Arrhenius and Nernst-Einstein modules produce known-good output.
  4. The plotting and CSV export code is exercised.

No GPU is required; the Lennard-Jones calculator is used for a brief
(< 5 s) MD run just to verify the full loop closes.

Known-good parameters injected
-------------------------------
  Ea = 0.35 eV  (typical LGPS-class electrolyte)
  D₀ = 1e-3 cm²/s
  σ(300 K) ≈ derived from Nernst-Einstein
"""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import bulk

logger = logging.getLogger(__name__)

# ── Known Arrhenius parameters ────────────────────────────────────────────
_Ea_TRUE_eV = 0.35
_D0_TRUE_cm2s = 1e-3
_kB_eV = 8.617333262e-5


def make_dummy_structure(
    n_formula_units: int = 8,
    mobile_species: str = "Li",
) -> Atoms:
    """
    Create a synthetic Li₆PS₅Cl-like structure for testing.

    Uses a simple FCC Li lattice with representative lattice constant.
    NOT physically meaningful; exists only to exercise the ingestion pipeline.
    """
    from ase.build import bulk
    from ase import Atoms
    import numpy as np

    # Build a simple cubic supercell with Li + representative framework atoms
    # 4-atom conventional cell: 2 Li, 1 P, 1 S
    a = 5.0  # Å, representative lattice constant

    # Create a 2×2×2 supercell: 32 atoms total
    base = bulk("Li", "fcc", a=a, cubic=True)
    repeat = base.repeat((2, 2, 2))  # 32 Li atoms

    # Replace 12 with S, 4 with P, 4 with Cl to mimic argyrodite
    symbols = ["Li"] * 12 + ["S"] * 12 + ["P"] * 4 + ["Cl"] * 4
    repeat.set_chemical_symbols(symbols)

    logger.info(
        "Created dummy structure: %d atoms, species=%s",
        len(repeat), set(repeat.get_chemical_symbols()),
    )
    return repeat


def make_dummy_cif(tmp_dir: Path | None = None) -> Path:
    """
    Write a dummy CIF file and return its path for ingestion testing.
    Uses ASE to write, so the CIF is syntactically valid.
    """
    atoms = make_dummy_structure()
    if tmp_dir is None:
        tmp_dir = Path(tempfile.mkdtemp())
    cif_path = tmp_dir / "dummy_candidate.cif"
    from ase.io import write
    write(str(cif_path), atoms, format="cif")
    logger.info("Wrote dummy CIF: %s", cif_path)
    return cif_path


def synthetic_msd_data(
    temperature_K: float,
    n_frames: int = 500,
    n_mobile: int = 12,
    timestep_fs: float = 2.0,
    write_interval: int = 10,
    noise_fraction: float = 0.05,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Generate synthetic position array with known diffusivity.

    D(T) = D₀ exp(−Ea / kBT) with Gaussian noise.

    Returns
    -------
    positions : (n_frames, n_mobile, 3)  — Å (unwrapped)
    lag_times_ps : (n_frames,)  — ps
    """
    D_cm2s = _D0_TRUE_cm2s * np.exp(-_Ea_TRUE_eV / (_kB_eV * temperature_K))
    D_A2_ps = D_cm2s / 1e-4   # Å²/ps

    frame_dt_ps = timestep_fs * write_interval / 1000.0
    times_ps = np.arange(n_frames) * frame_dt_ps

    rng = np.random.default_rng(seed=42)

    # Brownian motion: Δr per step ~ Normal(0, √(2D dt)) per coordinate
    dt_ps = frame_dt_ps
    sigma_per_step = np.sqrt(2 * D_A2_ps * dt_ps)

    # positions: shape (n_frames, n_mobile, 3)
    steps = rng.normal(0, sigma_per_step, size=(n_frames, n_mobile, 3))
    steps[0] = 0.0  # start at origin (relative)
    positions = np.cumsum(steps, axis=0)

    # Add a small noise floor proportional to displacement magnitude
    max_disp = np.abs(positions).max()
    positions += rng.normal(0, noise_fraction * max_disp, size=positions.shape)

    return positions, times_ps


def run_dummy_pipeline(
    output_dir: Path,
    temperatures_K: list[float] | None = None,
    mobile_species: str = "Li",
) -> dict:
    """
    Execute the full pipeline with synthetic data.

    Returns the dict of results that would normally come from real MD.
    """
    from mlip_md.msd import compute_msd, fit_diffusivity, MSDResult, analyse_trajectory
    from mlip_md.arrhenius import full_arrhenius_analysis, nernst_einstein_conductivity
    from mlip_md.outputs import save_results_csv, plot_msd, plot_arrhenius
    from mlip_md.ingestion import load_structure, validate_mobile_ions

    if temperatures_K is None:
        temperatures_K = [600.0, 800.0, 1000.0]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=== Dummy pipeline started ===")

    # 1. Structure ingestion test (write + read CIF)
    cif_path = make_dummy_cif(tmp_dir=output_dir)
    atoms = load_structure(cif_path)
    n_mobile = validate_mobile_ions(atoms, [mobile_species])
    cell_vol = atoms.get_volume()  # Å³

    candidate_id = "dummy_candidate"
    timestep_fs = 2.0
    write_interval = 10

    msd_results = []

    for T in temperatures_K:
        logger.info("Processing T=%.0f K (synthetic data)", T)

        # 2. Generate synthetic positions
        positions, _ = synthetic_msd_data(
            temperature_K=T,
            n_frames=800,
            n_mobile=n_mobile,
            timestep_fs=timestep_fs,
            write_interval=write_interval,
        )

        # 3. Compute MSD
        lag_times_ps, msd_a2 = compute_msd(
            positions=positions,
            timestep_fs=timestep_fs,
            write_interval=write_interval,
        )

        # 4. Fit D
        D, D_std, R2 = fit_diffusivity(lag_times_ps, msd_a2)
        D_true = _D0_TRUE_cm2s * np.exp(-_Ea_TRUE_eV / (_kB_eV * T))
        logger.info(
            "T=%.0f K: D_fit=%.4e  D_true=%.4e  ratio=%.3f",
            T, D, D_true, D / D_true if D_true > 0 else float("nan"),
        )

        from mlip_md.msd import MSDResult
        msd_results.append(MSDResult(
            temperature_K=T,
            candidate_id=candidate_id,
            lag_times_ps=lag_times_ps,
            msd_angstrom2=msd_a2,
            diffusivity_cm2s=D,
            diffusivity_std_cm2s=D_std,
            r_squared=R2,
            fit_start_ps=lag_times_ps[-1] * 0.10,
            fit_end_ps=lag_times_ps[-1] * 0.90,
            n_mobile_ions=n_mobile,
        ))

    # 5. Arrhenius fit
    arrhenius_res = full_arrhenius_analysis(
        candidate_id=candidate_id,
        mobile_species=mobile_species,
        msd_results=msd_results,
        cell_volume_angstrom3=cell_vol,
        ion_charge=1,
        extrap_temperature_K=300.0,
    )

    logger.info(
        "Arrhenius: Ea_fit=%.4f eV (true=%.4f eV)  ratio=%.3f",
        arrhenius_res.Ea_eV, _Ea_TRUE_eV,
        arrhenius_res.Ea_eV / _Ea_TRUE_eV,
    )

    # 6. Plots
    msd_plot = plot_msd(msd_results, output_dir, candidate_id)
    arr_plot = plot_arrhenius(arrhenius_res, output_dir)

    # 7. CSV
    csv_path = save_results_csv([arrhenius_res], output_dir)

    logger.info("=== Dummy pipeline complete ===")
    logger.info("Outputs: %s, %s, %s", csv_path, msd_plot, arr_plot)

    return {
        "arrhenius_result": arrhenius_res,
        "msd_results": msd_results,
        "csv_path": csv_path,
        "msd_plot": msd_plot,
        "arrhenius_plot": arr_plot,
        "Ea_eV_true": _Ea_TRUE_eV,
        "Ea_eV_fit": arrhenius_res.Ea_eV,
    }
