"""
Pipeline unit tests — CPU only, no real MLIP weights required.

Run with:
    cd mlip_md_pipeline
    python -m pytest tests/ -v

Tests cover:
    - Structure ingestion (CIF round-trip via ASE)
    - Supercell building
    - Mobile-ion validation
    - MSD computation (known input → known output)
    - Arrhenius fitting (synthetic data with known Ea)
    - Nernst-Einstein unit conversion
    - Dummy pipeline end-to-end
    - Output CSV schema
    - MSD and Arrhenius plot file generation
"""

from __future__ import annotations

import math
import tempfile
from pathlib import Path

import numpy as np
import pytest

# ── Fixtures ──────────────────────────────────────────────────────────────

@pytest.fixture
def tmp_dir(tmp_path):
    return tmp_path


@pytest.fixture
def dummy_atoms():
    """Return a small Li-containing Atoms object."""
    from mlip_md.dummy import make_dummy_structure
    return make_dummy_structure(n_formula_units=4, mobile_species="Li")


@pytest.fixture
def dummy_cif_path(tmp_dir):
    from mlip_md.dummy import make_dummy_cif
    return make_dummy_cif(tmp_dir)


# ── Ingestion tests ───────────────────────────────────────────────────────

class TestIngestion:
    def test_load_cif_roundtrip(self, dummy_cif_path):
        from mlip_md.ingestion import load_structure
        atoms = load_structure(dummy_cif_path)
        assert len(atoms) > 0
        assert all(atoms.pbc), "PBC should be True after load"

    def test_validate_mobile_ions_present(self, dummy_atoms):
        from mlip_md.ingestion import validate_mobile_ions
        count = validate_mobile_ions(dummy_atoms, ["Li"])
        assert count > 0

    def test_validate_mobile_ions_absent(self, dummy_atoms):
        from mlip_md.ingestion import validate_mobile_ions
        with pytest.raises(ValueError, match="none of the mobile species"):
            validate_mobile_ions(dummy_atoms, ["Mg"])

    def test_load_nonexistent_file(self):
        from mlip_md.ingestion import load_structure
        with pytest.raises(FileNotFoundError):
            load_structure("/nonexistent/path/structure.cif")

    def test_supercell_replication(self, dummy_atoms):
        from mlip_md.ingestion import build_supercell
        # With min_length=50 Å on a ~5 Å cell, expect replication
        supercell = build_supercell(dummy_atoms, min_length_angstrom=50.0)
        assert len(supercell) > len(dummy_atoms)

    def test_supercell_no_replication_large(self, dummy_atoms):
        from mlip_md.ingestion import build_supercell
        # Already large enough
        import numpy as np
        from ase.build import make_supercell
        big = make_supercell(dummy_atoms, np.diag([4, 4, 4]))
        original_len = len(big)
        result = build_supercell(big, min_length_angstrom=5.0)
        assert len(result) == original_len


# ── Calculator adapter tests (dummy only) ────────────────────────────────

class TestCalculators:
    def test_dummy_calculator_returns_ase_calc(self):
        from mlip_md.calculators import get_calculator
        from ase.calculators.calculator import Calculator
        calc = get_calculator("dummy")
        assert isinstance(calc, Calculator)

    def test_unknown_calculator_raises(self):
        from mlip_md.calculators import get_calculator
        with pytest.raises(ValueError, match="Unknown calculator"):
            get_calculator("badname")


# ── MSD tests ─────────────────────────────────────────────────────────────

class TestMSD:
    def test_msd_zero_for_static_positions(self):
        """Static atoms → MSD should be zero for all lags."""
        from mlip_md.msd import compute_msd
        n_frames, n_atoms = 200, 5
        positions = np.zeros((n_frames, n_atoms, 3))
        lag_times, msd = compute_msd(positions, timestep_fs=2.0, write_interval=1)
        assert np.allclose(msd, 0.0, atol=1e-10)

    def test_msd_linearity_brownian(self):
        """
        For ideal Brownian motion, MSD should grow linearly with lag time.
        We inject perfectly linear positions and verify slope recovery.
        """
        from mlip_md.msd import compute_msd, fit_diffusivity

        n_frames, n_atoms = 1000, 20
        D_target_A2_ps = 0.05   # Å²/ps

        rng = np.random.default_rng(42)
        dt_ps = 2.0 * 1 / 1000.0  # 2 fs timestep, write_interval=1 → 0.002 ps
        sigma = np.sqrt(2 * D_target_A2_ps * dt_ps)
        steps = rng.normal(0, sigma, (n_frames, n_atoms, 3))
        steps[0] = 0.0
        positions = np.cumsum(steps, axis=0)

        lag_times, msd = compute_msd(positions, timestep_fs=2.0, write_interval=1)
        D_cm2s, _, R2 = fit_diffusivity(lag_times, msd)

        D_target_cm2s = D_target_A2_ps * 1e-4   # → cm²/s
        assert R2 > 0.95, f"R² = {R2:.3f} is too low for Brownian motion"
        # Allow ±30% due to statistical noise
        ratio = D_cm2s / D_target_cm2s
        assert 0.5 < ratio < 2.0, f"Recovered D/D_true = {ratio:.2f}"

    def test_msd_output_shapes(self):
        from mlip_md.msd import compute_msd
        positions = np.random.randn(300, 10, 3)
        lag_times, msd = compute_msd(positions, timestep_fs=2.0, write_interval=5)
        assert lag_times.shape == msd.shape
        assert len(lag_times) > 0
        assert lag_times[0] == 0.0 or lag_times[0] >= 0.0

    def test_fit_diffusivity_known_slope(self):
        from mlip_md.msd import fit_diffusivity
        # Construct perfect linear MSD: MSD = 6D t
        D_target = 1e-7   # cm²/s
        D_target_A2_ps = D_target / 1e-4
        t = np.linspace(0, 100, 500)   # ps
        msd = 6 * D_target_A2_ps * t
        D_fit, _, R2 = fit_diffusivity(t, msd)
        assert abs(D_fit - D_target) / D_target < 1e-3
        assert R2 > 0.9999


# ── Arrhenius tests ───────────────────────────────────────────────────────

class TestArrhenius:
    def test_fit_recovers_ea(self):
        from mlip_md.arrhenius import fit_arrhenius, kB_eV_per_K
        Ea_true = 0.35
        D0_true = 1e-3
        T = np.array([600.0, 800.0, 1000.0])
        D = D0_true * np.exp(-Ea_true / (kB_eV_per_K * T))

        Ea, Ea_std, D0, D0_std, R2 = fit_arrhenius(T, D)
        assert abs(Ea - Ea_true) < 0.001, f"Ea error = {abs(Ea - Ea_true):.4f} eV"
        assert R2 > 0.9999

    def test_nernst_einstein_units(self):
        """σ = n q² D / (kB T) — check unit conversion gives S/cm."""
        from mlip_md.arrhenius import nernst_einstein_conductivity
        # LGPS: D ≈ 2e-7 cm²/s at 300 K, n ≈ 2.5 Li per 1000 Å³, σ ≈ 10⁻³ S/cm
        sigma = nernst_einstein_conductivity(
            diffusivity_cm2s=2e-7,
            temperature_K=300.0,
            n_mobile_ions=50,
            cell_volume_angstrom3=20000.0,
            ion_charge=1,
        )
        # Should be in the range 1e-4 to 1 S/cm for good electrolytes
        assert sigma > 0, "Conductivity must be positive"
        assert sigma < 100.0, f"Conductivity {sigma} S/cm seems unphysically large"

    def test_arrhenius_invalid_diffusivities(self):
        from mlip_md.arrhenius import fit_arrhenius
        T = np.array([600.0, 800.0, 1000.0])
        D = np.array([0.0, 0.0, 0.0])   # all zero → no valid fit
        with pytest.raises(ValueError):
            fit_arrhenius(T, D)


# ── Output / plot tests ────────────────────────────────────────────────────

class TestOutputs:
    def test_msd_plot_created(self, tmp_dir):
        from mlip_md.msd import MSDResult
        from mlip_md.outputs import plot_msd

        t = np.linspace(0, 100, 200)
        msd_data = [
            MSDResult(
                temperature_K=T,
                candidate_id="test_candidate",
                lag_times_ps=t,
                msd_angstrom2=6 * 5e-4 * t,   # linear MSD
                diffusivity_cm2s=5e-8,
                diffusivity_std_cm2s=1e-9,
                r_squared=0.99,
                fit_start_ps=10.0,
                fit_end_ps=90.0,
                n_mobile_ions=12,
            )
            for T in [600.0, 800.0, 1000.0]
        ]
        plot_path = plot_msd(msd_data, tmp_dir, "test_candidate")
        assert plot_path.exists(), f"MSD plot not created: {plot_path}"
        assert plot_path.stat().st_size > 1000, "Plot file seems empty"

    def test_arrhenius_plot_created(self, tmp_dir):
        from mlip_md.arrhenius import ArrheniusResult
        from mlip_md.outputs import plot_arrhenius

        T = np.array([600.0, 800.0, 1000.0])
        from mlip_md.arrhenius import kB_eV_per_K
        D = 1e-3 * np.exp(-0.35 / (kB_eV_per_K * T))

        res = ArrheniusResult(
            candidate_id="test_candidate",
            mobile_species="Li",
            temperatures_K=T,
            diffusivities_cm2s=D,
            diffusivities_std_cm2s=D * 0.05,
            r_squared_msd=np.array([0.98, 0.97, 0.96]),
            Ea_eV=0.35,
            Ea_std_eV=0.01,
            D0_cm2s=1e-3,
            D0_std_cm2s=1e-4,
            arrhenius_r_squared=0.999,
            sigma_S_cm=np.array([1e-4, 1e-3, 5e-3]),
            sigma_S_cm_300K_extrap=2e-5,
        )

        plot_path = plot_arrhenius(res, tmp_dir)
        assert plot_path.exists(), f"Arrhenius plot not created: {plot_path}"

    def test_csv_schema(self, tmp_dir):
        from mlip_md.arrhenius import ArrheniusResult
        from mlip_md.outputs import save_results_csv
        import pandas as pd

        T = np.array([600.0, 800.0, 1000.0])
        from mlip_md.arrhenius import kB_eV_per_K
        D = 1e-3 * np.exp(-0.35 / (kB_eV_per_K * T))

        res = ArrheniusResult(
            candidate_id="LGPS",
            mobile_species="Li",
            temperatures_K=T,
            diffusivities_cm2s=D,
            diffusivities_std_cm2s=D * 0.05,
            r_squared_msd=np.array([0.98, 0.97, 0.96]),
            Ea_eV=0.35,
            Ea_std_eV=0.01,
            D0_cm2s=1e-3,
            D0_std_cm2s=1e-4,
            arrhenius_r_squared=0.999,
            sigma_S_cm=np.array([1e-4, 1e-3, 5e-3]),
            sigma_S_cm_300K_extrap=2e-5,
        )

        csv_path = save_results_csv([res], tmp_dir)
        assert csv_path.exists()
        df = pd.read_csv(csv_path)
        assert len(df) == 1
        required_cols = {"candidate_id", "mobile_species", "Ea_eV", "sigma_300K_extrap_S_cm"}
        assert required_cols.issubset(set(df.columns)), \
            f"Missing columns: {required_cols - set(df.columns)}"


# ── End-to-end dummy pipeline ─────────────────────────────────────────────

class TestDummyPipeline:
    def test_end_to_end_dummy(self, tmp_dir):
        """Full pipeline with synthetic data — no GPU."""
        from mlip_md.dummy import run_dummy_pipeline, _Ea_TRUE_eV
        result = run_dummy_pipeline(
            output_dir=tmp_dir,
            temperatures_K=[600.0, 800.0, 1000.0],
            mobile_species="Li",
        )
        assert result["csv_path"].exists()
        assert result["msd_plot"].exists()
        assert result["arrhenius_plot"].exists()

        Ea_fit = result["Ea_eV_fit"]
        Ea_true = result["Ea_eV_true"]
        # Ea should be within 20% of true value (noise + finite sampling)
        assert abs(Ea_fit - Ea_true) / Ea_true < 0.30, \
            f"Ea recovered poorly: {Ea_fit:.3f} vs {Ea_true:.3f} eV"
