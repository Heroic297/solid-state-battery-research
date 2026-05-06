"""
NVT molecular dynamics runner using ASE Langevin thermostat.

Runs independent NVT trajectories at each requested temperature.
Writes ASE trajectory files (.traj) and optional XYZ checkpoints.

Key design choices
------------------
- Langevin thermostat (friction=0.01 fs⁻¹) is appropriate for diffusivity
  measurements; it does not bias free-energy differences.
- Trajectories are written at every `write_interval` steps; keeping a coarser
  interval (e.g. 10–50 fs effective frame spacing) avoids excessively large
  files while maintaining sufficient MSD sampling resolution.
- Equilibration steps use a higher friction (0.1 fs⁻¹) to accelerate
  thermalisation; production steps switch to lower friction.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
from ase import Atoms, units
from ase.io import Trajectory
from ase.io.trajectory import TrajectoryWriter
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation

logger = logging.getLogger(__name__)


@dataclass
class MDConfig:
    """
    Configuration for a single NVT MD run.

    All time quantities in femtoseconds unless noted.
    """
    temperatures: list[float] = field(default_factory=lambda: [600.0, 800.0, 1000.0])
    timestep_fs: float = 2.0
    total_time_ps: float = 200.0        # production time per temperature
    equilibration_time_ps: float = 10.0 # pre-production equilibration
    write_interval_steps: int = 10      # trajectory frame frequency
    friction_prod: float = 0.01         # fs⁻¹, production thermostat friction
    friction_equil: float = 0.10        # fs⁻¹, equilibration thermostat friction
    output_dir: Path = Path("outputs")
    logfile_interval: int = 100         # steps between log entries

    @property
    def production_steps(self) -> int:
        return int(self.total_time_ps * 1000 / self.timestep_fs)

    @property
    def equilibration_steps(self) -> int:
        return int(self.equilibration_time_ps * 1000 / self.timestep_fs)


def run_nvt_md(
    atoms: Atoms,
    calculator,
    config: MDConfig,
    candidate_id: str,
    dry_run: bool = False,
) -> dict[float, Path]:
    """
    Run NVT MD at each temperature in config.temperatures.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure (with PBC). A deep copy is made for each temperature.
    calculator : ASE-compatible calculator
    config : MDConfig
    candidate_id : str
        Label used in output filenames.
    dry_run : bool
        If True, run only 50 steps per temperature for fast pipeline testing.

    Returns
    -------
    dict mapping temperature (K) → Path of the .traj trajectory file.
    """
    import copy

    config.output_dir.mkdir(parents=True, exist_ok=True)
    traj_paths: dict[float, Path] = {}

    for T in config.temperatures:
        traj_path = _run_single_temperature(
            atoms=atoms,
            calculator=calculator,
            config=config,
            candidate_id=candidate_id,
            temperature=T,
            dry_run=dry_run,
        )
        traj_paths[T] = traj_path

    return traj_paths


def _run_single_temperature(
    atoms: Atoms,
    calculator,
    config: MDConfig,
    candidate_id: str,
    temperature: float,
    dry_run: bool = False,
) -> Path:
    """Run a single temperature NVT MD trajectory."""
    import copy

    atoms_copy = copy.deepcopy(atoms)
    atoms_copy.calc = calculator

    # Velocity initialisation
    MaxwellBoltzmannDistribution(atoms_copy, temperature_K=temperature)
    Stationary(atoms_copy)
    ZeroRotation(atoms_copy)

    traj_path = config.output_dir / f"{candidate_id}_T{int(temperature)}K.traj"
    log_path = config.output_dir / f"{candidate_id}_T{int(temperature)}K.log"

    prod_steps = 50 if dry_run else config.production_steps
    equil_steps = 20 if dry_run else config.equilibration_steps
    write_interval = max(1, config.write_interval_steps)

    logger.info(
        "MD: candidate=%s T=%.0f K  equil=%d steps  prod=%d steps  dt=%.2f fs",
        candidate_id, temperature, equil_steps, prod_steps, config.timestep_fs,
    )

    timestep_ase = config.timestep_fs * units.fs

    # ── Equilibration ─────────────────────────────────────────────────────
    dyn_equil = Langevin(
        atoms_copy,
        timestep=timestep_ase,
        temperature_K=temperature,
        friction=config.friction_equil / units.fs,
        logfile=str(log_path),
        loginterval=config.logfile_interval,
    )
    t0 = time.perf_counter()
    dyn_equil.run(equil_steps)
    logger.info("Equilibration done in %.1f s.", time.perf_counter() - t0)

    # ── Production ────────────────────────────────────────────────────────
    dyn_prod = Langevin(
        atoms_copy,
        timestep=timestep_ase,
        temperature_K=temperature,
        friction=config.friction_prod / units.fs,
        logfile=str(log_path),
        loginterval=config.logfile_interval,
    )

    traj_writer = TrajectoryWriter(str(traj_path), mode="w")

    def write_frame():
        traj_writer.write(atoms_copy)

    dyn_prod.attach(write_frame, interval=write_interval)

    t0 = time.perf_counter()
    dyn_prod.run(prod_steps)
    traj_writer.close()

    elapsed = time.perf_counter() - t0
    n_frames = prod_steps // write_interval
    logger.info(
        "Production done: %d steps in %.1f s → %s (%d frames expected)",
        prod_steps, elapsed, traj_path.name, n_frames,
    )

    return traj_path


def read_trajectory_positions(
    traj_path: Path,
    mobile_species: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    """
    Read trajectory and return unwrapped positions for mobile ions.

    Parameters
    ----------
    traj_path : Path
    mobile_species : list of str
        e.g. ["Li"] or ["Na"]

    Returns
    -------
    positions : np.ndarray, shape (n_frames, n_mobile, 3)  — Ångström
    times_ps : np.ndarray, shape (n_frames,)  — picoseconds
        Derived from frame index × timestep. Requires passing timestep separately
        if not stored in trajectory; caller should pass via the MDConfig.
    """
    from ase.io import Trajectory as TrajReader

    traj = TrajReader(str(traj_path), "r")
    frames = list(traj)
    traj.close()

    if len(frames) == 0:
        raise ValueError(f"Trajectory {traj_path} is empty.")

    # Identify mobile-ion indices from first frame
    symbols = frames[0].get_chemical_symbols()
    mobile_indices = [i for i, s in enumerate(symbols) if s in mobile_species]

    if not mobile_indices:
        raise ValueError(
            f"No mobile species {mobile_species} found in trajectory {traj_path}. "
            f"Present: {set(symbols)}"
        )

    # Extract scaled positions and convert to Cartesian for unwrapping
    n_frames = len(frames)
    n_mobile = len(mobile_indices)
    positions = np.zeros((n_frames, n_mobile, 3))

    for fi, frame in enumerate(frames):
        pos = frame.get_positions()
        positions[fi] = pos[mobile_indices]

    # Minimum-image unwrapping
    positions = _unwrap_positions(positions, frames[0].get_cell())

    return positions, np.arange(n_frames, dtype=float)


def _unwrap_positions(positions: np.ndarray, cell) -> np.ndarray:
    """
    Unwrap periodic positions using the minimum-image convention.

    positions : (n_frames, n_atoms, 3)
    cell : ASE Cell object or (3,3) array
    """
    from ase.geometry import find_mic

    cell_array = np.array(cell)
    unwrapped = positions.copy()

    for atom_idx in range(positions.shape[1]):
        for frame_idx in range(1, positions.shape[0]):
            dr = positions[frame_idx, atom_idx] - positions[frame_idx - 1, atom_idx]
            # Apply minimum image convention
            dr_mic, _ = find_mic(dr[np.newaxis], cell_array, pbc=True)
            unwrapped[frame_idx, atom_idx] = unwrapped[frame_idx - 1, atom_idx] + dr_mic[0]

    return unwrapped
