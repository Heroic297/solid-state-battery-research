"""
Multiprocessing orchestrator.

Distributes independent MD runs across candidates using Python's
multiprocessing.Pool (process-based parallelism).

GPU note
--------
When running multiple MLIP models on the same GPU, each process shares
the same GPU memory.  Depending on your GPU memory (e.g. 24 GB RTX 3090),
you may need to limit the pool to 1–2 workers per GPU.
Set `n_workers=1` if you run out of VRAM.

For multi-GPU systems, set CUDA_VISIBLE_DEVICES per worker via the
`gpu_ids` argument.

Architecture
------------
  main process        → spawns Pool of `n_workers` processes
  worker process      → loads calculator once, runs MD for one candidate
                         at all temperatures sequentially

Each worker is isolated: calculator loading, trajectory writing, and
analysis are fully self-contained.  Results are returned via the pool
return value.

Error handling
--------------
  Per-candidate failures are caught and logged; the pipeline continues
  with the remaining candidates and records "error" in the CSV.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
import os
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class CandidateConfig:
    """Per-candidate configuration passed to each worker."""
    candidate_id: str
    structure_path: str
    mobile_species: list[str]          # e.g. ["Li"]
    ion_charge: int                    # 1 for Li+/Na+
    calculator_name: str               # "mace-mp-0" | "chgnet" | "sevennet" | "dummy"
    device: str                        # "cuda" | "cuda:0" | "cpu"
    temperatures_K: list[float]
    timestep_fs: float
    total_time_ps: float
    equilibration_time_ps: float
    write_interval_steps: int
    supercell_min_length: float        # Å
    output_dir: Path
    dry_run: bool = False
    gpu_id: Optional[int] = None      # set CUDA_VISIBLE_DEVICES if not None


def _worker(cfg: CandidateConfig) -> dict[str, Any]:
    """
    Worker function: full pipeline for one candidate.

    Runs in a separate process.  Returns a result dict suitable for
    aggregation in the main process.
    """
    # Configure CUDA device isolation
    if cfg.gpu_id is not None:
        os.environ["CUDA_VISIBLE_DEVICES"] = str(cfg.gpu_id)

    # Set up per-worker logging to avoid interleaving
    logging.basicConfig(
        level=logging.INFO,
        format=f"[%(processName)s|{cfg.candidate_id}] %(levelname)s %(message)s",
    )
    log = logging.getLogger(__name__)

    try:
        from mlip_md.ingestion import load_structure, validate_mobile_ions, build_supercell
        from mlip_md.calculators import get_calculator
        from mlip_md.md_runner import MDConfig, run_nvt_md
        from mlip_md.msd import analyse_trajectory
        from mlip_md.arrhenius import full_arrhenius_analysis

        log.info("Worker started for '%s'", cfg.candidate_id)

        # ── 1. Load structure ──────────────────────────────────────────
        atoms = load_structure(cfg.structure_path)
        validate_mobile_ions(atoms, cfg.mobile_species)
        atoms = build_supercell(atoms, min_length_angstrom=cfg.supercell_min_length)
        cell_vol = atoms.get_volume()

        # ── 2. Initialise calculator ───────────────────────────────────
        calc = get_calculator(cfg.calculator_name, device=cfg.device)
        atoms.calc = calc

        # ── 3. Run NVT MD ──────────────────────────────────────────────
        md_cfg = MDConfig(
            temperatures=cfg.temperatures_K,
            timestep_fs=cfg.timestep_fs,
            total_time_ps=cfg.total_time_ps,
            equilibration_time_ps=cfg.equilibration_time_ps,
            write_interval_steps=cfg.write_interval_steps,
            output_dir=cfg.output_dir / cfg.candidate_id,
        )
        traj_paths = run_nvt_md(
            atoms=atoms,
            calculator=calc,
            config=md_cfg,
            candidate_id=cfg.candidate_id,
            dry_run=cfg.dry_run,
        )

        # ── 4. MSD analysis ────────────────────────────────────────────
        msd_results = []
        for T, traj_path in traj_paths.items():
            msd_res = analyse_trajectory(
                traj_path=traj_path,
                temperature_K=T,
                candidate_id=cfg.candidate_id,
                mobile_species=cfg.mobile_species,
                timestep_fs=cfg.timestep_fs,
                write_interval=cfg.write_interval_steps,
            )
            msd_results.append(msd_res)

        # ── 5. Arrhenius fit ───────────────────────────────────────────
        mobile_species_str = "+".join(cfg.mobile_species)
        arrhenius_res = full_arrhenius_analysis(
            candidate_id=cfg.candidate_id,
            mobile_species=mobile_species_str,
            msd_results=msd_results,
            cell_volume_angstrom3=cell_vol,
            ion_charge=cfg.ion_charge,
        )

        log.info(
            "'%s' done: Ea=%.4f eV  σ(300K)=%.4e S/cm",
            cfg.candidate_id, arrhenius_res.Ea_eV,
            arrhenius_res.sigma_S_cm_300K_extrap,
        )

        return {
            "status": "ok",
            "candidate_id": cfg.candidate_id,
            "arrhenius_result": arrhenius_res,
            "msd_results": msd_results,
        }

    except Exception as exc:
        tb = traceback.format_exc()
        log.error("FAILED for '%s': %s\n%s", cfg.candidate_id, exc, tb)
        return {
            "status": "error",
            "candidate_id": cfg.candidate_id,
            "error": str(exc),
            "traceback": tb,
            "arrhenius_result": None,
            "msd_results": [],
        }


def run_screening(
    candidate_configs: list[CandidateConfig],
    n_workers: int = 1,
    output_dir: Path = Path("outputs"),
) -> list[dict[str, Any]]:
    """
    Run the full screening pipeline across all candidates.

    Parameters
    ----------
    candidate_configs : list of CandidateConfig
    n_workers : int
        Number of parallel worker processes.  Set to 1 per GPU to avoid
        VRAM exhaustion.  Use 2+ only if you have multiple GPUs.
    output_dir : Path
        Root output directory; per-candidate subdirs are created automatically.

    Returns
    -------
    list of result dicts (one per candidate, including failures).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_candidates = len(candidate_configs)
    logger.info(
        "Starting screening: %d candidates, %d workers", n_candidates, n_workers
    )

    if n_workers == 1:
        # Serial execution — easier to debug
        results = [_worker(cfg) for cfg in candidate_configs]
    else:
        ctx = mp.get_context("spawn")   # spawn is safer with CUDA/PyTorch
        with ctx.Pool(processes=n_workers) as pool:
            results = pool.map(_worker, candidate_configs)

    n_ok = sum(1 for r in results if r["status"] == "ok")
    n_err = sum(1 for r in results if r["status"] != "ok")
    logger.info(
        "Screening complete: %d/%d candidates succeeded, %d failed.",
        n_ok, n_candidates, n_err,
    )

    return results


def load_candidates_from_csv(
    csv_path: str | Path,
    structures_dir: str | Path,
    calculator_name: str = "mace-mp-0",
    device: str = "cuda",
    temperatures_K: list[float] | None = None,
    timestep_fs: float = 2.0,
    total_time_ps: float = 200.0,
    equilibration_time_ps: float = 10.0,
    write_interval_steps: int = 10,
    supercell_min_length: float = 10.0,
    output_dir: Path = Path("outputs"),
    dry_run: bool = False,
) -> list[CandidateConfig]:
    """
    Build a list of CandidateConfig objects from the screening CSV.

    Expected CSV columns (minimum required)
    ----------------------------------------
    candidate_id     : str  — unique label
    structure_file   : str  — filename relative to structures_dir
    mobile_species   : str  — comma-separated, e.g. "Li" or "Li,Na"
    ion_charge       : int  — 1 for monovalent (optional, default 1)

    Optional columns that override global defaults
    ----------------------------------------------
    total_time_ps, equilibration_time_ps, timestep_fs,
    temperatures_K   : semicolon-separated, e.g. "600;800;1000"
    calculator       : override global calculator per candidate
    device           : override GPU device per candidate
    """
    if temperatures_K is None:
        temperatures_K = [600.0, 800.0, 1000.0]

    df = pd.read_csv(csv_path)
    logger.info("Loaded %d candidates from %s", len(df), csv_path)

    required = {"candidate_id", "structure_file", "mobile_species"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"CSV {csv_path} is missing required columns: {missing}. "
            f"Present columns: {list(df.columns)}"
        )

    configs = []
    for _, row in df.iterrows():
        cid = str(row["candidate_id"]).strip()
        struct_file = str(row["structure_file"]).strip()
        struct_path = Path(structures_dir) / struct_file

        mobile_str = str(row.get("mobile_species", "Li")).strip()
        mobile = [s.strip() for s in mobile_str.replace(";", ",").split(",")]

        charge = int(row.get("ion_charge", 1))

        # Per-candidate overrides
        calc_name = str(row.get("calculator", calculator_name)).strip()
        dev = str(row.get("device", device)).strip()
        dt = float(row.get("timestep_fs", timestep_fs))
        total_t = float(row.get("total_time_ps", total_time_ps))
        equil_t = float(row.get("equilibration_time_ps", equilibration_time_ps))

        if "temperatures_K" in df.columns and pd.notna(row["temperatures_K"]):
            T_list = [
                float(t.strip())
                for t in str(row["temperatures_K"]).split(";")
            ]
        else:
            T_list = temperatures_K

        configs.append(CandidateConfig(
            candidate_id=cid,
            structure_path=str(struct_path),
            mobile_species=mobile,
            ion_charge=charge,
            calculator_name=calc_name,
            device=dev,
            temperatures_K=T_list,
            timestep_fs=dt,
            total_time_ps=total_t,
            equilibration_time_ps=equil_t,
            write_interval_steps=write_interval_steps,
            supercell_min_length=supercell_min_length,
            output_dir=output_dir,
            dry_run=dry_run,
        ))

    return configs
