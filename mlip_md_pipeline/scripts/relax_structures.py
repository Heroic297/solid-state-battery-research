#!/usr/bin/env python
"""
Relax generated structures with CHGNet (default) or MACE-MP-0 (optional).

Reads one or more CIF/POSCAR files (or a directory glob), runs ASE-based
geometry optimization, and writes relaxed CIF/POSCAR plus a CSV row of
energies, max forces, and step counts.

Usage
-----
    # CHGNet, single structure
    python scripts/relax_structures.py \\
        --inputs structures/Li3Er0.5Sc0.5Cl6/*.vasp \\
        --calculator chgnet \\
        --device cpu \\
        --output-dir structures/Li3Er0.5Sc0.5Cl6/relaxed/

    # MACE-MP-0, multiple files (only if mace-torch is installed)
    python scripts/relax_structures.py \\
        --inputs structures/Li3In0.5Er0.5Cl6/*.cif \\
        --calculator mace-mp-0 \\
        --device cuda \\
        --output-dir structures/Li3In0.5Er0.5Cl6/relaxed/

Honesty notes
-------------
- CHGNet and MACE-MP-0 are universal MLIPs trained primarily on oxides; their
  accuracy on lithium halides (Li3MCl6) is variable. Treat relative energies
  as provisional. Confirm low-energy orderings with DFT (e.g. VASP/QE) before
  publishing or planning synthesis.
- If your dopant element is not in the model's training distribution (e.g.
  some lanthanides for some MLIPs), expect larger errors.
- This script does NOT enforce convergence on energies between calculators —
  consistent ordering across CHGNet and MACE is a useful sanity check.
"""

from __future__ import annotations

import csv
import glob
import logging
import sys
import time
from pathlib import Path

import click

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Calculator loading
# ---------------------------------------------------------------------------

def load_calculator(name: str, device: str):
    """Return an ASE Calculator for the requested backend.

    Reuses the project's mlip_md.calculators.get_calculator if available,
    otherwise falls back to local imports.
    """
    name = name.lower().strip()
    try:
        # Reuse the package-level adapter if importable from the script context.
        from mlip_md.calculators import get_calculator
        return get_calculator(name, device=device)
    except Exception as e:
        logger.warning("Falling back to local calculator import (mlip_md not on path): %s", e)

    if name == "chgnet":
        try:
            from chgnet.model.dynamics import CHGNetCalculator
            from chgnet.model import CHGNet
        except ImportError as imp_err:
            raise SystemExit(
                "chgnet not installed. On Miniforge:\n"
                "    pip install chgnet\n"
                f"Original error: {imp_err}"
            )
        return CHGNetCalculator(model=CHGNet.load(), use_device=device)

    if name in ("mace", "mace-mp-0", "mace_mp_0"):
        try:
            from mace.calculators import mace_mp
        except ImportError as imp_err:
            raise SystemExit(
                "mace-torch not installed. On Miniforge:\n"
                "    pip install mace-torch\n"
                f"Original error: {imp_err}"
            )
        return mace_mp(model="medium", default_dtype="float32", device=device)

    raise SystemExit(f"Unknown calculator '{name}'. Choose: chgnet, mace-mp-0.")


# ---------------------------------------------------------------------------
# Relaxation
# ---------------------------------------------------------------------------

def relax_one(
    input_path: Path,
    output_dir: Path,
    calculator,
    fmax: float,
    max_steps: int,
    optimizer_name: str,
    relax_cell: bool,
) -> dict:
    """Relax a single structure file. Returns a dict of results."""
    from ase.io import read, write

    atoms = read(str(input_path))
    atoms.calc = calculator

    target = atoms
    if relax_cell:
        from ase.constraints import ExpCellFilter
        target = ExpCellFilter(atoms)

    if optimizer_name.upper() == "FIRE":
        from ase.optimize import FIRE as Opt
    else:
        from ase.optimize import BFGS as Opt

    log_path = output_dir / f"{input_path.stem}_relax.log"
    output_dir.mkdir(parents=True, exist_ok=True)
    opt = Opt(target, logfile=str(log_path))

    t0 = time.time()
    try:
        opt.run(fmax=fmax, steps=max_steps)
        converged = bool(opt.converged())
        error = ""
    except Exception as exc:  # noqa: BLE001
        logger.error("Relaxation failed for %s: %s", input_path, exc)
        converged = False
        error = str(exc)
    elapsed = time.time() - t0

    # Force/energy metrics
    try:
        e_final = float(atoms.get_potential_energy())
    except Exception:  # noqa: BLE001
        e_final = float("nan")
    try:
        forces = atoms.get_forces()
        fmax_final = float(((forces ** 2).sum(axis=1) ** 0.5).max())
    except Exception:  # noqa: BLE001
        fmax_final = float("nan")

    # Write outputs
    out_cif = output_dir / f"{input_path.stem}_relaxed.cif"
    out_vasp = output_dir / f"{input_path.stem}_relaxed.vasp"
    try:
        write(str(out_cif), atoms)
        write(str(out_vasp), atoms, format="vasp")
    except Exception as exc:  # noqa: BLE001
        logger.error("Could not write relaxed outputs for %s: %s", input_path, exc)

    return {
        "input": str(input_path),
        "relaxed_cif": str(out_cif),
        "relaxed_poscar": str(out_vasp),
        "log": str(log_path),
        "converged": converged,
        "n_steps": int(getattr(opt, "nsteps", -1)),
        "fmax_target": fmax,
        "fmax_final": fmax_final,
        "energy_eV": e_final,
        "n_atoms": len(atoms),
        "elapsed_seconds": round(elapsed, 3),
        "relax_cell": relax_cell,
        "optimizer": optimizer_name.upper(),
        "error": error,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

@click.command(context_settings={"show_default": True})
@click.option("--inputs", "input_globs", multiple=True, required=True,
              help="One or more file paths or globs (e.g. 'structures/*/*.cif'). "
                   "Repeat the flag to provide multiple patterns.")
@click.option("--calculator", "calc_name", type=click.Choice(["chgnet", "mace-mp-0"]),
              default="chgnet", help="MLIP backend.")
@click.option("--device", default="cpu",
              help="Torch device: cpu, cuda, cuda:0, ...")
@click.option("--fmax", type=float, default=0.05,
              help="Target max force [eV/Å] for convergence.")
@click.option("--max-steps", type=int, default=200,
              help="Maximum optimizer steps.")
@click.option("--optimizer", "optimizer_name",
              type=click.Choice(["BFGS", "FIRE"], case_sensitive=False),
              default="FIRE", help="ASE optimizer.")
@click.option("--relax-cell/--no-relax-cell", default=True,
              help="Allow cell parameters to relax (variable-cell relaxation).")
@click.option("--output-dir", type=click.Path(file_okay=False), required=True,
              help="Directory for relaxed outputs and summary CSV.")
def main(input_globs, calc_name, device, fmax, max_steps, optimizer_name, relax_cell, output_dir):
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Expand globs
    input_files: list[Path] = []
    for pattern in input_globs:
        matches = [Path(p) for p in glob.glob(pattern)]
        if not matches and Path(pattern).exists():
            matches = [Path(pattern)]
        input_files.extend(matches)

    input_files = [p for p in input_files if p.is_file()]
    if not input_files:
        logger.error("No input files matched. Patterns: %s", list(input_globs))
        sys.exit(2)

    logger.info("Loading calculator '%s' on device='%s'.", calc_name, device)
    calc = load_calculator(calc_name, device)

    logger.warning(
        "Heads-up: %s on lithium halides may need DFT confirmation. "
        "Use these relative energies as a screen, not as ground truth.",
        calc_name,
    )

    summary_rows: list[dict] = []
    for idx, fpath in enumerate(input_files, start=1):
        logger.info("[%d/%d] Relaxing %s", idx, len(input_files), fpath)
        try:
            row = relax_one(
                input_path=fpath,
                output_dir=output_dir,
                calculator=calc,
                fmax=fmax,
                max_steps=max_steps,
                optimizer_name=optimizer_name,
                relax_cell=relax_cell,
            )
        except Exception as exc:  # noqa: BLE001
            logger.exception("Unhandled error on %s: %s", fpath, exc)
            row = {
                "input": str(fpath),
                "relaxed_cif": "",
                "relaxed_poscar": "",
                "log": "",
                "converged": False,
                "n_steps": -1,
                "fmax_target": fmax,
                "fmax_final": float("nan"),
                "energy_eV": float("nan"),
                "n_atoms": -1,
                "elapsed_seconds": 0.0,
                "relax_cell": relax_cell,
                "optimizer": optimizer_name.upper(),
                "error": str(exc),
            }
        summary_rows.append(row)
        logger.info(
            "  -> converged=%s, steps=%d, E=%.4f eV, fmax=%.4f",
            row["converged"], row["n_steps"], row["energy_eV"], row["fmax_final"],
        )

    # Write summary CSV
    summary_csv = output_dir / "relaxation_summary.csv"
    fieldnames = [
        "input", "relaxed_cif", "relaxed_poscar", "log",
        "converged", "n_steps", "fmax_target", "fmax_final",
        "energy_eV", "n_atoms", "elapsed_seconds",
        "relax_cell", "optimizer", "error",
    ]
    with open(summary_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    logger.info("Summary written to %s", summary_csv)
    logger.info("Done. Relaxed outputs in %s", output_dir.resolve())


if __name__ == "__main__":
    main()
