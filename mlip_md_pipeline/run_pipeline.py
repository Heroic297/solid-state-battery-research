#!/usr/bin/env python
"""
MLIP-MD Solid Electrolyte Screening Pipeline — Main Entry Point

Usage examples
--------------
# Validate pipeline (no GPU needed, ~5–30 s):
python run_pipeline.py --dummy --output-dir outputs/test

# Real screening from CSV:
python run_pipeline.py \
    --csv candidate_screening_initial.csv \
    --structures-dir structures/ \
    --calculator mace-mp-0 \
    --device cuda \
    --workers 1 \
    --output-dir outputs/

# Specific temperatures + shorter run for quick iteration:
python run_pipeline.py \
    --csv candidates.csv \
    --structures-dir structures/ \
    --calculator chgnet \
    --temperatures 600 800 1000 \
    --total-time 50 \
    --timestep 2.0 \
    --workers 2
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import click


def _configure_logging(verbose: bool):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


@click.command()
@click.option("--csv", "csv_path", default=None,
              help="Path to candidate_screening_initial.csv (or equivalent).")
@click.option("--structures-dir", "structures_dir", default="structures/",
              show_default=True,
              help="Directory containing CIF/POSCAR files referenced in CSV.")
@click.option("--calculator", "calculator_name", default="mace-mp-0",
              type=click.Choice(["mace-mp-0", "chgnet", "sevennet", "dummy"],
                                case_sensitive=False),
              show_default=True,
              help="MLIP calculator backend.")
@click.option("--device", default="cuda", show_default=True,
              help="PyTorch device string (e.g. 'cuda', 'cuda:0', 'cpu').")
@click.option("--temperatures", "-T", multiple=True, type=float,
              default=(600.0, 800.0, 1000.0), show_default=True,
              help="NVT temperatures in K (repeat flag for multiple).")
@click.option("--total-time", "total_time_ps", default=200.0,
              show_default=True, help="Production MD time per temperature (ps).")
@click.option("--equil-time", "equil_time_ps", default=10.0,
              show_default=True, help="Equilibration MD time (ps).")
@click.option("--timestep", "timestep_fs", default=2.0,
              show_default=True, help="MD timestep (fs).")
@click.option("--write-interval", "write_interval", default=10,
              show_default=True,
              help="Write trajectory frame every N MD steps.")
@click.option("--supercell-min", "supercell_min", default=10.0,
              show_default=True,
              help="Minimum supercell side length (Å); smaller cells are replicated.")
@click.option("--workers", "-n", default=1, show_default=True,
              help="Number of parallel worker processes.")
@click.option("--output-dir", "output_dir", default="outputs/",
              show_default=True, help="Root output directory.")
@click.option("--dummy", is_flag=True, default=False,
              help="Run synthetic test mode (no GPU, no real structures).")
@click.option("--verbose", "-v", is_flag=True, default=False,
              help="Enable DEBUG logging.")
def main(
    csv_path, structures_dir, calculator_name, device,
    temperatures, total_time_ps, equil_time_ps, timestep_fs,
    write_interval, supercell_min, workers, output_dir,
    dummy, verbose,
):
    _configure_logging(verbose)
    log = logging.getLogger("run_pipeline")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Dummy / test mode ───────────────────────────────────────────────
    if dummy:
        log.info("=== DUMMY / TEST MODE (no GPU required) ===")
        from mlip_md.dummy import run_dummy_pipeline
        result = run_dummy_pipeline(
            output_dir=output_dir,
            temperatures_K=list(temperatures),
            mobile_species="Li",
        )
        log.info("Dummy mode complete.")
        _print_summary([{"status": "ok",
                          "candidate_id": "dummy_candidate",
                          "arrhenius_result": result["arrhenius_result"],
                          "msd_results": result["msd_results"]}],
                        output_dir)
        return

    # ── Real screening mode ─────────────────────────────────────────────
    if csv_path is None:
        click.echo(
            "ERROR: Provide --csv <path> for real screening, or use --dummy for testing.",
            err=True,
        )
        sys.exit(1)

    log.info("Loading candidates from %s", csv_path)
    from mlip_md.orchestrator import load_candidates_from_csv, run_screening
    from mlip_md.outputs import save_results_csv, plot_msd, plot_arrhenius

    configs = load_candidates_from_csv(
        csv_path=csv_path,
        structures_dir=structures_dir,
        calculator_name=calculator_name,
        device=device,
        temperatures_K=list(temperatures),
        timestep_fs=timestep_fs,
        total_time_ps=total_time_ps,
        equilibration_time_ps=equil_time_ps,
        write_interval_steps=write_interval,
        supercell_min_length=supercell_min,
        output_dir=output_dir,
    )

    log.info("Running screening for %d candidates on %d worker(s)...", len(configs), workers)
    results = run_screening(configs, n_workers=workers, output_dir=output_dir)

    # ── Collect results and write outputs ───────────────────────────────
    ok_results = [r for r in results if r["status"] == "ok"]
    err_results = [r for r in results if r["status"] != "ok"]

    if ok_results:
        arrhenius_list = [r["arrhenius_result"] for r in ok_results]
        csv_out = save_results_csv(arrhenius_list, output_dir)
        log.info("Results CSV: %s", csv_out)

        for r in ok_results:
            cid = r["candidate_id"]
            plot_msd(r["msd_results"], output_dir / cid, cid)
            plot_arrhenius(r["arrhenius_result"], output_dir / cid)

    if err_results:
        import pandas as pd
        err_df = pd.DataFrame([
            {"candidate_id": r["candidate_id"], "error": r["error"]}
            for r in err_results
        ])
        err_path = output_dir / "failed_candidates.csv"
        err_df.to_csv(err_path, index=False)
        log.warning("%d candidates failed. See %s", len(err_results), err_path)

    _print_summary(results, output_dir)


def _print_summary(results: list, output_dir: Path):
    """Print a human-readable summary table."""
    log = logging.getLogger("run_pipeline")
    print("\n" + "="*70)
    print(f"{'Candidate':<25}  {'Eₐ (eV)':<10}  {'σ(300K) S/cm':<14}  {'Status'}")
    print("-"*70)
    for r in results:
        cid = r["candidate_id"]
        if r["status"] == "ok" and r["arrhenius_result"] is not None:
            arr = r["arrhenius_result"]
            print(f"{cid:<25}  {arr.Ea_eV:<10.4f}  {arr.sigma_S_cm_300K_extrap:<14.4e}  ok")
        else:
            err = str(r.get("error", "unknown"))[:40]
            print(f"{cid:<25}  {'N/A':<10}  {'N/A':<14}  FAILED: {err}")
    print("="*70)
    print(f"Outputs: {output_dir.resolve()}")
    print()


if __name__ == "__main__":
    main()
