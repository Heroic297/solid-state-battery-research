#!/usr/bin/env python
"""
Generate mixed-cation Li3MCl6 supercells from a parent CIF/POSCAR.

Designed for use cases such as:
    Li3Er0.5Sc0.5Cl6   (parent: Li3ErCl6 or Li3ScCl6)
    Li3In0.5Er0.5Cl6   (parent: Li3InCl6 or Li3ErCl6)

Approach
--------
1. Read parent structure (CIF or POSCAR) with pymatgen.
2. Identify the M cation sublattice as all sites whose species are NEITHER Li
   NOR an anion (Cl by default; configurable via --anion).
3. Build a supercell with the requested integer scaling.
4. Substitute a user-specified fraction of the M sites with the dopant element,
   producing a discrete, integer-count occupancy in the supercell.
5. For each requested random seed, generate a distinct random ordering and
   write CIF + POSCAR outputs.

IMPORTANT — honesty notes
-------------------------
- These orderings are *random* on the M sublattice. They do not represent
  ground-state orderings. Before drawing physical conclusions you must:
    * relax each ordering (see relax_structures.py),
    * generate enough orderings to estimate ordering-energy spread (e.g. 5-20),
    * preferably confirm the lowest-energy ordering with DFT.
- Universal MLIPs (CHGNet, MACE-MP-0) for halides such as Li3MCl6 are not
  always reliable; treat every relative energy as provisional until checked
  against DFT or experiment.

Usage
-----
    python scripts/generate_mixed_cation_supercell.py \\
        --parent structures/Li3ErCl6.cif \\
        --host-element Er \\
        --dopant Sc \\
        --dopant-fraction 0.5 \\
        --supercell 2 2 2 \\
        --seeds 0 1 2 3 4 \\
        --output-dir structures/Li3Er0.5Sc0.5Cl6/

If --host-element is omitted, the M sublattice is auto-detected as all sites
whose element is neither Li nor the configured anion (default: Cl).
"""

from __future__ import annotations

import json
import logging
import random
import sys
from pathlib import Path

import click

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Parent loading
# ---------------------------------------------------------------------------

def load_parent_structure(path: Path):
    """Load a CIF or POSCAR using pymatgen, with ASE fallback for unusual formats."""
    from pymatgen.core import Structure

    suffix = path.suffix.lower()
    try:
        if suffix in (".cif",):
            return Structure.from_file(str(path))
        if suffix in (".vasp", ".poscar") or path.name.upper() in ("POSCAR", "CONTCAR"):
            return Structure.from_file(str(path))
        # Generic try
        return Structure.from_file(str(path))
    except Exception as e:
        logger.warning("pymatgen.Structure.from_file failed (%s); trying ASE.", e)
        from ase.io import read
        from pymatgen.io.ase import AseAtomsAdaptor
        atoms = read(str(path))
        return AseAtomsAdaptor.get_structure(atoms)


# ---------------------------------------------------------------------------
# M-site identification
# ---------------------------------------------------------------------------

def identify_m_sites(
    structure,
    host_element: str | None,
    li_symbol: str = "Li",
    anion_symbol: str = "Cl",
) -> list[int]:
    """
    Return indices of M cation sites in the parent structure.

    If host_element is given, M sites are those whose dominant element matches.
    Otherwise, M sites are all sites whose element is neither Li nor anion.
    """
    indices: list[int] = []
    for i, site in enumerate(structure):
        # site.species is a Composition; we look at the element of the majority
        symbols = [str(sp.symbol) for sp in site.species.elements]
        if host_element is not None:
            if host_element in symbols:
                indices.append(i)
        else:
            if li_symbol not in symbols and anion_symbol not in symbols:
                indices.append(i)
    return indices


# ---------------------------------------------------------------------------
# Supercell + substitution
# ---------------------------------------------------------------------------

def build_supercell_and_substitute(
    parent,
    supercell_matrix: tuple[int, int, int],
    host_element: str,
    dopant_element: str,
    dopant_fraction: float,
    seed: int,
    li_symbol: str = "Li",
    anion_symbol: str = "Cl",
):
    """
    Return (substituted_structure, swap_indices, n_total_M, n_dopant) for a given seed.

    The substitution is performed on the supercell, on sites whose element is
    `host_element`. If `host_element` is auto (None), the most-common
    non-Li/non-anion species in the supercell is used as the host.
    """
    from pymatgen.core import Element

    sx, sy, sz = supercell_matrix
    super_struct = parent.copy()
    super_struct.make_supercell([sx, sy, sz])

    # Re-identify M sites in the supercell (by host element if given,
    # else by exclusion of Li and anion).
    if host_element is None:
        # auto-detect host as the most common cation that is not Li and not anion
        counts: dict[str, int] = {}
        for site in super_struct:
            for sp in site.species.elements:
                s = str(sp.symbol)
                if s == li_symbol or s == anion_symbol:
                    continue
                counts[s] = counts.get(s, 0) + 1
        if not counts:
            raise RuntimeError(
                f"Could not auto-detect M cation: every site is {li_symbol} or {anion_symbol}."
            )
        host_element = max(counts, key=counts.get)
        logger.info("Auto-detected host M element: %s (count=%d)", host_element, counts[host_element])

    m_indices = identify_m_sites(
        super_struct,
        host_element=host_element,
        li_symbol=li_symbol,
        anion_symbol=anion_symbol,
    )
    n_total = len(m_indices)
    if n_total == 0:
        raise RuntimeError(
            f"No host {host_element} sites found in the supercell. "
            "Check --host-element and --anion."
        )

    n_dopant = int(round(dopant_fraction * n_total))
    if n_dopant <= 0 or n_dopant >= n_total:
        logger.warning(
            "dopant fraction %.3f * %d sites => %d swap(s); "
            "this may not give a true mixed-cation cell. "
            "Increase --supercell or adjust --dopant-fraction.",
            dopant_fraction, n_total, n_dopant,
        )

    rng = random.Random(seed)
    swap_indices = rng.sample(m_indices, n_dopant)

    # Apply substitution. Use replace() per index so we keep a clean Element.
    dopant = Element(dopant_element)
    for idx in swap_indices:
        super_struct.replace(idx, dopant)

    return super_struct, swap_indices, n_total, n_dopant, host_element


# ---------------------------------------------------------------------------
# Output writing
# ---------------------------------------------------------------------------

def write_outputs(
    structure,
    out_dir: Path,
    base_name: str,
    seed: int,
    metadata: dict,
) -> dict[str, str]:
    """Write CIF + POSCAR + sidecar JSON. Return paths."""
    out_dir.mkdir(parents=True, exist_ok=True)

    cif_path = out_dir / f"{base_name}_seed{seed}.cif"
    poscar_path = out_dir / f"{base_name}_seed{seed}.vasp"
    meta_path = out_dir / f"{base_name}_seed{seed}.json"

    structure.to(fmt="cif", filename=str(cif_path))
    # pymatgen writes POSCAR via fmt="poscar"
    structure.to(fmt="poscar", filename=str(poscar_path))

    with open(meta_path, "w") as f:
        json.dump(metadata, f, indent=2, sort_keys=True)

    return {
        "cif": str(cif_path),
        "poscar": str(poscar_path),
        "metadata": str(meta_path),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

@click.command(context_settings={"show_default": True})
@click.option("--parent", "parent_path", type=click.Path(exists=True, dir_okay=False),
              required=True, help="Path to parent CIF or POSCAR.")
@click.option("--host-element", default=None,
              help="Element symbol of the host M cation (e.g. 'Er'). "
                   "If omitted, auto-detect as most-common non-Li/non-anion species.")
@click.option("--dopant", "dopant_element", required=True,
              help="Element symbol of the substituted (dopant) cation, e.g. 'Sc' or 'In'.")
@click.option("--dopant-fraction", type=float, default=0.5,
              help="Fraction of M sites to be replaced by dopant (0 < f < 1).")
@click.option("--supercell", nargs=3, type=int, default=(2, 2, 2),
              help="Integer supercell scaling (a b c).")
@click.option("--seeds", multiple=True, type=int, default=(0,),
              help="One or more random seeds (repeat the flag, or pass space-separated).")
@click.option("--anion", "anion_symbol", default="Cl",
              help="Anion species symbol used to exclude anion sites from M detection.")
@click.option("--li-symbol", default="Li",
              help="Mobile-ion symbol used to exclude Li sites from M detection.")
@click.option("--output-dir", type=click.Path(file_okay=False), required=True,
              help="Directory to write CIF/POSCAR/JSON outputs.")
@click.option("--name", "base_name", default=None,
              help="Optional base filename (default: derived from composition).")
def main(
    parent_path,
    host_element,
    dopant_element,
    dopant_fraction,
    supercell,
    seeds,
    anion_symbol,
    li_symbol,
    output_dir,
    base_name,
):
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")

    if not (0.0 < dopant_fraction < 1.0):
        logger.error("--dopant-fraction must be strictly between 0 and 1.")
        sys.exit(2)

    parent_path = Path(parent_path)
    output_dir = Path(output_dir)

    logger.info("Loading parent structure: %s", parent_path)
    parent = load_parent_structure(parent_path)
    logger.info("Parent formula: %s, sites: %d", parent.composition.reduced_formula, len(parent))

    if base_name is None:
        # Derive a base name from the eventual composition
        # e.g. Li3Er0.5Sc0.5Cl6
        host = host_element or "M"
        base_name = (
            f"Li3{host}{1.0 - dopant_fraction:g}"
            f"{dopant_element}{dopant_fraction:g}{anion_symbol}6"
        )
        # filesystem-safe (Windows): keep only alnum and dot
        base_name = "".join(c if (c.isalnum() or c in (".", "_", "-")) else "_" for c in base_name)

    seeds = list(seeds) if seeds else [0]
    logger.info("Will generate %d ordering(s) with seeds=%s", len(seeds), seeds)

    summary = []
    for seed in seeds:
        struct, swap_indices, n_total, n_dopant, used_host = build_supercell_and_substitute(
            parent=parent,
            supercell_matrix=tuple(supercell),
            host_element=host_element,
            dopant_element=dopant_element,
            dopant_fraction=dopant_fraction,
            seed=seed,
            li_symbol=li_symbol,
            anion_symbol=anion_symbol,
        )
        meta = {
            "parent_path": str(parent_path),
            "parent_formula": parent.composition.reduced_formula,
            "host_element": used_host,
            "dopant_element": dopant_element,
            "dopant_fraction_requested": dopant_fraction,
            "supercell": list(supercell),
            "seed": seed,
            "n_M_sites_total": n_total,
            "n_M_sites_replaced": n_dopant,
            "swap_indices_in_supercell": list(swap_indices),
            "anion_symbol": anion_symbol,
            "li_symbol": li_symbol,
            "warning": (
                "Random ordering on M sublattice; not a ground-state ordering. "
                "Validate energies with multiple seeds and DFT before drawing conclusions."
            ),
            "supercell_formula": struct.composition.formula,
        }
        paths = write_outputs(struct, output_dir, base_name, seed, meta)
        summary.append({"seed": seed, **paths, **meta})
        logger.info(
            "seed=%d  M=%s, %d/%d replaced -> %s",
            seed, used_host, n_dopant, n_total, paths["cif"],
        )

    # Write run-level summary
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "generation_summary.json", "w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)

    logger.info("Wrote %d ordering(s) to %s", len(seeds), out_dir.resolve())
    logger.warning(
        "Reminder: random orderings need validation. "
        "Relax each ordering and inspect the energy spread before claiming any single ordering is meaningful."
    )


if __name__ == "__main__":
    main()
