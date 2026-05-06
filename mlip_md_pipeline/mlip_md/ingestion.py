"""
Structure ingestion: CIF or POSCAR → ASE Atoms.

Supported formats (auto-detected by extension):
  .cif          → ASE read with format='cif'
  POSCAR/CONTCAR → ASE read with format='vasp'
  .extxyz       → ASE read with format='extxyz'

The atoms object is returned with periodic boundary conditions enforced.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

import ase.io
from ase import Atoms

logger = logging.getLogger(__name__)

_EXT_FORMAT_MAP: dict[str, str] = {
    ".cif": "cif",
    ".poscar": "vasp",
    ".vasp": "vasp",
    ".extxyz": "extxyz",
    ".xyz": "extxyz",
}


def load_structure(path: str | Path) -> Atoms:
    """
    Load a crystal structure from CIF, POSCAR, or extxyz file.

    Parameters
    ----------
    path : str or Path
        Path to the structure file.

    Returns
    -------
    ase.Atoms
        Atoms object with pbc=True enforced.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the format cannot be inferred.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Structure file not found: {path}")

    fmt = _infer_format(path)
    logger.info("Loading structure from %s (format=%s)", path, fmt)

    atoms = ase.io.read(str(path), format=fmt, index=0)
    atoms.pbc = True

    n_atoms = len(atoms)
    symbols = set(atoms.get_chemical_symbols())
    logger.info("Loaded %d atoms, species: %s", n_atoms, symbols)

    return atoms


def _infer_format(path: Path) -> str:
    """Infer ASE format string from file name."""
    name_lower = path.name.lower()

    # POSCAR / CONTCAR have no extension by convention
    if name_lower in ("poscar", "contcar") or name_lower.startswith("poscar"):
        return "vasp"
    if name_lower.startswith("contcar"):
        return "vasp"

    suffix = path.suffix.lower()
    if suffix in _EXT_FORMAT_MAP:
        return _EXT_FORMAT_MAP[suffix]

    # Fallback: let ASE guess
    logger.warning(
        "Cannot infer format for '%s'; letting ASE auto-detect.", path.name
    )
    return None  # ase.io.read accepts format=None for auto-detect


def validate_mobile_ions(atoms: Atoms, mobile_species: list[str]) -> int:
    """
    Check that the structure contains at least one mobile ion.

    Parameters
    ----------
    atoms : ase.Atoms
    mobile_species : list of str
        e.g. ["Li"] or ["Na"]

    Returns
    -------
    int
        Count of mobile ions found.

    Raises
    ------
    ValueError
        If no mobile ions are present.
    """
    symbols = atoms.get_chemical_symbols()
    count = sum(1 for s in symbols if s in mobile_species)
    if count == 0:
        raise ValueError(
            f"Structure contains none of the mobile species {mobile_species}. "
            f"Present species: {set(symbols)}"
        )
    logger.info(
        "Found %d mobile ions (%s) in structure.", count, mobile_species
    )
    return count


def build_supercell(atoms: Atoms, min_length_angstrom: float = 10.0) -> Atoms:
    """
    Replicate atoms to reach a minimum simulation-box side length.

    This ensures MSD statistics are not corrupted by PBC self-interaction
    within a single diffusion event.

    Parameters
    ----------
    atoms : ase.Atoms
    min_length_angstrom : float
        Target minimum lattice vector length in Angstrom.

    Returns
    -------
    ase.Atoms
        Possibly-replicated supercell.
    """
    import numpy as np
    from ase.build import make_supercell

    cell_lengths = np.linalg.norm(atoms.cell, axis=1)
    reps = [max(1, int(np.ceil(min_length_angstrom / L))) for L in cell_lengths]

    if reps == [1, 1, 1]:
        logger.info("Supercell not needed (cell already >= %.1f Å).", min_length_angstrom)
        return atoms

    logger.info("Building %dx%dx%d supercell.", *reps)
    P = np.diag(reps)
    supercell = make_supercell(atoms, P)
    supercell.pbc = True
    logger.info(
        "Supercell: %d atoms (was %d).", len(supercell), len(atoms)
    )
    return supercell
