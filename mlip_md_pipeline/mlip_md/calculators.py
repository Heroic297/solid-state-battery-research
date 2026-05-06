"""
MLIP calculator adapters.

Supported backends:
  mace-mp-0  — MACE-MP-0 universal potential (Batatia et al., 2023)
  chgnet     — CHGNet (Deng et al., 2023)
  sevennet   — SevenNet-0 (Park et al., 2024)
  dummy      — Lennard-Jones placeholder for CI / no-GPU testing

Usage
-----
>>> from mlip_md.calculators import get_calculator
>>> calc = get_calculator("mace-mp-0", device="cuda")
>>> atoms.calc = calc
"""

from __future__ import annotations

import logging
from typing import Literal

from ase import Atoms
from ase.calculators.calculator import Calculator

logger = logging.getLogger(__name__)

CalcName = Literal["mace-mp-0", "chgnet", "sevennet", "dummy"]


def get_calculator(
    name: CalcName,
    device: str = "cuda",
    dtype: str = "float32",
    **kwargs,
) -> Calculator:
    """
    Instantiate and return an ASE-compatible calculator.

    Parameters
    ----------
    name : str
        One of 'mace-mp-0', 'chgnet', 'sevennet', 'dummy'.
    device : str
        PyTorch device string, e.g. 'cuda', 'cuda:0', 'cpu'.
    dtype : str
        Floating-point precision: 'float32' or 'float64'.
    **kwargs
        Forwarded to the underlying calculator constructor.

    Returns
    -------
    ase.calculators.calculator.Calculator
    """
    name = name.lower().strip()
    logger.info("Initialising calculator '%s' on device='%s'", name, device)

    if name == "mace-mp-0":
        return _mace_mp0(device=device, dtype=dtype, **kwargs)
    elif name == "chgnet":
        return _chgnet(device=device, **kwargs)
    elif name in ("sevennet", "7net"):
        return _sevennet(device=device, **kwargs)
    elif name == "dummy":
        return _dummy_lj(**kwargs)
    else:
        raise ValueError(
            f"Unknown calculator '{name}'. "
            "Choose from: mace-mp-0, chgnet, sevennet, dummy."
        )


# ---------------------------------------------------------------------------
# MACE-MP-0
# ---------------------------------------------------------------------------

def _mace_mp0(device: str = "cuda", dtype: str = "float32", **kwargs) -> Calculator:
    """
    Load MACE-MP-0 from the mace-torch package.

    Models are downloaded automatically on first use to
    ~/.cache/mace/  (overridable via MACE_MODEL_PATH env var).

    Reference: Batatia et al., Science 2024.
    """
    try:
        from mace.calculators import mace_mp
    except ImportError as e:
        raise ImportError(
            "mace-torch is not installed. "
            "Run: pip install mace-torch\n"
            f"Original error: {e}"
        )

    float_dtype_map = {"float32": "float32", "float64": "float64"}
    dtype_str = float_dtype_map.get(dtype, "float32")

    calc = mace_mp(
        model="medium",          # 'small', 'medium', or 'large'
        dispersion=False,
        default_dtype=dtype_str,
        device=device,
        **kwargs,
    )
    logger.info("MACE-MP-0 (medium) loaded.")
    return calc


# ---------------------------------------------------------------------------
# CHGNet
# ---------------------------------------------------------------------------

def _chgnet(device: str = "cuda", **kwargs) -> Calculator:
    """
    Load CHGNet universal potential.

    Reference: Deng et al., Nature Machine Intelligence 2023.
    """
    try:
        from chgnet.model.dynamics import CHGNetCalculator
        from chgnet.model import CHGNet
    except ImportError as e:
        raise ImportError(
            "chgnet is not installed. "
            "Run: pip install chgnet\n"
            f"Original error: {e}"
        )

    model = CHGNet.load()
    calc = CHGNetCalculator(model=model, use_device=device, **kwargs)
    logger.info("CHGNet loaded (device=%s).", device)
    return calc


# ---------------------------------------------------------------------------
# SevenNet
# ---------------------------------------------------------------------------

def _sevennet(device: str = "cuda", **kwargs) -> Calculator:
    """
    Load SevenNet-0 (11July2023) universal potential.

    Reference: Park et al., J. Chem. Theory Comput. 2024.

    Install: pip install sevenn   (or from GitHub: https://github.com/MDIL-SNU/SevenNet)
    """
    try:
        from sevenn.sevennet_calculator import SevenNetCalculator
    except ImportError as e:
        raise ImportError(
            "sevenn is not installed. "
            "Run: pip install sevenn  (or install from https://github.com/MDIL-SNU/SevenNet)\n"
            f"Original error: {e}"
        )

    calc = SevenNetCalculator(
        model="7net-0",
        device=device,
        **kwargs,
    )
    logger.info("SevenNet-0 loaded (device=%s).", device)
    return calc


# ---------------------------------------------------------------------------
# Dummy / LJ (CPU-only, for testing)
# ---------------------------------------------------------------------------

def _dummy_lj(**kwargs) -> Calculator:
    """
    Lennard-Jones calculator — CPU only, zero GPU requirement.

    Used exclusively for pipeline validation (parsing, MSD, Arrhenius fitting)
    without access to a GPU or MLIP weights.  Forces and energies are
    unphysical; diffusivity values from dummy mode are meaningless.
    """
    from ase.calculators.lj import LennardJones

    logger.warning(
        "Using DUMMY (Lennard-Jones) calculator. "
        "Results are NOT physically meaningful — use only for pipeline testing."
    )
    return LennardJones(epsilon=0.01, sigma=2.5, rc=5.0)
