#!/usr/bin/env python
"""
Download CIF structures from Materials Project for all candidates in the screening CSV.

Requirements:
    pip install mp-api

Usage:
    export MP_API_KEY="your_key_here"   # from https://next-gen.materialsproject.org/api
    python scripts/fetch_mp_structures.py \
        --csv ../candidate_screening_initial.csv \
        --output-dir ../structures/

The script reads the `structure_file` column (e.g. "LGPS_mp-696128.cif"),
extracts the mp-XXXXXX Materials Project ID from the filename,
and downloads the conventional standard structure.

If a file already exists, the download is skipped.
"""

from __future__ import annotations

import logging
import os
import re
import sys
from pathlib import Path

import click
import pandas as pd

logger = logging.getLogger(__name__)


def extract_mp_id(filename: str) -> str | None:
    """Extract mp-XXXXXX from a filename like 'LGPS_mp-696128.cif'."""
    match = re.search(r"(mp-\d+)", filename)
    return match.group(1) if match else None


@click.command()
@click.option("--csv", "csv_path", required=True,
              help="Path to candidate_screening_initial.csv")
@click.option("--output-dir", "output_dir", default="structures/",
              show_default=True, help="Directory to save CIF files.")
@click.option("--api-key", "api_key", default=None,
              help="Materials Project API key (or set MP_API_KEY env var).")
def main(csv_path, output_dir, api_key):
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")

    api_key = api_key or os.environ.get("MP_API_KEY")
    if not api_key:
        logger.error(
            "No MP_API_KEY found. "
            "Set via --api-key or export MP_API_KEY=... "
            "Get your key at https://next-gen.materialsproject.org/api"
        )
        sys.exit(1)

    try:
        from mp_api.client import MPRester
    except ImportError:
        logger.error("mp-api not installed. Run: pip install mp-api")
        sys.exit(1)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(csv_path)
    logger.info("Loaded %d candidates from %s", len(df), csv_path)

    with MPRester(api_key) as mpr:
        for _, row in df.iterrows():
            filename = str(row["structure_file"]).strip()
            out_path = output_dir / filename

            if out_path.exists():
                logger.info("SKIP (exists): %s", out_path)
                continue

            mp_id = extract_mp_id(filename)
            if mp_id is None:
                logger.warning(
                    "Cannot extract MP ID from '%s' — skipping. "
                    "Rename the file to include mp-XXXXXX or download manually.",
                    filename,
                )
                continue

            logger.info("Downloading %s → %s", mp_id, out_path)
            try:
                structure = mpr.get_structure_by_material_id(
                    mp_id, conventional_unit_cell=True
                )
                structure.to(fmt="cif", filename=str(out_path))
                logger.info("Saved %s", out_path)
            except Exception as e:
                logger.error("Failed to download %s: %s", mp_id, e)

    logger.info("Done. Structures in: %s", output_dir.resolve())


if __name__ == "__main__":
    main()
