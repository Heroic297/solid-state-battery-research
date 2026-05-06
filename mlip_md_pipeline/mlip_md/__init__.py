"""
mlip_md — MLIP-MD screening pipeline for solid electrolyte candidates.

Modules:
  ingestion     — CIF/POSCAR → ASE Atoms
  calculators   — MACE-MP-0, CHGNet, SevenNet adapter factory
  md_runner     — NVT MD at multiple temperatures
  msd           — Li+/Na+ MSD extraction from trajectories
  arrhenius     — Arrhenius fit + Nernst-Einstein conversion to σ (S/cm)
  outputs       — CSV summary + MSD/Arrhenius plots
  orchestrator  — multiprocessing across candidates
  dummy         — synthetic test mode (no GPU required)
"""

__version__ = "0.1.0"
