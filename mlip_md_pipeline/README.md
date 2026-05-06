# MLIP-MD Solid Electrolyte Screening Pipeline

A reproducible, local-GPU pipeline for screening solid electrolyte candidates
using machine-learning interatomic potential molecular dynamics (MLIP-MD).
Extracts Li⁺/Na⁺ diffusivities, fits Arrhenius activation energies, and
converts to room-temperature conductivity estimates via the Nernst-Einstein
relation.

---

## Table of Contents

1. [Quick Start](#1-quick-start)
2. [Installation](#2-installation)
3. [Input Schema](#3-input-schema)
4. [Running the Pipeline](#4-running-the-pipeline)
5. [Expected Outputs](#5-expected-outputs)
6. [Output Schema for Downstream Agents](#6-output-schema-for-downstream-agents)
7. [Architecture Overview](#7-architecture-overview)
8. [Calculator Selection Guide](#8-calculator-selection-guide)
9. [Performance and GPU Notes](#9-performance-and-gpu-notes)
10. [Limitations and Scientific Caveats](#10-limitations-and-scientific-caveats)
11. [Troubleshooting](#11-troubleshooting)
12. [Mixed-Cation Structure Generation and Pre-MD Relaxation (Windows / Miniforge)](#12-mixed-cation-structure-generation-and-pre-md-relaxation-windows--miniforge)

---

## 1. Quick Start

```bash
# Clone / copy pipeline
cd mlip_md_pipeline

# Create conda environment
conda env create -f environment.yml
conda activate mlip_md

# Validate the full pipeline with synthetic data (no GPU, < 30 s):
python run_pipeline.py --dummy --output-dir outputs/test/

# Expected: screening_results.csv + MSD + Arrhenius plots in outputs/test/
```

---

## 2. Installation

### Recommended: conda environment

```bash
conda env create -f environment.yml
conda activate mlip_md
```

### Alternatively: pip only

```bash
pip install -r requirements.txt
```

### MLIP calculator installation

Each MLIP backend requires separate installation.  Install only what you need:

| Calculator | Install command | GPU required |
|------------|----------------|-------------|
| MACE-MP-0  | `pip install mace-torch` | Yes (CUDA) |
| CHGNet     | `pip install chgnet` | Recommended |
| SevenNet   | `pip install sevenn` | Yes (CUDA) |
| Dummy (LJ) | Built-in ASE | No |

MACE-MP-0 model weights (~100 MB) are downloaded automatically to
`~/.cache/mace/` on first use.

### PyTorch + CUDA

Install PyTorch with the CUDA version matching your driver:

```bash
# CUDA 12.1 (most common for RTX 30xx/40xx):
pip install torch --index-url https://download.pytorch.org/whl/cu121

# CUDA 11.8:
pip install torch --index-url https://download.pytorch.org/whl/cu118
```

Check your CUDA version: `nvidia-smi | head -4`

---

## 3. Input Schema

### candidate_screening_initial.csv

Located at `/home/user/workspace/candidate_screening_initial.csv`.
The pipeline reads this file directly.

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `candidate_id` | str | ✓ | Unique label (used in filenames) |
| `structure_file` | str | ✓ | CIF/POSCAR filename in `--structures-dir` |
| `mobile_species` | str | ✓ | Comma-separated species, e.g. `Li` or `Li,Na` |
| `ion_charge` | int | ✗ | Formal charge of mobile ion (default: 1) |
| `formula` | str | ✗ | Chemical formula (documentation only) |
| `space_group` | str | ✗ | Space group label (documentation only) |
| `calculator` | str | ✗ | Override global calculator per candidate |
| `device` | str | ✗ | Override GPU device per candidate |
| `total_time_ps` | float | ✗ | Override production time (ps) |
| `equilibration_time_ps` | float | ✗ | Override equilibration time (ps) |
| `timestep_fs` | float | ✗ | Override MD timestep (fs) |
| `temperatures_K` | str | ✗ | Semicolon-separated, e.g. `600;800;1000` |
| `notes` | str | ✗ | Free-text comments |

**Example row:**
```
LGPS,LGPS_mp-696128.cif,Li,1,Li10GeP2S12,P4_2/nmc,mace-mp-0,cuda,200,10,2.0,600;800;1000,Reference superionic
```

### Downloading structures from Materials Project

```bash
export MP_API_KEY="your_mp_api_key"   # from materialsproject.org/api
python scripts/fetch_mp_structures.py \
    --csv ../candidate_screening_initial.csv \
    --output-dir structures/
```

The script extracts `mp-XXXXXX` IDs from filenames and downloads
conventional standard cells.

### Structure file formats

Supported automatically (format inferred from filename):

| Extension | Format |
|-----------|--------|
| `.cif` | CIF (ASE) |
| `POSCAR`, `CONTCAR` | VASP format |
| `.vasp`, `.poscar` | VASP format |
| `.extxyz`, `.xyz` | Extended XYZ |

---

## 4. Running the Pipeline

### Test mode (no GPU required)

```bash
python run_pipeline.py --dummy --output-dir outputs/test/
```

Generates synthetic Brownian trajectories with known Ea = 0.35 eV and
validates every module end-to-end in < 30 s.

### Full screening run

```bash
python run_pipeline.py \
    --csv ../candidate_screening_initial.csv \
    --structures-dir structures/ \
    --calculator mace-mp-0 \
    --device cuda \
    --temperatures 600 800 1000 \
    --total-time 200 \
    --equil-time 10 \
    --timestep 2.0 \
    --write-interval 10 \
    --supercell-min 10.0 \
    --workers 1 \
    --output-dir outputs/run_001/
```

### Key flags

| Flag | Default | Notes |
|------|---------|-------|
| `--calculator` | `mace-mp-0` | MLIP backend |
| `--device` | `cuda` | PyTorch device |
| `--temperatures` | `600 800 1000` | Repeat for each T: `-T 600 -T 800 -T 1000` |
| `--total-time` | `200` | Production MD time per T (ps) |
| `--equil-time` | `10` | Equilibration (ps) — increase for sluggish systems |
| `--timestep` | `2.0` | fs; use 1.0 for light elements (Li-H compounds) |
| `--write-interval` | `10` | Frames written every N steps (10 × 2 fs = 20 fs/frame) |
| `--supercell-min` | `10.0` | Minimum box length (Å); smaller cells are replicated |
| `--workers` | `1` | Parallel processes; 1 per GPU recommended |
| `--dummy` | off | Synthetic test mode (no GPU) |

### Running unit tests

```bash
cd mlip_md_pipeline
pip install pytest
python -m pytest tests/ -v
```

All tests run on CPU; no MLIP weights required.

---

## 5. Expected Outputs

After a successful run, `--output-dir` contains:

```
outputs/run_001/
├── screening_results.csv          # Summary table (one row per candidate)
├── failed_candidates.csv          # Candidates that failed (if any)
├── LGPS/
│   ├── LGPS_T600K.traj            # ASE trajectory (binary)
│   ├── LGPS_T800K.traj
│   ├── LGPS_T1000K.traj
│   ├── LGPS_T600K.log             # Thermostat log (T, E, P per step)
│   ├── LGPS_T800K.log
│   ├── LGPS_T1000K.log
│   ├── LGPS_msd.png               # MSD curves at 600/800/1000 K
│   └── LGPS_arrhenius.png         # Arrhenius plot with fit
├── LLZO_cubic/
│   └── ...
└── ...
```

### Trajectory files

ASE binary `.traj` files.  Read with:
```python
from ase.io import Trajectory
traj = Trajectory("LGPS_T600K.traj")
for atoms in traj:
    print(atoms.get_positions())
```

Convert to extxyz: `ase convert LGPS_T600K.traj output.extxyz`

---

## 6. Output Schema for Downstream Agents

### screening_results.csv columns

| Column | Unit | Description |
|--------|------|-------------|
| `candidate_id` | — | Candidate label |
| `mobile_species` | — | e.g. `Li` |
| `Ea_eV` | eV | Arrhenius activation energy |
| `Ea_std_eV` | eV | 1σ uncertainty on Ea |
| `D0_cm2s` | cm²/s | Pre-exponential diffusivity |
| `D0_std_cm2s` | cm²/s | 1σ uncertainty on D0 |
| `arrhenius_R2` | — | R² of ln(D) vs 1/T linear fit |
| `sigma_300K_extrap_S_cm` | S/cm | Conductivity extrapolated to 300 K |
| `n_temperatures_fit` | — | Number of T points used in fit |
| `status` | — | `ok` or `error` |
| `T600_D_cm2s` | cm²/s | D at 600 K |
| `T800_D_cm2s` | cm²/s | D at 800 K |
| `T1000_D_cm2s` | cm²/s | D at 1000 K |
| `T600_D_std_cm2s` | cm²/s | σ(D) at 600 K |
| `T600_R2_msd` | — | R² of MSD linear fit at 600 K |
| `T600_sigma_S_cm` | S/cm | Nernst-Einstein σ at 600 K |

(Columns for T800, T1000 follow the same pattern.)

### Interpretation guide for downstream stability/report agents

- `sigma_300K_extrap_S_cm > 1e-3` → target range for solid electrolytes
- `Ea_eV < 0.3` → fast-ion conductor class (LGPS-like)
- `0.3 ≤ Ea_eV ≤ 0.6` → moderate; competitive with oxides
- `Ea_eV > 0.6` → poor ionic conductor at RT
- `arrhenius_R2 < 0.95` → Arrhenius fit unreliable; may need more temperatures
- `T600_R2_msd < 0.90` → MSD not in diffusive regime at 600 K; trajectory likely too short or material non-conducting at 600 K

---

## 7. Architecture Overview

```
run_pipeline.py   (CLI entry point via Click)
│
├── orchestrator.py
│     load_candidates_from_csv()   reads candidate_screening_initial.csv
│     run_screening()              multiprocessing.Pool dispatcher
│     _worker()                    per-candidate isolated process
│
├── ingestion.py
│     load_structure()             CIF / POSCAR → ASE Atoms
│     validate_mobile_ions()       checks Li/Na presence
│     build_supercell()            replicates to ≥ min_length Å
│
├── calculators.py
│     get_calculator()             factory → MACE / CHGNet / SevenNet / LJ
│
├── md_runner.py
│     MDConfig                     dataclass (temperatures, timestep, etc.)
│     run_nvt_md()                 Langevin NVT at each temperature
│     read_trajectory_positions()  .traj → unwrapped positions
│     _unwrap_positions()          minimum-image unwrapping
│
├── msd.py
│     compute_msd()                OTO MSD via FFT
│     fit_diffusivity()            linear fit → D [cm²/s]
│     analyse_trajectory()         combines above → MSDResult
│
├── arrhenius.py
│     fit_arrhenius()              weighted ln(D) vs 1/(kBT) regression → Ea, D0
│     nernst_einstein_conductivity()  D [cm²/s] → σ [S/cm]
│     full_arrhenius_analysis()    wraps both → ArrheniusResult
│
├── outputs.py
│     save_results_csv()           ArrheniusResult list → CSV
│     plot_msd()                   MSD curves + fit highlight
│     plot_arrhenius()             ln(D) vs 1000/T + secondary T axis
│
├── dummy.py
│     make_dummy_structure()       synthetic Li-argyrodite-like atoms
│     synthetic_msd_data()         Brownian motion with known Ea = 0.35 eV
│     run_dummy_pipeline()         full end-to-end test
│
└── tests/test_pipeline.py        pytest suite (CPU only)
```

---

## 8. Calculator Selection Guide

### MACE-MP-0

- **Paper:** Batatia et al., *Science* 384, 2024
- **Best for:** General screening; broad periodic-table coverage
- **Strengths:** Fast, well-validated across ~150,000 MP structures
- **Weaknesses:** Slightly lower accuracy on sulfides vs oxides; no dispersion
- **Model size:** ~5M parameters (medium); inference ~10–50 ms/step on RTX 3090
- **Download:** automatic via `mace_mp(model="medium")`

### CHGNet

- **Paper:** Deng et al., *Nature Machine Intelligence* 5, 2023
- **Best for:** Magnetic systems; explicit charge prediction
- **Strengths:** Magnetic moment prediction; good on oxides
- **Weaknesses:** Slower than MACE per step; narrower training set
- **Memory:** ~1.8 GB VRAM for 300-atom cell

### SevenNet (7net-0)

- **Paper:** Park et al., *J. Chem. Theory Comput.* 20, 2024
- **Best for:** High-throughput screening; speed-accuracy balance
- **Strengths:** SE(3)-equivariant; fast inference
- **Weaknesses:** Less community validation than MACE-MP-0 at time of writing
- **Install:** `pip install sevenn` or from GitHub

### Switching per candidate

Set the `calculator` column in your CSV to mix backends:
```
LGPS,LGPS.cif,Li,1,...,mace-mp-0,...
LLZO,LLZO.cif,Li,1,...,chgnet,...
```

---

## 9. Performance and GPU Notes

### Typical throughput (RTX 3090, 24 GB VRAM)

| System | Atoms | Steps/s (MACE medium) | 200 ps @ 2 fs/step |
|--------|-------|-----------------------|---------------------|
| Li6PS5Cl 1×1×1 (52 atoms) | 52 | ~800 | ~3 min |
| Li6PS5Cl 2×2×2 supercell | 416 | ~120 | ~19 min |
| LGPS 2×2×1 (256 atoms) | 256 | ~200 | ~11 min |

Full 7-candidate screening at 3 temperatures × 200 ps ≈ **2–5 hours on one GPU.**

### Multi-GPU parallelism

Use `--workers N` where N = number of available GPUs.
Per-candidate GPU assignment is controlled via `gpu_id` in `CandidateConfig`
(set manually or via CSV `device` column: `cuda:0`, `cuda:1`, etc.)

### VRAM requirements

| Calculator | 50 atoms | 200 atoms | 500 atoms |
|------------|----------|-----------|-----------|
| MACE medium | ~2 GB | ~4 GB | ~9 GB |
| CHGNet | ~1.5 GB | ~3 GB | ~7 GB |

If you get CUDA OOM errors: reduce `--supercell-min` or use `--calculator chgnet`.

### CPU fallback

All calculators work on CPU (`--device cpu`) but are 10–100× slower.
The dummy/LJ calculator always runs on CPU and is suitable for pipeline validation only.

---

## 10. Limitations and Scientific Caveats

### MLIP accuracy for ionic conductivity

1. **Force field transferability**: Universal MLIPs (MACE-MP-0, CHGNet, SevenNet)
   are trained on static DFT calculations (typically PBE/GGA).
   They may not accurately reproduce the flat potential-energy landscape
   required for fast ionic conduction, especially in novel compositions
   outside the training distribution.

2. **Short trajectories**: 200 ps NVT MD captures diffusion only if D > ~10⁻⁸ cm²/s
   at the simulation temperature.  For materials with Ea > 0.5 eV, the 600 K
   diffusivity may be near the detection limit; extrapolation to 300 K will have
   large uncertainty.

3. **Haven ratio**: The Nernst-Einstein equation assumes uncorrelated single-ion
   diffusion (Haven ratio H_R = 1).  Real fast-ion conductors often have
   H_R = 0.1–0.8 due to cooperative hopping and correlation effects.
   This pipeline does **not** apply Haven ratio corrections; conductivities
   may be overestimated by up to 10×.

4. **Finite-size effects**: Small supercells (< 200 atoms) can bias MSD
   via periodic image interactions.  The `--supercell-min 10.0` flag mitigates
   this; increase to 12–15 Å for materials with long correlation lengths.

5. **NVT vs NPT**: All MD here is at fixed volume (NVT).  Real electrolytes
   undergo thermal expansion; NPT would be more rigorous.  Lattice parameters
   from PBE-DFT are typically 1–3% overestimated, slightly inflating D.

6. **Arrhenius linearity**: The Arrhenius model assumes a single thermally
   activated hopping mechanism.  Many solid electrolytes show curved Arrhenius
   plots (multiple hop pathways, phase transitions).  Fitting over 600–1000 K
   may not represent the 300 K mechanism.

7. **Thermostat artefacts**: The Langevin thermostat introduces an artificial
   friction coupling to a thermal bath.  Low friction (0.01 fs⁻¹) minimises
   bias on D.  Avoid high friction (> 0.1 fs⁻¹) for production runs.

8. **Configurational sampling**: 200 ps is often insufficient for structures
   with slow framework relaxation.  Watch for MSD saturation (plateau) rather
   than linear growth — this indicates the system is not in the diffusive regime.

### Recommended validation workflow

1. Run `--dummy` to confirm pipeline integrity.
2. Test on LGPS (literature D ~10⁻⁷ cm²/s at 600 K) as a positive control.
3. Check MSD linearity at each temperature before trusting D values.
4. Flag any R² < 0.95 (MSD fit) or arrhenius_R2 < 0.90 for manual review.
5. Cross-reference predictions against NMR/AIMD benchmarks before experimental synthesis.

---

## 11. Troubleshooting

### `ImportError: mace-torch is not installed`
```bash
pip install mace-torch
```

### CUDA OOM (out-of-memory)
- Reduce `--supercell-min` (e.g. `--supercell-min 8.0`)
- Use CHGNet (lower VRAM): `--calculator chgnet`
- Add `--device cpu` as last resort (slow)

### `ValueError: No mobile species found`
The structure file lacks Li/Na atoms.  Check:
- Correct `mobile_species` column in CSV (case-sensitive: `Li` not `li`)
- Structure file is the correct phase (some polymorphs store Li as Li+)

### Flat MSD (MSD doesn't grow with time)
- Material may not be a conductor at this temperature
- Extend trajectory: `--total-time 500` and `--equil-time 50`
- Try higher temperatures (1200 K, 1500 K)
- Check structure validity: run `python -c "from mlip_md.ingestion import load_structure; print(load_structure('yourfile.cif'))"`

### Negative or zero diffusivity
- MSD linear fit failed; check R² in CSV
- Trajectory too short for the diffusive regime
- Material has Ea > 0.8 eV — not a conductor at simulation temperatures

### `mp.spawn` errors on Linux
If you see `RuntimeError: context has already been set`, try:
```bash
export PYTHONPATH=/path/to/mlip_md_pipeline:$PYTHONPATH
```
Or use `--workers 1` (serial mode, always stable).

---

## File Reference

```
mlip_md_pipeline/
├── README.md                     # This file
├── environment.yml               # conda environment spec
├── requirements.txt              # pip requirements
├── run_pipeline.py               # Main CLI entry point
├── mlip_md/
│   ├── __init__.py
│   ├── ingestion.py              # CIF/POSCAR → ASE Atoms
│   ├── calculators.py            # MACE/CHGNet/SevenNet/dummy adapters
│   ├── md_runner.py              # NVT MD + trajectory I/O
│   ├── msd.py                    # MSD extraction + D fitting
│   ├── arrhenius.py              # Arrhenius + Nernst-Einstein σ
│   ├── outputs.py                # CSV + MSD/Arrhenius plots
│   ├── orchestrator.py           # Multiprocessing dispatcher + CSV loader
│   └── dummy.py                  # Synthetic test mode
├── scripts/
│   ├── fetch_mp_structures.py    # Download CIFs from Materials Project
│   └── plot_convergence.py       # Post-hoc MSD convergence diagnostic
├── tests/
│   └── test_pipeline.py          # pytest unit tests (CPU only)
└── structures/
    └── example/                  # Place CIF/POSCAR files here
```

---

## 12. Mixed-Cation Structure Generation and Pre-MD Relaxation (Windows / Miniforge)

For mixed-cation Li3MCl6 candidates such as **Li3Er0.5Sc0.5Cl6** and
**Li3In0.5Er0.5Cl6**, two helper scripts live under `scripts/`:

- `scripts/generate_mixed_cation_supercell.py` — builds substituted supercells from a parent CIF/POSCAR, with multiple random seeds.
- `scripts/relax_structures.py` — relaxes the generated structures with CHGNet (default) or optional MACE-MP-0 before MD.

A small example manifest is at `manifests/mixed_cation_examples.csv`.

### 12.1 Honesty caveats (read before using results)

- The substitution ordering is **random** on the M sublattice. To approximate physical behavior, generate ≥ 5 seeds, relax all of them, and inspect the relative-energy spread. Treat the lowest-energy ordering as a **screening hypothesis**, not a ground state.
- Universal MLIPs (CHGNet, MACE-MP-0) for **lithium halides** are not guaranteed accurate. Confirm low-energy candidates with DFT (e.g. VASP or QE) before claiming a specific ordering or driving experimental synthesis.
- Convergence with `--fmax 0.05 --max-steps 200` is a reasonable default but may be insufficient for cells with strained dopant environments — increase `--max-steps` and re-check.

### 12.2 Windows Miniforge — exact commands

Open **Miniforge Prompt** and run, from the repository root:

```bat
:: 1) Activate the project env and install MLIP backends
conda activate mlip_md
pip install chgnet
pip install mace-torch    :: optional; only if you want MACE relaxations

cd mlip_md_pipeline

:: 2) Generate Li3Er0.5Sc0.5Cl6 supercells from a Li3ErCl6 parent (5 seeds)
python scripts\generate_mixed_cation_supercell.py ^
    --parent ..\structures\Li3ErCl6.cif ^
    --host-element Er ^
    --dopant Sc ^
    --dopant-fraction 0.5 ^
    --supercell 2 2 2 ^
    --seeds 0 --seeds 1 --seeds 2 --seeds 3 --seeds 4 ^
    --output-dir ..\structures\Li3Er0.5Sc0.5Cl6\

:: 3) Generate Li3In0.5Er0.5Cl6 supercells from a Li3InCl6 parent
python scripts\generate_mixed_cation_supercell.py ^
    --parent ..\structures\Li3InCl6.cif ^
    --host-element In ^
    --dopant Er ^
    --dopant-fraction 0.5 ^
    --supercell 2 2 2 ^
    --seeds 0 --seeds 1 --seeds 2 --seeds 3 --seeds 4 ^
    --output-dir ..\structures\Li3In0.5Er0.5Cl6\

:: 4) Relax all generated POSCARs with CHGNet (CPU-safe)
python scripts\relax_structures.py ^
    --inputs ..\structures\Li3Er0.5Sc0.5Cl6\*.vasp ^
    --calculator chgnet ^
    --device cpu ^
    --fmax 0.05 ^
    --max-steps 200 ^
    --output-dir ..\structures\Li3Er0.5Sc0.5Cl6\relaxed\

python scripts\relax_structures.py ^
    --inputs ..\structures\Li3In0.5Er0.5Cl6\*.vasp ^
    --calculator chgnet ^
    --device cpu ^
    --output-dir ..\structures\Li3In0.5Er0.5Cl6\relaxed\

:: 5) (Optional) Re-relax the lowest-energy seeds with MACE-MP-0 for a cross-check
python scripts\relax_structures.py ^
    --inputs ..\structures\Li3Er0.5Sc0.5Cl6\*_seed0.vasp ^
    --calculator mace-mp-0 ^
    --device cuda ^
    --output-dir ..\structures\Li3Er0.5Sc0.5Cl6\relaxed_mace\
```

The `^` is the Windows command-prompt line continuation. In PowerShell use a
backtick (`` ` ``) instead. On Linux/macOS use `\` and forward slashes.

### 12.3 What you get

For each seed, the generator writes:

```
<output-dir>/<base>_seed<N>.cif
<output-dir>/<base>_seed<N>.vasp
<output-dir>/<base>_seed<N>.json     :: metadata (swap indices, host, fraction, warnings)
<output-dir>/generation_summary.json
```

The relaxer writes, per input file:

```
<output-dir>/<stem>_relaxed.cif
<output-dir>/<stem>_relaxed.vasp
<output-dir>/<stem>_relax.log        :: ASE optimizer log
<output-dir>/relaxation_summary.csv  :: energies, fmax, n_steps, convergence flag
```

Use the energies in `relaxation_summary.csv` to pick the lowest-energy seed
per composition, then feed that relaxed POSCAR into the MD pipeline via
`tier1a_mlip_input_template.csv`.

---

*Pipeline version 0.1.0.  See `mlip_md/__init__.py` for changelog.*
