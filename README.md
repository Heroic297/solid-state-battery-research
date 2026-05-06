# Solid-State Battery Electrolyte Discovery

AI-assisted computational screening workspace for novel lithium and sodium solid electrolyte candidates, focused on high room-temperature ionic conductivity and wide electrochemical stability windows.

## Honest project status

This repository is **not yet a proven discovery**. It contains:

- a literature-backed and proxy-ranked candidate list;
- stability, interface, and prior-art triage tables;
- a runnable local-GPU MLIP-MD pipeline;
- pre-print-style notes and synthesis/validation plans.

The main unresolved step is production validation of the novel candidates by MLIP-MD, DFT relaxation/stability checks, and ultimately experimental synthesis/EIS.

## Top hypothesis

The strongest original hypothesis from the current screen is:

```text
Li3Er0.5Sc0.5Cl6
```

Rationale: a mixed-cation Li3MCl6 chloride designed to stabilize a cubic-close-packed chloride sublattice and promote 3D Li conduction while preserving the high oxidative stability typical of chloride halide electrolytes.

## Repository layout

```text
.
├── data/
│   ├── final_ranked_candidates.csv
│   ├── top15_candidate_summary.csv
│   ├── candidate_screening_initial.csv
│   ├── literature_matrix.csv
│   ├── stability_interface_analysis.csv
│   └── novelty_checks.csv
├── manifests/
│   ├── tier1a_mlip_input_template.csv
│   └── tier1a_structure_manifest.csv
├── mlip_md_pipeline/
│   ├── README.md
│   ├── environment.yml
│   ├── requirements.txt
│   ├── run_pipeline.py
│   ├── mlip_md/
│   ├── scripts/
│   └── tests/
├── plots/
└── reports/
```

## Quick start on a local GPU machine

```bash
git clone https://github.com/<your-user-or-org>/solid-state-battery-research.git
cd solid-state-battery-research/mlip_md_pipeline

conda env create -f environment.yml
conda activate mlip_md

pip install mace-torch chgnet

# Validate the full code path without GPU production MD.
python run_pipeline.py --dummy --output-dir outputs/test/
```

If the dummy run passes, the next step is to supply real CIF/POSCAR structures for the candidates in:

```text
../manifests/tier1a_mlip_input_template.csv
```

Then run:

```bash
python run_pipeline.py \
  --csv ../manifests/tier1a_mlip_input_template.csv \
  --structures-dir structures/ \
  --calculator mace-mp-0 \
  --device cuda \
  --total-time 200 \
  --workers 1 \
  --output-dir outputs/run_001/
```

## Minimum validation campaign

Do not screen all 44 candidates first. Start with:

| Role | Composition | Purpose |
|---|---|---|
| Novel candidate | Li3Er0.5Sc0.5Cl6 | Main original hypothesis |
| Novel candidate | Li3In0.5Er0.5Cl6 | Lower-risk mixed-cation chloride |
| Control | Li3InCl6 | Known chloride conductor |
| Control | Li3YCl6 | Known chloride conductor |

A result becomes scientifically interesting only if the controls are reproduced and a novel candidate shows:

- room-temperature \( \sigma > 10^{-3} \) S/cm;
- \(E_a < 0.30\) eV;
- diffusive MSD behavior rather than one-jump artifacts;
- stable relaxed structures across at least one trusted MLIP and preferably DFT;
- energy above hull below roughly 50 meV/atom;
- oxidative stability above 4 V by phase-diagram analysis.

## Important limitations

- Production MLIP-MD was not run before this repo was created.
- Several values are literature-derived or proxy-estimated, not newly computed.
- Universal MLIPs can be unreliable for unusual halide chemistries; validate against known controls.
- Nernst-Einstein conductivity can overestimate real conductivity if the Haven ratio is neglected.
- Direct Li-metal compatibility is unlikely for most chloride halides; these are usually catholyte or bilayer-electrolyte candidates.

## Next development tasks

1. Add structure-generation scripts for mixed-cation Li3MCl6 supercells.
2. Add CHGNet/MACE relaxation workflow before production MD.
3. Add seed-replicated MD runs and convergence diagnostics.
4. Add pymatgen phase-diagram notebooks for electrochemical window and interface analysis.
5. Add DFT/HPC workflow templates for final candidate validation.

