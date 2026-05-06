# Candidate Generation: Superionic Li/Na Solid Electrolyte Screening
## Methodological Notes

**Date generated:** 2026  
**Output file:** `/home/user/workspace/candidate_screening_initial.csv`  
**Total candidates:** 44 (10 known benchmarks, 34 novel/underexplored)

---

## 1. Scope and Families Targeted

The screening covers three families as specified:

| Family | Description | Sub-families |
|--------|-------------|--------------|
| Li-P-S-X | Argyrodite (Li6PS5X, X=Cl,Br,I), Thio-LISICON (LGPS-type), Li3PS4 | Ge/Si/Sn substitutions at P site; halide site mixing |
| Li-M-Cl/Br/I | Li3MCl6 rocksalt/layered; Li2M2/3Cl4 halospinel; LiMXCl4 non-close-packed | M = In, Y, Sc, Zr, Hf, Er; mixed-M; aliovalent substitution |
| Na-M-S/Se | Na3PS4/Se4, Na3SbS4/Se4, Na11Sn2PS12, Na7P3S/Se11 | W/Mo/As doping; Se/S mixed anion; halide doping |

Additional novel families explored:
- **Na-M-Cl/Br/I halides**: Na3MX6 (M = Gd, Y, In, Er; X = Br, I, mixed)
- **HE (high-entropy) halides**: Multi-cation Li3(In,Sc,Yb,Lu,Zr)Cl6

---

## 2. Data Sources Accessed

### 2.1 Successfully Accessed Sources

| Source | Access Method | Data Retrieved |
|--------|---------------|----------------|
| Materials Project (MP) | Web metadata via search/fetch (IP-blocked for direct API) | mp-IDs cited where found in literature; direct API access unavailable during this run due to rate-limiting block |
| Ceder Group Publications | PDF fetch via pplx content tool | LiMXCl4 family: E_hull (r2SCAN), conductivity (AIMD), activation energies |
| ACS/RSC/PMC literature | Web fetch + search | E_hull, band gaps (HSE06/optical), conductivities, space groups |
| USPEX team (Dallakyan 2026) | PDF fetch | Li3MX6 screening table: Eg, E_hull, σ, Ea for Li3InF6, Li3InI6, Li3YCl6, Li3YBr6, etc. |
| SandboxAQ arXiv 2604.02524 | Web search + fetch | AQVolt26 methodology, dataset scope (Li halides only, 322,656 DFT calcs) |
| Frontiers in Energy Research | Fetch | Na chalcogenide conductor summary table |
| UWO Nano Letters 2020 | Fetch | Li3Y1-xInxCl6 structural transition and conductivity data |
| OSTI.gov | Search | Na3SbS4, Na3PSe4 interface studies |

### 2.2 Materials Project API

The Materials Project public API was not accessible during this run—the web IP was rate-limited/blocked (returned `"error": "Your IP address has been (temporarily) blocked"`). Known MP IDs (mp-XXXXXX) are cited from literature cross-references where available. Direct E_hull and band gap retrieval via `mp-api` was not performed.

**Workaround**: E_hull and band gap values were extracted from peer-reviewed publications that explicitly report MP-derived values or their own DFT calculations (PBE/GGA or HSE06). All values are clearly flagged with source citations.

### 2.3 AFLOW/OQMD

AFLOW and OQMD were not directly queried via API endpoints in this run. However, the [OQMD grand canonical linear programming method](https://doi.org/10.1038/s41467-024-45258-3) is referenced in the high-entropy halide paper for chemical potential calculations, and AFLOW/MP/OQMD comparison methodology is cited from arXiv:2007.01988. No novel structural data was extracted from these sources that was not already covered by literature.

### 2.4 AQVolt26 / SandboxAQ Dataset

**Status: Available as public dataset (MLIP models on Hugging Face; arXiv paper 2604.02524); raw DFT data not individually queryable without institutional access.**

The AQVolt26 dataset ([SandboxAQ, 2026](https://arxiv.org/abs/2604.02524)) contains:
- **322,656 r²SCAN single-point DFT calculations** for lithium halide solid-state electrolytes
- Generated via surrogate-driven high-temperature (up to 2100 K) configurational sampling
- Covers ~5,000 unique Li halide structures; trained eSEN MLIP models available on Hugging Face

**Limitations documented here:**
1. The raw DFT dataset is not publicly downloadable as a structured CSV/JSON of individual compound properties (E_hull, band gap, etc. per composition).
2. The released artifacts are pretrained MLIP checkpoints, not a structured materials property database.
3. The paper explicitly covers **lithium halides only**—not sulfides, selenides, or sodium systems.
4. AQVolt26 is used here only as **contextual motivation** (demonstrating that ML-driven halide screening at scale is state-of-practice, and that high-temperature dynamics require targeted datasets) and **not as a data source** for individual candidate properties.

**Contextual relevance**: AQVolt26 validates that LiMXCl4 and Li3MX6 families screened here are the current high-priority halide targets in the field. The approach of coupling MLIP-accelerated MD with targeted high-T sampling directly motivates the selection of novel LiHfFCl4 and LiZrFCl4 candidates (from Ceder 2025 AIMD work) that would benefit from AQVolt26-type validation.

---

## 3. Filtering Criteria Applied

| Criterion | Threshold | Application |
|-----------|-----------|-------------|
| Thermodynamic stability | E_hull < 50 meV/atom (hard filter) | Applied where data available; compositions with E_hull > 50 meV/atom flagged or excluded. Two candidates (LiSbOCl4, not included) exceed 50 meV/atom in r2SCAN calculations |
| Band gap (electronic insulation) | > 3.5 eV preferred | Flagged in "band_gap_eV" column; low-gap compounds (e.g., Na3SbS4 family ~2.5 eV) retained as known benchmarks with limitation noted |
| Mobile ion connectivity | 3D > 2D > 1D | Encoded in connectivity_score sub-score |
| Anion framework | BCC-like/non-close-packed preferred | Encoded in anion_framework_score |
| Site disorder/vacancy | Disorder or vacancies required for fast conduction | Encoded in disorder_score |

Candidates with band gap below threshold (< 3.5 eV) are retained in the table with explicit rationale, as several champion conductors (Na3SbS4, argyrodites) fall in this range—their inclusion provides important benchmarks.

---

## 4. Proxy Scoring Model

**CRITICAL DISCLAIMER**: This is a heuristic proxy score inspired by bond valence sum (BVS) concepts and known fast-ion descriptors. **No actual BVS calculation was performed**—no crystal structure refinement, bond length data, or formal BVS summation (Σ_i [r_i/r_0]^{-N}) was executed. The scores are heuristic estimates anchored to literature-reported conductivities and structural analysis.

### 4.1 Score Formula

```
Proxy Score (0–100) = 10 × Σ_j [w_j × s_j]

Sub-scores s_j ∈ [0, 10], weights w_j:
  j=1: Anion framework type          w=0.20
  j=2: Site disorder / vacancy       w=0.20
  j=3: Bottleneck size proxy         w=0.15
  j=4: Mobile ion connectivity       w=0.15
  j=5: E_hull stability              w=0.15
  j=6: Band gap                      w=0.10
  j=7: Anion polarizability/softness w=0.05
```

### 4.2 Sub-score Definitions

**j=1: Anion framework type** (reflects BCC-like vs. close-packed lattice)
- BCC-like (I-43m, body-centered anion) = 10
- Non-close-packed (LiMXCl4 type) = 9
- ccp (cubic close-packed, Cl/S FCC-like sublattice) = 8
- hcp-to-ccp transition = 6–7
- hcp (hexagonal close-packed) = 5
- Distorted/1D layered = 3–4

**j=2: Site disorder / vacancy** (vacancy and chemical disorder are critical for fast conduction)
- Multiple site disorder + vacancies present = 10
- Single-type disorder OR vacancy only = 7–9
- Fully ordered, no vacancies = 3

**j=3: Bottleneck size proxy** (size of the saddle-point void through which mobile ions pass)
- Large cell + highly polarizable anion (Se2-, I-) + wide channel = 9–10
- Moderate cell + Br-/S2- = 7–8
- Small cell + Cl- = 5–7
- F- (very small) = 2–4

**j=4: Mobile ion connectivity** (dimensionality of diffusion network)
- 3D isotropic = 10
- 3D anisotropic = 7–9
- 2D layered (some 3D character) = 6–8
- 1D channels = 2–5

**j=5: E_hull stability score** (thermodynamic proximity to convex hull)
- < 5 meV/atom = 10
- 5–20 meV/atom = 9
- 20–35 meV/atom = 8
- 35–50 meV/atom = 7
- 50–100 meV/atom = 4
- > 100 meV/atom = 1
- Unknown = 5

**j=6: Band gap score** (electronic insulation)
- > 5 eV = 10
- 4–5 eV = 8
- 3.5–4 eV = 6
- 2.5–3.5 eV = 4
- < 2.5 eV = 1–2
- Unknown = 5

**j=7: Anion polarizability** (soft anion lowers migration barriers)
- Se2- or I- = 10
- Mixed Se/S or Br/I = 9
- S2- or Br- = 7
- Cl- = 5
- F- = 2

### 4.3 Inspiration from BVS Concepts

The scoring is *inspired by* bond valence sum analysis in the following ways:

1. **Anion framework type** ↔ BVS pathway energy landscape: BCC-like anion sublattices are known to provide low-barrier tetrahedral→octahedral site hopping. Non-BVS calculations but the physical reasoning is equivalent.
2. **Bottleneck size proxy** ↔ BVS bottleneck analysis: The critical bottleneck radius (ideally > 1.6 Å for Li+, > 1.8 Å for Na+) determines the migration barrier. Larger/more polarizable anions increase this.
3. **Disorder score** ↔ BVS site energy variance: Disordered structures have broader distributions of BVS minima and saddle points, generally lowering the maximum barrier.
4. **Connectivity** ↔ BVS pathway topology: 3D connected isoenergetic pathways in BVS maps indicate 3D fast-ion conduction.

True BVS calculations require: (a) experimental or DFT-relaxed crystal structure coordinates, (b) bond valence parameters r₀ and N for each ion pair, and (c) summation over all bonds. These were not performed.

---

## 5. Novelty Classification

| Flag | Meaning |
|------|---------|
| Known benchmark | Well-studied compound with >5 independent experimental reports |
| Known (lightly explored) | Synthesized and characterized but limited follow-up; limited optimization |
| Recently discovered benchmark | Experimental conductivity confirmed in last 2 years |
| Novel candidate (experimental) | Synthesized but only 1-2 reports; composition space underexplored |
| Novel candidate (DFT prediction) | Predicted computationally; no confirmed synthesis |
| Novel candidate (hypothetical) | Chemically plausible extrapolation; no literature precedent |

---

## 6. Structural Descriptors for Fast-Ion Conduction

Beyond the proxy score, the following structural features are used to qualitatively assess each candidate:

1. **BCC-like anion sublattice**: Body-centered arrangement of anions provides interconnected tetrahedral + octahedral sites with low-barrier hops (ΔE ~50–150 meV). Best exemplified by Na3PS4 (I-43m) and Na3SbS4.

2. **ccp (FCC-like) anion sublattice**: Face-centered halide packing (e.g., Li3InCl6 C2/m) provides 3D octahedral vacancy networks for Li+. Monoclinic C2/m is the dominant ccp variant; key for >1 mS/cm halide conductors.

3. **Non-close-packed framework**: LiMXCl4 family has 1D chains of MX₂Cl₄ octahedra bound by van der Waals forces. "Soft cradle effect" enables liquid-like Li+ transport. Projected 10–100 mS/cm (Jun et al., Matter 2025).

4. **Site disorder**: Chemical disorder (e.g., Cl-/S2- co-occupation in argyrodites; mixed M-cation in high-entropy halides) creates a flatter energy landscape for mobile ions.

5. **Aliovalent substitution**: Replacing M3+ with M4+ (Zr, Hf, Ti) in Li3MCl6 or Na3MCl6 creates Li/Na vacancies that activate fast-ion transport.

6. **Anion polarizability**: Polarizable anions (Se2- > S2- > Br- > Cl- > F-) reduce electrostatic barriers for mobile ion migration. Key parameter for chalcogenide systems.

7. **Connected mobile-ion sublattice**: Requires that Li+ or Na+ sites form a percolating network through the structure (as opposed to isolated cage occupancies).

---

## 7. Key Findings and Priority Candidates

### 7.1 Highest-Priority Novel Li Candidates (by score)

1. **Li6.6P0.4Ge0.6S5I** (score 87): Ge-substituted argyrodite; experimentally confirmed 2.9 mS/cm; near-hull stable; large-I argyrodite with Ge-driven Li excess. Priority for synthesis optimization at higher Ge loading.

2. **Li2.8In0.2Sc0.2Yb0.2Lu0.2Zr0.2Cl6** (HE-SE, score 86): High-entropy 5-cation halide; experimentally demonstrated 2.13 mS/cm; stacking faults enable 3D diffusion. Priority for cation composition sweeps.

3. **Li6PS5Cl0.5Br0.5** (score 86): Mixed-halide argyrodite; predicted ~2–3 mS/cm; near-hull (~9 meV/atom); fills gap between pure Cl and Br endpoint studies.

4. **LiHfFCl4 (Li-stuffed)** (score 84): Non-close-packed family; AIMD-predicted 137.5 mS/cm; E_hull 12–32 meV/atom depending on competing phases. Highest-novelty high-priority target; no experimental synthesis yet.

5. **Li3Er0.5Sc0.5Cl6** (score 81): Mixed Er-Sc chloride; In-inspired ccp phase engineering applied to underexplored Er system. Purely hypothetical; strong physical basis.

### 7.2 Highest-Priority Novel Na Candidates (by score)

1. **Na2.88Sb0.88Mo0.12S4** (score 85): Mo analogue of W-doped Na3SbS4 champion (32 mS/cm). Mo6+ mechanism identical to W6+; more abundant. No experimental report found.

2. **Na3GdBr3I3** (score 83): Mixed Br/I Na halide; 7.5 mS/cm AIMD-predicted; C2/m monoclinic ccp structure; highest-priority Na halide target.

3. **Na3SbSe3S** (score 84.5): Mixed Se/S cubic BCC; 3D isotropic; near-hull; anion disorder adds to BCC benefits.

4. **Na7P3Se11** (score 79): Full Se substitution of Na7P3S11 (>10 mS/cm predicted). Hypothetical; synthetic challenge is the polyanion stability.

### 7.3 Score Distribution

- Scores range from **56.5 to 87.0** across all 44 candidates
- Known benchmarks cluster at **65–87** (reflecting literature-validated high performance)
- Novel candidates span **56–87** depending on structural favorability
- Li-P-S-X sulfides score strongly on disorder/BCC but weakly on band gap
- Li-M-Cl halides score strongly on band gap but weakly on anion polarizability
- Na-M-S/Se sulfides score strongly on BCC framework and polarizability but weakly on band gap

---

## 8. Limitations and Caveats

1. **No API access to Materials Project**: E_hull values from MP are cited as reported in peer-reviewed papers, not directly retrieved. Some values may reflect older MP database states (pre-2024 r2SCAN corrections).

2. **Proxy scores are not predictive conductivities**: The score reflects a qualitative ranking of structural favorability. Candidates with similar scores may differ by orders of magnitude in conductivity due to subtle effects (e.g., synthesis-route-dependent disorder, grain boundary resistance, amorphous fraction).

3. **AIMD-predicted conductivities are generally overestimated**: High-temperature AIMD extrapolations typically overestimate room-temperature conductivity by 0.5–2× due to non-Arrhenius behavior at low temperature and structural relaxation effects. Values reported as "AIMD predicted" should be treated as upper bounds.

4. **E_hull values are functional-dependent**: PBE/GGA values for halides underestimate lattice energies; r2SCAN values (used in Ceder 2025 paper) are more accurate for halides. The ~12–30 meV/atom range for many halides is within the synthesizability threshold (~50 meV/atom) regardless of functional.

5. **Band gaps from GGA are underestimated**: Reported GGA band gaps are typically ~50% lower than experiment. HSE06 values are preferred where available. The ">3.5 eV" filter is applied to HSE06 or experimental optical values where stated; GGA values for the same material would be lower.

6. **AQVolt26 dataset**: The full DFT dataset is not publicly available as structured compound properties. Only MLIP model weights are released on Hugging Face. The 322,656 r2SCAN calculations cover Li halides at high temperature—not sulfides, selenides, or sodium systems. This dataset was not used as a data source for individual candidate properties in this screening.

7. **AFLOW/OQMD API**: Not directly queried in this run. Structural data from these databases is available via literature cross-references only.

8. **Hypothetical candidates**: 12 of the 44 candidates are purely hypothetical (no synthesis report). These are flagged as "Novel candidate (hypothetical)" and carry higher uncertainty in all estimated properties.

---

## 9. Recommended Next Steps

1. **DFT validation** of top-10 novel candidates: Run PBE/r2SCAN geometry optimization and phonon calculations for Li3Er0.5Sc0.5Cl6, LiHfFCl4, Na3GdBr3I3, and Li3In0.5Er0.5Cl6 to confirm E_hull < 50 meV/atom.

2. **AIMD screening**: Run 0.5–1 ns AIMD at 600–900 K for top-5 novel candidates to estimate diffusion coefficients and activation energies.

3. **AQVolt26 MLIP validation**: Apply SandboxAQ's published eSEN-AQVolt26 model (Hugging Face) to Li halide candidates for rapid high-temperature MD screening.

4. **Experimental synthesis**: LiHfFCl4 (mechanochemical ball-milling of LiF + HfCl4 + LiCl), Li3Er0.5Sc0.5Cl6 (solid-state, ErCl3 + ScCl3 + LiCl), Na3GdBr3I3 (GdBr3 + NaI + NaBr + Na).

5. **BVS mapping**: Run true bond valence sum pathway calculations (using BVEL/SoftBV tools) on DFT-relaxed structures of top candidates to validate migration pathways and bottleneck radii.

6. **Materials Cloud query**: Check Materials Cloud for DFT entries of hypothetical candidates via the AiiDA-based databases (https://www.materialscloud.org/).

---

## 10. Source References (Key)

1. Jun et al., *Matter* 2025 — LiMXCl4 soft cradle effect: https://ceder.berkeley.edu/publications/2025_KyuJung_Grace_LiMXCL4.pdf
2. Asano et al., *Nature* 2019 — Na2.88Sb0.88W0.12S4 champion: https://www.nature.com/articles/s41467-019-13178-2
3. Li et al., *Angew. Chem.* 2019 — Li3InCl6 water synthesis: https://www.eng.uwo.ca/nanoenergy/publications/2019/pdf/Sun_et_al-2019-Angewandte_Chemie_International_Edition1.pdf
4. Zhao et al., *Nano Lett.* 2020 — Li3Y1-xInxCl6: https://www.eng.uwo.ca/nanoenergy/publications/2020/pdf/Xiaona-2020-Nano-letter.pdf
5. Kim et al., arXiv 2025 — Li3YCl6-xBrx LYCB series: https://arxiv.org/abs/2510.09861
6. Hussain et al., *J. Phys. Chem. C* 2022 — Halospinel Li2M2/3Cl4: https://pubs.acs.org/doi/abs/10.1021/acs.jpcc.2c02511
7. Park et al., *J. Mater. Chem. A* 2022 — Na3MX6 design: https://pubs.rsc.org/en/content/getauthorversionpdf/D2TA05158C
8. Dong et al., *Nat. Commun.* 2024 — HE halide electrolyte: https://pmc.ncbi.nlm.nih.gov/articles/PMC10844219/
9. Schlenker et al., *RSC* 2025 — Argyrodite VdW E_hull: https://pubs.rsc.org/en/content/articlehtml/2025/ta/d4ta06603k
10. Dallakyan et al., USPEX 2026 — Li3MX6 screening: https://uspex-team.org/static/file/Dallakyan2026_JEC_Li3MX6.pdf
11. Kim et al., *PMC* 2018 — Na3SbSe4 cubic: https://pmc.ncbi.nlm.nih.gov/articles/PMC6002371/
12. SandboxAQ/Kim et al., arXiv 2026 — AQVolt26: https://arxiv.org/abs/2604.02524
13. ACS Appl. Mater. Interfaces 2023 — Li2.7In0.7Hf0.3Cl6: https://pubs.acs.org/doi/10.1021/acsami.2c21022
14. RSC J. Mater. Chem. A 2023 — Li3-xY1-xHfxCl6: https://pubs.rsc.org/en/content/articlelanding/2023/ta/d3ta02781c
15. Zhang et al., Pkusam 2025 — Bi-doped Li2ZrCl6: https://www.pkusam.com/uploads/upload/files/20251120/4c2aa6fc6acad72faac0f4679c57bda8.pdf
16. Frontiers Energy Res. 2020 — Na chalcogenide review: https://www.frontiersin.org/journals/energy-research/articles/10.3389/fenrg.2020.00097/full
