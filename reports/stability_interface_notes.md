# Stability & Interface Triage: Notes and Tier Criteria

## Scope and Inputs

This analysis covers all 44 candidates from `candidate_screening_initial.csv` with `migration_proxy_score >= 80` or reported/predicted conductivity > 1 mS/cm. Electrochemical stability windows and interface compatibility assessments are drawn from peer-reviewed publications, DFT grand-potential convex-hull calculations, and experimental post-mortem studies (XPS, AES, AIMD, in situ Raman). Where exact values for a candidate composition are unavailable, chemically justified family-level estimates are used and **explicitly flagged as ESTIMATED**.

---

## Tier Definitions

### Tier 1

A candidate is assigned **Tier 1** when all of the following criteria are satisfied:

1. **Conductivity**: Reported or AIMD-predicted ionic conductivity σ > 1 × 10⁻³ S/cm (1 mS/cm) at room temperature. Experimental confirmation is strongly preferred; AIMD predictions accepted provisionally with a "conditional" flag.
2. **Activation energy**: Reported or computed E_a < 0.30 eV where data exist. Where E_a is unavailable, strong proxy evidence (e.g., flat conductivity–temperature slope, very high AIMD σ) is accepted.
3. **Oxidative limit**: Practical or DFT oxidative stability limit **> 4.0 V vs Li/Li⁺** (or vs Na/Na⁺ for Na conductors), sufficient for 4 V-class layered oxide cathodes.
4. **Cathode compatibility**: Demonstrated or computationally confirmed compatibility with at least one 4 V-class cathode (LiCoO₂, LiNiO₂, or NMC/NCA equivalent). Compatibility may be direct or via a thin catholyte interlayer.
5. **Anode pathway**: A viable anode-facing strategy exists. This may be: (a) direct Li metal compatibility via passivating SEI, (b) graphite anode operation within the reductive window, (c) Li-In or In-Li alloy anode, or (d) a validated separate anode-side electrolyte interlayer (e.g., Li₆PS₅Cl).

**Tier 1 (catholyte role)** denotes materials meeting criteria 1–4 but relegated to the cathode-facing half of a bilayer or trilayer electrolyte architecture due to poor reductive stability. These are still Tier 1 because their σ, E_a, and oxidative window are outstanding and the catholyte architecture is experimentally validated.

### Tier 1A

A candidate is assigned **Tier 1A** when it meets the conductivity criterion (σ > 1 mS/cm or strong proxy evidence) AND either:

- **Tier 1A (oxidative)**: The practical or estimated oxidative stability window exceeds **4.0 V** (or 3.5 V for Na systems with NaCrO₂-class cathodes) AND a viable anode strategy exists (graphite, Li-In, or validated interlayer), but one or more Tier 1 requirements are not fully satisfied (e.g., E_a slightly above 0.30 eV, anode-facing compatibility requires coating, or the oxidative window is only practical/kinetic rather than thermodynamic).
- **Tier 1A (predicted)**: A candidate satisfying all numerical Tier 1 criteria but with conductivity/E_a data exclusively from AIMD (no experimental confirmation). These carry a "conditional" flag pending experimental synthesis.

**Tier 1A candidates require a protective strategy** at one or both interfaces to function in practical cells:
- Sulfide argyrodites (Li₆PS₅X): graphite anode to avoid direct Li metal contact; LiNbO₃ or LIC (Li₃InCl₆) coating on oxide cathodes; or kinetically stable SEI formation (LiF-rich from LiFSI additives, LPSC-LiCl composite).
- Halide Li₃MCl₆-type: LPSC protective interlayer between halide SE and Li metal anode (validated by Zeier group; negligible resistance contribution per Giessen thesis).

### Tier 2

Candidates failing to meet Tier 1A criteria. Common reasons:
- Conductivity below 1 mS/cm (e.g., Na₃SbSe₄ 0.85 mS/cm, Li₃InF₆ 0.55 mS/cm, Li₂ZrCl₆ 0.44 mS/cm).
- Oxidative window insufficient for any practical cathode chemistry (e.g., Na selenide/sulfide family: oxidative limits 2.15–3.0 V vs Na/Na⁺ — too narrow for oxide Na cathodes operating at 2.5–3.6 V vs Na).
- No experimental synthesis confirmed and AIMD σ unvalidated (e.g., Na₂.₈₈Sb₀.₈₈Mo₀.₁₂S₄, Na₇P₃Se₁₁, Na₂Zr₀.₇La₀.₃Cl₆, Na₂HfCl₆).
- Band gap below 2.5 eV creating electronic conductivity risk under applied voltage (e.g., Li₃InI₆ with band gap 2.5 eV and estimated oxidative limit 2.1 V).

---

## Electrochemical Stability Window Evidence Base

### Argyrodite Sulfide Family (Li₆PS₅X, X = Cl, Br, I)

**Thermodynamic windows** (DFT grand-potential phase diagram, Richards/Zhu approach):
- Li₆PS₅Cl: ~1.7–2.3 V vs Li/Li⁺ (thermodynamic). Oxidative decomposition at ~2.5 V (vs In/InLi, offset ~0.6 V → ~3.1 V vs Li) per Juelich Fz Manuscript ([juser.fz-juelich.de/record/907891](https://juser.fz-juelich.de/record/907891/files/Manuscript.pdf)). Reductive decomposition below ~0.6 V vs In/InLi.
- Li₆PS₅Br/I: Similar thermodynamic range 1.7–2.5 V ([repository.tudelft.nl](https://repository.tudelft.nl/file/File_09039068-ba33-48b5-a165-b352c3f2b700?preview=1)).
- General sulfide SE: thermodynamic range 1.7–2.4 V confirmed by multiple DFT studies ([pmc.ncbi.nlm.nih.gov/articles/PMC10910511](https://pmc.ncbi.nlm.nih.gov/articles/PMC10910511/)).

**Practical kinetic window**: Extended to ~4.0–4.2 V due to insulating decomposition products (Li₂S, LiCl, Li₃P at anode; S, LiCl, P₂Sx at cathode) forming passivating SEI/CEI layers. This kinetic extension is widely exploited in practical cells but is NOT thermodynamic stability — decomposition occurs continuously at slow rates.

**Li metal interface**: Li₆PS₅Cl reacts with Li metal forming Li₂S + LiCl + Li₃P (confirmed XPS: [pubs.rsc.org/doi/d5eb00101c](https://pubs.rsc.org/en/content/articlehtml/2025/eb/d5eb00101c)). Interphase is ionically conducting but electronically partially conducting, leading to continuous growth. Mitigation: LiF-rich SEI via LiFSI ([pubs.acs.org/doi/10.1021/acsami.3c14763](https://pubs.acs.org/doi/10.1021/acsami.3c14763)); LPSC-LiCl composite + PC traces; Li₃N coating; In-Li alloy anode.

**LiCoO₂ interface**: Li₆PS₅Cl oxidized at LCO interface into S, Li₂Sn, P₂Sx, phosphates, and LiCl (XPS/AES: [pubs.acs.org/doi/10.1021/acs.chemmater.6b04990](https://pubs.acs.org/doi/abs/10.1021/acs.chemmater.6b04990); in situ Raman: [livrepository.liverpool.ac.uk/3080776](https://livrepository.liverpool.ac.uk/3080776/)). Despite interface reactivity, good capacity retention (300 cycles) observed. LiNbO₃ coating on LCO reduces interfacial resistance.

**LiNiO₂ compatibility**: Worse than LiCoO₂ due to higher oxidation state driving more aggressive sulfide oxidation. Not recommended without coating ([hal.science](https://hal.science/hal-01481291/file/Redox%20activity%20of%20argyrodite%20Li6PS5Cl%20electrolyte.pdf)).

### Halide Li₃MCl₆ and Related Chloride Family

**Thermodynamic windows** (DFT grand potential, Richards approach):
- Li₃YCl₆: 0.62–4.21 V ([eng.uwo.ca Guofeng NE 2021](https://www.eng.uwo.ca/nanoenergy/publications/2021/pdf/Guofeng-NE-2021.pdf))
- Li₃ScCl₆: 0.87–4.21 V (same source)
- Li₃HoCl₆: 0.64–4.25 V; Li₃SmCl₆: 0.67–4.23 V (Western Engineering halide review)
- Li₃InCl₆: 0.62–4.21 V (family DFT); measured 2.68–4.22 V for Hf-doped Li₂.₇In₀.₇Hf₀.₃Cl₆ ([pubs.acs.org/doi/10.1021/acsami.2c21022](https://pubs.acs.org/doi/10.1021/acsami.2c21022))
- Li₂Sm₂/₃Cl₄ (halospinel): band gap = electrochemical window 4.26 V stated explicitly in JPCC 2022 ([pubs.acs.org/doi/10.1021/acs.jpcc.2c02511](https://pubs.acs.org/doi/abs/10.1021/acs.jpcc.2c02511))
- Li₃YBr₆: 0.59–3.15 V (Br anion reduces oxidative limit vs Cl)

**Chlorides exhibit highest oxidation potential (~4.3 V) among Li₃MX₆** — established principle from systematic DFT study ([osti.gov/1810664](https://www.osti.gov/pages/biblio/1810664)).

**Li metal interface**: All halide SEs are thermodynamically unstable vs Li metal. Reduction products: Y/In/Sc/Sm metal + LiCl (confirmed AIMD: [pmc.ncbi.nlm.nih.gov/articles/PMC11503608](https://pmc.ncbi.nlm.nih.gov/articles/PMC11503608/)). Exothermic reaction energy −1.83 eV/f.u. for Li₃YCl₅Br. Electronically conducting decomposition products → continuous interphase growth.

**LPSC interlayer strategy**: Validated solution — Li₆PS₅Cl protective layer between halide SE and Li metal anode. Li₆PS₅Cl is kinetically stable toward lithium and forms good heterocontact with Li₃YCl₆. Demonstrated: stable symmetric cycling >1000h vs 60h without ([pubs.rsc.org/doi/d1ta03042f](https://pubs.rsc.org/en/content/articlelanding/2021/ta/d1ta03042f)). Li₃InCl₆|Li₆PS₅Cl interface: negligible resistance contribution ([jlupub.ub.uni-giessen.de](https://jlupub.ub.uni-giessen.de/items/fc2123c9-f520-4516-91f0-68a88df33c47)).

**LiCoO₂ interface**: Li₃InCl₆ excellent LCO compatibility — 73.6% capacity retention at 5.2 V after 70 cycles; standard 4.2 V operation widely demonstrated ([osti.gov/purl/2572328](https://www.osti.gov/servlets/purl/2572328)). Li₃YCl₆/Li₃InCl₆ family: oxidative stability at 4.5 V confirmed CV ([eng.uwo.ca Guofeng](https://www.eng.uwo.ca/nanoenergy/publications/2021/pdf/Guofeng-NE-2021.pdf)).

**LiNiO₂ compatibility**: DFT chemical reaction energies between candidate halides and LiNiO₂ generally >100 meV/atom (unfavorable; highest among cathode series) per pkusam.com DFT screening ([pkusam.com](https://www.pkusam.com/uploads/upload/files/20251125/866fb79af99b4af17956be921af017dd.pdf)). LFePO₄ most compatible (<70 meV/atom). LiNiO₂ is problematic for all halide electrolytes studied.

### LiMXCl₄ Oxyhalide and Fluorochloride Family

**Thermodynamic windows** (DFT r2SCAN; Jun et al. 2025 Matter, [ceder.berkeley.edu](https://ceder.berkeley.edu/publications/2025_KyuJung_Grace_LiMXCL4.pdf)):
- LiNbOCl₄ (LNOC): total window 1.16 V; oxidative 4.06 V, reductive 2.24–2.90 V
- LiTaOCl₄ (LTOC): total window 1.82 V; oxidative 4.06 V, reductive ~2.24 V
- LiHfFCl₄: total window 2.54 V; oxidative 4.20 V, reductive 1.66 V (improved F-substitution)
- LiZrFCl₄: total window 2.20 V; oxidative 4.06 V, reductive ~1.86 V

All LiMXCl₄ systems oxidative limit ≥ 4.06 V → **compatible with 4 V layered oxide cathodes as catholytes**.
All LiMXCl₄ systems **poor reductive stability** → incompatible with Li metal; separate anode-side SE required.

**Fluoride substitution**: O→F widens window (LiHfFCl₄ achieves 4.20 V oxidative; 1.66 V reductive vs LiHfOCl₄ baseline).

**LCO and NMC experimental validation** (LiNbOCl₄ as catholyte):
- LCO at 4.6 V: 96.7% first-cycle Coulombic efficiency (cold-pressed ASSB; [xnergy.us/linbocl4](https://xnergy.us/linbocl4-lnoc-solid-electrolyte/))
- Ni83 NMC at 4.3 V: 87.8% first-cycle CE

**Decomposition analysis**:
- LTOC reductive products: LiCl + TaCl₄ + LiTa₄O₈. LiTa₄O₈ has E_g ~0.01 eV → electronically conducting → continuous SEI growth risk.
- LNOC reductive products: LiCl + NbCl₂O. NbCl₂O has larger band gap (not insulating, but better than LTOC). LiCl is insulating — partial passivation possible.

### Na Sulfide/Selenide Family (Na₃SbS₄, Na₃PSe₄, Na₃SbSe₄)

**Thermodynamic windows** (DFT grand potential; Ceder group):
- Na₃PS₄: 0.9–2.5 V vs Na/Na⁺
- Na₃PSe₄: 1.25–2.35 V vs Na/Na⁺ (explicitly calculated; [ceder.berkeley.edu compatible](https://ceder.berkeley.edu/publications/2017_yaosen_electrode_electrolyte_compat.pdf))
- Na₃SbS₄: ~1.5–3.0 V vs Na/Na⁺ (wider than Na₃PS₄ due to Sb; cathodic ~1.5 V, anodic ~3.0 V; [sciencedirect.com S0013468624001634](https://www.sciencedirect.com/science/article/abs/pii/S0013468624001634))

**Na halide analogue window** (Na₃₋ₓY₁₋ₓZrₓCl₆ DFT): oxidative limit ~3.8 V vs Na/Na⁺ — significantly wider than Na sulfides ([pmc.ncbi.nlm.nih.gov/articles/PMC7902639](https://pmc.ncbi.nlm.nih.gov/articles/PMC7902639/)).

**Na metal interface** (Na₃SbS₄): Most thermodynamically favorable reaction: 8 Na + Na₃SbS₄ → 4 Na₂S + Na₃Sb. Forms Na₂S + Na₃Sb decomposition products. Kinetic passivation via hydrate coating (Na₃SbS₄·8H₂O surface) partially stabilizes interface per Ceder group ([ceder.berkeley.edu](https://ceder.berkeley.edu/publications/2018_Yaosen_Na_hydrate_coating.pdf)).

**Na₃SbSe₄**: Reports 'good compatibility with metallic Na' and 'excellent electrochemical stability' in RT Na cells ([pmc.ncbi.nlm.nih.gov/articles/PMC6002371](https://pmc.ncbi.nlm.nih.gov/articles/PMC6002371/)). However, DFT window estimated at 1.25–2.35 V from Se-framework analogy — likely kinetically stabilized.

**Cathode implications**: All Na sulfide/selenide candidates have oxidative limits of 2.1–3.0 V vs Na/Na⁺. NaCrO₂ operates at 2.0–3.6 V vs Na/Na⁺ — only partially within window for S-based systems. High-voltage Na oxide cathodes (>3.8 V) incompatible with sulfide/selenide electrolytes.

---

## Interface Compatibility Summary Table

| Family | Li metal | Graphite | LiNiO₂ | LiCoO₂ | Key Reference |
|---|---|---|---|---|---|
| Argyrodite Li₆PS₅X (sulfide) | Poor; SEI (Li₂S+LiCl+Li₃P); kinetic passivation possible with LiF-rich coating | Good (standard) | Poor without coating | Moderate (kinetic ~4.2 V; partial oxidation) | [pubs.rsc.org d5eb00101c](https://pubs.rsc.org/en/content/articlehtml/2025/eb/d5eb00101c) |
| Halide Li₃MCl₆ chloride | Poor direct; LPSC interlayer validated | Moderate with interlayer | Poor (>100 meV/atom rxn energy) | Good (4.2 V confirmed; LIC even at 5.2 V) | [jlupub.ub.uni-giessen.de](https://jlupub.ub.uni-giessen.de/items/fc2123c9-f520-4516-91f0-68a88df33c47) |
| Halide Li₃MBr₆/mixed | Poor direct; same interlayer strategy | Moderate with interlayer | Poor | Moderate (oxidative limit 3.15 V limits LCO) | [pmc.ncbi.nlm.nih.gov PMC11503608](https://pmc.ncbi.nlm.nih.gov/articles/PMC11503608/) |
| LiMXCl₄ oxyhalide (LNOC/LTOC) | Poor (DFT confirmed); LPSC+Li-In anode | Moderate with anode-side SE | Poor (DFT rxn energy) | Excellent (4.6 V LCO confirmed experimental) | [xnergy.us/linbocl4](https://xnergy.us/linbocl4-lnoc-solid-electrolyte/) |
| LiMFCl₄ fluorochloride | Poor (DFT; improved vs LNOC but still reductive limit 1.66 V) | Moderate with interlayer | Moderate (4.20 V oxidative borderline) | Good (4.20 V oxidative within window) | [ceder.berkeley.edu LiMXCl4](https://ceder.berkeley.edu/publications/2025_KyuJung_Grace_LiMXCL4.pdf) |
| Na₃SbS₄ family | Poor (Na₂S+Na₃Sb passivation; hydrate coating improves) | N/A | Poor (not Na cathode) | N/A (Na system) | [ceder.berkeley.edu Na hydrate](https://ceder.berkeley.edu/publications/2018_Yaosen_Na_hydrate_coating.pdf) |
| Na selenide (Na₃SbSe₄, Na₃PSe₄) | Poor-Moderate (some kinetic compatibility reported) | N/A | Poor | N/A (Na system) | [pmc.ncbi.nlm.nih.gov PMC6002371](https://pmc.ncbi.nlm.nih.gov/articles/PMC6002371/) |
| LGPS Li₁₀GeP₂S₁₂ | Poor (severe Ge reduction; Li₁₅Ge₄+Li₂S+Li₃P; 950h possible with LiH₂PO₄ coating) | Poor (below reductive limit; Li-In anode standard) | Poor (thermodynamically) | Poor (degradation at interface confirmed STEM/XPS) | [ceder.berkeley.edu LGPS](https://ceder.berkeley.edu/publications/YiFei_Lithium_Conductor.pdf); [pubs.acs.org doi 10.1021/acsami.8b05132](https://pubs.acs.org/doi/10.1021/acsami.8b05132) |

---

## Tier Assignments Rationale: Key Decisions

### Tier 1 Assignments

**LiNbOCl₄** and **LiTaOCl₄** (catholyte Tier 1):
- σ: 10.7–12.4 mS/cm experimental >> 1 mS/cm
- E_a: 0.18–0.256 eV < 0.30 eV
- Oxidative limit: 4.06 V (DFT, confirmed experimentally at 4.6 V with LCO)
- LCO compatibility: 96.7% CE at 4.6 V — direct experimental validation
- Anode: catholyte role mitigates reductive instability; LPSC+Li-In anode-side strategy validated by Tanaka et al.

**Li₂Sm₂/₃Cl₄** (Tier 1 conditional — no experimental synthesis):
- σ: 15.3 mS/cm (AIMD) >> 1 mS/cm
- E_a: 0.195 eV < 0.20 eV (lowest in halospinel family)
- Oxidative limit: 4.26 V (DFT; band gap as proxy)
- Halospinel 3D isotropic Li diffusion — structural advantage
- LCO: compatible per chloride family
- **Flagged**: No experimental synthesis; all σ/E_a from AIMD only

**Li₃InCl₆** (catholyte Tier 1):
- σ: 1.0–4.0 mS/cm experimental (standard to wet-chemistry routes)
- Oxidative limit: 4.21–4.22 V (DFT + measured for Hf-doped variant)
- LCO: 73.6% retention at 5.2 V; standard 4.2 V widely validated
- Anode: LPSC interlayer validated with negligible resistance addition
- Li metal: poor direct, but LPSC interlayer solution industrially applicable

### Tier 1A Decisions (Key Cases)

**Li₆.₆P₀.₄Ge₀.₆S₅I and Li₆.₂₅P₀.₇₅Si₀.₂₅S₅I** (argyrodite Ge/Si-substituted, I-anion):
- σ: 2.9 mS/cm and 1.7 mS/cm experimental → meet conductivity threshold
- Oxidative: practical kinetic window ~4.2 V (not thermodynamic; flagged)
- Anode: graphite standard strategy; Li-In alternative
- Limitation: thermodynamic window only 1.7–2.5 V; requires ongoing SEI management

**Li₂.₈In₀.₂Sc₀.₂Yb₀.₂Lu₀.₂Zr₀.₂Cl₆** (high-entropy chloride):
- σ: 2.13 mS/cm experimental → above threshold
- Oxidative: ~4.3 V estimated from Li₃MCl₆ family (not measured for this composition)
- Anode: LPSC interlayer strategy from halide family
- Limitation: No direct electrochemical window measurement; E_a not reported

**Li₂.₇In₀.₇Hf₀.₃Cl₆**:
- σ: 1.28 mS/cm experimental → above threshold
- Oxidative: 4.22 V measured by CV — **direct experimental evidence**
- Anode: halide family LPSC interlayer
- Strongest halide candidate below Tier 1 (higher confidence than purely predicted materials)

**LiHfFCl₄ and LiZrFCl₄** (Tier 1A; fluorochloride predictions):
- σ: 10–137.5 mS/cm AIMD only; no experimental synthesis
- Oxidative: 4.06–4.20 V DFT
- Reductive: 1.66–2.20 V (wider than LNOC/LTOC; Zr/Hf d⁰ stability)
- Tier 1A rather than Tier 1 because: experimental synthesis not confirmed; all σ from AIMD

### Tier 2 Exclusions (Key Cases)

**Na₂.₈₈Sb₀.₈₈Mo₀.₁₂S₄** and all Na sulfide/selenide family except high-conductivity S variants:
- Even where σ > 1 mS/cm (Na₃SbS₃.₇₅Se₀.₂₅: 4.03 mS/cm), the oxidative limit ~3.0 V vs Na/Na⁺ is insufficient for high-voltage Na oxide cathodes
- Suitable only for NaCrO₂-class (2–3.6 V) or conversion cathodes — narrow application window
- Practical difficulty: Na metal anode instability not resolved

**Li₃InI₆** (Tier 2 despite AIMD σ 2.18 mS/cm):
- Band gap 2.5 eV (below 3.5 eV threshold) → electronic leakage risk
- Estimated oxidative limit ~2.1 V (iodide family same window as sulfide)
- Cannot be paired with any 4 V cathode

**Li₂ZrCl₆, Li₂HfCl₆** (Tier 2):
- σ below 1 mS/cm threshold (0.44–0.64 mS/cm)
- Bi-doped variant (Li₂ZrBiCl₆) raises to Tier 1A based on AIMD prediction, but synthesis unconfirmed

---

## Evidence Hierarchy

1. **Direct experimental measurement** (preferred): CV, EIS, XPS post-mortem, STEM, in situ Raman — highest confidence.
2. **DFT grand-potential phase diagram** (Richards/Zhu approach): Thermodynamic limits only; kinetic extensions are not predicted. r2SCAN meta-GGA preferred over PBE.
3. **AIMD conductivity/E_a**: Reliable at high temperatures extrapolated to RT; ~factor 2–5 uncertainty in absolute conductivity. Accepted for provisional Tier 1/1A with "conditional" flag.
4. **Family-level chemical analogy**: Applied where composition has no direct data; flagged ESTIMATED. Justified by systematic trends within sulfide, halide chloride, and halide bromide/iodide families.
5. **No DFT calculations were fabricated**: All DFT results cited here reference specific publications with accessible DOIs/URLs. Pymatgen calculations were not run for this triage.

---

## Source Citations

Key references underpinning this analysis:

- Argyrodite family stability: [pubs.rsc.org/en/content/articlehtml/2025/eb/d5eb00101c](https://pubs.rsc.org/en/content/articlehtml/2025/eb/d5eb00101c)
- Argyrodite LCO/NMC/LMO interface XPS/AES: [pubs.acs.org/doi/10.1021/acs.chemmater.6b04990](https://pubs.acs.org/doi/abs/10.1021/acs.chemmater.6b04990)
- Argyrodite in situ Raman vs LCO and Li metal: [livrepository.liverpool.ac.uk/3080776](https://livrepository.liverpool.ac.uk/3080776/)
- Li₆PS₅Cl oxidative stability (Juelich; Ge/Si substitutions): [juser.fz-juelich.de/record/907891](https://juser.fz-juelich.de/record/907891/files/Manuscript.pdf)
- Li halide family DFT windows (Li₃YCl₆/Li₃ScCl₆/Li₃InCl₆): [eng.uwo.ca Guofeng NE 2021](https://www.eng.uwo.ca/nanoenergy/publications/2021/pdf/Guofeng-NE-2021.pdf)
- Li halide material design strategy (OSTI/ACS): [pubs.acs.org/doi/10.1021/acs.chemmater.1c00555](https://pubs.acs.org/doi/10.1021/acs.chemmater.1c00555)
- Li halide vs Li metal (LPSC interlayer): [pubs.rsc.org/doi/d1ta03042f](https://pubs.rsc.org/en/content/articlelanding/2021/ta/d1ta03042f)
- Li₃InCl₆ Li metal instability (XPS); LPSC interlayer zero-resistance: [jlupub.ub.uni-giessen.de](https://jlupub.ub.uni-giessen.de/items/fc2123c9-f520-4516-91f0-68a88df33c47)
- Li₃InCl₆ high-voltage LCO 5.2 V: [osti.gov/purl/2572328](https://www.osti.gov/servlets/purl/2572328)
- Li₂.₇In₀.₇Hf₀.₃Cl₆ window 2.68–4.22 V measured: [pubs.acs.org/doi/10.1021/acsami.2c21022](https://pubs.acs.org/doi/10.1021/acsami.2c21022)
- LiMXCl₄ family electrochemical windows (Jun 2025 Matter): [ceder.berkeley.edu LiMXCl4](https://ceder.berkeley.edu/publications/2025_KyuJung_Grace_LiMXCL4.pdf)
- LiNbOCl₄ LCO 4.6 V validation: [xnergy.us/linbocl4](https://xnergy.us/linbocl4-lnoc-solid-electrolyte/)
- Li₂Sm₂/₃Cl₄ halospinel DFT: [pubs.acs.org/doi/10.1021/acs.jpcc.2c02511](https://pubs.acs.org/doi/abs/10.1021/acs.jpcc.2c02511)
- Li₃YCl₅Br AIMD interface (NMC and Li metal): [pmc.ncbi.nlm.nih.gov/PMC11503608](https://pmc.ncbi.nlm.nih.gov/articles/PMC11503608/)
- Na halide family (Na₃₋ₓY₁₋ₓZrₓCl₆ DFT window ~3.8 V; NaCrO₂ compatibility): [pmc.ncbi.nlm.nih.gov/PMC7902639](https://pmc.ncbi.nlm.nih.gov/articles/PMC7902639/)
- Na₃PSe₄/Na₃PS₄ DFT windows 1.25–2.35 V / 0.9–2.5 V: [ceder.berkeley.edu compatible](https://ceder.berkeley.edu/publications/2017_yaosen_electrode_electrolyte_compat.pdf)
- Na₃SbS₄ Na metal interface (DFT; hydrate coating strategy): [ceder.berkeley.edu Na hydrate](https://ceder.berkeley.edu/publications/2018_Yaosen_Na_hydrate_coating.pdf)
- Na₃SbSe₄ cells cycling and 'excellent electrochemical stability': [pmc.ncbi.nlm.nih.gov/PMC6002371](https://pmc.ncbi.nlm.nih.gov/articles/PMC6002371/)
- LGPS electrochemical window DFT 1.19–2.38 V: [ceder.berkeley.edu LGPS](https://ceder.berkeley.edu/publications/YiFei_Lithium_Conductor.pdf)
- LGPS LCO interface degradation (STEM/XPS): [pubs.acs.org/doi/10.1021/acsami.8b05132](https://pubs.acs.org/doi/10.1021/acsami.8b05132)
- Halide DFT chemical stability vs LiNiO₂ (>100 meV/atom): [pkusam.com DFT screening](https://www.pkusam.com/uploads/upload/files/20251125/866fb79af99b4af17956be921af017dd.pdf)
- Grand potential approach (Nature Review): [ceder.berkeley.edu interface review](https://ceder.berkeley.edu/publications/2019_xiao_nature_review.pdf)
- Li₃InI₆ computed (Dallakyan 2026 USPEX): [uspex-team.org](https://uspex-team.org/static/file/Dallakyan2026_JEC_Li3MX6.pdf)
- Li₃YCl₁.₅Br₄.₅ peak conductivity 5.36 mS/cm experimental: [arxiv.org/abs/2510.09861](https://arxiv.org/abs/2510.09861)
- Na₃GdBr₃I₃ / Na₃YBr₃I₃ AIMD predictions: [pubs.rsc.org D2TA05158C](https://pubs.rsc.org/en/content/getauthorversionpdf/D2TA05158C)
