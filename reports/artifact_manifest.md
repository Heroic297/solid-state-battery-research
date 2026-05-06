# Solid-State Electrolyte Discovery Artifact Manifest

## Main report

- `solid-state-electrolyte-discovery-brief.pplx.md`: cited scientific brief with ranked candidates, mechanistic interpretation, synthesis briefs, limitations, and pre-print outline.

## Candidate and evidence tables

- `final_ranked_candidates.csv`: merged final ranking with conductivity basis, stability windows, novelty score, tier, interface notes, and URLs.
- `top15_candidate_summary.csv`: compact top-15 subset for quick review.
- `candidate_screening_initial.csv`: original 44-candidate descriptor-screening table.
- `literature_matrix.csv`: literature matrix by compound/family.
- `stability_interface_analysis.csv`: stability and interface triage for all 44 candidates.
- `novelty_checks.csv`: prior-art check results.

## Notes

- `candidate_generation_notes.md`: candidate-generation methodology and data-source limitations.
- `literature_synthesis.md`: literature synthesis.
- `stability_interface_notes.md`: stability-tier definitions and family-level evidence.
- `novelty_notes.md`: novelty assessment summary.

## MLIP-MD workflow

- `mlip_md_pipeline/`: runnable ASE + MACE/CHGNet/SevenNet-compatible pipeline for local-GPU MD screening.
- `tier1a_mlip_input_template.csv`: input template for Tier 1A novel candidates once structure files are supplied.
- `tier1a_structure_manifest.csv`: structure-file status and recommended next action for Tier 1A novel candidates.

## Plots

- `candidate_conductivity_screen.png`: top candidates versus the \(10^{-3}\) S/cm conductivity target.
- `novelty_vs_migration_proxy.png`: novelty versus fast-ion structural proxy.

## Important caveat

Production MLIP-MD trajectories were not run in the cloud environment. The MD workflow was validated in dummy mode only; candidate conductivities in the final ranking are labeled as experimental/literature, AIMD/literature prediction, or family/proxy estimate.
