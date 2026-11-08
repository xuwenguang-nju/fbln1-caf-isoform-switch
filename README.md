# Code and figure provenance — FBLN1 1C→1D isoform switching in CAFs

This repository contains the analysis scripts that generate every figure in the
manuscript, together with the derived tables those scripts read. Raw sequencing data
are **not** included; they are public and listed in `DATA_SOURCES.md`.

---

## Quick start

```bash
pip install -r requirements.txt
cd figures_release            # the folder containing this README
python analysis/scripts/63_fig1_depth_gate.py     # reproduces Figure 2
```

Each figure script writes its PNG into `figures/`.

**Path convention.** Every script computes its locations as

```python
ROOT = Path(__file__).resolve().parents[2]
```

so the scripts must stay at `<repository-root>/analysis/scripts/<name>.py`.
The directory layout is therefore not cosmetic — it is what makes the scripts run
unmodified. Do not flatten it.

---

## Which script draws which figure

| Figure | Script | Output |
|---|---|---|
| Figure 1 | `87_fig1_design_and_power.py` | `fig_three_arm_design.png` |
| Figure 2 | `63_fig1_depth_gate.py` | `fig_depth_gate.png` |
| Figure 3 | `67_fig2_fbln1_isoforms.py` | `fig_fbln1_paired.png` |
| Figure 4 | `126_fig_s1_singlecell.py` | `supp_fig_s1_singlecell_fbln1.png` |
| Figure 5 | `71_fig3_domains.py` | `fig_fbln1_domains.png` |
| Figure 6 | `69_fig4_tcga.py` | `fig_fbln1_tcga.png` |
| Figure 7 | `204_fig_pancancer.py` | `fig_pancancer.png` |
| Figure 8 | `70_fig5_cross_system.py` | `fig_fbln1_cross_system.png` |
| Figure 9 | `72_fig6_singlecell.py` | `fig_singlecell_locus.png` |
| Figure S1 | `138_fig_s2_ledger.py` | `supp_fig_s2_ledger.png` |
| Figure S2 | `136_fig_s3_controls.py` | `supp_fig_s3_controls.png` |
| Figure S3 | `140_fig_s5_positioning.py` | `supp_fig_s5_positioning.png` |

---

## Figure → input files → producing script

Paths are relative to the repository root. "Produced by" names the script in this
repository that writes the input; where that step ran on a compute server, the server
pipeline is in `analysis/scripts/remote_step1/`, `c_arm/` and `replication_caf/`.

### Figure 1 — three-arm design and power
**Script** `87_fig1_design_and_power.py`
| Input | Produced by |
|---|---|
| `data/processed/fig1_arm_b_numbers.json` | `186_fig1_arm_b_numbers.py` |
| `results/tables/caf_power_measured_fixed_referent.csv` | `191_power_from_fixed_referent.py` |

### Figure 2 — the depth gate
**Script** `63_fig1_depth_gate.py`
| Input | Produced by |
|---|---|
| `results/tables/switch_test_full_fixed_referent.tsv` | `181_gate_under_fixed_referent.py` |
| `results/tables/switch_test_sub10_fixed_referent.tsv` | `181_gate_under_fixed_referent.py` |
| `results/tables/candidates_gated_fixed_referent.tsv` | `181_gate_under_fixed_referent.py` |
| `results/tables/caf_switch_threshold_sweep_fixed_referent.csv` | `183_sweep_fixed_referent.py` |
| `results/tables/caf_panel_depth_consistency_fixed_referent.csv` | `187_panel_depth_consistency.py` |
| `results/tables/caf_isoform_noise_by_depth.csv` | `38_isoform_noise.py` |
| `results/tables/caf_depth_requirement.csv` | `39_power_from_measured.py` |
| `results/tables/isoform_identifiability_genes.tsv` | `92_isoform_identifiability.py` |
| `data/processed/fig1_arm_b_numbers.json` | `186_fig1_arm_b_numbers.py` |

### Figure 3 — the paired FBLN1 switch
**Script** `67_fig2_fbln1_isoforms.py`
| Input | Produced by |
|---|---|
| `results/tables/fbln1_C_D_ratio_by_patient.csv` | server pipeline (`dIF` per patient) |
| `results/tables/caf_panel_test_deep.csv` | server pipeline (panel-deep tests) |
| `results/tables/fbln1_isoform_shares.tsv` | server pipeline (per-library shares) |
| `data/processed/isoform_fracs_panel_deep/*.tsv` (36 files) | server pipeline (k-mer quantitation) |

### Figure 4 — single-cell validation
**Script** `126_fig_s1_singlecell.py`
| Input | Produced by |
|---|---|
| `remote/easyai/sc_percell.tsv` | server pipeline (per-cell isoform counts) |
| `remote/easyai/sc_qc.tsv` | server pipeline |
| `remote/easyai/sc_marker_dot.tsv` | `129_marker_dot.py` |
| `data/processed/fbln1_singlecell_validation.json` | `125_singlecell_fbln1_validation.py` |

### Figure 5 — protein consequence
**Script** `71_fig3_domains.py`
| Input | Produced by |
|---|---|
| `data/processed/fbln1_structure.json` | `120_physchem.py` |
| `data/processed/fbln1_sequences.json` | `120_physchem.py` |
| `data/processed/fbln1_nmd.json` | `119_nmd_prediction.py` |
| `data/processed/fbln1_physchem.json` | `120_physchem.py` |

### Figure 6 — TCGA-HNSC
**Script** `69_fig4_tcga.py`
| Input | Produced by |
|---|---|
| `results/tables/fbln1_tcga_events.csv` | `47_fbln1_tcga.py` |
| `results/tables/tcga_caf_events_tvn.csv` | `45_tcga_bulk_validation.py` |
| `results/tables/tcga_caf_events_cox.csv` | `42_tcga_caf_events.py` / `44_tcga_caf_events_full.py` |
| `results/tables/fbln1_tcga_composition_by_cohort.csv` | `151_fbln1_composition_by_cohort.py` |
| `results/tables/c_arm_scoring_variants.tsv` | `c_arm/02_scoring_variants.py` |
| `results/tables/c_arm_literature_panels.tsv` | `c_arm/03_literature_markers.py` |
| `results/tables/c_arm_splicing_factors_response.tsv` | `c_arm/01_c_arm_signature.py` |
| `results/tables/tcga_hnsc_site_breakdown.csv` | `27_arm_sample_ledger.py` |
| `results/tables/tcga_hnsc_oral_power.csv` | `27_arm_sample_ledger.py` |
| `data/processed/fbln1_tcga.json` | `47_fbln1_tcga.py` |

### Figure 7 — pan-cancer
**Script** `204_fig_pancancer.py`
| Input | Produced by |
|---|---|
| `results/tables/pancancer_fbln1_events.csv` | `150b_pancancer_fbln1.py` |
| `results/tables/pancancer_caf_fraction.csv` | `150b_pancancer_fbln1.py` |
| `data/processed/pancancer_fbln1.json` | `150b_pancancer_fbln1.py` |

### Figure 8 — sorted-CAF cross-system test
**Script** `70_fig5_cross_system.py`
| Input | Produced by |
|---|---|
| `results/tables/caf_fbln1_isoform.csv` | `replication_caf/analyze_caf.py` |

### Figure 9 — contrast specificity
**Script** `72_fig6_singlecell.py`
| Input | Produced by |
|---|---|
| `remote/easyai/stringent_events.tsv` | server pipeline (`psi_diff`) |
| `remote/easyai/global_tests.json` | server pipeline |
| `remote/easyai/nd_donor_events.tsv` | server pipeline |
| `remote/easyai/nd_donor_summary.json` | server pipeline |
| `remote/easyai/caf_q0_assignment.csv` | server pipeline |

### Figure S1 — cohort and power ledger
**Script** `138_fig_s2_ledger.py`
| Input | Produced by |
|---|---|
| `results/tables/arm_sample_ledger.csv` | `27_arm_sample_ledger.py` |
| `results/tables/tcga_hnsc_site_ledger.csv` | `27_arm_sample_ledger.py` |
| `results/tables/sc_dataset_feasibility.csv` | `11_sc_feasibility_tables.py` |
| `results/tables/sc_competition_tiers.csv` | `11_sc_feasibility_tables.py` |
| `results/tables/sc_competition_sweep.csv` | `09_sc_competition_sweep.py` |
| `results/tables/as_archetype_overview.csv` | `07_archetype_table.py` |
| `results/tables/as_bioinformatics_exemplars.csv` | `08_citation_labels.py` |
| `results/tables/icb_power_9v19.csv` | `30_icb_power.py` |
| `results/tables/tcga_hnsc_psi_event_depth.csv` | server pipeline |
| `results/tables/caf_power_measured_fixed_referent.csv` | `191_power_from_fixed_referent.py` |

### Figure S2 — controls
**Script** `136_fig_s3_controls.py`
| Input | Produced by |
|---|---|
| `remote/easyai/global_tests.json` | server pipeline |
| `remote/easyai/nd_donor_events.tsv` | server pipeline |
| `remote/easyai/stringent_events.tsv` | server pipeline |

### Figure S3 — positioning
**Script** `140_fig_s5_positioning.py`
| Input | Produced by |
|---|---|
| `results/tables/sc_competition_tiers.csv` | `11_sc_feasibility_tables.py` |
| `results/tables/sc_competition_sweep.csv` | `09_sc_competition_sweep.py` |
| `results/tables/as_archetype_overview.csv` | `07_archetype_table.py` |
| `results/tables/as_bioinformatics_exemplars.csv` | `08_citation_labels.py` |

---

## One deliberate deviation from the original scripts

All twelve figure scripts loaded their plotting style from an **absolute path on the
machine that produced the figures** — a path inside a separate skills directory that
resolves nowhere else:

```python
STYLE = Path(r"<absolute-path-on-the-original-machine>/runtime.py")
exec(compile(STYLE.read_text(encoding="utf-8"), "runtime.py", "exec"), globals())
```

That path resolves nowhere else. The runtime is vendored here as
`analysis/scripts/_style_runtime.py` and each script's single `STYLE = ...` line is
rewritten to:

```python
STYLE = Path(__file__).resolve().parent / "_style_runtime.py"
```

This is the **only** change to the released scripts. Nothing else was edited, so line
numbers and logic match the versions that produced the figures.

---

## Verification status — please read

The release was tested by running every figure script from the release directory and
checking that it exits cleanly and writes its PNG:

**9 of 12 reproduce without intervention.** Figures 1, 2, 4, 5, 6, 7, 8, S1 and S3
run end to end and regenerate their panels from the included tables.

**3 do not, and the reason is known:**

| Figure | Missing input | Nature of the gap |
|---|---|---|
| 3 | `results/tables/candidates_gated.tsv` | an intermediate table from the depth-gate chain that was not carried into this release |
| 9 | `remote/easyai/nd_donor_summary.json` | server-side output of the `psi_diff` pipeline |
| S2 | `remote/easyai/global_tests.json` reads empty | one of its inputs is present but yields no rows, so an axis limit cannot be computed |

These three are **not** silently broken: the scripts fail with a clear
`FileNotFoundError` or `ValueError` at the first missing piece. To close them, either
add the named files, or re-run the server pipeline in `analysis/scripts/remote_step1/`
against the raw data listed in `DATA_SOURCES.md`.

A helper is included to check the current state at any time:

```bash
python analysis/scripts/482_smoke_test_release.py
```

---

## What is deliberately not here

- **Raw sequencing data** — public; see `DATA_SOURCES.md`. The two full-depth
  quantification matrices are ~90 MB each and are regenerated by the server pipeline.
- **One-off helper scripts** — the working directory they came from contains several
  hundred scripts for figure iteration, patching and verification. Only the 12 figure
  scripts and the 29 scripts that produce their inputs are included.
- **Intermediate tables with no consumer** — including superseded quantitation
  (dominant-isoform `dIF`), which is documented in the manuscript as a corrected
  definition. Including it would invite recomputation of numbers that differ from the
  published ones.

---

## Licence

MIT Licence — see `LICENSE`. Copyright (c) 2026 Wei Han.

You are free to use, modify and redistribute this code, including commercially,
provided the copyright notice is retained. If you use these scripts or the
accompanying tables in work of your own, please cite the paper.
