# Code and figure provenance — FBLN1 1C→1D isoform switching in CAFs

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24055128.svg)](https://doi.org/10.5281/zenodo.24055128)

Archived on Zenodo. The concept DOI [10.5281/zenodo.24055128](https://doi.org/10.5281/zenodo.24055128) always resolves to the newest version. **This is v1.1.0**, [10.5281/zenodo.24084490](https://doi.org/10.5281/zenodo.24084490), which corrects the two statistical errors described under "Changes in v1.1.0" below; the superseded v1.0.0 snapshot is [10.5281/zenodo.24055129](https://doi.org/10.5281/zenodo.24055129). Cite the version DOI of the release you actually used rather than a branch, so the snapshot stays fixed.

This repository contains the analysis scripts that generate every figure in the
manuscript, together with the derived tables those scripts read. Raw sequencing data
are **not** included; they are public and listed in `DATA_SOURCES.md`.

---

## Quick start

```bash
git clone https://github.com/xuwenguang-nju/fbln1-caf-isoform-switch.git
cd fbln1-caf-isoform-switch       # the folder containing this README
pip install -r requirements.txt
python analysis/scripts/63_fig1_depth_gate.py       # reproduces Figure 2
python analysis/scripts/482_smoke_test_release.py   # reproduces all 12, about a minute
```

Each figure script writes its PNG into `figures/`; no step contacts a remote host. Two
by-products also appear and are ignored by git: `69_` and `204_` drop panel crops into
`.cache/figure-style/`, and `87_` writes `data/processed/three_arm_power.json`. Nothing
reads them back, and no script overwrites any shipped input.

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
| `results/tables/candidates_gated.tsv` | server pipeline (depth gate, pre-fixed-referent; panel f) |
| `results/tables/switch_test_full.tsv` | server pipeline (genome-wide `dIF`, pre-fixed-referent; panel f note) |

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
| `remote/easyai/donor_OCSCC1_Post_summary.json` | server pipeline (lineage composition, panel c) |
| `remote/easyai/donor_OCSCC3_Post_summary.json` | server pipeline (lineage composition, panel c) |
| `remote/easyai/donor_OCSCC13_Post_summary.json` | server pipeline (lineage composition, panel c) |
| `remote/easyai/donor_OCSCC36_Post_summary.json` | server pipeline (lineage composition, panel c) |

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
| `remote/easyai/cross_donor_control.json` | server pipeline (panel b) |
| `remote/easyai/mapq_compare.json` | server pipeline (panel d) |
| `remote/easyai/doublets.json` | server pipeline (panel e) |
| `remote/easyai/scrublet_effect.json` | server pipeline (panel e annotation) |
| `remote/easyai/epi_heterogeneity.json` | server pipeline (panel f) |

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

The plotting helpers every script imports from `_figstyle.py` (`panel`,
`fit_panel_titles`) are vendored next to it for the same reason: the release has to run
without the original machine's skills directory.

That single line is the **only** change to the released scripts. Diffing every released
script against the working copy it came from shows one differing line per figure script
and no difference at all in the 29 upstream scripts, so line numbers and logic match the
versions that produced the figures.

---

## Verification status

**12 of 12 figure scripts reproduce.** Run

```bash
python analysis/scripts/482_smoke_test_release.py
```

from the release root. It runs each figure script in a subprocess with a fresh output
path — any stale PNG is deleted first — and passes it only on a clean exit plus a PNG in
`figures/`, so nothing can pass on a leftover file. A subset can be checked by label,
for example `... 482_smoke_test_release.py 3 9 S2`.

Verified on Python 3.11.17 with numpy 2.4.6, pandas 3.0.6, matplotlib 3.11.2,
scipy 1.17.1 and pillow 13.0.0 (Windows). Every regenerated PNG was **byte-identical**
(SHA-256) to the file the manuscript figure was exported from, so the shipped tables and
scripts reproduce the published panels exactly, not merely approximately. Exact PNG bytes
depend on the matplotlib version, so treat the exit status as the portable criterion and
byte-identity as evidence that the shipped inputs are the published ones.

Two things to keep in mind when reading the code:

- **Figure 3 panel f reads the pre-fixed-referent depth gate**
  (`results/tables/candidates_gated.tsv`, `results/tables/switch_test_full.tsv`), because
  that is the layer the published panel was drawn from. The manuscript's quantitative
  statements use the fixed-referent tables (`*_fixed_referent.*`), which Figures 1, 2 and
  S1 read. Both layers ship here; `181_gate_under_fixed_referent.py` is the script that
  bridges them and prints how far the two definitions agree at the gene-call level.
- **Figures 9 and S2 read small JSON summaries** written by the server-side `psi_diff`
  pipeline (`remote/easyai/global_tests.json`, `nd_donor_summary.json`,
  `donor_*_summary.json` and the five control summaries). They ship as-is so those two
  figures run without a compute host.

### Changes since the first upload (commit `42021da`)

- Added the inputs three figure scripts were missing — 12 files, all listed in
  `_inputs_manifest.json`: `results/tables/candidates_gated.tsv`,
  `results/tables/switch_test_full.tsv`, `remote/easyai/nd_donor_summary.json`, the four
  `remote/easyai/donor_*_summary.json` lineage summaries, and the five control summaries
  `mapq_compare.json`, `doublets.json`, `scrublet_effect.json`, `epi_heterogeneity.json`,
  `cross_donor_control.json`. Figures 1, 2, 4–8, S1 and S3 had been reproducing already.
- Added `analysis/scripts/482_smoke_test_release.py`. The previous README pointed at it,
  but it was not in the repository.
- `requirements.txt` now lists `pillow`, which every figure script imports to re-open its
  own PNG, and names `statsmodels`, `scanpy` and `anndata` as upstream-only extras.
- `_inputs_manifest.json` byte counts now match the shipped files; two of them had
  recorded the size of the un-redacted local copy rather than the released one.

No figure script was edited: the single `STYLE = ...` line described above is still the
only difference from the scripts that produced the figures.

### Changes in v1.1.0 (scientific review corrections)

A review of the analysis found two statistical errors. Both are corrected here, and the
figures and the manuscript were regenerated from the corrected code. The numbers that
change are listed so a reader can see exactly what moved.

- **Intercept inference used the slope's uncertainty.** `150b_pancancer_fbln1.py` and its
  TCGA counterpart (not shipped here) took `linregress(...).pvalue` and `.stderr` for the
  *intercept* p-value and standard error. Those attributes belong to the **slope**, so every
  reported composition-adjusted intercept carried the wrong interval and p-value. The bare
  intercept point estimate was unaffected. `150b_pancancer_fbln1.py` now uses
  `intercept_stderr` and its own t statistic on `n - 2` df.
- **The 50-nt NMD rule was applied without direction.** `119_nmd_prediction.py` tested
  `abs(stop - last_junction) > 50`, which cannot separate an exon-junction-complex substrate
  (stop *upstream* of the last junction) from the normal configuration (stop *downstream*,
  inside the terminal exon). Both FBLN1 stops lie downstream by 137 nt (1D) and 352 nt (1C),
  so **neither isoform is an NMD substrate**; the earlier flag marked both as substrates.
  The corrected output ships as `data/processed/fbln1_nmd.json`, which Figure 5 panel E
  reads, so the panel and the JSON now agree.

Figure 5 panel E was redrawn from the corrected classification, and the regenerated PNG is
byte-identical to the figure in the revised manuscript. `119_nmd_prediction.py` needs the
GENCODE GTF and runs on the analysis server; the JSON it produced there was corrected and
is shipped so Figure 5 reproduces without a compute host.

Two things were checked and found already correct, so nothing changed: the depth-gate
tolerance is 0.04 in both the Methods and the Figure 2 legend, and the terminal-exon 3'UTR
lengths printed in panel E (692 nt for 1D, 99 nt for 1C) match the exon coordinates.

---

## What is deliberately not here

- **Raw sequencing data** — public; see `DATA_SOURCES.md`. The two full-depth
  quantification matrices are ~90 MB each and are regenerated by the server pipeline.
- **One-off helper scripts** — the working directory they came from contains several
  hundred scripts for figure iteration, patching and verification. Only the 12 figure
  scripts and the 29 scripts that produce their inputs are included.
- **Intermediate tables with no consumer.** The one superseded layer that *is* here —
  the pre-fixed-referent gate, `results/tables/candidates_gated.tsv` and
  `switch_test_full.tsv` — is included only because Figure 3 panel f reads it. The
  manuscript documents the fixed-referent definition as the corrected one, and every
  reported number comes from those tables. Nothing else from the several hundred
  intermediate files in the working directory is included.

---

## Licence

MIT Licence — see `LICENSE`. Copyright (c) 2026 Wei Han.

You are free to use, modify and redistribute this code, including commercially,
provided the copyright notice is retained. If you use these scripts or the
accompanying tables in work of your own, please cite the paper.
