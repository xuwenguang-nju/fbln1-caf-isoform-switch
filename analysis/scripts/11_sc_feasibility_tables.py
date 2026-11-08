"""Step 11: build the three decision tables for the oral-cancer single-cell splicing plan.

  results/tables/sc_competition_tiers.csv  — evidence for each specificity tier
  results/tables/sc_dataset_feasibility.csv — public data assets and what they support
  results/tables/sc_tool_stack.csv          — candidate methods and their input requirement

Every count comes from data/processed/sc_competition_counts.json or the Europe PMC
title/abstract checks recorded in the report; dataset facts were read from the GEO
record and, where noted, from the original paper's full text.
"""
import csv, json

counts = json.load(open("data/processed/sc_competition_counts.json", encoding="utf-8"))
trend = json.load(open("data/processed/sc_trends.json", encoding="utf-8"))

TIERS = [
    ("T1", "single-cell + alternative splicing (any disease)", counts["sc_as_all"]["count"],
     "PubMed TIAB", "field exists and is growing ~2x/yr"),
    ("T2", "... restricted to cancer", counts["sc_as_cancer"]["count"],
     "PubMed TIAB", "glioma atlas (Genome Biol 2025) is the closest template"),
    ("T3", "... in the tumour microenvironment", counts["tme_splicing"]["count"],
     "PubMed TIAB", "NAR 2025 ICB splicing map; HCC T-cell splicing (Cancer Sci 2026)"),
    ("T4", "... isoform switch / usage at single-cell level", counts["isoform_switch_sc"]["count"],
     "PubMed TIAB", "tool-driven, few tumour applications"),
    ("T5", "... in oral / head-and-neck tissue", counts["sc_as_hnsc_or_oral"]["count"],
     "PubMed TIAB", "HPV-integration splicing (CCR 2026); no cell-type splicing map"),
    ("T6", "oral/HNSCC single-cell + alternative splicing (title/abstract)", 1,
     "Europe PMC TITLE_ABS", "one 2026 Sci Bull review-style omics portrait; no atlas"),
    ("T7", "CAF + alternative splicing / isoform (any cancer, all designs)", 4,
     "Europe PMC TITLE_ABS", "FGFR2c-CAF (Biology 2023) is the only mechanistic isoform paper"),
    ("T8", "macrophage + splicing / isoform + tumour", 54,
     "Europe PMC TITLE_ABS", "MS4A7, CD47, neuropilin isoforms already published"),
]

DATASETS = [
    # gse, study, tissue/scope, n_cells_or_samples, chemistry, access, splice_fit, role
    ("GSE103322", "Puram 2017 Cell", "OSCC primary + matched lymph-node metastasis, 18 patients",
     "5901 cells", "Smart-seq2, full-length (verified in paper methods)", "GEO open",
     "HIGH: full-length, median 1.34M reads/cell",
     "primary discovery set for CAF / macrophage / malignant states"),
    ("HRA001006", "JCI 2022 (TDO2+ myofibroblasts)", "13 OSCC + 3 oral leukoplakia + 8 adjacent normal",
     "131702 cells", "10x 3' (14 samples) + 10x 5' with TCR (10 samples)", "GSA-Human, controlled access",
     "MEDIUM: droplet; 5' subset suits 5'-aware callers; only 3 OLK",
     "premalignant→cancer axis, CAF-myeloid crosstalk"),
    ("GSE181919", "Nat Commun 2023", "normal → leukoplakia → primary HNSCC → metastasis",
     "37 samples", "droplet scRNA-seq (per study)", "GEO open",
     "MEDIUM-LOW: droplet", "premalignant axis; bulk-style pseudobulk splicing only"),
    ("GSE195832", "Clin Cancer Res 2022", "HNSCC CAF subsets, neoadjuvant ICB (pre/post)",
     "86 samples", "droplet scRNA-seq (per study)", "GEO open",
     "MEDIUM-LOW: droplet", "CAF phenotype ↔ ICB response mapping"),
    ("GSE200996", "2022 (tissue-resident memory T cells)", "oral cancer, pre-surgical ICB",
     "204 samples", "mixed platforms (GPL18573/20795/20301)", "GEO open",
     "MEDIUM: mixed", "ICB context, T-cell states"),
    ("GSE301741", "Cell Rep Med 2026", "16 HNSCC, pre/post neoadjuvant pembrolizumab",
     "58 samples", "droplet scRNA-seq (per study)", "GEO open",
     "LOW-MEDIUM", "treatment-response context"),
    ("GSE310797", "Nat Commun 2026", "OSCC stage series + spatial (MDSC/ANXA1-FPRs)",
     "19 samples", "single-cell + spatial", "GEO open",
     "LOW-MEDIUM", "competitor paper; spatial validation route"),
    ("GSE117872", "2018", "OSCC cell lines: cisplatin sensitive/resistant/holiday",
     "1431 cells", "droplet scRNA-seq (per study)", "GEO open",
     "LOW-MEDIUM", "drug-resistance hypothesis testing"),
    ("TCGA-HNSC", "TCGA", "bulk RNA-seq HNSCC (with clinical + survival)",
     "~500 tumours", "bulk RNA-seq", "GDC / OncoSplicing 3.0 PSI",
     "HIGH for bulk PSI (precomputed)", "cell-type events → bulk PSI → survival"),
    ("GSE270047", "2024 (OncoPrism-HNSCC)", "HNSCC ICB-treated bulk cohort",
     "434 samples", "bulk RNA-seq", "GEO open",
     "HIGH for bulk PSI", "splicing signature ↔ ICB response"),
    ("GSE288199", "2025", "neoadjuvant immunotherapy CD8 dynamics (bulk)",
     "180 samples", "bulk RNA-seq", "GEO open", "HIGH for bulk PSI", "independent ICB validation"),
    ("GSE41613 / GSE42743", "2012-2013", "OSCC cohorts with survival",
     "97 / 103 samples", "microarray", "GEO open",
     "NOT USABLE for PSI (array)", "gene-level only; not a splicing validation set"),
    ("GSE319719", "2026", "oral submucous fibrosis and malignant transformation",
     "13 samples", "bulk RNA-seq", "GEO open", "MEDIUM (small n)", "OSF fibrosis axis"),
]

TOOLS = [
    ("BRIE2", "Bayesian single-cell PSI with covariate regression", "BAM + annotation",
     "per-cell PSI + differential test", "10.1186/s13059-021-02461-5"),
    ("scQuint", "annotation-free, intron-centric global PSI", "BAM/FASTQ",
     "global PSI, cell clustering", "10.7554/elife.73520"),
    ("psix", "cell-state-associated splicing events + co-regulation", "junction counts + cell states",
     "event-level statistics", "10.1101/gr.276109.121"),
    ("SpliZ", "annotation-free splicing Z score, droplet-tolerant", "BAM",
     "per-cell splicing score", "10.1038/s41592-022-01400-x"),
    ("SCASL", "splicing-landscape clustering / PSI imputation", "PSI matrix",
     "cell identity from splicing", "10.1038/s41467-024-46480-9"),
    ("MARVEL", "integrated scRNA splicing platform", "BAM + annotation",
     "PSI distribution, differential events", "10.1093/nar/gkac1260"),
    ("DOLPHIN", "exon + junction read model for single-cell", "BAM",
     "isoform-level readout", "10.1038/s41467-025-61580-w"),
    ("pseudobulk + rMATS", "aggregate cells per state, then bulk differential splicing",
     "count matrix + cell labels", "robust event-level PSI (recommended backbone)",
     "10.1073/pnas.1419161111"),
    ("pseudobulk + LeafCutter", "intron-cluster quantification without transcript annotation",
     "BAM + cell labels", "annotation-free event discovery", "10.1038/s41588-017-0004-9"),
    ("5'-aware caller (NAR 2025)", "AS detection designed for 10x 5' libraries", "5' BAM",
     "AS events linked to ICB response", "10.1093/nar/gkaf1171"),
    ("IsoformSwitchAnalyzeR", "isoform switch + domain/NMD consequences", "transcript quantification",
     "functional consequence of switches", "10.1093/bioinformatics/btz247"),
    ("benchmark (iScience 2026)", "6-method comparison on full-length scRNA-seq",
     "—", "scQuint most robust; droplet coverage is the binding constraint",
     "10.1016/j.isci.2026.116090"),
]


def main():
    with open("results/tables/sc_competition_tiers.csv", "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tier", "definition", "hits", "source", "closest_prior_art_or_note", "trend_2026"])
        for tier, definition, hits, src, note in TIERS:
            w.writerow([tier, definition, hits, src, note, ""])

    with open("results/tables/sc_dataset_feasibility.csv", "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["accession", "study", "scope", "n", "chemistry", "access", "splice_fit", "role"])
        w.writerows(DATASETS)

    with open("results/tables/sc_tool_stack.csv", "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tool", "what_it_does", "input", "output", "doi"])
        w.writerows(TOOLS)

    print("tiers:", len(TIERS), "datasets:", len(DATASETS), "tools:", len(TOOLS))
    print("trend sc_as_all 2026:", trend["sc_as_all"]["2026"], "oral_sc 2026:", trend["oral_sc"]["2026"])


if __name__ == "__main__":
    main()
