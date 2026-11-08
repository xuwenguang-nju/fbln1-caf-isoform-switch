# Data sources

All data are public. No new sequencing was generated for this study.

| Dataset | Accession | Content | Used for |
|---|---|---|---|
| Paired NF/CAF RNA-seq | **GEO GSE303690** | 6 head and neck cancer patients, normal-fibroblast / cancer-associated-fibroblast pairs, 36 libraries (3 replicates each), NextSeq 2000 | Discovery layer — Figures 1–3 |
| Single-cell RNA-seq | **GEO GSE301741** | Head and neck tumours before/after neoadjuvant pembrolizumab, 58 samples, 15 donors. This study analysed 7 libraries from 5 donors (29,295 cells) | Figure 4, and the CAF-versus-epithelial contrast arm |
| Nivolumab cohort | **GEO GSE195832** | Neoadjuvant nivolumab HNSCC series (Obradovic et al., *Clin Cancer Res* 2022). Supplies two distinct subsets: 11 sorted-CAF libraries (GSM5851573–GSM5851588) and a 28-patient bulk expression cohort | Figure 8 (7 libraries after the expression filter); Figure 6F (immunotherapy) |
| TCGA-HNSC | via **UCSC Xena** and **TCGA SpliceSeq** | 501 tumours, 43 adjacent normal; PSI values plus expression with clinical annotation | Figure 6 |
| TCGA-PAAD, BRCA, LUAD, ESCA | via **UCSC Xena** and **TCGA SpliceSeq** | Same two layers | Figure 7 |

## Annotation and reference sets

| Resource | Version | Used for |
|---|---|---|
| GENCODE | v44 | Transcript models, CDS coordinates, effective lengths |
| MSigDB | 7.5.1, `NABA_MATRISOME` | Matrisome gene definition |
| UniProt | P23142 (FBLN1) | Domain boundaries |
| Pfam | PF22914 | C-terminal Ig-like domain |
| Cancer-associated-fibroblast marker signature | 559-marker set | NNLS deconvolution of tumour composition |

## Notes on reuse

- **GSE195832 is used twice.** The sorted-CAF libraries (Figure 8) and the bulk
  immunotherapy cohort (Figure 6F) come from the same GEO series and the same clinical
  trial. They are different library preparations from the same study, not independent
  cohorts. The manuscript states this in the Data Availability section and in the
  Figure 8 legend.
- **GSE303690** carries no linked publication in GEO. It is cited as
  "dataset, accession GSE303690".
- **GSE301741** contains 58 samples; this study used a 7-library subset. Downloading
  the full series will therefore yield more libraries than the manuscript reports.
- **GSE195832** contains 86 libraries in total; 7 sorted-CAF libraries were analysed.

## Regenerating the derived tables

The per-library quantification that feeds every table is produced by the pipeline in
`analysis/scripts/remote_step1/` (alignment and quantitation) and
`analysis/scripts/c_arm/`, `replication_caf/` (per-arm scoring). Those steps expect the
raw FASTQ files from the accessions above and are not run as part of the figure scripts.
