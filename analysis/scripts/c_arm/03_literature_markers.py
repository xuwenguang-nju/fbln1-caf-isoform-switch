#!/usr/bin/env python3
"""Two defensible second tests before concluding the C arm cannot do this.

Test 1  Literature-defined CAF subtype markers, independent of the authors'
        ex vivo files. myCAF / iCAF / apCAF panels are the standard ones used
        across the CAF literature (Elyada 2019, Dominguez 2020, Kerdidani 2021
        style panels). If these do not separate R from NR either, the limitation
        is the data, not my signature construction.

Test 2  Rebuild the ex vivo contrast with housekeeping genes removed. The first
        signature was topped by RPS10, RPS17, MRPS24, TOMM6, NDUFA13, EEF1G,
        NME2 - genes that track cellularity/purity rather than subtype. Removing
        ribosomal, mitochondrial and a small housekeeping list is a standard
        sensitivity analysis, and it directly addresses the observed rho = -0.69
        between score and fibroblast content.
"""
import os
import sys
import gzip
import numpy as np
import pandas as pd
from scipy import stats

WORK = sys.argv[1] if len(sys.argv) > 1 else "."
OUT = os.path.join(WORK, "results")
FASTA = os.environ.get("GENCODE_FA", "/root/hnscc_caf/ref/gencode.v44.transcripts.fa.gz")

MYCAF = ["ACTA2", "TAGLN", "POSTN", "COMP", "COL11A1", "COL10A1", "COL1A1",
         "COL3A1", "FN1", "HGF", "TGFB1", "CTGF", "CCN2", "MYL9", "TPM2"]
ICAF = ["IL6", "CXCL12", "CXCL1", "CXCL2", "PDGFRA", "HAS1", "CFD", "DPT",
        "LMNA", "IL11", "LIF", "VEGFA", "CCL2", "NFKB1"]
APCAF = ["HLA-DRA", "HLA-DRB1", "CD74", "HLA-DPA1", "HLA-DPB1"]
FIBRO = ["COL1A2", "DCN", "LUM", "FAP", "PDGFRB", "THY1", "COL6A1", "COL6A2"]
HOUSEKEEPING = {"EEF1G", "NME2", "TOMM6", "WDR83OS", "RNASEK", "ADIRF", "MIF",
                "BCYRN1", "RN7SL1", "RN7SK", "RPPH1"}


def ensg2sym_map():
    m = {}
    with gzip.open(FASTA, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                p = line[1:].split("|")
                if len(p) >= 6:
                    m.setdefault(p[1].split(".")[0], p[5].strip())
    return m


def read_ex_vivo(path, cols, m):
    df = pd.read_csv(path)
    df = df.set_index("gene_id") if "gene_id" in df.columns else df.set_index(df.columns[0])
    df = df[cols]
    df.index = [m.get(str(i).split(".")[0]) for i in df.index]
    df = df[[i is not None for i in df.index]]
    return df.groupby(level=0).mean()


def is_junk(g):
    if g in HOUSEKEEPING:
        return True
    if g.startswith(("RPS", "RPL", "MRPS", "MRPL", "MT-", "RNA5", "SCARNA")):
        return True
    return False


def main():
    m = ensg2sym_map()
    counts = pd.read_csv(os.path.join(WORK, "GSE195832_TJU_featurecounts_symbol.tsv.gz"),
                         sep="\t", index_col=0)
    meta = pd.read_csv(os.path.join(WORK, "metadata_TJcohort.csv"),
                       index_col=0, encoding="utf-8-sig")
    lib2pat = dict(zip(meta["Sample"], meta["Patient"]))
    rcol = "Clinical Response(Primary_Site_TE)"
    pat2resp = meta[meta["Treatment"] == "Pre-Treated"].set_index("Patient")[rcol].to_dict()
    resp = {p: ("R" if v in {"Complete Responder", "Responder"} else "NR")
            for p, v in pat2resp.items()}

    cpm = counts / counts.sum(axis=0) * 1e6
    lc = np.log2(cpm + 1)
    pre = [s for s in lc.columns if meta.set_index("Sample").loc[s, "Treatment"] == "Pre-Treated"]
    grp = {s: resp[lib2pat[s]] for s in pre}
    rlib = [s for s in pre if grp[s] == "R"]
    nlib = [s for s in pre if grp[s] == "NR"]

    zs = lc[pre].sub(lc[pre].mean(axis=1), axis=0).div(lc[pre].std(axis=1).replace(0, np.nan), axis=0)

    def report(label, score, extra=""):
        v = pd.Series(np.asarray(score, float), index=pre)
        r, n = v[rlib], v[nlib]
        u = stats.mannwhitneyu(r, n, alternative="two-sided")
        auc = u.statistic / (len(r) * len(n))
        print("  %-46s AUC=%.3f  p=%.3f  R=%+.3f NR=%+.3f %s"
              % (label, auc, u.pvalue, r.median(), n.median(), extra))
        return dict(panel=label, auc=auc, p_two_sided=u.pvalue,
                    median_R=r.median(), median_NR=n.median())

    rows = []

    print("\n=== Test 1: literature CAF subtype panels (gene-wise z, two-sided) ===")
    for name, panel in (("myCAF", MYCAF), ("iCAF", ICAF), ("apCAF", APCAF),
                        ("any CAF (fibroblast content)", FIBRO)):
        genes = [g for g in panel if g in zs.index]
        if len(genes) < 3:
            print("  %-46s skipped (only %d genes present)" % (name, len(genes)))
            continue
        s = zs.loc[genes].mean()
        rows.append(report("%s (n=%d)" % (name, len(genes)), s))
        if name in ("myCAF", "iCAF"):
            globals()["s_" + name] = s

    if "s_myCAF" in globals() and "s_iCAF" in globals():
        ratio = s_myCAF - s_iCAF
        rows.append(report("myCAF minus iCAF", ratio))

    print("\n=== Test 2: rebuild ex vivo contrast without housekeeping genes ===")
    a = read_ex_vivo(os.path.join(WORK, "sorted_HNCAF0.3_for_cell_coculture.csv"),
                     ["Replicate.1", "Replicate.2", "Replicate.3"], m)
    b = read_ex_vivo(os.path.join(WORK, "sorted_HNCAF1_for_cell_coculture.csv"),
                     ["FIBRO01", "FIBR02", "FIBR03"], m)
    common = a.index.intersection(b.index)
    common = common[common.isin(lc.index)]
    clean = [g for g in common if not is_junk(g)]
    print("  genes: %d common, %d after removing housekeeping/ribosomal/mito"
          % (len(common), len(clean)))

    acpm = a.loc[clean].sum(axis=1) / a.loc[clean].sum(axis=1).sum() * 1e6
    bcpm = b.loc[clean].sum(axis=1) / b.loc[clean].sum(axis=1).sum() * 1e6
    keep = (acpm >= 10) | (bcpm >= 10)
    delta = (np.log2(acpm + 1) - np.log2(bcpm + 1))[keep].sort_values(ascending=False)
    up03, up1 = delta.head(200).index, delta.tail(200).index
    print("  top 15 up in HNCAF-0/3:", list(up03[:15]))
    print("  top 15 up in HNCAF-1  :", list(up1[:15]))

    s2 = zs.loc[up03].mean() - zs.loc[up1].mean()
    rows.append(report("cleaned signature (0/3 minus 1)", s2))

    fib_genes = [g for g in FIBRO if g in zs.index]
    sF = zs.loc[fib_genes].mean()
    rho, pp = stats.spearmanr(np.asarray(s2, float), np.asarray(sF, float))
    print("\n  corr(cleaned score, fibroblast content): rho=%+.3f p=%.3g  (was -0.59 to -0.69 when dirty)"
          % (rho, pp))

    pd.DataFrame(rows).to_csv(os.path.join(OUT, "c_arm_literature_panels.tsv"),
                              sep="\t", index=False)
    pd.DataFrame({"gene": delta.index, "log2FC_clean": delta.values}).to_csv(
        os.path.join(OUT, "c_arm_signature_logfc_clean.tsv"), sep="\t", index=False)
    print("\nwritten:", OUT)


if __name__ == "__main__":
    main()
