#!/usr/bin/env python3
"""Why did the positive control fail? Alternative scorings + a content confounder.

Two flaws in script 01 are fixed/tested here:
  1. Signature construction. 01 z-scored each profile ACROSS GENES, so the ranking
     was dominated by expression level rather than subtype specificity. Here both
     ex vivo profiles are CPM-scaled, log2'd, and differenced, with an expression
     floor, so the signature is a real fold change.
  2. Scoring. 01 averaged within-sample z-scores, which couples the score to
     overall expression level. The standard score z-scores each gene ACROSS
     SAMPLES. Both are reported, plus a rank-based variant, so the choice is
     visible instead of hidden.

It also tests the obvious confounder: if fibroblast CONTENT differs between
responders and non-responders, every CAF-subtype score is confounded by purity.
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

FIBRO_MARKERS = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP",
                 "PDGFRB", "THY1", "COL6A1", "COL6A2", "FN1", "POSTN"]
N_TOP = 200
MIN_PROFILE_CPM = 10.0   # gene must reach this in at least one subtype


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


def cpm_log(v):
    return np.log2(v / v.sum() * 1e6 + 1)


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

    a = read_ex_vivo(os.path.join(WORK, "sorted_HNCAF0.3_for_cell_coculture.csv"),
                     ["Replicate.1", "Replicate.2", "Replicate.3"], m)
    b = read_ex_vivo(os.path.join(WORK, "sorted_HNCAF1_for_cell_coculture.csv"),
                     ["FIBRO01", "FIBR02", "FIBR03"], m)
    common = a.index.intersection(b.index)
    common = common[common.isin(lc.index)]
    print("genes in both ex vivo files and the bulk matrix: %d" % len(common))

    acpm = a.sum(axis=1) / a.sum(axis=1).sum() * 1e6     # per-replicate mean CPM
    acpm = a.loc[common].sum(axis=1) / a.loc[common].sum(axis=1).sum() * 1e6
    bcpm = b.loc[common].sum(axis=1) / b.loc[common].sum(axis=1).sum() * 1e6
    la, lb = np.log2(acpm + 1), np.log2(bcpm + 1)
    keep = (acpm >= MIN_PROFILE_CPM) | (bcpm >= MIN_PROFILE_CPM)
    delta = (la - lb)[keep].sort_values(ascending=False)
    up03, up1 = delta.head(N_TOP).index, delta.tail(N_TOP).index
    print("contrast genes after expression floor: %d" % keep.sum())
    print("\ntop 15 up in HNCAF-0/3:", list(up03[:15]))
    print("top 15 up in HNCAF-1 :", list(up1[:15]))

    pre = [s for s in lc.columns if meta.set_index("Sample").loc[s, "Treatment"] == "Pre-Treated"]
    grp = {s: resp[lib2pat[s]] for s in pre}
    rlib = [s for s in pre if grp[s] == "R"]
    nlib = [s for s in pre if grp[s] == "NR"]
    print("\nR libraries: %d  NR libraries: %d" % (len(rlib), len(nlib)))

    def stats_of(score):
        v = pd.Series(np.asarray(score), index=pre)
        r, n = v[rlib], v[nlib]
        u = stats.mannwhitneyu(r, n, alternative="greater")
        return {"auc": u.statistic / (len(r) * len(n)), "p_one_sided": u.pvalue,
                "median_R": r.median(), "median_NR": n.median()}

    rows = []
    zpre = lc[pre].apply(stats.zscore, axis=0)
    sA = zpre.loc[up03].mean() - zpre.loc[up1].mean()
    rows.append(dict(method="A within-sample z across genes (as in 01)", **stats_of(sA)))

    zs = lc[pre].sub(lc[pre].mean(axis=1), axis=0).div(lc[pre].std(axis=1).replace(0, np.nan), axis=0)
    sB = zs.loc[up03].mean() - zs.loc[up1].mean()
    rows.append(dict(method="B gene-wise z across samples (standard)", **stats_of(sB)))

    ranks = lc[pre].rank(axis=0, pct=True)
    sC = ranks.loc[up03].mean() - ranks.loc[up1].mean()
    rows.append(dict(method="C within-sample rank (ssGSEA-like)", **stats_of(sC)))

    fib = [g for g in FIBRO_MARKERS if g in lc.index]
    sF = zs.loc[fib].mean()
    rows.append(dict(method="FIBROBLAST CONTENT (confounder)", **stats_of(sF)))

    res = pd.DataFrame(rows)
    print("\n=== positive control under four scorings (one-sided R > NR) ===")
    print(res.to_string(index=False, float_format=lambda x: "%.4g" % x))
    print("\n  fibroblast markers used: %s" % fib)

    for nm, s in (("A", sA), ("B", sB), ("C", sC)):
        rho, pp = stats.spearmanr(np.asarray(s, float), np.asarray(sF, float))
        print("  corr(score%s, fibroblast content): rho=%+.3f p=%.3g" % (nm, rho, pp))

    # is the subtype ratio itself just tracking fibroblast content?
    res.to_csv(os.path.join(OUT, "c_arm_scoring_variants.tsv"), sep="\t", index=False)
    pd.DataFrame({"gene": delta.index, "log2FC_03_vs_1": delta.values}).to_csv(
        os.path.join(OUT, "c_arm_signature_logfc.tsv"), sep="\t", index=False)
    pd.DataFrame({"lib": pre, "patient": [lib2pat[s] for s in pre],
                  "group": [grp[s] for s in pre],
                  "scoreA": np.asarray(sA), "scoreB": np.asarray(sB),
                  "scoreC": np.asarray(sC), "fibroblast": np.asarray(sF)}).to_csv(
        os.path.join(OUT, "c_arm_scores_per_sample.tsv"), sep="\t", index=False)
    print("\nwritten:", OUT)


if __name__ == "__main__":
    main()
