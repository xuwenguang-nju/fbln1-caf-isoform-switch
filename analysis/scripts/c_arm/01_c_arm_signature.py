#!/usr/bin/env python3
"""S8 (C-arm) first pass: CAF-subtype score vs nivolumab response in the TJU cohort.

Q1 (positive control): does the authors' own HNCAF-0/3 vs HNCAF-1 contrast
separate responders from non-responders in the 28 pre-treatment bulk samples?
Q2 (AS link): do splicing factors differ between responders and non-responders?

Verified design facts:
  * 56 bulk libraries = 28 patients x (pre, post), fully paired, no gaps.
  * Response CR 3 + R 6 = 9 R; MR 7 + NR 12 = 19 NR (paper: R=9, NR=19).
  * All 86 GSE195832 samples are oral-cavity SCC (TJU cohort).
  * HNCAF0.3: 60,609 versioned ENSG, normalised. HNCAF1: 58,736 unversioned, counts.

Scoring is z-based because the two ex vivo files were processed differently;
each profile is z-scored across genes before differencing. Same class of trick
VIPER uses. This is a positive control, NOT a replication of their VIPER run.
"""
import os
import sys
import json
import numpy as np
import pandas as pd
from scipy import stats

MIN_EXPR_CPM = 1.0
N_TOP = 200


def strip_version(g):
    return str(g).split(".")[0]


def load_ensg2sym(fasta):
    """ENSG -> gene symbol from GENCODE FASTA headers.

    Header layout: >ENST....|ENSG....|OTTHUMG....|OTTHUMT....|TXNAME|SYMBOL|...
    No download needed: this file is already on disk because the Salmon index
    was built from it.
    """
    import gzip
    m = {}
    with gzip.open(fasta, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = line[1:].split("|")
            if len(parts) >= 6:
                ensg = parts[1].split(".")[0]
                sym = parts[5].strip()
                if ensg and sym and ensg not in m:
                    m[ensg] = sym
    return m


def to_symbol(df, ensg2sym):
    """Re-index an ENSG-keyed frame onto gene symbols, collapsing duplicates."""
    mapped = df.copy()
    mapped.index = [ensg2sym.get(g) for g in df.index]
    mapped = mapped[[i is not None for i in mapped.index]]
    return mapped.groupby(level=0).mean()


def read_ex_vivo(path, value_cols=None):
    """Read an ex vivo subtype file.

    The two files do not share a layout: HNCAF0.3 puts the Ensembl id in column 0,
    HNCAF1 puts a row number there and the id in a 'gene_id' column. Detect it
    rather than assume, because assuming silently yields an empty intersection.
    """
    df = pd.read_csv(path)
    if "gene_id" in df.columns:
        df = df.set_index("gene_id")
    else:
        df = df.set_index(df.columns[0])
    if value_cols:
        df = df[value_cols]
    df.index = [strip_version(i) for i in df.index]
    df.index.name = "gene"
    return df.groupby(level=0).mean()


def main():
    work = sys.argv[1] if len(sys.argv) > 1 else "."
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(work, "results")
    os.makedirs(out, exist_ok=True)

    counts = pd.read_csv(os.path.join(work, "GSE195832_TJU_featurecounts_symbol.tsv.gz"),
                         sep="\t", index_col=0)
    meta = pd.read_csv(os.path.join(work, "metadata_TJcohort.csv"),
                       index_col=0, encoding="utf-8-sig")
    print("bulk: %d genes x %d libraries" % counts.shape)
    print("meta: %d rows, columns: %s" % (meta.shape[0], list(meta.columns)))

    lib2pat = dict(zip(meta["Sample"], meta["Patient"]))
    resp_col = "Clinical Response(Primary_Site_TE)"
    pat2resp = (meta[meta["Treatment"] == "Pre-Treated"]
                .set_index("Patient")[resp_col].to_dict())
    pat2therapy = (meta[meta["Treatment"] == "Pre-Treated"]
                   .set_index("Patient")["Therapy"].to_dict())

    responders = {"Complete Responder", "Responder"}
    resp = {p: ("R" if r in responders else "NR") for p, r in pat2resp.items()}
    n_r = sum(1 for v in resp.values() if v == "R")
    n_nr = sum(1 for v in resp.values() if v == "NR")
    print("groups: R=%d  NR=%d  (paper: R=9, NR=19)" % (n_r, n_nr))

    cpm = counts / counts.sum(axis=0) * 1e6
    logcpm = np.log2(cpm + 1)
    expressed = (cpm >= MIN_EXPR_CPM).sum(axis=1) >= 3
    print("expressed genes (CPM>=%s in >=3 libs): %d" % (MIN_EXPR_CPM, expressed.sum()))

    # The bulk matrix is symbol-keyed while both ex vivo files are ENSG-keyed, so an
    # unmapped intersection is empty. Map before intersecting.
    fasta = os.environ.get("GENCODE_FA", "/root/hnscc_caf/ref/gencode.v44.transcripts.fa.gz")
    ensg2sym = load_ensg2sym(fasta)
    print("ENSG->symbol map: %d entries from %s" % (len(ensg2sym), fasta))

    a = to_symbol(read_ex_vivo(os.path.join(work, "sorted_HNCAF0.3_for_cell_coculture.csv")), ensg2sym)
    b = to_symbol(read_ex_vivo(os.path.join(work, "sorted_HNCAF1_for_cell_coculture.csv"),
                               value_cols=["FIBRO01", "FIBR02", "FIBR03"]), ensg2sym)
    a_mean = np.log2(a.mean(axis=1) + 1)
    b_mean = np.log2(b.mean(axis=1) + 1)
    print("ex vivo HNCAF-0/3: %s  HNCAF-1: %s" % (a.shape, b.shape))

    cc = float("nan")
    try:
        rc = pd.read_csv(os.path.join(work, "GSE195832_readcount.txt.gz"), sep="\t")
        rc["gene"] = rc.iloc[:, 0].map(strip_version)
        rc = rc.set_index("gene")[["FIBRO01", "FIBR02", "FIBR03"]].groupby(level=0).mean()
        shared = b.index.intersection(rc.index)
        cc = np.corrcoef(np.log2(b.loc[shared].mean(axis=1) + 1),
                         np.log2(rc.loc[shared].mean(axis=1) + 1))[0, 1]
        print("cross-check HNCAF1 vs readcount.txt.gz: r=%.4f on %d genes" % (cc, len(shared)))
    except Exception as e:
        print("cross-check skipped:", e)

    common = a_mean.index.intersection(b_mean.index)
    za = stats.zscore(a_mean.loc[common])
    zb = stats.zscore(b_mean.loc[common])
    delta = pd.Series(za - zb, index=common).sort_values(ascending=False)
    # mask must come from delta.index, not from `common`: sort_values reordered the
    # rows, so a mask built in `common` order would be applied to the wrong rows.
    print("common genes between the two ex vivo files: %d" % len(common))
    delta = delta[delta.index.isin(logcpm.index[expressed])]
    print("of those, expressed in the bulk cohort: %d" % len(delta))
    up03 = delta.head(N_TOP).index
    up1 = delta.tail(N_TOP).index
    print("signature: %d up in HNCAF-0/3, %d up in HNCAF-1" % (len(up03), len(up1)))

    treat = meta.set_index("Sample")["Treatment"]
    pre = [s for s in logcpm.columns if treat.loc[s] == "Pre-Treated"]
    post = [s for s in logcpm.columns if treat.loc[s] == "Post-Treated"]
    print("pre-treatment libraries: %d  post: %d" % (len(pre), len(post)))

    z_pre = logcpm[pre].apply(stats.zscore, axis=0)
    score = pd.DataFrame({
        "lib": pre,
        "patient": [lib2pat[s] for s in pre],
        "group": [resp.get(lib2pat[s], "?") for s in pre],
        "therapy": [pat2therapy.get(lib2pat[s], "?") for s in pre],
        "score_0_3": z_pre.loc[up03].mean().values,
        "score_1": z_pre.loc[up1].mean().values,
    })
    score["diff"] = score["score_0_3"] - score["score_1"]

    r = score.loc[score.group == "R", "diff"]
    nr = score.loc[score.group == "NR", "diff"]
    u, p_mw = stats.mannwhitneyu(r, nr, alternative="greater")
    auc = u / (len(r) * len(nr))
    rb = 2 * auc - 1
    print("\n=== Q1 positive control: HNCAF-0/3 relative score, R vs NR ===")
    print("  R  n=%d  median=%+.3f" % (len(r), r.median()))
    print("  NR n=%d  median=%+.3f" % (len(nr), nr.median()))
    print("  Mann-Whitney (one-sided R>NR): U=%.1f  p=%.4g" % (u, p_mw))
    print("  AUC=%.3f  rank-biserial=%+.3f" % (auc, rb))
    level_med = score.groupby(score.patient.map(pat2resp))["diff"].agg(["median", "size"])
    print("\n  by response level:")
    print(level_med.to_string())

    sf = ["SRSF1", "SRSF2", "SRSF3", "SRSF4", "SRSF5", "SRSF6", "SRSF7", "SRSF9",
          "SRSF10", "SRSF11", "HNRNPA1", "HNRNPA2B1", "HNRNPC", "HNRNPH1",
          "HNRNPK", "HNRNPL", "HNRNPM", "HNRNPU", "PTBP1", "PTBP3",
          "RBFOX1", "RBFOX2", "ESRP1", "ESRP2", "QKI", "RBM5", "RBM10",
          "RBM39", "SF1", "SF3B1", "U2AF1", "U2AF2", "SRRM1", "SRRM2",
          "TRA2B", "NOVA1", "NOVA2", "MBNL1", "MBNL2", "CELF1", "CELF2"]
    present = [g for g in sf if g in logcpm.index]
    missing = [g for g in sf if g not in logcpm.index]
    print("\n=== Q2 splicing factors: %d present ===" % len(present))
    if missing:
        print("  absent from matrix:", missing)

    grp = dict(zip(score.lib, score.group))
    rows = []
    for g in present:
        vr = logcpm.loc[g, [l for l in pre if grp[l] == "R"]]
        vn = logcpm.loc[g, [l for l in pre if grp[l] == "NR"]]
        if len(vr) < 3 or len(vn) < 3:
            continue
        uu, pp = stats.mannwhitneyu(vr, vn, alternative="two-sided")
        rows.append({"gene": g, "median_R": vr.median(), "median_NR": vn.median(),
                     "delta_logcpm": vr.median() - vn.median(),
                     "auc": uu / (len(vr) * len(vn)), "p": pp})
    sft = pd.DataFrame(rows).sort_values("p").reset_index(drop=True)
    m = len(sft)
    q = sft["p"].values * m / np.arange(1, m + 1)
    sft["q"] = np.clip(np.minimum.accumulate(q[::-1])[::-1], 0, 1)
    print(sft.head(12).to_string(index=False, float_format=lambda x: "%.4g" % x))
    print("\n  splicing factors with q<0.05: %d of %d" % ((sft.q < 0.05).sum(), m))

    pairs = []
    for s in pre:
        cand = [t for t in post if lib2pat[t] == lib2pat[s]]
        if len(cand) == 1:
            pairs.append((s, cand[0]))
    print("\n=== paired pre/post: %d patients ===" % len(pairs))
    d03 = []
    if pairs:
        z_all = logcpm.apply(stats.zscore, axis=0)
        for u2, v2 in pairs:
            d03.append(float(z_all.loc[up03, v2].mean() - z_all.loc[up03, u2].mean()))
        d03 = np.array(d03)
        print("  HNCAF-0/3 score change post-pre: median=%+.3f, Wilcoxon p=%.4g, n=%d"
              % (np.median(d03), stats.wilcoxon(d03).pvalue, len(d03)))

    score.to_csv(os.path.join(out, "c_arm_subtype_score_per_sample.tsv"), sep="\t", index=False)
    delta.rename("delta_z").to_csv(os.path.join(out, "c_arm_hncaf03_vs_hncaf1_signature.tsv"), sep="\t")
    sft.to_csv(os.path.join(out, "c_arm_splicing_factors_response.tsv"), sep="\t", index=False)
    summary = {
        "n_patients": int(score.patient.nunique()),
        "n_responders": n_r, "n_nonresponders": n_nr,
        "expressed_genes": int(expressed.sum()),
        "signature_size_each_direction": N_TOP,
        "hncaf1_vs_readcount_correlation": (None if np.isnan(cc) else round(float(cc), 4)),
        "positive_control": {
            "median_R": round(float(r.median()), 4),
            "median_NR": round(float(nr.median()), 4),
            "mannwhitney_U": round(float(u), 2),
            "p_one_sided": float(p_mw),
            "auc": round(float(auc), 4),
            "rank_biserial": round(float(rb), 4)},
        "paired_pre_post_n": len(pairs),
        "paired_pre_post_median_change": (round(float(np.median(d03)), 4) if len(d03) else None),
        "splicing_factors_tested": int(len(sft)),
        "splicing_factors_q_lt_0.05": int((sft.q < 0.05).sum()),
    }
    with open(os.path.join(out, "c_arm_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print("\nwritten:", out)


if __name__ == "__main__":
    main()
