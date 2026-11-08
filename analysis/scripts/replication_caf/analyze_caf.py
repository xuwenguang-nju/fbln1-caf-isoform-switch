#!/usr/bin/env python3
"""FBLN1 C/D in pure sorted CAF: replication + subtype-resolved test.

Three questions, in order of how much they cost to answer:

  Q0  Is the hypothesised Rep->subtype mapping real?
      GEO titles say only "Rep 1..16"; the Mendeley matrix labels columns by
      subtype. The mapping is inferred from column order and MUST be verified,
      not assumed. Verification: correlate our gene-level TPMs against that
      matrix (Spearman on expressed genes). A correct mapping gives rho ~0.9;
      a wrong one gives markedly lower. Reported per library.

  Q1  Is C:D in pure sorted CAF near the GSE303690 CAF value (~0.19) rather
      than the normal-fibroblast value (~0.50)? This cohort has no paired
      normal fibroblasts, so the comparison is against those established
      ranges, not a paired fold change.

  Q2  Does C:D differ between CAF subtypes S1-S4? Only answerable for subtypes
      with >=3 libraries; with the pilot (1 per subtype) this is descriptive.

Isoform identity is fixed by UniProt/1997-cloning evidence:
  ENST00000327858 = FBLN1-202 = fibulin-1D (703 aa)
  ENST00000262722 = FBLN1-201 = fibulin-1C (683 aa)
"""
import os
import re
import sys
import gzip
import numpy as np
import pandas as pd
from scipy import stats

CAF = "/root/hnscc_caf/gse195832_caf"
WORK = "/root/hnscc_caf"
IDX_FA = os.path.join(WORK, "ref", "gencode.v44.transcripts.fa.gz")
MENDELEY = os.path.join(WORK, "c_arm", "sorted_CAF_normalizedGeneCounts.csv")

TID_C = "ENST00000262722"   # fibulin-1C, 683 aa
TID_D = "ENST00000327858"   # fibulin-1D, 703 aa
FBLN1 = "ENSG00000077942"
MIN_READS = 10

# established ranges from the GSE303690 paired analysis
GSE303690_CAF_MEDIAN_CD = 0.189
GSE303690_NF_MEDIAN_CD = 0.497


def tx_to_gene(fa):
    m = {}
    with gzip.open(fa, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                p = line[1:].split("|")
                if len(p) >= 6:
                    m[p[0].split(".")[0]] = (p[1].split(".")[0], p[5].strip(), p[4])
    return m


def ensg_to_symbol(fa):
    m = {}
    with gzip.open(fa, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                p = line[1:].split("|")
                if len(p) >= 6:
                    m.setdefault(p[1].split(".")[0], p[5].strip())
    return m


def to_tpm(sf):
    q = pd.read_csv(sf, sep="\t")
    q = q[q["NumReads"] > 0].copy()
    eff = q["EffectiveLength"].clip(lower=1)
    rate = q["NumReads"] / eff
    q["TPM"] = rate / rate.sum() * 1e6
    return q


def main():
    out = os.path.join(CAF, "results")
    os.makedirs(out, exist_ok=True)

    man = pd.read_csv(os.path.join(CAF, "sorted_caf_manifest.tsv"), sep="\t")
    ann = tx_to_gene(IDX_FA)
    sym = ensg_to_symbol(IDX_FA)

    quant_dirs = sorted(d for d in os.listdir(os.path.join(CAF, "quant"))
                        if os.path.exists(os.path.join(CAF, "quant", d, "quant.sf")))
    print("libraries quantified: %d" % len(quant_dirs))
    if not quant_dirs:
        sys.exit("no quant.sf yet")

    # ---------------- per-library FBLN1 isoform fractions ----------------
    rows = []
    tpm_by_lib = {}
    for run in quant_dirs:
        q = to_tpm(os.path.join(CAF, "quant", run, "quant.sf"))
        q["gene"] = q["Name"].map(lambda t: ann.get(t.split(".")[0], (None,))[0])
        # gene-level TPM for the mapping check
        g = q.dropna(subset=["gene"]).groupby("gene")["TPM"].sum()
        g.index = [sym.get(i, i) for i in g.index]
        tpm_by_lib[run] = g.groupby(level=0).sum()

        f = q[q["gene"] == FBLN1].copy()
        if f.empty:
            continue
        f["tid"] = f["Name"].str.split(".").str[0]
        rec = {"run": run, "FBLN1_total_reads": f["NumReads"].sum(),
               "FBLN1_TPM": f["TPM"].sum()}
        for label, tid in (("C", TID_C), ("D", TID_D)):
            sub = f[f["tid"] == tid]
            rec[f"{label}_reads"] = float(sub["NumReads"].sum()) if len(sub) else 0.0
            rec[f"{label}_TPM"] = float(sub["TPM"].sum()) if len(sub) else 0.0
        tot = rec["C_reads"] + rec["D_reads"]
        rec["C_plus_D_reads"] = tot
        rec["C_fraction"] = rec["C_reads"] / tot if tot else np.nan
        rec["D_fraction"] = rec["D_reads"] / tot if tot else np.nan
        rec["C_over_D"] = rec["C_reads"] / rec["D_reads"] if rec["D_reads"] else np.nan
        rows.append(rec)

    iso = pd.DataFrame(rows).merge(
        man[["run", "rep", "subtype_hypothesis"]], on="run", how="left")
    iso = iso.sort_values("rep")
    iso.to_csv(os.path.join(out, "caf_fbln1_isoform.csv"), index=False, encoding="utf-8-sig")

    pd.set_option("display.width", 220)
    print("\n=== Q1: FBLN1 C:D in pure sorted CAF ===")
    show = iso[["run", "rep", "subtype_hypothesis", "FBLN1_total_reads",
                "C_reads", "D_reads", "C_fraction", "D_fraction", "C_over_D"]]
    print(show.round(4).to_string(index=False))

    ok = iso[(iso["C_plus_D_reads"] >= MIN_READS) & iso["C_over_D"].notna()]
    print("\n  libraries with usable FBLN1 reads (>=%d): %d / %d" % (MIN_READS, len(ok), len(iso)))
    if len(ok):
        print("  median C:D = %.3f   range %.3f-%.3f" % (
            ok["C_over_D"].median(), ok["C_over_D"].min(), ok["C_over_D"].max()))
        print("  GSE303690 reference: CAF %.3f | NF %.3f"
              % (GSE303690_CAF_MEDIAN_CD, GSE303690_NF_MEDIAN_CD))
        # one-sample test against the NF value: does pure CAF sit below it?
        t, p = stats.ttest_1samp(np.log2(ok["C_over_D"]), np.log2(GSE303690_NF_MEDIAN_CD))
        print("  one-sample t vs NF median (log2): t=%.3f p=%.4g" % (t, p))
        t2, p2 = stats.ttest_1samp(np.log2(ok["C_over_D"]), np.log2(GSE303690_CAF_MEDIAN_CD))
        print("  one-sample t vs CAF median (log2): t=%.3f p=%.4g" % (t2, p2))

    # ---------------- Q2: subtype comparison ----------------
    print("\n=== Q2: C:D by CAF subtype (hypothesised mapping) ===")
    if len(ok) >= 3:
        g = ok.groupby("subtype_hypothesis")["C_over_D"].agg(["size", "median", "mean"])
        print(g.round(3).to_string())
        groups = [v["C_over_D"].values for _, v in ok.groupby("subtype_hypothesis")
                  if len(v) >= 2]
        if len(groups) >= 2:
            H, pk = stats.kruskal(*groups)
            print("  Kruskal-Wallis across subtypes: H=%.3f p=%.4g (n per group %s)"
                  % (H, pk, [len(x) for x in groups]))
        else:
            print("  too few replicates per subtype for a test (pilot has 1 each)")
    else:
        print("  fewer than 3 usable libraries; descriptive only")

    # ---------------- Q0: verify Rep -> subtype mapping ----------------
    print("\n=== Q0: is the Rep->subtype mapping real? (vs authors' matrix) ===")
    if not os.path.exists(MENDELEY):
        print("  authors' matrix not found at %s; skipping" % MENDELEY)
    else:
        mdf = pd.read_csv(MENDELEY, index_col=0)
        mdf.index = [str(i).split(".")[0] for i in mdf.index]
        mdf = mdf.groupby(level=0).mean()
        mdf.index = [sym.get(i, i) for i in mdf.index]
        mdf = mdf.groupby(level=0).mean()
        print("  matrix: %d genes x %d libraries" % mdf.shape)

        corr_rows = []
        for _, r in iso.iterrows():
            run = r["run"]
            if run not in tpm_by_lib:
                continue
            ours = tpm_by_lib[run]
            best = []
            for col in mdf.columns:
                theirs = mdf[col]
                common = ours.index.intersection(theirs.index)
                if len(common) < 2000:
                    continue
                rho, _ = stats.spearmanr(np.log1p(ours.loc[common]),
                                         np.log1p(theirs.loc[common]))
                best.append((rho, col))
            if not best:
                continue
            best.sort(reverse=True)
            top_rho, top_col = best[0]
            corr_rows.append(dict(run=run, rep=r["rep"],
                                  hyp=r["subtype_hypothesis"],
                                  best_column=top_col, best_rho=round(top_rho, 4),
                                  second_rho=round(best[1][0], 4) if len(best) > 1 else np.nan))
        cc = pd.DataFrame(corr_rows).sort_values("rep")
        print(cc.to_string(index=False))
        if len(cc):
            cc.to_csv(os.path.join(out, "caf_mapping_check.csv"),
                      index=False, encoding="utf-8-sig")
            # the column name encodes the subtype the authors assigned
            cc["author_subtype"] = cc["best_column"].str.split("_").str[-1].str.split(".").str[0]
            print("\n  agreement between hypothesised and author-inferred subtype: %d / %d"
                  % ((cc["author_subtype"] == cc["hyp"]).sum(), len(cc)))
            print("  median best rho = %.3f (a correct library-to-column match should be high "
                  "and clearly above the second best)" % cc["best_rho"].median())
            print("  median gap to second best = %.3f" % (cc["best_rho"] - cc["second_rho"]).median())

    print("\nwritten:", out)


if __name__ == "__main__":
    main()
