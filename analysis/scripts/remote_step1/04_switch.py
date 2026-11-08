#!/usr/bin/env python3
"""Step 4: isoform fractions and the paired switch test from Salmon quants.

Replaces the R version: this node has no R/Bioconductor and scipy here is 1.8 (no
false_discovery_control), so the statistics that IsoformSwitchAnalyzeR's switchTest performs
for the fraction step are reproduced directly. The quantity tested is IF1, the dominant
transcript's share of its gene, and dIF is the per-patient CAF-minus-NF difference of that
share, averaged over the three replicate libraries of each condition.

Design fixed before results were seen (from the local validation):
  * unit of inference = patient (n = 6), three replicate libraries averaged
  * one-sample t-test on the six dIF values per gene
  * two FDR families reported separately: the whole transcriptome, and the pre-registered panel
  * effect size tiers: switch >= 0.15, small shift 0.05-0.15, negligible < 0.05

Usage:
  python3 04_switch.py <workdir> <full|sub10> <panel.tsv> <outdir>
Output:
  <outdir>/fractions_<mode>.tsv         per-library isoform fractions
  <outdir>/switch_test_<mode>.tsv       per-gene paired test
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

MIN_PATIENTS = 5
MIN_READS = 10


def bh(p):
    """Benjamini-Hochberg q-values (no scipy dependency)."""
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order]
    qv = ranked * n / (np.arange(n) + 1)
    qv = np.minimum.accumulate(qv[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(qv, 0, 1)
    q[ok] = out
    return q


def tx_to_gene(fasta_gz):
    """Map transcript id -> gene id from the GENCODE FASTA header (no GTF needed)."""
    import gzip
    m = {}
    opener = gzip.open if fasta_gz.endswith(".gz") else open
    with opener(fasta_gz, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                f = line[1:].split("|")
                tid = f[0].strip()
                gid = f[1].split(".")[0] if len(f) > 1 else tid
                m[tid] = gid
    return m


def main():
    if len(sys.argv) < 5:
        sys.exit("usage: 04_switch.py <workdir> <full|sub10> <panel.tsv> <outdir>")
    work, mode, panel_path, out = sys.argv[1:5]
    os.makedirs(out, exist_ok=True)

    quant_dir = os.path.join(work, "quant_" + mode)
    if not os.path.isdir(quant_dir):
        sys.exit("missing " + quant_dir)

    man = pd.read_csv(os.path.join(work, "gse303690_ena_manifest.tsv"), sep="\t",
                      dtype={"patient": str})
    man["quant"] = [os.path.join(quant_dir, r, "quant.sf") for r in man["run"]]
    man = man[[os.path.exists(q) for q in man["quant"]]]
    man = man.sort_values(["patient", "condition", "replicate"]).reset_index(drop=True)
    print("libraries with quant.sf: %d / 36" % len(man))
    if len(man) < 30:
        sys.exit("too few libraries quantified for mode %s" % mode)

    ann = tx_to_gene(os.path.join(work, "ref", "gencode.v44.transcripts.fa.gz"))
    print("transcripts in annotation: %d" % len(ann))

    frames = []
    for row in man.itertuples():
        q = pd.read_csv(row.quant, sep="\t")
        q = q[q["NumReads"] > 0].copy()
        q["gene"] = q["Name"].map(ann)
        q = q[q["gene"].notna()]
        tot = q.groupby("gene")["NumReads"].transform("sum")
        q["frac"] = q["NumReads"] / tot
        frames.append(pd.DataFrame({
            "run": row.run, "patient": row.patient, "condition": row.condition,
            "replicate": row.replicate, "gene": q["gene"].values,
            "tx": q["Name"].values, "numreads": q["NumReads"].values,
            "frac": q["frac"].values}))
    fr = pd.concat(frames, ignore_index=True)
    fr.to_csv(os.path.join(out, "fractions_%s.tsv" % mode), sep="\t", index=False)
    print("library-gene rows: %d" % len(fr))

    # dominant transcript per gene per library
    dom = fr.loc[fr.groupby(["run", "gene"])["frac"].idxmax(),
                 ["run", "patient", "condition", "gene", "tx", "frac", "numreads"]]
    dom = dom[dom["numreads"] >= MIN_READS]

    # replicate mean per (gene, patient, condition)
    per_cell = (dom.groupby(["gene", "patient", "condition"])
                .agg(IF1=("frac", "mean"), dom_tx=("tx", lambda s: s.mode().iat[0]),
                     n_libs=("run", "nunique"))
                .reset_index())

    caf = per_cell[per_cell["condition"] == "CAF"].set_index(["gene", "patient"])
    nf = per_cell[per_cell["condition"] == "NF"].set_index(["gene", "patient"])
    j = caf[["IF1", "dom_tx"]].join(nf[["IF1", "dom_tx"]], lsuffix="_CAF", rsuffix="_NF",
                                    how="inner").reset_index()
    j["dIF"] = j["IF1_CAF"] - j["IF1_NF"]
    print("genes with paired dIF: %d | patients per gene median %.0f"
          % (j["gene"].nunique(), j.groupby("gene")["patient"].nunique().median()))

    rows = []
    for g, d in j.groupby("gene"):
        x = d["dIF"].to_numpy()
        if len(x) < MIN_PATIENTS:
            continue
        t, p = stats.ttest_1samp(x, 0.0)
        rows.append({"gene": g, "n_patients": len(x), "mean_dIF": x.mean(),
                     "sd_dIF": x.std(ddof=1), "p": p, "n_pos": int((x > 0).sum()),
                     "dom_tx_CAF": d["dom_tx_CAF"].iat[0]})
    res = pd.DataFrame(rows)
    res = res[np.isfinite(res["p"])].copy()
    res["q_all"] = bh(res["p"].values)

    panel = pd.read_csv(panel_path, sep="\t")
    res = res.merge(panel, on="gene", how="left")
    res["in_panel"] = res["symbol"].notna()
    res["q_panel"] = np.nan
    if res["in_panel"].any():
        res.loc[res["in_panel"], "q_panel"] = bh(res.loc[res["in_panel"], "p"].values)
    res["same_direction"] = (res["n_pos"] == 0) | (res["n_pos"] == res["n_patients"])
    res = res.sort_values("p")
    res.to_csv(os.path.join(out, "switch_test_%s.tsv" % mode), sep="\t", index=False)

    print("\n== genes tested: %d | panel genes tested: %d"
          % (len(res), int(res["in_panel"].sum())))
    if res["in_panel"].any():
        pp = res[res["in_panel"] & (res["q_panel"] < 0.05)]
        print("== panel genes with q_panel < 0.05: %d" % len(pp))
        if len(pp):
            print(pp.sort_values("q_panel")[
                ["symbol", "n_patients", "mean_dIF", "sd_dIF", "q_panel", "n_pos"]]
                .round(4).to_string(index=False))
    oo = res[(res["q_all"] < 0.05) & (res["mean_dIF"].abs() >= 0.15)]
    print("\n== transcriptome-wide: q_all<0.05 and |dIF|>=0.15: %d" % len(oo))
    if len(oo):
        print(oo.head(20)[["gene", "n_patients", "mean_dIF", "q_all", "n_pos"]]
              .round(4).to_string(index=False))
    print("\nDONE %s" % mode)


if __name__ == "__main__":
    main()
