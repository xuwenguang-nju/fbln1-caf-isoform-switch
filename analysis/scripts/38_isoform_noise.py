"""Validation step 7: measure, on real data, how reproducible an isoform fraction is,
and how that reproducibility depends on read support.

Design question: with 6 matched patients, is a reported isoform switch (dIF) separable
from noise? The between-patient SD of dIF was previously a scenario assumption; here it is
measured, and the depth dependence of the noise is measured too, because at 600k reads per
library most genes carry only a few dozen reads and counting noise dominates.

Restrictions: gene must have >= 2 transcripts a 75 bp read can distinguish; the reported
quantity is the fraction of the gene's dominant transcript (IF1), the quantity that
isoform-switch analyses test.

Inputs : data/processed/isoform_fracs/*.tsv
         data/processed/kmer_index/{tx_identifiability.tsv,gene_identifiability.tsv}
Outputs: results/tables/caf_isoform_reproducibility.csv
         results/tables/caf_isoform_noise_by_depth.csv
         data/processed/caf_isoform_noise.json
"""
import glob, json
import numpy as np
import pandas as pd

FRACS = "data/processed/isoform_fracs"
IDX = "data/processed/kmer_index"
MIN_READS = 5
DEPTH_BINS = [5, 10, 20, 40, 80, 160, 320, 10 ** 6]


def main():
    frames = [pd.read_csv(p, sep="\t") for p in sorted(glob.glob(f"{FRACS}/*.tsv"))]
    est = pd.concat([f for f in frames if not f.empty], ignore_index=True)
    print("rows:", len(est), flush=True)

    txi = pd.read_csv(f"{IDX}/tx_identifiability.tsv", sep="\t")
    ident = set(pd.read_csv(f"{IDX}/gene_identifiability.tsv", sep="\t")
                .query("gene_is_identifiable")["gene"])
    est = est.merge(txi[["gene", "tx_id", "is_distinguishable"]], on=["gene", "tx_id"], how="left")
    est = est[(est["gene"].isin(ident)) & (est["is_distinguishable"] == True)]

    # dominant transcript per gene per library
    dom = est.loc[est.groupby(["gene", "run"])["frac"].idxmax(),
                  ["gene", "run", "condition", "patient", "replicate", "tx_id", "frac", "n_reads"]] \
        .rename(columns={"tx_id": "dom_tx", "frac": "IF1"})
    dom = dom[dom["n_reads"] >= MIN_READS]
    print("gene-library observations:", len(dom), "genes:", dom["gene"].nunique(), flush=True)

    # ---- replicate reproducibility within (gene, condition, patient), 3 replicates
    rep = []
    for (g, c, p), sub in dom.groupby(["gene", "condition", "patient"]):
        if len(sub) < 3:
            continue
        rep.append({"gene": g, "condition": c, "patient": p,
                    "IF1_mean": sub["IF1"].mean(),
                    "IF1_sd": sub["IF1"].std(ddof=1),
                    "reads_mean": sub["n_reads"].mean()})
    rep = pd.DataFrame(rep)
    print("gene x condition x patient groups with triplicates:", len(rep), flush=True)

    rep["depth_bin"] = pd.cut(rep["reads_mean"], DEPTH_BINS, right=False)
    bydepth = (rep.groupby("depth_bin", observed=True)
               .agg(n_groups=("IF1_sd", "size"),
                    median_reads=("reads_mean", "median"),
                    median_replicate_sd=("IF1_sd", "median"),
                    p90_replicate_sd=("IF1_sd", lambda s: s.quantile(0.9)))
               .reset_index())
    bydepth["depth_bin"] = bydepth["depth_bin"].astype(str)
    bydepth.to_csv("results/tables/caf_isoform_noise_by_depth.csv", index=False,
                   encoding="utf-8-sig")
    print(bydepth.to_string(index=False), flush=True)

    # ---- between-patient spread of dIF, using replicate means (noise already averaged down)
    tmp = rep.copy()
    tmp["patient"] = tmp["patient"].astype(str)
    pooled = tmp.groupby(["gene", "condition", "patient"])["IF1_mean"].first() \
        .unstack(["condition", "patient"])
    if ("CAF" in pooled.columns.get_level_values(0)):
        # every patient that has both conditions, not just the two with triplicates
        caf_pats = {p for c, p in pooled.columns if c == "CAF"}
        nf_pats = {p for c, p in pooled.columns if c == "NF"}
        diffs = {}
        for pat in sorted(caf_pats & nf_pats):
            diffs[pat] = pooled[("CAF", pat)] - pooled[("NF", pat)]
        d = pd.DataFrame(diffs).dropna()
        out = {"n_genes_paired": int(len(d))}
        if len(d) >= 20:
            spread = d.std(axis=1, ddof=1)
            out.update({
                "n_patients": int(d.shape[1]),
                "between_patient_sd_of_dIF_median": float(spread.median()),
                "between_patient_sd_of_dIF_p75": float(spread.quantile(0.75)),
                "abs_dIF_median": float(d.abs().stack().median()),
                "abs_dIF_p90": float(d.abs().stack().quantile(0.90)),
                "frac_dIF_above_0.10": float((d.abs().stack() > 0.10).mean()),
                "frac_dIF_above_0.20": float((d.abs().stack() > 0.20).mean()),
                "note": "replicate noise is averaged over the available replicates before "
                        "differencing; patients contributing a dIF are listed in n_patients",
            })
            d.to_csv("results/tables/caf_dIF_per_patient.csv", encoding="utf-8-sig")
    else:
        out = {"error": "no CAF/NF pairs"}

    noise = {
        "n_genes_quantified": int(dom["gene"].nunique()),
        "n_groups_triplicate": int(len(rep)),
        "median_replicate_sd_overall": float(rep["IF1_sd"].median()),
        "median_IF1_mean": float(rep["IF1_mean"].median()),
        "depth_dependence": bydepth.to_dict(orient="records"),
        "dif": out,
    }
    json.dump(noise, open("data/processed/caf_isoform_noise.json", "w", encoding="utf-8"), indent=1)
    rep.to_csv("results/tables/caf_isoform_reproducibility.csv", index=False, encoding="utf-8-sig")
    print(json.dumps(noise["dif"], indent=1))
    print("median replicate SD overall: %.3f" % noise["median_replicate_sd_overall"])


if __name__ == "__main__":
    main()
