"""Authoritative rescan: fixed-referent dIF on the SAME input the old scan used.

Why not data/processed/isoform_fracs
------------------------------------
`switch_test_full.tsv` was built from `fractions_full.tsv` (12 libraries, rep1
only, 40 934 genes, per-transcript reads). `data/processed/isoform_fracs` is a
*separate* quantification of the same 12 runs: FBLN1's 1D share differs by up to
0.1214 per library and the gene universe differs (13 658 vs 40 934). Rescoring
isoform_fracs and comparing it with the old scan therefore confounds two changes.
This script uses `fractions_full.tsv`, so the only difference from the original
is the definition.

The defect being corrected
--------------------------
`04_switch.py` set IF1 = the fraction of each library's own dominant transcript.
When the dominant transcript differs between conditions, dIF subtracts one
isoform's fraction from another's. For FBLN1 the dominant call is not even stable
inside a condition (10 of 18 NF libraries 1D, 8 of 18 1C).

The corrected definition
------------------------
    ref(g)  = transcript with the most reads in the pooled NF libraries
    share(g, lib) = reads(ref) / reads(gene)              [read-weighted]
    dIF(g)  = mean over patients of ( share_CAF - share_NF )

Reference is chosen on normal fibroblasts only, so it is not selected for
changing in cancer. Read weighting is legitimate here because `numreads` is a
per-transcript count, unlike the gene-level depth column in isoform_fracs.

Outputs
-------
results/tables/switch_test_fixed_referent.tsv         fixed-referent scan
results/tables/switch_test_olddefinition_replayed.tsv replay of 04_switch.py
data/processed/dif_fixed_referent.json
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
FF = (ROOT / "remote/easyai/hnscc_essential_20261020T1000/"
      "results/fractions_full.tsv")
T = ROOT / "results/tables"
P = ROOT / "data/processed"

MIN_PATIENTS = 5
MIN_READS = 10
DIF_SWITCH, DIF_SMALL = 0.15, 0.05
FBL = "ENSG00000077942"


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return q
    order = np.argsort(pv)
    qv = pv[order] * n / (np.arange(n) + 1)
    qv = np.minimum.accumulate(qv[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(qv, 0, 1)
    q[ok] = out
    return q


def load() -> pd.DataFrame:
    fr = pd.read_csv(FF, sep="\t", low_memory=False)
    fr["patient"] = fr["patient"].astype(str)
    fr["tx"] = fr["tx"].map(lambda t: re.sub(r"\.\d+$", "", str(t)))
    fr["gene"] = fr["gene"].astype(str)
    print(f"fractions_full: {len(fr):,} rows, {fr['run'].nunique()} libraries, "
          f"{fr['gene'].nunique():,} genes, {fr['tx'].nunique():,} transcripts")
    return fr


def replay_old(fr: pd.DataFrame) -> pd.DataFrame:
    """Exactly 04_switch.py: per-library dominant transcript, averaged per cell."""
    dom = fr.loc[fr.groupby(["run", "gene"])["frac"].idxmax(),
                 ["run", "patient", "condition", "gene", "tx", "frac", "numreads"]]
    dom = dom[dom["numreads"] >= MIN_READS]
    per_cell = (dom.groupby(["gene", "patient", "condition"])
                   .agg(IF1=("frac", "mean")).reset_index())
    wide = per_cell.pivot_table(index=["gene", "patient"], columns="condition",
                                values="IF1", aggfunc="mean")
    wide = wide.reindex(columns=["NF", "CAF"]).dropna()
    wide["dIF"] = wide["CAF"] - wide["NF"]
    rows = []
    for g, d in wide.groupby(level="gene"):
        x = d["dIF"].to_numpy()
        if len(x) < MIN_PATIENTS:
            continue
        t, p = stats.ttest_1samp(x, 0.0)
        rows.append(dict(gene=g, n_patients=len(x), mean_dIF=float(x.mean()),
                         p=float(p), t=float(t), n_pos=int((x > 0).sum())))
    r = pd.DataFrame(rows)
    r["q_all"] = bh(r["p"].to_numpy())
    return r.sort_values("p").reset_index(drop=True)


def scan_fixed(fr: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    gene_tot = (fr.groupby(["gene", "run"])["numreads"].sum()
                  .rename("gene_reads").reset_index())
    nf = fr[fr["condition"] == "NF"]
    tot_nf = nf.groupby(["gene", "tx"])["numreads"].sum().reset_index()
    ref = (tot_nf.loc[tot_nf.groupby("gene")["numreads"].idxmax()]
                 .set_index("gene")["tx"])
    sub = fr.merge(ref.rename("ref"), left_on="gene", right_index=True)
    obs = (sub[sub["tx"] == sub["ref"]]
             .groupby(["gene", "run"])["numreads"].sum().rename("ref_reads")
             .reset_index())
    per_lib = gene_tot.merge(obs, on=["gene", "run"], how="left")
    per_lib["ref_reads"] = per_lib["ref_reads"].fillna(0.0)
    per_lib["share"] = per_lib["ref_reads"] / per_lib["gene_reads"].replace(0, np.nan)
    meta = fr.drop_duplicates("run").set_index("run")[["patient", "condition"]]
    per_lib = per_lib.join(meta, on="run")
    wide = per_lib.pivot_table(index=["gene", "patient"], columns="condition",
                               values="share", aggfunc="mean")
    wide = wide.reindex(columns=["NF", "CAF"]).dropna()
    wide["dIF"] = wide["CAF"] - wide["NF"]
    rows = []
    for g, d in wide.groupby(level="gene"):
        x = d["dIF"].to_numpy()
        if len(x) < MIN_PATIENTS:
            continue
        t, p = stats.ttest_1samp(x, 0.0)
        rows.append(dict(gene=g, n_patients=len(x), mean_dIF=float(x.mean()),
                         sd_dIF=float(x.std(ddof=1)), p=float(p), t=float(t),
                         n_pos=int((x > 0).sum()),
                         share_NF=float(d["NF"].mean()),
                         share_CAF=float(d["CAF"].mean()),
                         ref_tx=str(ref.loc[g])))
    r = pd.DataFrame(rows)
    r = r[np.isfinite(r["p"])].copy()
    r["q_all"] = bh(r["p"].to_numpy())
    r["tier"] = np.where(r["mean_dIF"].abs() >= DIF_SWITCH, "switch",
                         np.where(r["mean_dIF"].abs() >= DIF_SMALL,
                                  "small_shift", "negligible"))
    r["same_direction"] = (r["n_pos"] == 0) | (r["n_pos"] == r["n_patients"])
    return r.sort_values("p").reset_index(drop=True), ref


def main() -> None:
    fr = load()

    print("\n=== replay of the original definition on its own input ===")
    old = replay_old(fr)
    filed = pd.read_csv(T / "switch_test_full.tsv", sep="\t")
    j = old.merge(filed[["gene", "mean_dIF", "q_all"]], on="gene",
                  suffixes=("_replay", "_filed"))
    r = stats.spearmanr(j["mean_dIF_replay"], j["mean_dIF_filed"])
    print(f"  genes: replay {len(old):,} vs filed {len(filed):,} vs common {len(j):,}")
    print(f"  Spearman(replay, filed) = {r.statistic:+.4f}")
    print(f"  max |mean_dIF difference| = "
          f"{(j['mean_dIF_replay'] - j['mean_dIF_filed']).abs().max():.5f}")
    f = j[j["gene"] == FBL]
    if len(f):
        print(f"  FBLN1: replay {float(f['mean_dIF_replay'].iloc[0]):+.4f} "
              f"vs filed {float(f['mean_dIF_filed'].iloc[0]):+.4f}")

    print("\n=== fixed-referent scan, same input ===")
    new, ref = scan_fixed(fr)
    print(f"  genes tested: {len(new):,}")
    fn = new[new["gene"] == FBL].iloc[0]
    fo = filed[filed["gene"] == FBL].iloc[0]
    print(f"  FBLN1 old (dominant)        dIF {fo['mean_dIF']:+.4f}  q {fo['q_all']:.4g}")
    print(f"  FBLN1 new (fixed referent)  dIF {fn['mean_dIF']:+.4f}  "
          f"q {fn['q_all']:.4g}  t {fn['t']:+.2f}  n_pos {int(fn['n_pos'])}  "
          f"ref {fn['ref_tx']}")
    print(f"  ref share: NF {fn['share_NF']:.4f} -> CAF {fn['share_CAF']:.4f}")

    jn = new.merge(old[["gene", "mean_dIF"]], on="gene", suffixes=("", "_old"))
    rr = stats.spearmanr(jn["mean_dIF"], jn["mean_dIF_old"])
    sf = int((np.sign(jn["mean_dIF"]) != np.sign(jn["mean_dIF_old"])).sum())
    print(f"\n  genes in both definitions: {len(jn):,}")
    print(f"  Spearman(new, old) = {rr.statistic:+.3f}  sign flips "
          f"{sf:,} ({sf / len(jn):.1%})")
    print(f"  q<0.05  old {int((old['q_all'] < 0.05).sum()):,}  "
          f"new {int((new['q_all'] < 0.05).sum()):,}")
    print(f"  |dIF|>=0.15  old {int((old['mean_dIF'].abs() >= 0.15).sum()):,}  "
          f"new {int((new['mean_dIF'].abs() >= 0.15).sum()):,}")
    print(f"  q<0.05 and |dIF|>=0.15  old "
          f"{int(((old['q_all'] < 0.05) & (old['mean_dIF'].abs() >= 0.15)).sum()):,}  "
          f"new {int(((new['q_all'] < 0.05) & (new['mean_dIF'].abs() >= 0.15)).sum()):,}")

    # stratified: does the definition only matter where the dominant call flips?
    dom = (fr.loc[fr.groupby(["run", "gene"])["frac"].idxmax(),
                  ["gene", "condition", "tx"]]
             .groupby(["gene", "condition"])["tx"]
             .agg(lambda s: s.value_counts().idxmax()).unstack("condition").dropna())
    dom["flips"] = dom["NF"] != dom["CAF"]
    js = jn.merge(dom["flips"].rename("flips"), left_on="gene", right_index=True,
                  how="left")
    print("\n  stratified by whether the old dominant call flips:")
    for label, grp in (("stable", js[js["flips"] == False]),
                       ("flips", js[js["flips"] == True])):
        if len(grp) < 10:
            continue
        r2 = stats.spearmanr(grp["mean_dIF"], grp["mean_dIF_old"])
        s2 = int((np.sign(grp["mean_dIF"]) != np.sign(grp["mean_dIF_old"])).sum())
        print(f"    {label:7s} n={len(grp):,}  rho={r2.statistic:+.3f}  "
              f"sign flips {s2:,} ({s2 / len(grp):.1%})  "
              f"median|old|={grp['mean_dIF_old'].abs().median():.4f}  "
              f"median|new|={grp['mean_dIF'].abs().median():.4f}")

    new.to_csv(T / "switch_test_fixed_referent.tsv", sep="\t", index=False)
    old.to_csv(T / "switch_test_olddefinition_replayed.tsv", sep="\t", index=False)

    info = dict(
        input="remote/easyai/hnscc_essential_20261020T1000/results/fractions_full.tsv "
              "(12 libraries, rep1, 40,934 genes, per-transcript reads)",
        definition="ref(g)=transcript with most reads in pooled NF libraries; "
                   "share=reads(ref)/reads(gene) per library; "
                   "dIF=mean over patients of (CAF-NF)",
        defect_corrected="04_switch.py used each library's own dominant transcript, "
                         "subtracting different isoforms when the dominant call "
                         "changes (results/reports/dif_definition_defect.md)",
        replay_spearman_vs_filed=float(r.statistic),
        replay_max_abs_diff=float((j["mean_dIF_replay"] - j["mean_dIF_filed"]).abs().max()),
        n_genes_old=int(len(filed)), n_genes_new=int(len(new)),
        fbln1=dict(old_dif=float(fo["mean_dIF"]), old_q=float(fo["q_all"]),
                   new_dif=float(fn["mean_dIF"]), new_q=float(fn["q_all"]),
                   new_t=float(fn["t"]), n_pos=int(fn["n_pos"]),
                   ref_tx=str(fn["ref_tx"]),
                   share_NF=float(fn["share_NF"]), share_CAF=float(fn["share_CAF"])),
        q05_old=int((old["q_all"] < 0.05).sum()),
        q05_new=int((new["q_all"] < 0.05).sum()),
        abs_015_old=int((old["mean_dIF"].abs() >= 0.15).sum()),
        abs_015_new=int((new["mean_dIF"].abs() >= 0.15).sum()),
        spearman_vs_old=float(rr.statistic), sign_flips_vs_old=sf,
        min_patients=MIN_PATIENTS, min_reads=MIN_READS,
        tiers=dict(switch=DIF_SWITCH, small=DIF_SMALL),
    )
    (P / "dif_fixed_referent.json").write_text(
        json.dumps(info, indent=1, ensure_ascii=False), encoding="utf-8")

    print("\ntop 12 by q_all (fixed referent):")
    print(new.head(12)[["gene", "n_patients", "mean_dIF", "q_all", "n_pos", "tier"]]
          .round(5).to_string(index=False))
    print("\nwritten:", T / "switch_test_fixed_referent.tsv")


if __name__ == "__main__":
    main()
