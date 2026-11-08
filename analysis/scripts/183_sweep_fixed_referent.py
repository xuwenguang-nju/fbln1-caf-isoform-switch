"""Recompute the threshold sweep with a fixed-referent dIF.

`41_switch_detection.py` defines IF1 as the *dominant* isoform's fraction per
library, so a gene whose dominant isoform differs between conditions has
IF1(CAF) - IF1(NF) computed between two different transcripts. This mirrors that
script exactly but subtracts the *same* reference transcript's share in both
conditions.

Two facts make the share easy to get here: `isoform_fracs` carries a per-transcript
`frac`, and share-of-transcript = reads_t / gene_reads = frac_t. The `n_reads`
column is the gene-level depth repeated on every isoform row, so it is used only
as a read-support filter, never as a denominator.

Runs the old definition first and asserts it reproduces
`caf_switch_threshold_sweep.csv` exactly, so the comparison is a single-variable
substitution.
"""
from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
FRACS = ROOT / "data/processed/isoform_fracs"
IDX = ROOT / "data/processed/kmer_index"
T = ROOT / "results/tables"

MIN_READS = 30
MIN_PATIENTS = 5
SWEEP_D = (0.05, 0.10, 0.15, 0.20, 0.25)
SWEEP_Q = (0.10, 0.05, 0.01)


def load_filtered() -> pd.DataFrame:
    est = pd.concat([pd.read_csv(p, sep="\t")
                     for p in sorted(glob.glob(str(FRACS / "*.tsv")))],
                    ignore_index=True)
    txi = pd.read_csv(IDX / "tx_identifiability.tsv", sep="\t")
    ident = set(pd.read_csv(IDX / "gene_identifiability.tsv", sep="\t")
                .query("gene_is_identifiable")["gene"])
    est = est.merge(txi[["gene", "tx_id", "is_distinguishable"]],
                    on=["gene", "tx_id"], how="left")
    est = est[(est["gene"].isin(ident)) & (est["is_distinguishable"] == True)
              & (est["n_reads"] >= MIN_READS)]
    est["patient"] = est["patient"].astype(str)
    return est


def summarise(difs: pd.Series) -> dict:
    d = difs.to_numpy()
    n = len(d)
    if n < MIN_PATIENTS:
        return dict(n_patients=n, mean_dIF=np.nan, sd_dIF=np.nan, p=np.nan, n_pos=0)
    t, p = stats.ttest_1samp(d, 0.0)
    return dict(n_patients=n, mean_dIF=float(d.mean()), sd_dIF=float(d.std(ddof=1)),
                p=float(p), n_pos=int((d > 0).sum()))


def per_gene(est: pd.DataFrame, values: pd.Series) -> pd.DataFrame:
    """values indexed like est.index: the per-row quantity to average per cell."""
    tmp = est.assign(v=values)
    cell = (tmp.groupby(["gene", "patient", "condition"])["v"].mean()
               .unstack("condition").dropna())
    cell["dIF"] = cell["CAF"] - cell["NF"]
    rows = []
    for g, s in cell.groupby(level="gene"):
        r = summarise(s["dIF"])
        if not np.isfinite(r["mean_dIF"]) or not np.isfinite(r["p"]):
            continue
        r["gene"] = g
        rows.append(r)
    out = pd.DataFrame(rows)
    out["q"] = stats.false_discovery_control(out["p"].to_numpy(), method="bh")
    out["same_direction"] = (out["n_pos"] == 0) | (out["n_pos"] == out["n_patients"])
    return out


def sweep(res: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for dmin in SWEEP_D:
        for q in SWEEP_Q:
            hit = res[(res["mean_dIF"].abs() >= dmin) & (res["q"] < q)
                      & res["same_direction"]]
            nod = res[(res["mean_dIF"].abs() >= dmin) & (res["q"] < q)]
            rows.append(dict(min_abs_dIF=dmin, fdr=q,
                             n_candidates=int(len(hit)),
                             n_without_direction_filter=int(len(nod))))
    return pd.DataFrame(rows)


def main() -> None:
    est = load_filtered()
    print(f"libraries {est['run'].nunique()} | genes {est['gene'].nunique():,} | "
          f"rows {len(est):,}")

    # ---- A. replay the old definition and check it against the filed sweep
    dom = est.loc[est.groupby(["gene", "run"])["frac"].idxmax()].copy()
    old = per_gene(dom.assign(patient=dom["patient"].astype(str)),
                   dom["frac"])
    sw_old = sweep(old)
    filed = pd.read_csv(T / "caf_switch_threshold_sweep.csv")
    key = ["min_abs_dIF", "fdr"]
    cmp = sw_old.merge(filed, on=key, suffixes=("_replay", "_filed"))
    same = (cmp["n_candidates_replay"] == cmp["n_candidates_filed"]).all()
    print(f"\n=== A. old definition replay: {len(old):,} genes tested ===")
    print(cmp[["min_abs_dIF", "fdr", "n_candidates_replay", "n_candidates_filed"]]
          .to_string(index=False))
    print("  sweep reproduced exactly:", bool(same))
    assert same, "replay does not match the filed sweep"

    # ---- B. fixed referent: one transcript per gene, the tallest in NF
    nf = est[est["condition"] == "NF"]
    nf_mean = nf.groupby(["gene", "tx_id"])["frac"].sum().reset_index()
    ref = nf_mean.loc[nf_mean.groupby("gene")["frac"].idxmax()].set_index("gene")["tx_id"]
    print(f"\n=== B. fixed referent ===\n  references chosen: {len(ref):,}")

    is_ref = est["tx_id"].eq(est["gene"].map(ref))
    # keep only the reference rows, then average over replicates -- averaging over
    # all transcript rows would divide the share by the isoform count
    ref_rows = est[is_ref].copy()
    fixed = per_gene(ref_rows, ref_rows["frac"])
    sw_new = sweep(fixed)
    print(f"  genes tested: {len(fixed):,}")
    print(sw_new.to_string(index=False))

    print("\n=== what changes at the pre-registered thresholds (0.15, FDR<0.05) ===")
    for lbl, r in (("old", old), ("fixed", fixed)):
        s = r[(r["mean_dIF"].abs() >= 0.15) & (r["q"] < 0.05) & r["same_direction"]]
        print(f"  {lbl:6s} n={len(s):2d}  {s['gene'].tolist()}")

    print("\n  genes passing |dIF|>=0.15 in each definition:")
    o = set(old.loc[old["mean_dIF"].abs() >= 0.15, "gene"])
    f = set(fixed.loc[fixed["mean_dIF"].abs() >= 0.15, "gene"])
    print(f"    old {len(o):,}  fixed {len(f):,}  shared {len(o & f):,}")

    fbl = "ENSG00000077942"
    for lbl, r in (("old", old), ("fixed", fixed)):
        row = r[r["gene"] == fbl]
        if len(row):
            x = row.iloc[0]
            print(f"  FBLN1 {lbl:6s}: dIF {x['mean_dIF']:+.4f} p {x['p']:.4g} "
                  f"q {x['q']:.4f} n_pos {int(x['n_pos'])}/{int(x['n_patients'])}")

    sw_new.to_csv(T / "caf_switch_threshold_sweep_fixed_referent.csv", index=False)
    fixed.to_csv(T / "switch_detection_fixed_referent.tsv", sep="\t", index=False)
    json.dump(dict(
        definition="fixed referent (tallest NF transcript), share = frac",
        libraries=int(est["run"].nunique()),
        genes_tested=int(len(fixed)),
        genes_tested_old=int(len(old)),
        sweep_replay_matches_filed=bool(same),
        sweep=sw_new.to_dict(orient="records"),
        n_candidates_015_q05=int(len(fixed[(fixed["mean_dIF"].abs() >= 0.15)
                                           & (fixed["q"] < 0.05)
                                           & fixed["same_direction"]])),
        n_genes_abs_dIF_ge_015=int(len(f)),
        n_genes_abs_dIF_ge_015_old=int(len(o)),
        fbln1=fixed.loc[fixed["gene"] == fbl,
                        ["mean_dIF", "p", "q", "n_pos", "n_patients"]].to_dict("records"),
    ), open(ROOT / "data/processed/sweep_fixed_referent.json", "w", encoding="utf-8"),
        indent=1)
    print("\nwritten:", T / "caf_switch_threshold_sweep_fixed_referent.csv")


if __name__ == "__main__":
    main()
