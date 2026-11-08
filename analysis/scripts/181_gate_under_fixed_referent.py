"""Does the depth gate still keep exactly one switch under a fixed reference?

This is the question that decides whether Figures 1-2 must be redrawn. It runs
the gate (05_gate.sh: sign_same AND |delta dIF| <= 0.04, then effect/q tiers)
twice from the same inputs:

  A. old definition, both depths  -> must reproduce candidates_gated.tsv (validation)
  B. fixed referent, both depths  -> what the gate would report if the statistic
                                     were corrected

The reference transcript is chosen ONCE, from the full-depth NF libraries, and
applied to both depths. Choosing it per depth would let a reference swap be
mistaken for a measurement change, which is exactly what the gate is meant to
test.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "remote/easyai/hnscc_essential_20261020T1000/results"
T = ROOT / "results/tables"
PANEL = ROOT / "analysis/scripts/remote_step1/panel_genes.tsv"

MIN_PATIENTS = 5
MIN_READS = 10
GATE_DELTA = 0.04
DIF_SWITCH, DIF_SMALL = 0.15, 0.05


def bh(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return q
    o = np.argsort(pv)
    qv = pv[o] * n / (np.arange(n) + 1)
    qv = np.minimum.accumulate(qv[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(qv, 0, 1)
    q[ok] = out
    return q


def load(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df["patient"] = df["patient"].astype(str)
    df["gene"] = df["gene"].astype(str)
    df["tx"] = df["tx"].map(lambda t: re.sub(r"\.\d+$", "", str(t)))
    return df


def test_gene(x: np.ndarray) -> tuple[float, float, float, int]:
    if len(x) < MIN_PATIENTS:
        return np.nan, np.nan, np.nan, 0
    t, p = stats.ttest_1samp(x, 0.0)
    return float(x.mean()), float(p), float(t), int((x > 0).sum())


def _per_gene(difs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for g, d in difs.groupby(level="gene"):
        x = d["dIF"].to_numpy()
        m, p, t, np_ = test_gene(x)
        if not np.isfinite(m) or not np.isfinite(p):
            continue
        rows.append(dict(gene=g, n_patients=len(x), mean_dIF=m,
                         sd_dIF=float(np.std(x, ddof=1)) if len(x) > 1 else np.nan,
                         p=p, n_pos=np_))
    r = pd.DataFrame(rows)
    r["q_all"] = bh(r["p"].to_numpy())
    r["same_direction"] = (r["n_pos"] == 0) | (r["n_pos"] == r["n_patients"])
    return r


def old_definition(df: pd.DataFrame) -> pd.DataFrame:
    dom = df.loc[df.groupby(["run", "gene"])["frac"].idxmax(),
                 ["run", "patient", "condition", "gene", "frac", "numreads"]]
    dom = dom[dom["numreads"] >= MIN_READS]
    per_cell = (dom.groupby(["gene", "patient", "condition"])
                   .agg(IF1=("frac", "mean")).reset_index())
    wide = per_cell.pivot_table(index=["gene", "patient"], columns="condition",
                                values="IF1", aggfunc="mean")
    wide = wide.reindex(columns=["NF", "CAF"]).dropna()
    wide["dIF"] = wide["CAF"] - wide["NF"]
    return _per_gene(wide)


def fixed_definition(df: pd.DataFrame, ref: pd.Series) -> pd.DataFrame:
    gene_tot = df.groupby(["gene", "run"])["numreads"].sum().rename("gene_reads").reset_index()
    sub = df.merge(ref.rename("ref"), left_on="gene", right_index=True)
    obs = (sub[sub["tx"] == sub["ref"]]
             .groupby(["gene", "run"])["numreads"].sum().rename("ref_reads").reset_index())
    pl = gene_tot.merge(obs, on=["gene", "run"], how="left")
    pl["ref_reads"] = pl["ref_reads"].fillna(0.0)
    pl["share"] = pl["ref_reads"] / pl["gene_reads"].replace(0, np.nan)
    meta = df.drop_duplicates("run").set_index("run")[["patient", "condition"]]
    pl = pl.join(meta, on="run")
    wide = pl.pivot_table(index=["gene", "patient"], columns="condition",
                          values="share", aggfunc="mean")
    wide = wide.reindex(columns=["NF", "CAF"]).dropna()
    wide["dIF"] = wide["CAF"] - wide["NF"]
    return _per_gene(wide)


def gate(full: pd.DataFrame, sub: pd.DataFrame, panel: pd.DataFrame) -> pd.DataFrame:
    a = full.rename(columns={"mean_dIF": "dIF_full", "p": "p_full"})
    b = sub.rename(columns={"mean_dIF": "dIF_sub", "p": "p_sub"})
    m = a.merge(b[["gene", "dIF_sub"]], on="gene", how="inner")
    m = m.merge(panel, on="gene", how="left")
    m["abs_change"] = (m["dIF_full"] - m["dIF_sub"]).abs()
    m["sign_same"] = np.sign(m["dIF_full"]) == np.sign(m["dIF_sub"])
    m["reproducible"] = m["sign_same"] & (m["abs_change"] <= GATE_DELTA)
    # q_panel is BH within the panel family, as 04_switch.py computed it
    in_p = m["symbol"].notna()
    m["q_panel"] = np.nan
    if in_p.any():
        m.loc[in_p, "q_panel"] = bh(m.loc[in_p, "p_full"].to_numpy())

    def tier(r):
        if not r["reproducible"]:
            return "excluded_depth_unstable"
        a_ = abs(r["dIF_full"])
        q = r["q_panel"]
        if a_ >= DIF_SWITCH and np.isfinite(q) and q < 0.05:
            return "switch"
        if a_ >= DIF_SMALL and np.isfinite(q) and q < 0.05:
            return "small_shift"
        if a_ < DIF_SMALL and np.isfinite(q) and q < 0.05:
            return "negligible_but_reproducible"
        return "not_significant"

    m["call"] = m.apply(tier, axis=1)
    return m


def main() -> None:
    full = load(BASE / "fractions_full.tsv")
    sub10 = load(BASE / "fractions_sub10.tsv")
    panel = pd.read_csv(PANEL, sep="\t")
    panel["gene"] = panel["gene"].astype(str)
    print(f"full: {len(full):,} rows {full['run'].nunique()} libs | "
          f"sub10: {len(sub10):,} rows {sub10['run'].nunique()} libs")

    # reference fixed on full-depth NF, applied to both depths
    nf = full[full["condition"] == "NF"]
    tot = nf.groupby(["gene", "tx"])["numreads"].sum().reset_index()
    ref = tot.loc[tot.groupby("gene")["numreads"].idxmax()].set_index("gene")["tx"]
    print(f"reference transcripts: {len(ref):,}")

    print("\n=== A. old definition (validation against candidates_gated.tsv) ===")
    of, os_ = old_definition(full), old_definition(sub10)
    ga = gate(of, os_, panel)
    print(ga["call"].value_counts().to_string())
    filed = pd.read_csv(T / "candidates_gated.tsv", sep="\t")
    cmp = ga.merge(filed[["gene", "call"]], on="gene", suffixes=("_replay", "_filed"))
    agree = int((cmp["call_replay"] == cmp["call_filed"]).sum())
    print(f"  calls reproduced: {agree:,} of {len(cmp):,} "
          f"({agree / len(cmp):.2%})")

    print("\n=== B. fixed referent (both depths, same reference) ===")
    nf_, ns_ = fixed_definition(full, ref), fixed_definition(sub10, ref)
    gb = gate(nf_, ns_, panel)
    print(gb["call"].value_counts().to_string())

    print("\n=== comparison of what survives ===")
    for label, g in (("old", ga), ("fixed", gb)):
        sw = g[g["call"] == "switch"]
        print(f"  {label:6s} switch n={len(sw)}  "
              f"genes={sw['gene'].tolist()[:6]}")
        fbl = g[g["gene"] == "ENSG00000077942"]
        if len(fbl):
            r = fbl.iloc[0]
            print(f"         FBLN1: dIF_full {r['dIF_full']:+.4f} "
                  f"dIF_sub {r['dIF_sub']:+.4f} "
                  f"|delta| {r['abs_change']:.4f} "
                  f"q_panel {r['q_panel']:.4f} call={r['call']}")

    print("\n  excluded depth-unstable: "
          f"old {int((ga['call'] == 'excluded_depth_unstable').sum()):,}  "
          f"fixed {int((gb['call'] == 'excluded_depth_unstable').sum()):,}")
    print(f"  genes compared: old {len(ga):,}  fixed {len(gb):,}")
    both = ga.merge(gb[["gene", "call"]], on="gene", suffixes=("_old", "_fixed"))
    print(f"  call agreement between definitions: "
          f"{int((both['call_old'] == both['call_fixed']).sum()):,} "
          f"of {len(both):,} ({(both['call_old'] == both['call_fixed']).mean():.1%})")

    ga.to_csv(T / "gate_olddefinition_replayed.tsv", sep="\t", index=False)
    gb.to_csv(T / "gate_fixed_referent.tsv", sep="\t", index=False)

    # ---- emit the fixed-referent tables in the schema Figure 2 already reads,
    # so the figure follows its inputs instead of hard-coded numbers
    panel = panel.copy()
    sym = panel.drop_duplicates("gene").set_index("gene")["symbol"]
    for name, r in (("switch_test_full_fixed_referent.tsv", nf_),
                    ("switch_test_sub10_fixed_referent.tsv", ns_)):
        out = r.copy()
        out["symbol"] = out["gene"].map(sym)
        out["in_panel"] = out["symbol"].notna()
        out["q_panel"] = np.nan
        if out["in_panel"].any():
            out.loc[out["in_panel"], "q_panel"] = bh(
                out.loc[out["in_panel"], "p"].to_numpy())
        out.to_csv(T / name, sep="\t", index=False)
        print(f"  wrote {name}: {len(out):,} genes")
    gate_cols = ["gene", "mean_dIF_full", "mean_dIF_sub", "abs_change", "sign_same",
                 "reproducible", "q_panel_full", "n_pos_full", "q_all_full", "call"]
    out = gb.rename(columns={"dIF_full": "mean_dIF_full", "dIF_sub": "mean_dIF_sub",
                             "q_panel": "q_panel_full", "q_all": "q_all_full",
                             "n_pos": "n_pos_full"})[gate_cols]
    out.to_csv(T / "candidates_gated_fixed_referent.tsv", sep="\t", index=False)
    print(f"  wrote candidates_gated_fixed_referent.tsv: {len(out):,} genes")

    print("\nwritten:", T / "gate_fixed_referent.tsv")


if __name__ == "__main__":
    main()
