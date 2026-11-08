"""Exact numbers Figure 1 needs, read from the fixed-referent tables.

Fig 1a's arm-B card states the observed dIF, the paired SD (which feeds the power
curve), the q, and the gate's funnel. All four change when the statistic changes,
so they are read here instead of being typed into the plotting script.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
T = ROOT / "results/tables"
FBLN1 = "ENSG00000077942"

full = pd.read_csv(T / "switch_test_full_fixed_referent.tsv", sep="\t")
gated = pd.read_csv(T / "candidates_gated_fixed_referent.tsv", sep="\t")
det = json.load(open(ROOT / "data/processed/sweep_fixed_referent.json", encoding="utf-8"))

row = full[full["gene"] == FBLN1].iloc[0]
g = gated[gated["gene"] == FBLN1].iloc[0]

print("=== Figure 1 arm B (fixed referent) ===")
print(f"  quantified (full depth)      : {len(full):,}")
print(f"  tested at both depths        : {len(gated):,}")
print(f"  depth-unstable (excluded)    : {int((gated['call'] == 'excluded_depth_unstable').sum()):,}")
print(f"  switch                       : {int((gated['call'] == 'switch').sum())}")
print()
print(f"  FBLN1 mean dIF               : {row['mean_dIF']:+.4f}")
print(f"  FBLN1 paired SD              : {row['sd_dIF']:.4f}   <- MEAS_SD for the power curve")
print(f"  FBLN1 q_all / q_panel        : {row['q_all']:.4f} / {row['q_panel']:.4f}")
print(f"  FBLN1 n_pos / n_patients     : {int(row['n_pos'])}/{int(row['n_patients'])}")
print(f"  FBLN1 |delta dIF| across depths: {g['abs_change']:.4f}")
print()
print(f"  genes with |dIF| >= 0.15     : {int((full['mean_dIF'].abs() >= 0.15).sum()):,}")
print(f"  of those, reproducible       : "
      f"{int(((full['mean_dIF'].abs() >= 0.15) & full['gene'].isin(gated.loc[gated['reproducible'], 'gene'])).sum()):,}")

# effect size for 80% power at n = 6
from scipy import stats


def power_of_d(d, n=6, alpha=0.05):
    df = n - 1
    tc = stats.t.ppf(1 - alpha / 2, df)
    return stats.nct.sf(tc, df, d * np.sqrt(n)) + stats.nct.cdf(-tc, df, d * np.sqrt(n))


lo, hi = 0.0, 6.0
for _ in range(200):
    mid = (lo + hi) / 2
    if power_of_d(mid) < 0.80:
        lo = mid
    else:
        hi = mid
print(f"\n  d needed for 80% power, n=6  : {hi:.3f}")
print(f"  observed d = dIF/SD          : {row['mean_dIF'] / row['sd_dIF']:.3f}")
print(f"  power at observed d          : {power_of_d(row['mean_dIF'] / row['sd_dIF']):.3f}")

json.dump(dict(
    quantified=int(len(full)),
    tested_both_depths=int(len(gated)),
    depth_unstable=int((gated["call"] == "excluded_depth_unstable").sum()),
    n_switch=int((gated["call"] == "switch").sum()),
    fbln1_mean_dIF=float(row["mean_dIF"]),
    fbln1_sd=float(row["sd_dIF"]),
    fbln1_q=float(row["q_panel"]),
    fbln1_q_all=float(row["q_all"]),
    fbln1_n_pos=int(row["n_pos"]),
    fbln1_n=int(row["n_patients"]),
    measured_sd=float(row["sd_dIF"]),
    d_for_80pct=float(hi),
), open(ROOT / "data/processed/fig1_arm_b_numbers.json", "w", encoding="utf-8"), indent=1)
print("\nwritten: data/processed/fig1_arm_b_numbers.json")
