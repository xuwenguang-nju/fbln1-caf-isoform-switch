"""`caf_power_measured.csv` under the fixed referent.

The filed table's SDs (0.074 median, 0.117 p75) were measured on the dominant-isoform
dIF. Under the fixed referent the cross-gene distribution is wider, so both the
minimum detectable effect and the patient count needed for a 0.10 effect change.
Figure S1 panels d and e read this table, so it has to follow the statistic.

The SDs are the between-patient SD of each gene's paired dIF, summarised across
genes; the power arithmetic mirrors 39_power_from_measured.py (noncentral t,
n = 6 pairs, simulation with the same seed).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
T = ROOT / "results/tables"
RNG = np.random.default_rng(20261017)
N_PATIENTS = 6


def power_curve(sd_diff, alphas=(0.05, 0.01), n=N_PATIENTS, grid=None, n_iter=6000):
    grid = grid if grid is not None else np.round(np.arange(0.02, 0.42, 0.01), 3)
    out = {}
    for a in alphas:
        ps = []
        for d in grid:
            hits = 0
            for _ in range(n_iter):
                eff = RNG.normal(d, sd_diff, n)
                if stats.ttest_1samp(eff, 0.0).pvalue < a:
                    hits += 1
            ps.append(hits / n_iter)
        out[a] = ps
    return grid, out


def mde(sd_diff, alpha, power=0.8, n=N_PATIENTS):
    lo, hi = 0.0, 0.8
    for _ in range(16):
        mid = (lo + hi) / 2
        _, res = power_curve(sd_diff, alphas=(alpha,), n=n,
                             grid=np.array([mid]), n_iter=2500)
        if res[alpha][0] < power:
            lo = mid
        else:
            hi = mid
    return round(hi, 3)


def n_for(dif, sd_diff, alpha=0.05, power=0.8):
    for n in range(3, 41):
        _, res = power_curve(sd_diff, alphas=(alpha,), n=n,
                             grid=np.array([dif]), n_iter=2500)
        if res[alpha][0] >= power:
            return n
    return 40


def main() -> None:
    g = pd.read_csv(T / "gate_fixed_referent.tsv", sep="\t")
    s = g["sd_dIF"].dropna()
    med, p75 = float(s.median()), float(s.quantile(0.75))
    print(f"fixed-referent between-patient SD: median {med:.4f}  p75 {p75:.4f}  "
          f"(n = {len(s):,} genes)")

    old = pd.read_csv(T / "caf_power_measured.csv")
    old_med = float(old.loc[old["sd_source"] == "measured_median_dIF_SD", "sd"].iloc[0])
    print(f"filed table (dominant-isoform dIF): median {old_med:.3f}")

    sds = {
        "measured_median_dIF_SD": med,
        "measured_p75_dIF_SD": p75,
        "scenario_previous_assumption": 0.10,
    }
    rows = []
    for label, sd in sds.items():
        for alpha in (0.05, 0.01):
            rows.append({"sd_source": label, "sd": round(sd, 3), "alpha": alpha,
                         "min_detectable_dIF_80pct_power_n6": mde(sd, alpha),
                         "n_needed_for_dIF_0.10": n_for(0.10, sd, alpha)})
        print(f"  {label}: SD={sd:.3f} done")
    out = pd.DataFrame(rows)
    print("\n" + out.to_string(index=False))
    out.to_csv(T / "caf_power_measured_fixed_referent.csv", index=False,
               encoding="utf-8-sig")

    json.dump(dict(definition="fixed referent", n_genes=int(len(s)),
                   median_sd=med, p75_sd=p75,
                   filed_median_sd_old_definition=old_med,
                   table=out.to_dict(orient="records")),
              open(ROOT / "data/processed/caf_power_measured_fixed_referent.json",
                   "w", encoding="utf-8"), indent=1)
    print(f"\nwritten: {T / 'caf_power_measured_fixed_referent.csv'}")


if __name__ == "__main__":
    main()
