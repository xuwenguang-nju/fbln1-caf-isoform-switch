"""Validation step 8: turn the measured isoform noise into a design decision.

Measured (script 27, real GSE303690 data):
  - replicate SD of the dominant-isoform fraction falls from 0.138 at ~8 reads/gene to
    0.010 at ~460 reads/gene; it crosses 0.05 at about 50 reads per gene per library;
  - paired CAF-NF dIF with replicates averaged: median |dIF| 0.049, 27% above 0.10,
    8% above 0.20; between-patient SD of dIF median 0.049, p75 0.103 (n = 2, lower bound).

This script converts those into (i) the reads-per-library needed for a target noise and
(ii) the power of the 6-patient paired design at the measured and conservative SDs.

Inputs : data/processed/caf_isoform_noise.json
         data/processed/caf_power_arithmetic.json (scenario version, for comparison)
Outputs: results/tables/caf_power_measured.csv
         results/tables/caf_depth_requirement.csv
         data/processed/caf_power_measured.json
"""
import json
import numpy as np
import pandas as pd
from scipy import stats

NOISE = "data/processed/caf_isoform_noise.json"
OUT_JSON = "data/processed/caf_power_measured.json"
RNG = np.random.default_rng(20261017)
N_PATIENTS = 6
ASSIGN_RATE = 0.467          # measured: 280k of 600k reads assigned in script 26
N_GENES = 14176              # genes in the index


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
        grid, res = power_curve(sd_diff, alphas=(alpha,), n=n, grid=np.array([mid]), n_iter=2500)
        if res[alpha][0] < power:
            lo = mid
        else:
            hi = mid
    return round(hi, 3)


def n_for(dif, sd_diff, alpha=0.05, power=0.8):
    for n in range(3, 41):
        _, res = power_curve(sd_diff, alphas=(alpha,), n=n, grid=np.array([dif]), n_iter=2500)
        if res[alpha][0] >= power:
            return n
    return None


def main():
    noise = json.load(open(NOISE, encoding="utf-8"))
    depth = pd.DataFrame(noise["depth_dependence"])
    # fit log-linear model: SD ~ a * reads^-b, to interpolate the reads needed for a target SD
    x = np.log(depth["median_reads"].astype(float))
    y = np.log(depth["median_replicate_sd"].astype(float))
    b, a = np.polyfit(x, y, 1)
    pred = np.exp(a + b * np.log(100.0))
    print(f"noise model: SD = {np.exp(a):.3f} * reads^({b:.3f})  (SD at 100 reads = {pred:.3f})")

    rows = []
    for target in (0.05, 0.04, 0.03, 0.02):
        reads = float(np.exp((np.log(target) - a) / b))
        rows.append({"target_replicate_sd": target,
                     "reads_per_gene_needed": round(reads),
                     "reads_per_library_needed_median_gene": round(reads * N_GENES / ASSIGN_RATE),
                     "note": "per library; the median gene needs this, low-expressed genes more"})
    depth_req = pd.DataFrame(rows)
    depth_req.to_csv("results/tables/caf_depth_requirement.csv", index=False, encoding="utf-8-sig")
    print("\n" + depth_req.to_string(index=False))

    sds = {
        "measured_median_dIF_SD": noise["dif"]["between_patient_sd_of_dIF_median"],
        "measured_p75_dIF_SD": noise["dif"]["between_patient_sd_of_dIF_p75"],
        "scenario_previous_assumption": 0.10,
    }
    prows, out = [], {}
    for label, sd in sds.items():
        for alpha in (0.05, 0.01):
            e = mde(sd, alpha)
            n100 = n_for(0.10, sd, alpha)
            prows.append({"sd_source": label, "sd": round(sd, 3), "alpha": alpha,
                          "min_detectable_dIF_80pct_power_n6": e,
                          "n_needed_for_dIF_0.10": n100})
            out[f"{label}_alpha{alpha}"] = {"mde_n6": e, "n_for_0.10": n100}
        print(f"{label}: SD={sd:.3f} -> MDE(n=6, a=0.05)={out[f'{label}_alpha0.05']['mde_n6']:.3f}, "
              f"a=0.01  {out[f'{label}_alpha0.01']['mde_n6']:.3f}")
    pd.DataFrame(prows).to_csv("results/tables/caf_power_measured.csv", index=False,
                               encoding="utf-8-sig")
    out["noise_model"] = {"form": "SD = c * reads^b", "c": float(np.exp(a)), "b": float(b)}
    out["depth_requirement"] = depth_req.to_dict(orient="records")
    out["measured_dIF"] = noise["dif"]
    json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), indent=1)


if __name__ == "__main__":
    main()
