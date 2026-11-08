"""Stage 4 step 4: what the 28-patient ICB cohort (9 R / 19 NR) can actually detect.

The C-arm treatment association rests on GSE195832 bulk RNA-seq: 28 patients,
9 responders, 19 non-responders, pre/post paired. This script converts the
qualitative caveat "only direction consistency" into numbers.

Two testing strategies are compared, because they have wildly different power at
this sample size:

  per-event   Test each splicing event separately (t-test on PSI), then correct
              for multiplicity. This is what an rMATS/MAJIQ genome-wide run does.
  set-level   Aggregate a PRE-SPECIFIED event set (the CAF-subtype events found in
              arm A) into one score per patient, then test that score once. Power
              depends on how correlated PSI is across events within a sample.

Scenarios (all explicitly assumptions, not measurements):
  PSI between-sample SD  0.10 / 0.15 / 0.20
  PSI shape              Beta with mean 0.5 (best case: no boundary pile-up)
  inter-event correlation rho 0.2 / 0.6 for the set-level aggregate

Power is computed with the noncentral t (and cross-checked by Monte Carlo on the
Beta scale, which is what reports the boundary effect).

Inputs : none (geometry of the design + stated scenarios)
Outputs: results/tables/icb_power_9v19.csv
         data/processed/icb_power_9v19.json
         figures/icb_power_9v19.png   (via analysis/scripts/31_fig_icb_power.py)
"""
import json

import numpy as np
import pandas as pd
from scipy.stats import nct, t

N1, N2 = 9, 19          # responders / non-responders
ALPHA = 0.05
POWER_TARGET = 0.80
M_EVENTS = 1000         # events tested in a genome-wide per-event run
SET_M = 100             # size of a pre-specified CAF-subtype event set
RHOS = (0.2, 0.6)
SDS = (0.10, 0.15, 0.20)
DELTAS = np.round(np.arange(0.02, 0.36, 0.01), 3)


def power_two_sample(delta, sd_eff, n1=N1, n2=N2, alpha=ALPHA):
    """Power of a two-sample t-test (equal variance) for a mean difference."""
    df = n1 + n2 - 2
    ncp = delta / (sd_eff * np.sqrt(1 / n1 + 1 / n2))
    crit = t.ppf(1 - alpha / 2, df)
    return nct.sf(crit, df, ncp) + nct.cdf(-crit, df, ncp)


def mdd(sd_eff, alpha=ALPHA, target=POWER_TARGET, n1=N1, n2=N2):
    """Smallest detectable mean difference at the target power (bisection)."""
    lo, hi = 0.0, 10 * sd_eff
    for _ in range(200):
        mid = (lo + hi) / 2
        if power_two_sample(mid, sd_eff, n1, n2, alpha) < target:
            lo = mid
        else:
            hi = mid
    return hi


def sd_of_set_mean(sd, rho, m):
    """SD of a per-patient mean over m events with within-sample correlation rho."""
    return sd * np.sqrt(rho + (1 - rho) / m)


def bonferroni_alpha(m):
    return ALPHA / m


def mc_check(delta, sd, alpha=ALPHA, reps=4000, seed=20261018):
    """Monte Carlo on the Beta scale: does the bounded PSI change the geometry?"""
    rng = np.random.default_rng(seed)
    mu_r, mu_n = 0.5 + delta / 2, 0.5 - delta / 2
    out = {}
    for lab, mu in (("R", mu_r), ("NR", mu_n)):
        nu = mu * (1 - mu) / sd ** 2 - 1
        if nu <= 0:
            raise ValueError(f"Beta cannot have mean {mu} and SD {sd}")
        a, b = mu * nu, (1 - mu) * nu
        out[lab] = rng.beta(a, b, size=(reps, N1 if lab == "R" else N2))
    reject = 0
    for i in range(reps):
        x, y = out["R"][i], out["NR"][i]
        se = np.sqrt(x.var(ddof=1) / N1 + y.var(ddof=1) / N2)
        if se > 0 and abs((x.mean() - y.mean()) / se) > t.ppf(1 - alpha / 2, N1 + N2 - 2):
            reject += 1
    return reject / reps


def main():
    rows = []
    for sd in SDS:
        for alpha_lab, alpha in (("per-event, alpha=0.05",
                                  ALPHA),
                                 (f"per-event, Bonferroni m={M_EVENTS}",
                                  bonferroni_alpha(M_EVENTS))):
            d = mdd(sd, alpha)
            rows.append(dict(strategy=alpha_lab, assumption=f"PSI SD={sd}",
                             psi_sd=sd, events_tested=M_EVENTS if "Bonferroni" in alpha_lab else 1,
                             rho=np.nan, min_detectable_psi_delta=round(d, 3),
                             standardized_d=round(d / sd, 2)))
        for rho in RHOS:
            s_eff = sd_of_set_mean(sd, rho, SET_M)
            d = mdd(s_eff)
            rows.append(dict(strategy=f"set-level (m={SET_M}, rho={rho})",
                             assumption=f"PSI SD={sd}", psi_sd=sd, events_tested=SET_M,
                             rho=rho, min_detectable_psi_delta=round(d, 3),
                             standardized_d=round(d / sd, 2)))
    tab = pd.DataFrame(rows)
    tab.to_csv("results/tables/icb_power_9v19.csv", index=False)
    print(tab.to_string(index=False))

    # required n per group for a realistic effect, at the set-level geometry
    req = []
    for delta in (0.10, 0.15, 0.20):
        for rho in RHOS:
            for n in range(5, 201):
                if power_two_sample(delta, sd_of_set_mean(0.15, rho, SET_M), n, n) >= POWER_TARGET:
                    req.append(dict(target_psi_delta=delta, rho=rho,
                                    psi_sd=0.15, balanced_n_per_group=n))
                    break
    reqt = pd.DataFrame(req)
    print("\nrequired balanced n per group (PSI SD=0.15, set-level m=100):")
    print(reqt.to_string(index=False))

    # MC cross-check of the normal geometry on the bounded PSI scale
    checks = []
    for sd in (0.15,):
        for delta in (0.10, 0.20, 0.30):
            mc = mc_check(delta, sd)
            an = float(power_two_sample(delta, sd))
            checks.append(dict(psi_sd=sd, psi_delta=delta,
                               power_noncentral_t=round(an, 3), power_monte_carlo=round(mc, 3),
                               abs_diff=round(abs(an - mc), 3)))
    checkt = pd.DataFrame(checks)
    print("\nMonte Carlo cross-check (Beta-distributed PSI, 4000 reps):")
    print(checkt.to_string(index=False))

    json.dump(dict(n_responders=N1, n_nonresponders=N2, alpha=ALPHA,
                   m_events_tested=M_EVENTS, set_size=SET_M,
                   power_target=POWER_TARGET, scenarios=tab.to_dict("records"),
                   required_n=reqt.to_dict("records"),
                   monte_carlo_check=checkt.to_dict("records"),
                   caveat="PSI SD, Beta shape, and inter-event rho are scenarios, not "
                          "measurements. Read depth is not the binding constraint: ENA "
                          "reports 6.62e9 reads over 86 GSE195832 runs (~77M per library)."),
          open("data/processed/icb_power_9v19.json", "w", encoding="utf-8"),
          indent=1, ensure_ascii=False)
    print("\nwritten: results/tables/icb_power_9v19.csv, data/processed/icb_power_9v19.json")


if __name__ == "__main__":
    main()
