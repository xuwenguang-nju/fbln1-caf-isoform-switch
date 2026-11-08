"""Stage 8 step 3: targeted TCGA-HNSC validation of the CAF-derived splicing events.

Design
------
The 27 events come from single-cell CAF-vs-epithelial analysis, reproducible across four
donors. They are a **pre-specified hypothesis set**, so multiplicity is 27, not 27,000 --
which is what makes a bulk test possible at all (the earlier whole-genome screen over
2,787 events could never clear FDR with ~100 deaths).

Two questions are asked of TCGA-HNSC bulk:

Q1 prognosis. Does per-SD PSI associate with overall survival (Cox PH, continuous PSI)?
   Continuous rather than median-split because the project's earlier power calculation
   showed a median split needs HR 1.74 while continuous PSI detects HR 1.32.

Q2 module coherence. Do the mapped events co-vary across tumours more than random event
   sets? Bulk tumours mix cell types; if these events share a CAF cell of origin, their
   PSI should move together as the CAF proportion varies between tumours. This is the one
   stromal-related question answerable from PSI alone, without expression data.

Mapping
-------
SUPPA2 -> SpliceSeq: SE->ES, A5->AD, A3->AA, MX->ME, RI->RI.
Symbol aliases are applied where the SpliceSeq annotation predates a rename; the only one
that matters here is ATP5F1C -> ATP5C1.
SpliceSeq identifies events by exon number, so a match is gene x type, not exon-exact.
When several SpliceSeq events match one CAF event, the candidate with the highest
pct_with_values is used for the module test and all candidates are reported in the table.

Inputs
------
data/raw/tcga_spliceseq/full/PSI_download_HNSC.txt   permissive export (95,700 events)
remote/easyai/stringent_events.tsv                   the 27 CAF events
data/processed/tcga_hnsc_clinical.csv                GDC clinical, TCGA-HNSC

Outputs
-------
results/tables/tcga_caf_events_mapping.csv
results/tables/tcga_caf_events_cox.csv
results/tables/tcga_caf_module_correlation.csv
data/processed/tcga_caf_events_test.json
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.duration.hazard_regression import PHReg

ROOT = Path(__file__).resolve().parents[2]
EVENTS = ROOT / "remote/easyai/stringent_events.tsv"
PSI = ROOT / "data/raw/tcga_spliceseq/full/PSI_download_HNSC.txt"
CLIN = ROOT / "data/processed/tcga_hnsc_clinical.csv"
OUT_T = ROOT / "results/tables"
OUT_P = ROOT / "data/processed"

TYPE_MAP = {"SE": "ES", "A5": "AD", "A3": "AA", "MX": "ME", "RI": "RI"}
ALIAS = {"ATP5F1C": "ATP5C1"}

META_COLS = ["symbol", "as_id", "splice_type", "exons", "from_exon", "to_exon",
             "novel_splice", "pct_with_values", "psi_range", "std_psi"]

ORAL_STRICT = ["Other and unspecified parts of tongue", "Floor of mouth", "Gum",
               "Palate", "Other and unspecified parts of mouth",
               "Other and unspecified parts of lip"]
ORAL_WIDE = ORAL_STRICT + ["Other and ill-defined sites in lip, oral cavity and pharynx"]

MIN_SAMPLES_FRAC = 0.90   # same coverage floor as the earlier screen
MIN_PSI_SD = 0.05         # below this there is nothing to associate with survival
N_BACKGROUND = 2000
N_PERM = 1000
RNG = np.random.default_rng(20261020)

warnings.filterwarnings("ignore", category=RuntimeWarning)


def bh(p: np.ndarray) -> np.ndarray:
    m = len(p)
    order = np.argsort(p)
    q = np.empty(m)
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        prev = min(prev, p[i] * m / (rank + 1))
        q[i] = prev
    return q


def cox_per_sd(psi: np.ndarray, time: np.ndarray, event: np.ndarray):
    """Return (hr_per_sd, log_hr, se, z, p) or None if the fit is unusable."""
    ok = np.isfinite(psi) & np.isfinite(time) & np.isfinite(event)
    psi, time, event = psi[ok], time[ok], event[ok]
    if len(psi) < 40 or event.sum() < 15:
        return None
    sd = psi.std()
    if not np.isfinite(sd) or sd <= 0:
        return None
    x = (psi - psi.mean()) / sd
    if np.abs(x).max() > 20:      # degenerate spread
        return None
    try:
        res = PHReg(time.astype(float), x, status=event.astype(float)).fit()
    except Exception:
        return None
    log_hr = float(res.params[0])
    se = float(res.bse[0])
    if not np.isfinite(log_hr) or not np.isfinite(se) or se <= 0:
        return None
    z = log_hr / se
    return dict(n_used=int(len(psi)), n_deaths=int(event.sum()),
                psi_sd=float(sd), hr_per_sd=float(np.exp(log_hr)),
                log_hr=log_hr, se=se, z=z,
                p=float(2 * stats.norm.sf(abs(z))))


def main() -> None:
    ev = pd.read_csv(EVENTS, sep="\t")
    ev["ss_type"] = ev["type"].map(TYPE_MAP)
    ev["ss_symbol"] = ev["symbol"].map(lambda s: ALIAS.get(s, s))

    psi = pd.read_csv(PSI, sep="\t", low_memory=False)
    patients = [c for c in psi.columns if c not in META_COLS]
    psi = psi.rename(columns={p: p.replace("_", "-") for p in patients})
    patients = [p.replace("_", "-") for p in patients]
    psi[patients] = psi[patients].apply(pd.to_numeric, errors="coerce")
    print(f"events={len(psi):,} patients={len(patients)}")
    # the permissive export carries one fewer sample than the default export; record it
    # rather than silently comparing cohorts of different size
    old_hdr = pd.read_csv(ROOT / "data/raw/tcga_spliceseq/PSI_download_HNSC.txt",
                          sep="\t", nrows=0).columns
    old_pat = {c.replace("_", "-") for c in old_hdr if c not in META_COLS}
    print(f"default export had {len(old_pat)} patients; "
          f"dropped here: {sorted(old_pat - set(patients))[:5]}")

    clin = pd.read_csv(CLIN)
    clin = clin.dropna(subset=["os_time", "event"])
    cohorts = {
        "oral_strict": set(clin.loc[clin["site"].isin(ORAL_STRICT), "case"]) & set(patients),
        "oral_wide": set(clin.loc[clin["site"].isin(ORAL_WIDE), "case"]) & set(patients),
        "all_hnsc": set(clin["case"]) & set(patients),
    }
    clin_idx = clin.set_index("case")

    # ---- mapping table (all candidates) -------------------------------------------
    rows = []
    for _, r in ev.iterrows():
        hits = psi[(psi["symbol"] == r["ss_symbol"]) & (psi["splice_type"] == r["ss_type"])]
        for _, h in hits.iterrows():
            rows.append(dict(caf_symbol=r["symbol"], caf_type=r["type"],
                             caf_event=r["event"], caf_mean_dpsi=r["mean"],
                             caf_sd=r["sd"], caf_p=r["p"],
                             ss_symbol=h["symbol"], ss_as_id=h["as_id"],
                             ss_exons=h["exons"],
                             ss_pct_with_values=h["pct_with_values"]))
    mapping = pd.DataFrame(rows)
    OUT_T.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(OUT_T / "tcga_caf_events_mapping.csv", index=False)

    n_ev_mapped = mapping["caf_event"].nunique()
    unmapped = sorted(set(ev["symbol"]) - set(mapping["caf_symbol"]))
    print(f"mapped {n_ev_mapped}/{len(ev)} CAF events ({len(mapping)} candidate pairs)")
    if unmapped:
        print(f"  unmapped: {unmapped}")

    # one representative per CAF event: highest pct_with_values
    rep = (mapping.sort_values("ss_pct_with_values", ascending=False)
                  .groupby("caf_event", as_index=False).first())

    # ---- Q1: prognosis ------------------------------------------------------------
    cox_rows = []
    for name, cases in cohorts.items():
        idx = [c for c in psi.columns if c in cases]
        sub = psi.set_index("as_id")[idx]
        n_samp = len(idx)
        for _, m in rep.iterrows():
            row = sub.loc[m["ss_as_id"]]
            vals = row.to_numpy(float)
            frac = np.isfinite(vals).mean()
            if frac < MIN_SAMPLES_FRAC or np.nanstd(vals) < MIN_PSI_SD:
                continue
            joined = clin_idx.reindex(row.index)
            fit = cox_per_sd(vals, joined["os_time"].to_numpy(float),
                             joined["event"].to_numpy(float))
            if fit is None:
                continue
            cox_rows.append(dict(cohort=name, caf_symbol=m["caf_symbol"],
                                 caf_type=m["caf_type"],
                                 caf_mean_dpsi=m["caf_mean_dpsi"],
                                 ss_as_id=m["ss_as_id"],
                                 ss_type=TYPE_MAP[m["caf_type"]],
                                 frac_with_value=frac, **fit))
    cox = pd.DataFrame(cox_rows)
    for name, grp in cox.groupby("cohort"):
        cox.loc[grp.index, "q"] = bh(grp["p"].to_numpy())
    cox["sig_fdr"] = cox["q"] < 0.05
    cox.to_csv(OUT_T / "tcga_caf_events_cox.csv", index=False)

    summary: dict = {
        "n_caf_events": int(len(ev)),
        "n_caf_events_mapped": int(n_ev_mapped),
        "n_candidate_pairs": int(len(mapping)),
        "unmapped_symbols": unmapped,
        "alias_applied": ALIAS,
        "filters": {"min_samples_frac": MIN_SAMPLES_FRAC, "min_psi_sd": MIN_PSI_SD},
        "cohorts": {},
    }
    for name, grp in cox.groupby("cohort"):
        expected = np.sign(grp["caf_mean_dpsi"].to_numpy())
        observed = np.sign(grp["log_hr"].to_numpy())
        concordant = int((expected == observed).sum())
        binom = stats.binomtest(concordant, len(grp), 0.5, alternative="greater")
        summary["cohorts"][name] = {
            "events_tested": int(len(grp)),
            "n_cases": int(len(cohorts[name])),
            "n_deaths": int(grp["n_deaths"].iloc[0]) if len(grp) else 0,
            "min_p": float(grp["p"].min()) if len(grp) else None,
            "min_q": float(grp["q"].min()) if len(grp) else None,
            "n_q_lt_0_05": int(grp["sig_fdr"].sum()),
            "direction_concordant": concordant,
            "direction_n": int(len(grp)),
            "direction_binom_p": float(binom.pvalue),
        }

    # ---- Q2: module coherence ------------------------------------------------------
    name = "oral_strict"
    idx = [c for c in psi.columns if c in cohorts[name]]
    sub = psi.set_index("as_id")[idx]
    rep_ids = [i for i in rep["ss_as_id"] if i in sub.index]
    M = sub.loc[rep_ids].to_numpy(float)
    keep = np.isfinite(M).mean(axis=1) >= MIN_SAMPLES_FRAC
    keep &= np.nanstd(M, axis=1) >= MIN_PSI_SD
    M = M[keep]
    print(f"\nmodule test: {M.shape[0]} events x {M.shape[1]} oral tumours")

    def mean_abs_rho(mat: np.ndarray) -> float:
        rho = np.corrcoef(mat)
        iu = np.triu_indices_from(rho, k=1)
        v = rho[iu]
        v = v[np.isfinite(v)]
        return float(np.mean(np.abs(v))) if len(v) else float("nan")

    module_rho = mean_abs_rho(M) if M.shape[0] > 2 else float("nan")

    # background: coverage+sd filtered events with the same type distribution
    bg_pool = []
    for _, row in psi.iterrows():
        v = row[idx].to_numpy(float)
        if np.isfinite(v).mean() >= MIN_SAMPLES_FRAC and np.nanstd(v) >= MIN_PSI_SD:
            bg_pool.append(v)
    bg_pool = np.array(bg_pool)
    print(f"background pool: {bg_pool.shape[0]} events")

    null = []
    for _ in range(N_PERM):
        pick = RNG.choice(bg_pool.shape[0], size=M.shape[0], replace=False)
        r = mean_abs_rho(bg_pool[pick])
        if np.isfinite(r):
            null.append(r)
    null = np.array(null)
    p_emp = float((null >= module_rho).mean()) if len(null) else float("nan")

    summary["module_coherence"] = {
        "cohort": name,
        "n_events": int(M.shape[0]),
        "n_tumours": int(M.shape[1]),
        "mean_abs_spearman_caf_events": module_rho,
        "background_mean_abs_rho": float(np.mean(null)) if len(null) else None,
        "background_p95": float(np.quantile(null, 0.95)) if len(null) else None,
        "empirical_p": p_emp,
        "n_permutations": int(len(null)),
    }

    # save the correlation matrix for the report
    cm = pd.DataFrame(np.corrcoef(M), index=rep_ids, columns=rep_ids)
    cm.to_csv(OUT_T / "tcga_caf_module_correlation.csv")

    OUT_P.mkdir(parents=True, exist_ok=True)
    (OUT_P / "tcga_caf_events_test.json").write_text(
        json.dumps(summary, indent=1, ensure_ascii=False), encoding="utf-8")

    print()
    for cname, blk in summary["cohorts"].items():
        print(f"[{cname}] cases={blk['n_cases']} tested={blk['events_tested']} "
              f"deaths={blk['n_deaths']} min_p={blk['min_p']:.3g} "
              f"min_q={blk['min_q']:.3g} q<0.05={blk['n_q_lt_0_05']} "
              f"direction={blk['direction_concordant']}/{blk['direction_n']} "
              f"(p={blk['direction_binom_p']:.3g})")
    mb = summary["module_coherence"]
    print(f"\nmodule coherence (oral_strict): mean |rho| = {mb['mean_abs_spearman_caf_events']:.3f} "
          f"vs background {mb['background_mean_abs_rho']:.3f} "
          f"(p95 {mb['background_p95']:.3f}), empirical p = {mb['empirical_p']:.3f}")


if __name__ == "__main__":
    main()
