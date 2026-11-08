"""Stage 8 step 4: TCGA-HNSC bulk validation of the CAF-derived splicing events.

Correction to earlier work
-------------------------
`39_psi_pilot_spliceseq.py` recorded "no normal-tissue contrast available in this file".
That was wrong: the SpliceSeq export carries normal samples with a `_Norm` suffix --
43 of them, 42 with a matched primary tumour. This changes the primary test, because
tumour-vs-normal is a far more direct probe of the CAF axis than survival: normal oral
mucosa is epithelium with essentially no CAF, so for an event where CAF PSI > epithelial
PSI the tumour should sit higher than the normal.

Three questions, in order of how directly they test the single-cell finding
----------------------------------------------------------------------------
Q1 tumour vs normal. Paired Wilcoxon where a patient has both, else Mann-Whitney.
   Directional prediction: sign(tumour - normal) should equal sign(CAF - epithelial).
Q2 prognosis. Cox PH on continuous PSI per SD (median split needs HR 1.74, continuous 1.32).
Q3 module coherence. Do the mapped events co-vary across tumours more than random sets?
   Bulk mixes cell types, so a shared CAF origin predicts co-variation.

The event set is pre-specified from single-cell data (27 events, 4 donors). That is the
whole point: 18 tests can clear FDR; the earlier 2,787-event screen could not.

Mapping: SUPPA2 SE/A5/A3/MX/RI -> SpliceSeq ES/AD/AA/ME/RI, plus the alias
ATP5F1C -> ATP5C1 (SpliceSeq's annotation predates the rename).
SpliceSeq indexes events by exon number, so matches are gene x type, not exon-exact.

Outputs
-------
results/tables/tcga_caf_events_mapping.csv
results/tables/tcga_caf_events_tvn.csv        tumour-vs-normal per event
results/tables/tcga_caf_events_cox.csv        prognosis per event
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

MIN_FRAC = 0.50          # coverage floor on both sides of the tumour/normal test
MIN_SD = 0.03            # SD floor, lower than the screen because the export is permissive
N_PERM = 1000
RNG = np.random.default_rng(20261020)

warnings.filterwarnings("ignore")


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


def read_header() -> list[str]:
    with open(PSI, encoding="utf-8", errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
    return [c for c in hdr if c != ""]      # the export carries a trailing empty field


def load_patients(hdr: list[str], wanted_ids: set[str]) -> pd.DataFrame:
    """Stream the 274 MB export, keeping only the mapped events."""
    idx = {c: i for i, c in enumerate(hdr)}
    pcols = hdr[10:]
    keep_pos = [idx[c] for c in pcols]
    rows = []
    with open(PSI, encoding="utf-8", errors="replace") as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= keep_pos[-1]:
                continue
            if f[idx["as_id"]] not in wanted_ids:
                continue
            vals = []
            for p in keep_pos:
                x = f[p]
                vals.append(np.nan if x in ("", "null", "NA", "nan") else float(x))
            rows.append([f[idx["symbol"]], f[idx["as_id"]], f[idx["splice_type"]],
                         f[idx["pct_with_values"]]] + vals)
    return pd.DataFrame(rows, columns=["symbol", "as_id", "splice_type",
                                       "pct_with_values"] + pcols)


def main() -> None:
    ev = pd.read_csv(EVENTS, sep="\t")
    ev["ss_type"] = ev["type"].map(TYPE_MAP)
    ev["ss_symbol"] = ev["symbol"].map(lambda s: ALIAS.get(s, s))

    hdr = read_header()
    pcols = hdr[10:]
    tumours = [c for c in pcols if not c.endswith("_Norm")]
    normals = [c for c in pcols if c.endswith("_Norm")]
    print(f"export: {len(pcols)} sample columns -> {len(tumours)} tumour, "
          f"{len(normals)} normal")

    clin = pd.read_csv(CLIN).dropna(subset=["os_time", "event"])
    clin_idx = clin.set_index("case")
    # PSI columns use underscores (TCGA_CV_6933); the clinical index uses hyphens.
    # Passing the underscored names to reindex silently yields all-NaN.
    as_case = lambda c: c.replace("_Norm", "").replace("_", "-")

    # ---- mapping (needs symbol/type lookup over the whole export, not just mapped ids)
    # a two-pass read would double the 274 MB scan; instead read meta columns only
    meta = pd.read_csv(PSI, sep="\t", usecols=["symbol", "as_id", "splice_type",
                                               "exons", "pct_with_values"],
                       low_memory=False, dtype={"as_id": str})
    print(f"export events: {len(meta):,} across {meta['symbol'].nunique():,} symbols")

    rows = []
    for _, r in ev.iterrows():
        hits = meta[(meta["symbol"] == r["ss_symbol"]) & (meta["splice_type"] == r["ss_type"])]
        for _, h in hits.iterrows():
            rows.append(dict(caf_symbol=r["symbol"], caf_type=r["type"], caf_event=r["event"],
                             caf_mean_dpsi=r["mean"], caf_sd=r["sd"], caf_p=r["p"],
                             ss_symbol=h["symbol"], ss_as_id=str(h["as_id"]),
                             ss_type=h["splice_type"],
                             ss_exons=h["exons"], ss_pct_with_values=h["pct_with_values"]))
    mapping = pd.DataFrame(rows)
    OUT_T.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(OUT_T / "tcga_caf_events_mapping.csv", index=False)
    n_mapped = mapping["caf_event"].nunique()
    unmapped = sorted(set(ev["symbol"]) - set(mapping["caf_symbol"]))
    print(f"mapped {n_mapped}/{len(ev)} CAF events ({len(mapping)} candidate pairs)")
    print(f"  unmapped symbols: {unmapped}")

    # ---- carry the PSI values for candidate events
    sub = load_patients(hdr, set(mapping["ss_as_id"]))
    sub = sub.set_index("as_id")
    print(f"PSI loaded for {len(sub)} SpliceSeq events")

    # cohort masks over tumour columns
    site_of = {c: clin_idx["site"].get(as_case(c)) for c in tumours}
    cohort_cases = {
        "oral_strict": {c for c, s in site_of.items() if s in ORAL_STRICT},
        "oral_wide": {c for c, s in site_of.items() if s in ORAL_WIDE},
        "all_hnsc": set(tumours),
    }

    # ---- Q1: tumour vs normal -----------------------------------------------------
    tvn_rows, cox_rows = [], []
    for _, m in mapping.iterrows():
        aid = m["ss_as_id"]
        if aid not in sub.index:
            continue
        row = sub.loc[aid]
        tvals_all = row[tumours].to_numpy(float)
        nvals_all = row[normals].to_numpy(float)
        if np.isfinite(tvals_all).mean() < MIN_FRAC or np.isfinite(nvals_all).mean() < MIN_FRAC:
            continue

        for cname, cases in cohort_cases.items():
            cidx = [c for c in tumours if c in cases]
            t = row[cidx].to_numpy(float)
            if np.isfinite(t).mean() < MIN_FRAC:
                continue
            # paired where the same patient has a normal
            pairs = [(c, c + "_Norm") for c in cidx if (c + "_Norm") in normals]
            pairs = [(a, b) for a, b in pairs
                     if np.isfinite(row[a]) and np.isfinite(row[b])]
            if len(pairs) >= 8:
                d = np.array([row[a] - row[b] for a, b in pairs])
                stat, p = stats.wilcoxon(d)
                test = "wilcoxon_paired"
                delta = float(d.mean())
                n_t, n_n = len(pairs), len(pairs)
            else:
                n = row[normals].to_numpy(float)
                n = n[np.isfinite(n)]
                tt = t[np.isfinite(t)]
                if len(n) < 5 or len(tt) < 20:
                    continue
                stat, p = stats.mannwhitneyu(tt, n, alternative="two-sided")
                test = "mannwhitney"
                delta = float(np.mean(tt) - np.mean(n))
                n_t, n_n = len(tt), len(n)
            tvn_rows.append(dict(cohort=cname, caf_symbol=m["caf_symbol"],
                                 caf_type=m["caf_type"], caf_mean_dpsi=m["caf_mean_dpsi"],
                                 ss_as_id=aid, ss_type=m["ss_type"],
                                 test=test, n_tumour=n_t, n_normal=n_n,
                                 mean_delta=delta, stat=float(stat), p=float(p)))
    tvn = pd.DataFrame(tvn_rows)

    # ---- Q2: prognosis ------------------------------------------------------------
    for cname, cases in cohort_cases.items():
        cidx = [c for c in tumours if c in cases]
        for _, m in mapping.iterrows():
            aid = m["ss_as_id"]
            if aid not in sub.index:
                continue
            v = sub.loc[aid, cidx].to_numpy(float)
            if np.isfinite(v).mean() < MIN_FRAC or np.nanstd(v) < MIN_SD:
                continue
            j = clin_idx.reindex([as_case(c) for c in cidx])
            ok = np.isfinite(v) & j["os_time"].notna().to_numpy() & j["event"].notna().to_numpy()
            if ok.sum() < 40 or j["event"].to_numpy(float)[ok].sum() < 15:
                continue
            x = v[ok]
            x = (x - x.mean()) / x.std()
            try:
                res = PHReg(j["os_time"].to_numpy(float)[ok], x,
                            status=j["event"].to_numpy(float)[ok]).fit()
            except Exception:
                continue
            lh, se = float(res.params[0]), float(res.bse[0])
            if not (np.isfinite(lh) and np.isfinite(se) and se > 0):
                continue
            z = lh / se
            cox_rows.append(dict(cohort=cname, caf_symbol=m["caf_symbol"],
                                 caf_type=m["caf_type"], caf_mean_dpsi=m["caf_mean_dpsi"],
                                 ss_as_id=aid, ss_type=m["ss_type"],
                                 n_used=int(ok.sum()),
                                 n_deaths=int(j["event"].to_numpy(float)[ok].sum()),
                                 hr_per_sd=float(np.exp(lh)), log_hr=lh, se=se, z=z,
                                 p=float(2 * stats.norm.sf(abs(z)))))
    cox = pd.DataFrame(cox_rows)

    for df, name in ((tvn, "tvn"), (cox, "cox")):
        if len(df):
            for cname, grp in df.groupby("cohort"):
                df.loc[grp.index, "q"] = bh(grp["p"].to_numpy())
            df["sig_fdr"] = df["q"] < 0.05
    if len(tvn):
        tvn.sort_values(["cohort", "p"]).to_csv(OUT_T / "tcga_caf_events_tvn.csv", index=False)
    if len(cox):
        cox.sort_values(["cohort", "p"]).to_csv(OUT_T / "tcga_caf_events_cox.csv", index=False)

    summary: dict = {
        "n_caf_events": int(len(ev)),
        "n_caf_events_mapped": int(n_mapped),
        "n_candidate_pairs": int(len(mapping)),
        "unmapped_symbols": unmapped,
        "alias_applied": ALIAS,
        "samples": {"tumour": len(tumours), "normal": len(normals)},
        "filters": {"min_frac": MIN_FRAC, "min_sd": MIN_SD},
        "tumour_vs_normal": {},
        "prognosis": {},
    }

    for cname, grp in (tvn.groupby("cohort") if len(tvn) else []):
        # One vote per CAF event. The candidate pairs are not independent -- a CAF event
        # with 7 matching SpliceSeq events would otherwise contribute 7 correlated votes.
        per_event = grp.groupby("caf_symbol").agg(
            delta=("mean_delta", "median"), dpsi=("caf_mean_dpsi", "first"))
        exp = np.sign(per_event["dpsi"].to_numpy())
        obs = np.sign(per_event["delta"].to_numpy())
        conc = int((exp == obs).sum())
        bt = stats.binomtest(conc, len(per_event), 0.5, alternative="greater")
        # effect size floor: a significant p on |delta| = 0.0005 is not a biological finding
        big = per_event[per_event["delta"].abs() >= 0.05]
        big_conc = int((np.sign(big["dpsi"]) == np.sign(big["delta"])).sum())
        summary["tumour_vs_normal"][cname] = {
            "pairs_tested": int(len(grp)),
            "events_tested": int(len(per_event)),
            "min_p": float(grp["p"].min()), "min_q": float(grp["q"].min()),
            "n_pairs_q_lt_0_05": int(grp["sig_fdr"].sum()),
            "n_p_lt_0_05": int((grp["p"] < 0.05).sum()),
            "direction_concordant": conc, "direction_n": int(len(per_event)),
            "direction_binom_p": float(bt.pvalue),
            "abs_delta_ge_0_05": int(len(big)),
            "abs_delta_ge_0_05_concordant": big_conc,
            "max_abs_delta": float(per_event["delta"].abs().max()),
            "median_abs_delta": float(grp["mean_delta"].abs().median()),
        }
    for cname, grp in (cox.groupby("cohort") if len(cox) else []):
        exp = np.sign(grp["caf_mean_dpsi"].to_numpy())
        obs = np.sign(grp["log_hr"].to_numpy())
        conc = int((exp == obs).sum())
        bt = stats.binomtest(conc, len(grp), 0.5, alternative="greater")
        summary["prognosis"][cname] = {
            "events_tested": int(len(grp)),
            "min_p": float(grp["p"].min()), "min_q": float(grp["q"].min()),
            "n_q_lt_0_05": int(grp["sig_fdr"].sum()),
            "direction_concordant": conc, "direction_n": int(len(grp)),
            "direction_binom_p": float(bt.pvalue),
        }

    # ---- Q3: module coherence -----------------------------------------------------
    cname = "oral_strict"
    cidx = [c for c in tumours if c in cohort_cases[cname]]
    # Pick one SpliceSeq event per CAF event by ORAL-cohort coverage. Ranking by the
    # whole-export pct_with_values (the previous rule) left only 2 usable events,
    # which made the module test vacuous.
    ids, mats = [], []
    for ev_id, grp in mapping.groupby("caf_event"):
        best, best_frac = None, -1.0
        for aid in grp["ss_as_id"]:
            if aid not in sub.index:
                continue
            v = sub.loc[aid, cidx].to_numpy(float)
            frac = np.isfinite(v).mean()
            if frac > best_frac and np.nanstd(v) >= MIN_SD:
                best, best_frac = v, frac
        if best is not None and best_frac >= MIN_FRAC:
            ids.append(ev_id); mats.append(best)
    if mats:
        M = np.array(mats)
        keep = np.isfinite(M).all(axis=1)
        M = M[keep]
        ids_k = [i for i, k in zip(ids, keep) if k]
        print(f"\nmodule test: {M.shape[0]} events x {M.shape[1]} oral tumours")

        def mean_abs_rho(mat):
            r = np.corrcoef(mat)
            iu = np.triu_indices_from(r, k=1)
            v = r[iu]; v = v[np.isfinite(v)]
            return float(np.mean(np.abs(v))) if len(v) else float("nan")

        if M.shape[0] > 2:
            obs = mean_abs_rho(M)
            pool = []
            for aid in sub.index:
                v = sub.loc[aid, cidx].to_numpy(float)
                if np.isfinite(v).all() and np.nanstd(v) >= MIN_SD:
                    pool.append(v)
            pool = np.array(pool)
            null = []
            for _ in range(N_PERM):
                pick = RNG.choice(pool.shape[0], size=M.shape[0], replace=False)
                r = mean_abs_rho(pool[pick])
                if np.isfinite(r):
                    null.append(r)
            null = np.array(null)
            summary["module_coherence"] = {
                "cohort": cname, "n_events": int(M.shape[0]), "n_tumours": int(M.shape[1]),
                "pool_events": int(pool.shape[0]),
                "mean_abs_rho_observed": obs,
                "null_mean": float(null.mean()) if len(null) else None,
                "null_p95": float(np.quantile(null, 0.95)) if len(null) else None,
                "empirical_p": float((null >= obs).mean()) if len(null) else None,
            }
            pd.DataFrame(np.corrcoef(M), index=ids_k, columns=ids_k).to_csv(
                OUT_T / "tcga_caf_module_correlation.csv")

    OUT_P.mkdir(parents=True, exist_ok=True)
    (OUT_P / "tcga_caf_events_test.json").write_text(
        json.dumps(summary, indent=1, ensure_ascii=False), encoding="utf-8")

    print()
    for cname, b in summary["tumour_vs_normal"].items():
        print(f"[TvN {cname}] pairs={b['pairs_tested']} events={b['events_tested']} min_p={b['min_p']:.3g} "
              f"min_q={b['min_q']:.3g} pairs_q<0.05={b['n_pairs_q_lt_0_05']} "
              f"p<0.05={b['n_p_lt_0_05']} direction={b['direction_concordant']}/{b['direction_n']} "
              f"(p={b['direction_binom_p']:.3g}) max|delta|={b['max_abs_delta']:.3f} "
              f"|delta|>=0.05: {b['abs_delta_ge_0_05']} ({b['abs_delta_ge_0_05_concordant']} concordant)")
    for cname, b in summary["prognosis"].items():
        print(f"[Cox {cname}] tested={b['events_tested']} min_p={b['min_p']:.3g} "
              f"min_q={b['min_q']:.3g} q<0.05={b['n_q_lt_0_05']} "
              f"direction={b['direction_concordant']}/{b['direction_n']} "
              f"(p={b['direction_binom_p']:.3g})")
    if "module_coherence" in summary:
        mb = summary["module_coherence"]
        print(f"\nmodule: mean|rho|={mb['mean_abs_rho_observed']:.3f} vs null "
              f"{mb['null_mean']:.3f} (p95 {mb['null_p95']:.3f}), "
              f"empirical p={mb['empirical_p']:.3f}")


if __name__ == "__main__":
    main()
