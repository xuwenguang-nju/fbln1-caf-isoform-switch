"""Stage 8 step 1: map the 27 CAF-derived splicing events onto TCGA-HNSC bulk and test them.

Why this differs from the earlier TCGA work
-------------------------------------------
`39_psi_pilot_spliceseq.py` screened all 2,787 testable SpliceSeq events against survival
in the oral cohort and found nothing at FDR (min q = 0.56). That was a whole-genome screen
with 101 deaths: the multiplicity penalty alone made FDR unreachable, not an absence of
biology.

This script instead tests a **pre-specified set of 27 events** derived independently from
single-cell data (CAF vs epithelial, reproducible across four donors). Twenty-seven tests
can clear FDR; 2,787 cannot. The hypothesis was fixed before looking at TCGA survival, so
this is a targeted validation, not a data-driven hunt.

Mapping
-------
SUPPA2 event types map onto SpliceSeq types as SE->ES, A5->AD, A3->AA, MX->ME, RI->RI.
SpliceSeq identifies events by exon number, not genomic coordinate, so a match is
gene x splice-type rather than exon-exact. Where a gene has several SpliceSeq events of
the same type, all are retained and the multiplicity is reported.

Outputs
-------
results/tables/tcga_caf_events_mapping.csv     one row per (CAF event, SpliceSeq event)
results/tables/tcga_caf_events_cox.csv         Cox result per mapped SpliceSeq event
data/processed/tcga_caf_events_test.json       summary + test statistics
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.duration.hazard_regression import PHReg

ROOT = Path(__file__).resolve().parents[2]
EVENTS = ROOT / "remote/easyai/stringent_events.tsv"
PSI = ROOT / "data/raw/tcga_spliceseq/PSI_download_HNSC.txt"
CLIN = ROOT / "data/processed/tcga_hnsc_clinical.csv"
OUT_T = ROOT / "results/tables"
OUT_P = ROOT / "data/processed"

# SUPPA2 -> SpliceSeq splice-type family
TYPE_MAP = {"SE": "ES", "A5": "AD", "A3": "AA", "MX": "ME", "RI": "RI"}

ORAL_STRICT = [
    "Other and unspecified parts of tongue",
    "Floor of mouth",
    "Gum",
    "Palate",
    "Other and unspecified parts of mouth",
    "Other and unspecified parts of lip",
]
ORAL_WIDE = ORAL_STRICT + ["Other and ill-defined sites in lip, oral cavity and pharynx"]

META_COLS = [
    "symbol", "as_id", "splice_type", "exons", "from_exon", "to_exon",
    "novel_splice", "pct_with_values", "psi_range", "std_psi",
]


def load_caf_events() -> pd.DataFrame:
    ev = pd.read_csv(EVENTS, sep="\t")
    ev["ss_type"] = ev["type"].map(TYPE_MAP)
    missing = ev.loc[ev["ss_type"].isna(), "type"].unique()
    if len(missing):
        raise ValueError(f"unmapped SUPPA2 types: {missing}")
    return ev


def load_psi() -> tuple[pd.DataFrame, list[str]]:
    psi = pd.read_csv(PSI, sep="\t", low_memory=False)
    patients = [c for c in psi.columns if c not in META_COLS]
    psi = psi[META_COLS + patients].copy()
    # SpliceSeq barcodes look like TCGA_4P_AA8J; clinical case ids use dashes.
    norm = {p: p.replace("_", "-") for p in patients}
    psi = psi.rename(columns=norm)
    return psi, [norm[p] for p in patients]


def map_events(ev: pd.DataFrame, psi: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in ev.iterrows():
        hits = psi[(psi["symbol"] == r["symbol"]) & (psi["splice_type"] == r["ss_type"])]
        for _, h in hits.iterrows():
            rows.append(dict(
                caf_symbol=r["symbol"], caf_type=r["type"], caf_event=r["event"],
                caf_mean_dpsi=r["mean"], caf_sd=r["sd"], caf_p=r["p"],
                ss_as_id=h["as_id"], ss_exons=h["exons"],
                ss_pct_with_values=h["pct_with_values"],
            ))
    return pd.DataFrame(rows)


def fit_cox(psi_row: pd.Series, clin: pd.DataFrame, cases: set[str]) -> dict | None:
    vals = psi_row[[c for c in psi_row.index if c in cases]]
    vals = pd.to_numeric(vals, errors="coerce")
    df = pd.DataFrame({"psi": vals}).join(clin.set_index("case"), how="inner")
    df = df.dropna(subset=["psi", "os_time", "event"])
    if len(df) < 40 or df["event"].sum() < 15:
        return None
    sd = df["psi"].std()
    if not np.isfinite(sd) or sd <= 0:
        return None
    x = ((df["psi"] - df["psi"].mean()) / sd).to_numpy()
    model = PHReg(df["os_time"].to_numpy(float), x, status=df["event"].to_numpy(float))
    res = model.fit()
    log_hr = float(res.params[0])
    se = float(res.bse[0])
    z = log_hr / se
    return dict(n_used=int(len(df)), n_deaths=int(df["event"].sum()),
                psi_sd=float(sd), hr_per_sd=float(np.exp(log_hr)),
                log_hr=log_hr, se=se, z=z,
                p=float(2 * stats.norm.sf(abs(z))))


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


def main() -> None:
    ev = load_caf_events()
    psi, patients = load_psi()
    clin = pd.read_csv(CLIN)

    mapped = map_events(ev, psi)
    print(f"CAF events: {len(ev)}")
    print(f"mapped to >=1 SpliceSeq event: {mapped['caf_event'].nunique()} "
          f"({len(mapped)} pairs)")
    unmapped = sorted(set(ev["symbol"]) - set(mapped["caf_symbol"]))
    if unmapped:
        print(f"symbols with no SpliceSeq entry: {unmapped}")

    OUT_T.mkdir(parents=True, exist_ok=True)
    mapped.to_csv(OUT_T / "tcga_caf_events_mapping.csv", index=False)

    cohorts = {"oral_strict": set(clin.loc[clin["site"].isin(ORAL_STRICT), "case"]),
               "oral_wide": set(clin.loc[clin["site"].isin(ORAL_WIDE), "case"]),
               "all_hnsc": set(clin["case"])}

    results = []
    for name, cases in cohorts.items():
        for _, m in mapped.iterrows():
            row = psi[psi["as_id"] == m["ss_as_id"]].iloc[0]
            fit = fit_cox(row, clin, cases)
            if fit is None:
                continue
            results.append(dict(cohort=name, caf_symbol=m["caf_symbol"],
                                caf_type=m["caf_type"], caf_mean_dpsi=m["caf_mean_dpsi"],
                                ss_as_id=m["ss_as_id"], ss_type=row["splice_type"],
                                **fit))
    cox = pd.DataFrame(results)

    # FDR within each cohort, over the pre-specified mapped set only
    for name, grp in cox.groupby("cohort"):
        cox.loc[grp.index, "q"] = bh(grp["p"].to_numpy())
    cox["sig_fdr"] = cox["q"] < 0.05
    cox = cox.sort_values(["cohort", "p"])
    cox.to_csv(OUT_T / "tcga_caf_events_cox.csv", index=False)

    summary = {"n_caf_events": int(len(ev)),
               "n_caf_events_mapped": int(mapped["caf_event"].nunique()),
               "n_pairs": int(len(mapped)),
               "unmapped_symbols": unmapped,
               "cohorts": {}}

    for name, grp in cox.groupby("cohort"):
        block = {
            "events_tested": int(len(grp)),
            "n_deaths": int(grp["n_deaths"].iloc[0]),
            "min_p": float(grp["p"].min()),
            "min_q": float(grp["q"].min()),
            "n_q_lt_0_05": int(grp["sig_fdr"].sum()),
        }

        # Directional prediction, fixed before looking at the result:
        # a positive scRNA dPSI means CAF PSI > epithelial PSI, so in a CAF-rich tumour the
        # bulk PSI should be higher. If CAF content is adverse in HNSCC, that predicts HR>1.
        expected = np.sign(grp["caf_mean_dpsi"].to_numpy())
        observed = np.sign(grp["log_hr"].to_numpy())
        concordant = int((expected == observed).sum())
        binom = stats.binomtest(concordant, len(grp), 0.5, alternative="greater")
        block["direction_concordant"] = concordant
        block["direction_n"] = int(len(grp))
        block["direction_binom_p"] = float(binom.pvalue)

        # Enrichment vs the full SpliceSeq screen: are the pre-specified events more
        # strongly prognostic than the background of all testable events?
        block["median_abs_z_caf"] = float(grp["z"].abs().median())
        summary["cohorts"][name] = block

    # background comparison using every SpliceSeq event, same filters as the earlier screen
    bg_rows = []
    for _, row in psi.iterrows():
        fit = fit_cox(row, clin, cohorts["oral_strict"])
        if fit is not None:
            bg_rows.append(fit["z"])
    bg_z = np.abs(np.array(bg_rows))
    caf_z = np.abs(cox.loc[cox["cohort"] == "oral_strict", "z"].to_numpy())
    if len(caf_z) > 5 and len(bg_z) > 100:
        mw = stats.mannwhitneyu(caf_z, bg_z, alternative="greater")
        summary["background_vs_caf_oral"] = {
            "n_background_events": int(len(bg_z)),
            "background_median_abs_z": float(np.median(bg_z)),
            "caf_median_abs_z": float(np.median(caf_z)),
            "mannwhitney_p": float(mw.pvalue),
        }

    OUT_P.mkdir(parents=True, exist_ok=True)
    (OUT_P / "tcga_caf_events_test.json").write_text(
        json.dumps(summary, indent=1, ensure_ascii=False), encoding="utf-8")

    print()
    for name, block in summary["cohorts"].items():
        print(f"[{name}] tested={block['events_tested']} deaths={block['n_deaths']} "
              f"min_p={block['min_p']:.3g} min_q={block['min_q']:.3g} "
              f"q<0.05={block['n_q_lt_0_05']} "
              f"direction={block['direction_concordant']}/{block['direction_n']} "
              f"(p={block['direction_binom_p']:.3g})")
    if "background_vs_caf_oral" in summary:
        b = summary["background_vs_caf_oral"]
        print(f"\nvs background (oral_strict): CAF median |z|={b['caf_median_abs_z']:.3f} "
              f"vs {b['background_median_abs_z']:.3f} over {b['n_background_events']} events, "
              f"Mann-Whitney p={b['mannwhitney_p']:.3g}")


if __name__ == "__main__":
    main()
