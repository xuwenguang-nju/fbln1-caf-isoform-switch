"""Composition correlation for the FBLN1 AT events, per TCGA-HNSC cohort.

Replaces the legacy `rho_fibro` / `rho_epi` columns of fbln1_tcga_events.csv,
which correlated event PSI against marker-gene splice-event PSI rather than
marker expression and came out with inverted signs (see
results/reports/dtu_and_composition.md, section 5).

Here the composition axis is the NNLS CAF / epithelial fraction produced by
148_tcga_deconvolution.py, and the correlation is reported *signed* because the
two AT events are mirrors of one 3'-end choice.

Cohorts are the same three used elsewhere in the figure set: oral_strict
(tongue, floor of mouth, gum, palate, lip, other mouth), oral_wide (that plus
"other and ill-defined sites in lip, oral cavity and pharynx"), and all_hnsc.

Outputs
-------
results/tables/fbln1_tcga_composition_by_cohort.csv
data/processed/fbln1_composition_by_cohort.json
"""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
PSI_ZIP = ROOT / "data/raw/tcga_spliceseq/full/HNSC_full.zip"
CLIN = ROOT / "data/processed/tcga_hnsc_clinical.csv"
FRAC = ROOT / "results/tables/tcga_caf_fraction.csv"
T = ROOT / "results/tables"
P = ROOT / "data/processed"

EVENTS = {"62665": "exon 17", "62666": "exon 22"}
ORAL_STRICT = ["Other and unspecified parts of tongue", "Floor of mouth", "Gum",
               "Palate", "Lip", "Other and unspecified parts of mouth"]
ORAL_WIDE = ORAL_STRICT + ["Other and ill-defined sites in lip, oral cavity and pharynx"]


def fbln1_psi() -> pd.DataFrame:
    """Tumour-sample FBLN1 PSI, one row per (event, sample)."""
    with zipfile.ZipFile(PSI_ZIP) as z:
        with z.open(z.namelist()[0]) as fh:
            hdr = (fh.readline().decode("utf-8", "replace").replace("\r", "")
                   .rstrip("\n").split("\t"))
            smp = [c for c in hdr[10:] if c and c != "avg_exp_pct"]
            keep = [i for i, c in enumerate(hdr) if c in smp]
            recs = []
            for line in fh:
                f = line.decode("utf-8", "replace").replace("\r", "").rstrip("\n").split("\t")
                if not f or f[0] != "FBLN1":
                    continue
                if len(f) < len(hdr):
                    f += [""] * (len(hdr) - len(f))
                vals = [np.nan if f[i] in ("", "null", "NA", "nan") else float(f[i])
                        for i in keep]
                for c, v in zip(smp, vals):
                    recs.append((f[1], c, v))
    d = pd.DataFrame(recs, columns=["as_id", "psi_sample", "psi"]).dropna(subset=["psi"])
    d = d[d["as_id"].isin(EVENTS)]
    d["case"] = d["psi_sample"].str.replace("_", "-", regex=False)
    d = d[~d["case"].str.endswith("-Norm")]
    return d


def main() -> None:
    psi = fbln1_psi()
    fr = pd.read_csv(FRAC, index_col=0)
    fr["case"] = ["-".join(i.split(".")[:3]) for i in fr.index]
    fr["stype"] = [i.split(".")[-1] for i in fr.index]
    fr = fr[fr["stype"] == "01"].drop_duplicates("case")
    clin = pd.read_csv(CLIN).set_index("case")

    cohorts = {
        "oral_strict": [c for c in fr["case"] if clin["site"].get(c) in ORAL_STRICT],
        "oral_wide": [c for c in fr["case"] if clin["site"].get(c) in ORAL_WIDE],
        "all_hnsc": list(fr["case"]),
    }

    rows = []
    for name, cases in cohorts.items():
        cs = set(cases)
        sub = psi[psi["case"].isin(cs)].merge(
            fr[["case", "CAF", "epithelial"]], on="case", how="inner")
        for aid, lab in EVENTS.items():
            g = sub[sub["as_id"] == aid]
            rc = stats.spearmanr(g["psi"], g["CAF"])
            re_ = stats.spearmanr(g["psi"], g["epithelial"])
            rows.append(dict(
                cohort=name, as_id=aid, exon=lab, n=len(g),
                rho_caf=float(rc.statistic), p_caf=float(rc.pvalue),
                r2_caf=float(rc.statistic ** 2),
                rho_epithelial=float(re_.statistic), p_epithelial=float(re_.pvalue)))
    out = pd.DataFrame(rows)
    out.to_csv(T / "fbln1_tcga_composition_by_cohort.csv", index=False)
    print(out.round(4).to_string(index=False))

    info = dict(
        composition_axis="NNLS CAF / epithelial fraction (148_tcga_deconvolution.py)",
        n_tumour_per_cohort={k: len(v) for k, v in cohorts.items()},
        rho_by_cohort={
            f"{r['cohort']}|{r['as_id']}": {
                "n": int(r["n"]), "rho_caf": round(r["rho_caf"], 4),
                "p_caf": float(f"{r['p_caf']:.3g}"),
                "rho_epithelial": round(r["rho_epithelial"], 4),
                "p_epithelial": float(f"{r['p_epithelial']:.3g}")}
            for r in out.to_dict("records")},
        supersedes="rho_fibro / rho_epi in fbln1_tcga_events.csv (inverted sign; "
                   "see results/reports/dtu_and_composition.md section 5)",
    )
    (P / "fbln1_composition_by_cohort.json").write_text(
        json.dumps(info, indent=1, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
