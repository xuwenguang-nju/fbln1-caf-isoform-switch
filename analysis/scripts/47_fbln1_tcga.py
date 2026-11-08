"""Stage 10: does FBLN1 alternative splicing track the CAF axis in TCGA-HNSC bulk?

Why this is a separate analysis
------------------------------
`45_tcga_bulk_validation.py` tested only the 27 single-cell events. FBLN1 is the
main line of the project (B arm: 1C->1D switch, 6/6 patients, q_panel=0.007) yet
it was never tested in bulk. Its SpliceSeq events are known:

  62663 AT exons=16  1.000      62667 RI exons=18.2  0.4945
  62664 AT exons=18.4 1.000     62668 ES exons=3     0.9743
  62665 AT exons=17  1.000      96227 ES exons=13    0.9651
  62666 AT exons=22  1.000

Design
------
Q1 tumour vs normal (43 normals, 42 paired). Normal oral mucosa is mostly
   epithelium, so if FBLN1 splicing marks fibroblasts/CAF the tumour should shift
   toward the CAF direction. (Caveat: normal mucosa still has lamina-propria
   fibroblasts, so this is a diluted contrast, not CAF-vs-epithelium.)
Q2 prognosis (Cox on continuous PSI per SD).
Q3 CAF-content proxy: correlate each FBLN1 event's PSI across tumours with
   (a) a fibroblast marker score and (b) an epithelial marker score. If the
   splicing change is a stromal-composition readout, PSI should track the
   fibroblast score. This is the closest bulk proxy for "is this CAF-driven"
   available without deconvolution.

All three cohorts (oral_strict / oral_wide / all_hnsc) as before.
"""
from __future__ import annotations
import json, warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.duration.hazard_regression import PHReg

ROOT = Path(__file__).resolve().parents[2]
PSI = ROOT / "data/raw/tcga_spliceseq/full/PSI_download_HNSC.txt"
CLIN = ROOT / "data/processed/tcga_hnsc_clinical.csv"
OUT_T = ROOT / "results/tables"
OUT_P = ROOT / "data/processed"

GENE = "FBLN1"
META_N = 10
MIN_FRAC, MIN_SD = 0.50, 0.03
N_PERM = 2000

ORAL_STRICT = ["Other and unspecified parts of tongue", "Floor of mouth", "Gum",
               "Palate", "Other and unspecified parts of mouth",
               "Other and unspecified parts of lip"]
ORAL_WIDE = ORAL_STRICT + ["Other and ill-defined sites in lip, oral cavity and pharynx"]

# marker panels for the composition check
FIB_MARKERS = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "PDGFRB",
               "FAP", "THY1", "POSTN", "FN1", "MMP2"]
EPI_MARKERS = ["KRT5", "KRT14", "KRT13", "KRT4", "CDH1", "EPCAM", "SFN", "KRT6A"]

RNG = np.random.default_rng(20261021)
warnings.filterwarnings("ignore")


def bh(p):
    m = len(p); order = np.argsort(p); q = np.empty(m); prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]; prev = min(prev, p[i] * m / (rank + 1)); q[i] = prev
    return q


def main():
    with open(PSI, encoding="utf-8", errors="replace") as fh:
        hdr = [c for c in fh.readline().rstrip("\n").split("\t") if c]
    pcols = hdr[10:]
    tumours = [c for c in pcols if not c.endswith("_Norm")]
    normals = [c for c in pcols if c.endswith("_Norm")]
    as_case = lambda c: c.replace("_Norm", "").replace("_", "-")
    clin = pd.read_csv(CLIN).set_index("case")

    # --- FBLN1 events -----------------------------------------------------------------
    want = set()
    fbln1 = []
    with open(PSI, encoding="utf-8", errors="replace") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(h) if c}
        for line in fh:
            if not line.startswith(GENE + "\t"):
                continue
            f = line.rstrip("\n").split("\t")
            fbln1.append(dict(as_id=f[idx["as_id"]], type=f[idx["splice_type"]],
                              exons=f[idx["exons"]], pct=float(f[idx["pct_with_values"]])))
            want.add(f[idx["as_id"]])
    print("FBLN1 events: %d -> %s" % (len(fbln1), [(e["as_id"], e["type"], e["exons"]) for e in fbln1]))

    # --- FBLN1 PSI + marker gene expression -------------------------------------------
    psi_rows, marker_rows = {}, {}
    genes_of_interest = {GENE} | set(FIB_MARKERS) | set(EPI_MARKERS)
    with open(PSI, encoding="utf-8", errors="replace") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(h) if c}
        keep = [idx[c] for c in pcols]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= keep[-1]:
                continue
            sym = f[idx["symbol"]]
            if sym != GENE and sym not in genes_of_interest:
                continue
            vals = np.array([np.nan if f[p] in ("", "null", "NA", "nan") else float(f[p]) for p in keep])
            if sym == GENE:
                psi_rows[f[idx["as_id"]]] = vals
            else:
                marker_rows[f[idx["symbol"]]] = vals
    print("FBLN1 PSI rows: %d | marker genes found: %d" % (len(psi_rows), len(marker_rows)))

    pm = pd.DataFrame(psi_rows, index=pcols).T          # events x samples
    mk = pd.DataFrame(marker_rows, index=pcols).T

    cohort_cases = {
        "oral_strict": [c for c in tumours if clin["site"].get(as_case(c)) in ORAL_STRICT],
        "oral_wide": [c for c in tumours if clin["site"].get(as_case(c)) in ORAL_WIDE],
        "all_hnsc": list(tumours),
    }

    out = {"gene": GENE, "events": fbln1, "cohorts": {}}

    for cname, cidx in cohort_cases.items():
        print("\n=== %s: %d tumours / %d normals ===" % (cname, len(cidx), len(normals)))
        recs = []
        for aid, row in pm.iterrows():
            t = row[cidx].to_numpy(float)
            if np.isfinite(t).mean() < MIN_FRAC or np.nanstd(t) < MIN_SD:
                continue
            pairs = [(c, c + "_Norm") for c in cidx if (c + "_Norm") in normals]
            pairs = [(a, b) for a, b in pairs if np.isfinite(row[a]) and np.isfinite(row[b])]
            if len(pairs) >= 8:
                d = np.array([row[a] - row[b] for a, b in pairs])
                stat, p = stats.wilcoxon(d); test, delta = "wilcoxon_paired", float(d.mean())
                nt = nn = len(pairs)
            else:
                nv = row[normals].to_numpy(float); nv = nv[np.isfinite(nv)]
                tv = t[np.isfinite(t)]
                if len(nv) < 5 or len(tv) < 20:
                    continue
                stat, p = stats.mannwhitneyu(tv, nv, alternative="two-sided")
                test, delta = "mannwhitney", float(tv.mean() - nv.mean())
                nt, nn = len(tv), len(nv)

            # prognosis
            j = clin.reindex([as_case(c) for c in cidx])
            ok = np.isfinite(t) & j["os_time"].notna().to_numpy() & j["event"].notna().to_numpy()
            hr = pcox = np.nan
            if ok.sum() >= 40 and j["event"].to_numpy(float)[ok].sum() >= 15:
                x = t[ok]; x = (x - x.mean()) / x.std()
                try:
                    r = PHReg(j["os_time"].to_numpy(float)[ok], x,
                              status=j["event"].to_numpy(float)[ok]).fit()
                    hr = float(np.exp(r.params[0]))
                    pcox = float(2 * stats.norm.sf(abs(r.params[0] / r.bse[0])))
                except Exception:
                    pass

            # composition proxy: Spearman of PSI with fibroblast / epithelial score.
            # Use pairwise-complete observations; a blanket isfinite().all() guard
            # silently returns NaN whenever one sample lacks a value.
            def score(panel):
                g = [x for x in panel if x in mk.index]
                if not g:
                    return None
                return np.nanmean(np.log1p(mk.loc[g, cidx].to_numpy(float)), axis=0)
            fs, es = score(FIB_MARKERS), score(EPI_MARKERS)

            def sp(a, b):
                ok = np.isfinite(a) & np.isfinite(b)
                if ok.sum() < 20:
                    return np.nan
                return float(stats.spearmanr(a[ok], b[ok]).correlation)
            rf = sp(t, fs) if fs is not None else np.nan
            re_ = sp(t, es) if es is not None else np.nan

            ev = next(e for e in fbln1 if e["as_id"] == aid)
            recs.append(dict(cohort=cname, as_id=aid, type=ev["type"], exons=ev["exons"],
                             test=test, n_tumour=nt, n_normal=nn, mean_delta=delta,
                             p=float(p), hr_per_sd=hr, p_cox=pcox,
                             rho_fibro=rf, rho_epi=re_))

        df = pd.DataFrame(recs)
        if len(df):
            df["q"] = bh(df["p"].to_numpy())
            df["q_cox"] = bh(df["p_cox"].fillna(1).to_numpy())
            print(df.round(4).to_string(index=False))
            out["cohorts"][cname] = df.to_dict("records")

    if out["cohorts"]:
        allrows = pd.concat([pd.DataFrame(v) for v in out["cohorts"].values()])
        allrows.to_csv(OUT_T / "fbln1_tcga_events.csv", index=False)

        # random-gene null for the composition correlation, so rho_fibro is interpretable
        # (a pure-depth artefact would also correlate; a CAF-specific one should exceed it)
        print("\n=== composition-correlation context (oral_strict) ===")
        sub = allrows[allrows["cohort"] == "oral_strict"]
        if len(sub):
            print("  FBLN1 rho_fibro: %s" % sub["rho_fibro"].round(3).tolist())
            print("  FBLN1 rho_epi:   %s" % sub["rho_epi"].round(3).tolist())

    OUT_P.mkdir(parents=True, exist_ok=True)
    (OUT_P / "fbln1_tcga.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
