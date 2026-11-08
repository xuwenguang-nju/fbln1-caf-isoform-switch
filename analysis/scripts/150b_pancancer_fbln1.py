"""Task 5: does the FBLN1 3'-end (1C/1D) usage-stroma relationship generalise?

Why the metric is the AT event and not "1D fraction"
----------------------------------------------------
GENCODE v44 places FBLN1 on the plus strand. The two isoforms that differ by
20 aa (ENST00000262722 = 1C, 683 aa; ENST00000327858 = 1D, 703 aa) share their
first exon (chr22:45 502 883-45 503 064) and diverge only at the 3' end: 1C
terminates in a 451-nt exon (45 562 912-45 563 362) while 1D runs on through
45 574 511 / 45 576 977 into an 829-nt terminal exon (45 600 307-45 601 135).
The choice is therefore an alternative-terminal-exon / alternative-
polyadenylation event, i.e. exactly SpliceSeq's `AT` class -- which is why
`AT` PSI is the TCGA-measurable proxy for the 1C/1D choice. Nothing here claims
to measure the 1C:1D ratio itself; SpliceSeq cannot, because it only reports
predefined event PSI and the two isoforms are not an annotated exon-skip event.

Design (fixed before running)
-----------------------------
For each of HNSC / PAAD / BRCA / LUAD / ESCA:
  1. FBLN1 AT PSI (62665 primary, 62666 secondary) from the unfiltered SpliceSeq
     export (`full/*_full.zip`, pctwithval=0 -- the filtered export drops 62665);
  2. CAF fraction per sample by the same NNLS signature used in
     `148_tcga_deconvolution.py` (no re-tuning per cancer);
  3. three quantities: within-tumour Spearman rho(psi, CAF); tumour-vs-normal
     median delta psi; and, where >= 8 matched patients exist, the delta psi vs
     delta CAF regression intercept (the composition-adjusted effect).
  4. Fisher-z meta-analysis of rho across cancers.

Outputs
-------
results/tables/pancancer_fbln1_events.csv    per cancer x event statistics
results/tables/pancancer_caf_fraction.csv    per sample CAF fraction, all cancers
data/processed/pancancer_fbln1.json
"""

from __future__ import annotations

import importlib.util
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
T = ROOT / "results/tables"
P = ROOT / "data/processed"
PSI_DIR = ROOT / "data/raw/tcga_spliceseq/full"
XENA_DIR = ROOT / "data/raw/tcga_xena"

CANCERS = ["HNSC", "PAAD", "BRCA", "LUAD", "ESCA"]
PRIMARY = "62665"      # AT exon 17
SECONDARY = "62666"    # AT exon 22


def _load_module():
    """Reuse the signature/deconvolution code from task 2 verbatim."""
    path = Path(__file__).with_name("148_tcga_deconvolution.py")
    spec = importlib.util.spec_from_file_location("deconv148", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


d148 = _load_module()


def read_fbln1_stream(cancer: str) -> pd.DataFrame:
    """Stream the export and keep only FBLN1 rows (files reach 1 GB+)."""
    zpath = PSI_DIR / f"{cancer}_full.zip"
    with zipfile.ZipFile(zpath) as z:
        name = z.namelist()[0]
        with z.open(name) as fh:
            hdr = (fh.readline().decode("utf-8", "replace").replace("\r", "")
                   .rstrip("\n").split("\t"))
            smp = [c for c in hdr[10:] if c and c != "avg_exp_pct"]
            keep = [i for i, c in enumerate(hdr) if c in smp]
            recs = []
            for line in fh:
                # these exports are CRLF; a bare rstrip("\n") leaves "\r" on the
                # last field and every "null\r" then fails float()
                f = line.decode("utf-8", "replace").replace("\r", "").rstrip("\n").split("\t")
                if not f or f[0] != "FBLN1":
                    continue
                if len(f) < len(hdr):          # short rows occur in these exports
                    f = f + [""] * (len(hdr) - len(f))
                vals = [np.nan if f[i] in ("", "null", "NA", "nan") else float(f[i])
                        for i in keep]
                for c, v in zip(smp, vals):
                    recs.append((f[1], f[2], f[3], f[7], f[9], c, v))
    out = pd.DataFrame(recs, columns=["as_id", "type", "exons", "pct_with_values",
                                      "std_psi", "psi_sample", "psi"])
    return out.dropna(subset=["psi"])


def sample_keys(s: pd.Series) -> pd.DataFrame:
    """TCGA_XX_YYYY / TCGA_XX_YYYY_Norm -> (patient, '01'/'11')."""
    t = s.str.replace("_", "-", regex=False)
    is_norm = t.str.endswith("-Norm")
    pat = t.str.replace("-Norm", "", regex=False)
    return pd.DataFrame({"patient": pat, "sample_type": np.where(is_norm, "11", "01")},
                        index=s.index)


def caf_fraction(cancer: str, gene: pd.DataFrame, sig: dict[str, list[str]]):
    x = pd.read_csv(XENA_DIR / f"{cancer}_HiSeqV2.gz", sep="\t", index_col=0)
    x.index = x.index.astype(str)
    x = np.power(2.0, x) - 1.0
    x.columns = [c.replace("-", ".") for c in x.columns]
    common = [g for g in sorted({m for ct in d148.CELL_TYPES for m in sig[ct]})
              if g in x.index and g in gene.index]
    A = gene.loc[common, d148.CELL_TYPES].to_numpy()
    frac = np.array([d148.nnls_fraction(A, x.loc[common, s].to_numpy())
                     for s in x.columns])
    fr = pd.DataFrame(frac, index=x.columns, columns=d148.CELL_TYPES)
    keys = pd.DataFrame({"patient": ["-".join(c.split(".")[:3]) for c in fr.index],
                         "sample_type": [c.split(".")[-1] for c in fr.index]},
                        index=fr.index)
    return fr.join(keys), len(common)


def main() -> None:
    gene, sig = d148.load_signature()
    rows, per_sample = [], []
    for cancer in CANCERS:
        psi = read_fbln1_stream(cancer)
        keys = sample_keys(psi["psi_sample"])
        psi = psi.join(keys)
        fr, n_markers = caf_fraction(cancer, gene, sig)
        fr = fr.drop_duplicates(["patient", "sample_type"])
        d = psi.merge(fr.reset_index().rename(columns={"index": "xena_id"}),
                      on=["patient", "sample_type"], how="inner")
        print(f"\n{cancer}: markers={n_markers}  psi rows={len(psi):,}  "
              f"joined={len(d):,}  samples-with-CAF={len(fr):,}")
        per_sample.append(d.assign(cancer=cancer))

        for aid in [PRIMARY, SECONDARY, "62663", "62664", "62668", "96227"]:
            g = d[d["as_id"] == aid]
            if g.empty:
                continue
            tum = g[g["sample_type"] == "01"]
            nor = g[g["sample_type"] == "11"]
            if len(tum) < 20:
                continue
            rho = stats.spearmanr(tum["psi"], tum["CAF"])
            # paired: patients with both a tumour and a normal
            t = tum.set_index("patient")
            n = nor.set_index("patient")
            both = t.index.intersection(n.index)
            dpsi = t.loc[both, "psi"] - n.loc[both, "psi"]
            dcaf = t.loc[both, "CAF"] - n.loc[both, "CAF"]
            inter = p_int = slope = np.nan
            if len(both) >= 8:
                lr = stats.linregress(dcaf.to_numpy(), dpsi.to_numpy())
                inter, p_int, slope = lr.intercept, lr.pvalue, lr.slope
            wil = np.nan
            if len(both) >= 6:
                try:
                    wil = stats.wilcoxon(dpsi).pvalue
                except ValueError:
                    wil = np.nan
            rows.append(dict(
                cancer=cancer, as_id=aid, type=g["type"].iloc[0],
                exons=g["exons"].iloc[0], std_psi_all=float(g["std_psi"].iloc[0]),
                n_tumour=len(tum), n_normal=len(nor), n_paired=len(both),
                median_psi_tumour=float(tum["psi"].median()),
                median_psi_normal=(float(nor["psi"].median()) if len(nor) else np.nan),
                median_delta_psi=(float(dpsi.median()) if len(both) else np.nan),
                n_pos_delta=int((dpsi > 0).sum()), wilcoxon_p=wil,
                rho_psi_caf_tumour=float(rho.statistic), p_rho=float(rho.pvalue),
                median_caf_tumour=float(tum["CAF"].median()),
                median_delta_caf=(float(dcaf.median()) if len(both) else np.nan),
                slope_dpsi_dcaf=slope, intercept_dpsi=inter, p_intercept=p_int))

    res = pd.DataFrame(rows)
    res.to_csv(T / "pancancer_fbln1_events.csv", index=False)
    pd.concat(per_sample, ignore_index=True).to_csv(
        T / "pancancer_caf_fraction.csv", index=False)
    print("\n" + res[res["as_id"] == PRIMARY].to_string(index=False))
    print("\n" + res[res["as_id"] == SECONDARY].to_string(index=False))

    # Fisher-z meta-analysis of rho(psi, CAF) across cancers, per event
    meta = {}
    for aid in [PRIMARY, SECONDARY]:
        sub = res[(res["as_id"] == aid) & res["rho_psi_caf_tumour"].notna()]
        if len(sub) < 2:
            continue
        z = np.arctanh(sub["rho_psi_caf_tumour"].clip(-0.999, 0.999))
        # inverse-variance weights for a Spearman rho: n - 3 per study
        w = np.clip(sub["n_tumour"].to_numpy(dtype=float) - 3.0, 1.0, None)
        zbar = float(np.sum(z * w) / w.sum())
        se = float(np.sqrt(1 / w.sum()))
        meta[aid] = dict(n_cancers=len(sub), rho_pooled=float(np.tanh(zbar)),
                         ci_lo=float(np.tanh(zbar - 1.96 * se)),
                         ci_hi=float(np.tanh(zbar + 1.96 * se)),
                         p=float(2 * stats.norm.sf(abs(zbar / se))),
                         n_same_sign=int((sub["rho_psi_caf_tumour"] > 0).sum()))
        print(f"\nmeta {aid}: pooled rho={meta[aid]['rho_pooled']:+.3f} "
              f"({meta[aid]['ci_lo']:+.3f} to {meta[aid]['ci_hi']:+.3f}) "
              f"p={meta[aid]['p']:.2e}; positive in "
              f"{meta[aid]['n_same_sign']}/{meta[aid]['n_cancers']} cancers")

    info = dict(
        metric_note=(
            "SpliceSeq reports event PSI only; the 1C/1D choice is an alternative "
            "terminal exon (plus strand: 1C ends at chr22:45 563 362, 1D continues "
            "into an 829-nt terminal exon at chr22:45 600 307-45 601 135), so AT "
            "PSI is a proxy, not the 1C:1D ratio"),
        events=dict(primary=PRIMARY, secondary=SECONDARY),
        cancers=CANCERS,
        source_psi="data/raw/tcga_spliceseq/full/<cancer>_full.zip (pctwithval=0)",
        source_expr="data/raw/tcga_xena/<cancer>_HiSeqV2.gz",
        signature="data/processed/server_restore/expr/celltype_tpm.tsv (same as task 2)",
        meta_rho=meta,
        caveats=[
            "PAAD (4) and ESCA (13) have too few normals for a paired design",
            "CAF fractions are relative (no immune reference), so cross-cancer "
            "comparisons are rank-based only",
            "SpliceSeq PSI is precomputed by its own pipeline; it is not "
            "reproducible from open junction files in this project",
        ],
    )
    (P / "pancancer_fbln1.json").write_text(
        json.dumps(info, indent=1, ensure_ascii=False), encoding="utf-8")
    print("\nwritten:", T / "pancancer_fbln1_events.csv")


if __name__ == "__main__":
    main()
