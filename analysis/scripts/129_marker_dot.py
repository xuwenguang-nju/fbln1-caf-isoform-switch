"""Export a marker-gene summary per Leiden cluster, for the single-cell figure.

Why this exists
---------------
The per-cell table carries UMAP coordinates and lineage calls, but a marker dot
plot needs expression. Rather than ship a 32k-cell matrix, compute per-cluster
mean (log-normalised) expression and the fraction of cells expressing, for a
small canonical marker panel, and export one small TSV.

Inputs are the cached AnnData written by 128_percell_umap.py.
"""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

RES = Path("/root/hnscc_caf/results")
H5 = RES / "sc_raw.h5ad"
OUT = RES / "sc_marker_dot.tsv"

# canonical markers per lineage; only those present in GENCODE symbol space
PANEL = {
    "CAF": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "PDGFRB",
            "FAP", "THY1", "ACTA2", "POSTN", "FBLN1"],
    "epithelial": ["EPCAM", "KRT5", "KRT14", "KRT13", "KRT17", "CDH1", "SFN"],
    "immune": ["PTPRC", "CD3D", "CD3E", "CD2", "MS4A1", "CD79A", "LYZ",
               "CD68", "AIF1", "CD14"],
    "endothelial": ["PECAM1", "VWF", "CDH5", "PLVAP", "ENG", "KDR"],
}


def main() -> None:
    A = ad.read_h5ad(H5)
    print(f"loaded {A.shape}", flush=True)

    A.layers["counts"] = A.X.copy()
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)

    # The cached h5ad was written before Leiden ran, so join the cluster labels
    # from the per-cell table instead of expecting them in obs.
    lc = pd.read_csv(RES / "sc_percell.tsv", sep="\t")
    # The cached h5ad keeps only a global integer index, so align by row order
    # and verify with the donor column rather than trusting position blindly.
    assert len(lc) == A.n_obs, f"row count mismatch {len(lc)} vs {A.n_obs}"
    assert (lc["donor"].to_numpy() == A.obs["donor"].to_numpy()).all(), \
        "per-cell table is not in the same order as the h5ad"
    clusters = pd.Series(lc["leiden"].astype(int).astype(str).to_numpy(),
                         index=A.obs_names)

    rows = []
    for lineage, genes in PANEL.items():
        present = [g for g in genes if g in A.var_names]
        for cl in sorted(clusters.unique(), key=lambda s: int(s)):
            mask = (clusters == cl).to_numpy()
            if mask.sum() == 0:
                continue
            for g in present:
                v = np.asarray(A[mask, g].X.todense()).ravel()
                rows.append(dict(cluster=cl, lineage=lineage, gene=g,
                                 n_cells=int(mask.sum()),
                                 mean_expr=float(v.mean()),
                                 frac_expr=float((v > 0).mean())))

    df = pd.DataFrame(rows)
    df.to_csv(OUT, sep="\t", index=False)
    print(f"written {OUT}  rows={len(df)}  genes={df['gene'].nunique()}")

    # which cluster is which lineage: majority vote from the per-cell calls
    tab = (lc.groupby(["leiden", "lineage"]).size().unstack(fill_value=0))
    tab.to_csv(RES / "sc_cluster_lineage.tsv", sep="\t")
    print("\ncluster x lineage (majority):")
    print(tab.idxmax(axis=1).to_string())


if __name__ == "__main__":
    main()
