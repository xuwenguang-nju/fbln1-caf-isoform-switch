"""Per-cell extraction + UMAP for the GSE301741 single-cell validation figure.

Why a second pass
-----------------
`bin/build_celltype_tpm.py` computed the lineage of every cell but only wrote the
*aggregated* per-cell-type TPM; the per-cell labels were never saved. The typical
single-cell panels (UMAP, marker dot plot, QC, FBLN1 detection rate) all need
per-cell information, so this script rebuilds it with the identical lineage rule
and additionally runs the standard scanpy workflow.

Lineage rule is copied verbatim from build_celltype_tpm.py so the per-cell labels
here are guaranteed consistent with the pseudobulk used in the paper:
    best = argmax(marker score); best >= 3 and best > 2 * second  -> assigned

Outputs (under /root/hnscc_caf/results/):
  sc_percell.tsv   donor, cell, lineage, umi, genes, n_1C, n_1D, cluster, UMAP1..2
  sc_markers.tsv   donor, cell, log-normalised expression of the marker panel
  sc_umap_summary.json
"""
from __future__ import annotations
import collections
import gzip
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

DONORS = ["OCSCC1_Post", "OCSCC1_Pre", "OCSCC3_Post", "OCSCC13_Post",
          "OCSCC36_Post", "OCSCC36_Pre", "LSCC20_Pre"]
BASE = Path("/root/as_psi/donors")
REF = Path("/root/as_psi/ref")
OUT = Path("/root/hnscc_caf/results")

MARK = dict(
    CAF=["COL1A1", "COL1A2", "DCN", "LUM", "COL3A1", "COL6A1", "COL6A2",
         "PDGFRA", "PDGFRB", "FAP", "THY1"],
    epithelial=["EPCAM", "KRT5", "KRT13", "KRT14", "SFN", "KRT17"],
    immune=["PTPRC", "CD3D", "CD3E", "CD2", "MS4A1", "CD68", "LYZ"],
    endothelial=["PECAM1", "VWF", "CDH5", "PLVAP"],
)
EXTRA = ["FBLN1", "ACTA2", "TAGLN", "MKI67", "COL1A1"]
TX_C, TX_D = "ENST00000262722", "ENST00000327858"


def load_maps():
    t2g, g2s = {}, {}
    for line in open(REF / "tx2gene.tsv"):
        a, b = line.rstrip("\n").split("\t")[:2]
        t2g[a.split(".")[0]] = b.split(".")[0]
    for line in open(REF / "pc_transcripts.fa"):
        if line[0] != ">":
            continue
        p = line[1:].split("|")
        if len(p) >= 6:
            g2s[p[1].split(".")[0]] = p[5]
    return t2g, g2s


def donor_matrix(donor, t2g, g2s, gene_index, n_genes):
    """Return (gene-level CSR for this donor, obs dict of lists)."""
    adir = BASE / donor / "alevin/alevin"
    txs = [l.strip() for l in open(adir / "quants_mat_cols.txt")]
    cells = [l.strip() for l in open(adir / "quants_mat_rows.txt")]
    ti = {t: i for i, t in enumerate(txs)}
    # transcript -> gene row index, or -1
    tx_gene = np.array([gene_index.get(g2s.get(t2g.get(t.split(".")[0], ""), ""), -1)
                        for t in txs], dtype=np.int32)
    tx_is_c = np.array([t.startswith(TX_C) for t in txs])
    tx_is_d = np.array([t.startswith(TX_D) for t in txs])

    # marker gene indices we need per-cell scores for
    sym_of = {}
    for gi, sym in enumerate(gene_index):
        sym_of.setdefault(sym, gi)
    mark_idx = {k: [sym_of[g] for g in gs if g in sym_of] for k, gs in MARK.items()}
    score_idx = np.zeros(n_genes, dtype=np.int8)   # 0 none, else bitmask-ish id
    flat = {}
    for k, idxs in mark_idx.items():
        for gi in idxs:
            flat[gi] = k

    n_cells = len(cells)
    rows, cols, vals = [], [], []
    umi = np.zeros(n_cells, dtype=np.int64)
    genes = np.zeros(n_cells, dtype=np.int32)
    n_c = np.zeros(n_cells, dtype=np.int32)
    n_d = np.zeros(n_cells, dtype=np.int32)
    sc = {k: np.zeros(n_cells, dtype=np.float64) for k in MARK}
    per_cell_nz = np.zeros(n_cells, dtype=np.int32)

    with gzip.open(adir / "quants_mat.mtx.gz", "rt") as fh:
        fh.readline(); fh.readline()
        for line in fh:
            i, j, v = line.split()
            ci, tj = int(i) - 1, int(j) - 1
            val = float(v)
            umi[ci] += val
            per_cell_nz[ci] += 1
            if tx_is_c[tj]:
                n_c[ci] += val
            elif tx_is_d[tj]:
                n_d[ci] += val
            gi = tx_gene[tj]
            if gi >= 0:
                rows.append(ci); cols.append(gi); vals.append(val)
                k = flat.get(gi)
                if k:
                    sc[k][ci] += val

    X = sp.coo_matrix((vals, (rows, cols)), shape=(n_cells, n_genes),
                      dtype=np.float32).tocsr()
    genes = np.asarray(X.getnnz(axis=1)).ravel().astype(np.int32)
    lineage = []
    for c in range(n_cells):
        s = {k: sc[k][c] for k in MARK}
        best = max(s, key=s.get)
        second = sorted(s.values())[-2]
        lineage.append(best if s[best] >= 3 and s[best] > 2 * max(second, 1) else "unassigned")
    return X, dict(cells=cells, umi=umi, genes=genes, n_c=n_c, n_d=n_d,
                   lineage=np.array(lineage))


def main():
    t2g, g2s = load_maps()
    # one metadata pass to fix the gene axis
    gsym = {}
    for line in open(REF / "pc_transcripts.fa"):
        if line[0] == ">":
            p = line[1:].split("|")
            if len(p) >= 6:
                gsym[p[1].split(".")[0]] = p[5]
    # gene axis = union of genes reachable from the quantified transcripts
    gene_index, symbols = {}, []
    for line in open(REF / "tx2gene.tsv"):
        a, b = line.rstrip("\n").split("\t")[:2]
        gid = b.split(".")[0]
        sym = gsym.get(gid)
        if not sym or sym in gene_index:
            continue
        gene_index[sym] = len(symbols)
        symbols.append(sym)
    n_genes = len(symbols)
    print(f"gene axis: {n_genes:,} symbols", flush=True)

    blocks, obs_all, donors_col = [], [], []
    for d in DONORS:
        X, obs = donor_matrix(d, t2g, g2s, gene_index, n_genes)
        print(f"{d}: {X.shape[0]:,} cells, median UMI {np.median(obs['umi']):.0f}, "
              f"{dict(collections.Counter(obs['lineage']))}", flush=True)
        blocks.append(X)
        obs_all.append(obs)
        donors_col.extend([d] * X.shape[0])

    import anndata as ad
    import scanpy as sc
    X = sp.vstack(blocks, format="csr")
    var = [s for s in symbols]
    A = ad.AnnData(X=X)
    # var_names must BE the symbols: passing var={"symbol": ...} alone leaves
    # var_names as 0..n-1 strings, which silently breaks every later gene lookup.
    A.var_names = pd.Index(var, name="symbol")
    A.obs["donor"] = donors_col
    A.obs["lineage"] = np.concatenate([o["lineage"] for o in obs_all])
    A.obs["umi"] = np.concatenate([o["umi"] for o in obs_all])
    A.obs["genes"] = np.concatenate([o["genes"] for o in obs_all])
    A.obs["n_1C"] = np.concatenate([o["n_c"] for o in obs_all])
    A.obs["n_1D"] = np.concatenate([o["n_d"] for o in obs_all])
    assert "COL1A1" in A.var_names and "EPCAM" in A.var_names, \
        "gene symbols missing from var_names"
    print(f"merged: {A.shape}", flush=True)

    raw_h5 = OUT / "sc_raw.h5ad"
    if raw_h5.exists():
        print(f"reusing cached {raw_h5.name}", flush=True)
        A = ad.read_h5ad(raw_h5)
    else:
        sc.pp.filter_cells(A, min_genes=200)
        sc.pp.filter_genes(A, min_cells=10)
        print(f"after QC filter: {A.shape}", flush=True)
        A.write_h5ad(raw_h5)

    A.layers["counts"] = A.X.copy()
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    sc.pp.highly_variable_genes(A, n_top_genes=2000, batch_key="donor")
    sc.pp.scale(A, max_value=10)
    sc.tl.pca(A, n_comps=30, svd_solver="arpack")
    sc.pp.neighbors(A, n_neighbors=15, n_pcs=30)
    sc.tl.umap(A, min_dist=0.3)
    sc.tl.leiden(A, resolution=0.8, flavor="igraph", n_iterations=2, directed=False)
    print(f"clusters: {A.obs['leiden'].nunique()}", flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    markers = [g for gs in MARK.values() for g in gs] + EXTRA
    markers = list(dict.fromkeys(markers))
    L = A[:, [g for g in markers if g in A.var_names]].X
    if sp.issparse(L):
        L = L.toarray()
    md = A.obs.copy()
    md["cell"] = A.obs_names
    md[["donor", "cell", "lineage", "umi", "genes", "n_1C", "n_1D", "leiden",
        "UMAP1", "UMAP2"]] = np.column_stack([
        md["donor"].values, md["cell"].values, md["lineage"].values,
        md["umi"].values, md["genes"].values, md["n_1C"].values, md["n_1D"].values,
        md["leiden"].values, A.obsm["X_umap"][:, 0], A.obsm["X_umap"][:, 1]])
    md[["donor", "cell", "lineage", "leiden", "UMAP1", "UMAP2"]].to_csv(
        OUT / "sc_percell.tsv", sep="\t", index=False)
    mdf = A.obs[["donor", "lineage"]].copy()
    mdf["cell"] = A.obs_names
    for k, g in enumerate([x for x in markers if x in A.var_names]):
        mdf[g] = L[:, k]
    mdf.to_csv(OUT / "sc_markers.tsv", sep="\t", index=False)

    summ = dict(n_cells=int(A.n_obs), n_genes=int(A.n_vars),
                clusters=sorted(A.obs["leiden"].unique().tolist()),
                per_donor=A.obs["donor"].value_counts().to_dict(),
                per_lineage=A.obs["lineage"].value_counts().to_dict(),
                markers=[x for x in markers if x in A.var_names])
    (OUT / "sc_umap_summary.json").write_text(json.dumps(summ, indent=1))
    print("written", OUT / "sc_percell.tsv")
    print(json.dumps(summ, indent=1)[:600])


if __name__ == "__main__":
    main()
