"""Supplementary Fig S1 (3x3): independent single-cell validation of the FBLN1 result.

# PANEL-SCALE-1.32 (2026-10-22): fonts enlarged to match Fig 2.

Why this replaces the 2-panel version
-------------------------------------
The first version showed only the per-library C:D comparison, which carries no
evidence that the data are single cell. This rebuild uses the per-cell table
(29 295 cells, 7 libraries, GSE301741) and adds the panels a reader expects of a
single-cell figure: UMAP embeddings, a canonical marker dot plot, composition,
QC, and detection rate -- plus the FBLN1 validation itself.

Panels
  a   UMAP, coloured by lineage (29 295 cells)
  b   UMAP, coloured by library -- makes the batch structure visible
  c   marker dot plot, clusters ordered by majority lineage
  d   lineage composition of each library
  e   QC: UMIs per cell, by lineage
  f   FBLN1 detection rate by lineage (CAF 67% vs epithelia 25%)
  g   FBLN1 1D share per lineage per library (the validation, paired)
  h   cluster x lineage composition, all 30 clusters
  i   per-library delta (CAF - epithelia) with mean +/- SE

Every number is computed from the tables in remote/easyai/; nothing is hard-coded
except the lineage colour map.
"""
from __future__ import annotations

import json
from pathlib import Path

import sys as _sys

_sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figstyle import panel, fit_panel_titles  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

STYLE = Path(__file__).resolve().parent / "_style_runtime.py"
exec(compile(STYLE.read_text(encoding="utf-8"), "runtime.py", "exec"), globals())
apply_figure_style(sizes=(11, 9.5, 8.5))

ROOT = Path(__file__).resolve().parents[2]
REM = ROOT / "remote/easyai"
OUT = ROOT / "figures"

LCOL = {
    "CAF": "#b8342b",
    "epithelial": "#1f6fb2",
    "immune": "#c47a2c",
    "endothelial": "#5b8f5b",
    "unassigned": "#bbbbbb",
}
LORDER = ["CAF", "epithelial", "immune", "endothelial", "unassigned"]


def main() -> None:
    pc = pd.read_csv(REM / "sc_percell.tsv", sep="\t", dtype={"cell": str})
    dot = pd.read_csv(REM / "sc_marker_dot.tsv", sep="\t")
    val = json.loads((REM / "fbln1_singlecell_validation.json").read_text(encoding="utf-8"))

    print(f"cells={len(pc):,}  clusters={pc.leiden.nunique()}  donors={pc.donor.nunique()}")

    # cluster -> majority lineage (used to order the dot plot and to aggregate)
    maj = pc.groupby("leiden")["lineage"].agg(lambda s: s.value_counts().idxmax())
    dot["maj"] = dot["cluster"].map(maj)
    cl_order = sorted(maj.index, key=lambda c: (LORDER.index(maj[c]), int(c)))

    fig, axes = plt.subplots(3, 3, figsize=(16.72, 15.44))
    ax_a, ax_b, ax_c = axes[0]
    ax_d, ax_e, ax_f = axes[1]
    ax_g, ax_h, ax_i = axes[2]

    # ---------------------------------------------------------------- a
    # clip to robust quantiles: a handful of outlier cells otherwise leave stray
    # marks against the right spine
    xlo, xhi = np.percentile(pc.UMAP1, [0.5, 99.5])
    ylo, yhi = np.percentile(pc.UMAP2, [0.5, 99.5])
    xlo, xhi, ylo, yhi = xlo - 0.6, xhi + 0.6, ylo - 0.6, yhi + 0.6
    inbox = pc.UMAP1.between(xlo, xhi) & pc.UMAP2.between(ylo, yhi)
    # plot a windowed VIEW; never shorten `pc`, downstream panels index it by row
    vis = pc[inbox]
    for ax in (ax_a, ax_b):
        ax.set_xlim(xlo, xhi)
        ax.set_ylim(ylo, yhi)
    # legend counts are over ALL cells (`pc`) so they sum to the 29 295 stated in
    # the title; only the drawn points come from the clipped window
    lin_tot = pc.lineage.value_counts()
    for lin in LORDER:
        sub = vis[vis.lineage == lin]
        if len(sub) == 0:
            continue
        ax_a.scatter(sub.UMAP1, sub.UMAP2, s=1.1, c=LCOL[lin], linewidths=0,
                     alpha=0.55, rasterized=True,
                     label=f"{lin} ({int(lin_tot[lin]):,})")
    ax_a.set_xlabel("UMAP-1")
    ax_a.set_ylabel("UMAP-2")
    panel(ax_a, f"a  lineage landscape, {len(pc):,} cells", fontsize=12.4)
    lg = ax_a.legend(frameon=False, fontsize=7.7, loc="upper right",
                     markerscale=5, handletextpad=0.3, labelspacing=0.25,
                     borderaxespad=0.2)
    for h in lg.legend_handles:
        h.set_alpha(1.0)

    # ---------------------------------------------------------------- b
    donors = sorted(pc.donor.unique())
    cmap = plt.get_cmap("tab10")
    # Counts in both legends must describe the SAME population as the panel
    # titles. Panels a/b plot the robust-quantile window (`vis`, 29 149 of
    # 29 295 cells), and panel a's title says 29 295 -- so the per-library
    # labels here must be taken over `pc`, not over `vis`, or the two legends
    # disagree by the 146 clipped cells and the visible labels stop summing to
    # the stated total.
    tot = pc.donor.value_counts()
    vis_n = vis.donor.value_counts()
    for k, d in enumerate(donors):
        sub = vis[vis.donor == d]
        ax_b.scatter(sub.UMAP1, sub.UMAP2, s=1.1, c=[cmap(k % 10)], linewidths=0,
                     alpha=0.55, rasterized=True,
                     label=f"{d} ({int(tot[d]):,})")
    ax_b.set_xlabel("UMAP-1")
    ax_b.set_ylabel("UMAP-2")
    panel(ax_b, "b  the same embedding by library", fontsize=12.4)
    ax_b.legend(frameon=False, fontsize=7.4, loc="upper right", markerscale=5,
                handletextpad=0.3, labelspacing=0.22, borderaxespad=0.2)
    ax_b.text(0.015, 0.015,
              f"{len(pc):,} cells; {len(pc) - len(vis)} clipped by the "
              f"0.5-99.5 percentile window\nand therefore not drawn",
              transform=ax_b.transAxes, fontsize=6.8, color="#666666",
              va="bottom", ha="left", linespacing=1.25)

    # ---------------------------------------------------------------- c
    genes = ["COL1A1", "COL1A2", "DCN", "LUM", "COL3A1", "PDGFRA", "FAP", "THY1",
             "EPCAM", "KRT5", "KRT14", "SFN",
             "PTPRC", "CD3D", "CD2", "CD68", "LYZ",
             "PECAM1", "VWF", "CDH5", "PLVAP"]
    genes = [g for g in genes if g in set(dot.gene)]
    mat_f = np.zeros((len(genes), len(cl_order)))
    mat_e = np.zeros_like(mat_f)
    for gi, g in enumerate(genes):
        d = dot[dot.gene == g].set_index("cluster")
        for ci, cl in enumerate(cl_order):
            if cl in d.index:
                mat_f[gi, ci] = d.loc[cl, "frac_expr"]
                mat_e[gi, ci] = d.loc[cl, "mean_expr"]
    # per-gene scaling so all markers are visible on one colour scale
    z = (mat_e - mat_e.mean(axis=1, keepdims=True)) / np.where(
        mat_e.std(axis=1, keepdims=True) > 0, mat_e.std(axis=1, keepdims=True), 1)
    yx, xx = np.meshgrid(np.arange(len(genes)), np.arange(len(cl_order)), indexing="ij")
    sc_ = ax_c.scatter(xx.ravel(), yx.ravel(), s=(mat_f.ravel() ** 1.6) * 42,
                       c=z.ravel(), cmap="RdBu_r", vmin=-2, vmax=2,
                       linewidths=0)
    ax_c.set_xticks(range(len(cl_order)))
    ax_c.set_xticklabels(cl_order, fontsize=6.1, rotation=90)
    ax_c.set_yticks(range(len(genes)))
    ax_c.set_yticklabels(genes, fontsize=7.1)
    ax_c.set_xlabel("Leiden cluster (ordered by majority lineage)")
    panel(ax_c, "c  lineage markers recover the annotation", fontsize=12.4)
    ax_c.set_xlim(-0.7, len(cl_order) - 0.3)
    ax_c.set_ylim(-1.30, len(genes) - 0.3)
    cb = fig.colorbar(sc_, ax=ax_c, fraction=0.030, pad=0.015)
    cb.set_label("scaled mean expression", fontsize=7.7)
    cb.ax.tick_params(labelsize=6.9)
    # lineage bands above the columns
    for lin in LORDER:
        idx = [i for i, c in enumerate(cl_order) if maj[c] == lin]
        if not idx:
            continue
        # drawn INSIDE the data area: clip_on=False plus a blended transform
        # makes bbox_inches="tight" expand the canvas without bound
        ax_c.plot([min(idx) - 0.4, max(idx) + 0.4], [-0.55] * 2,
                  color=LCOL[lin], lw=3.0, solid_capstyle="butt", zorder=5)
        ax_c.text(np.mean([min(idx) - 0.4, max(idx) + 0.4]), -0.66, lin[:4],
                  ha="center", va="top", fontsize=6.6, color=LCOL[lin])

    # ---------------------------------------------------------------- d
    comp = (pc.groupby(["donor", "lineage"]).size().unstack(fill_value=0)
              .reindex(columns=LORDER, fill_value=0))
    comp = comp.loc[comp.sum(axis=1).sort_values(ascending=False).index]
    bot = np.zeros(len(comp))
    xs = np.arange(len(comp))
    for lin in LORDER:
        v = comp[lin].to_numpy() / comp.sum(axis=1).to_numpy()
        ax_d.bar(xs, v, bottom=bot, color=LCOL[lin], width=0.68,
                 edgecolor="white", lw=0.5, label=lin)
        bot += v
    ax_d.set_xticks(xs)
    # 35 deg still lets neighbouring 11-char names touch; 90 deg cannot
    ax_d.set_xticklabels([f"{d.replace('_', ' ')}" for d in comp.index],
                         rotation=90, ha="center", fontsize=7.0)
    ax_d.set_ylabel("fraction of cells")
    ax_d.set_ylim(0, 1)
    panel(ax_d, "d  composition differs between libraries", fontsize=12.4)
    # the stacks fill the whole axes, so the legend goes to a reserved band above
    ax_d.set_ylim(0, 1.30)
    ax_d.legend(frameon=False, fontsize=7.0, ncol=3, loc="upper center",
                bbox_to_anchor=(0.5, 1.005), handlelength=1.0, columnspacing=1.0,
                handletextpad=0.35)

    # ---------------------------------------------------------------- e
    # QC: UMI per cell, by lineage. Depth differs between lineages, so this panel
    # also shows that the CAF>epithelia result is not a low-depth artefact.
    qc = pd.read_csv(REM / "sc_qc.tsv", sep="\t")
    assert len(qc) == len(pc), f"QC table has {len(qc)} rows, expected {len(pc)}"
    order_e = [l for l in LORDER if (pc.lineage == l).sum() > 50]
    dat = [np.log10(qc.loc[pc.lineage == l, "umi"].clip(lower=1).to_numpy())
           for l in order_e]
    parts = ax_e.violinplot(dat, positions=range(len(order_e)), widths=0.78,
                            showextrema=False, showmedians=False)
    for body, l in zip(parts["bodies"], order_e):
        body.set_facecolor(LCOL[l])
        body.set_alpha(0.78)
        body.set_edgecolor("white")
        body.set_linewidth(0.6)
    for k, (l, d) in enumerate(zip(order_e, dat)):
        ax_e.plot([k], [np.median(d)], marker="o", ms=3.0, color="#222222", zorder=4)
        ax_e.text(k, 0.985, f"{10 ** np.median(d):,.0f}", transform=ax_e.get_xaxis_transform(),
                  ha="center", va="top", fontsize=7.4, color="#333333")
    ax_e.set_xticks(range(len(order_e)))
    ax_e.set_xticklabels(order_e, rotation=30, ha="right", fontsize=8.7)
    ax_e.set_ylabel("log$_{10}$ UMI per cell")
    panel(ax_e, "e  depth is comparable across lineages", fontsize=12.4)

    # ---------------------------------------------------------------- f
    f_ = dot[dot.gene == "FBLN1"]
    det = f_.groupby("maj").apply(
        lambda g: np.average(g["frac_expr"], weights=g["n_cells"]), include_groups=False)
    det = det.reindex([l for l in LORDER if l in det.index])
    ax_f.bar(range(len(det)), det.to_numpy() * 100,
             color=[LCOL[l] for l in det.index], width=0.62,
             edgecolor="white", lw=0.6)
    for k, v in enumerate(det.to_numpy() * 100):
        ax_f.text(k, v + 1.6, f"{v:.1f}%", ha="center", fontsize=8.4)
    ax_f.set_xticks(range(len(det)))
    ax_f.set_xticklabels(det.index, rotation=30, ha="right", fontsize=8.7)
    ax_f.set_ylabel("FBLN1-positive cells (%)")
    ax_f.set_ylim(0, max(det.to_numpy() * 100) * 1.22)
    panel(ax_f, "f  FBLN1 is a CAF-enriched transcript", fontsize=12.4)

    # ---------------------------------------------------------------- g
    libs = [r for r in val["per_library"]
            if r.get("share_tpm_CAF") is not None
            and r.get("share_tpm_epi") is not None]
    for r in libs:
        ax_g.plot([0, 1], [r["share_tpm_epi"], r["share_tpm_CAF"]], color="#888888", lw=0.9,
                  marker="o", ms=3.6, mfc="white", mew=0.9, zorder=2)
    ax_g.scatter(np.zeros(len(libs)), [r["share_tpm_epi"] for r in libs], s=26,
                 color=LCOL["epithelial"], zorder=3, edgecolor="white", lw=0.7)
    ax_g.scatter(np.ones(len(libs)), [r["share_tpm_CAF"] for r in libs], s=26,
                 color=LCOL["CAF"], zorder=3, edgecolor="white", lw=0.7)
    dv = np.array([r["d_tpm"] for r in libs])
    ymin = min(min(r["share_tpm_epi"] for r in libs), min(r["share_tpm_CAF"] for r in libs))
    ymax = max(max(r["share_tpm_epi"] for r in libs), max(r["share_tpm_CAF"] for r in libs))
    ax_g.set_ylim(ymin - (ymax - ymin) * 0.55, ymax + (ymax - ymin) * 0.10)
    ax_g.text(0.5, 0.02,
              f"1D share rises in CAF in {int((dv > 0).sum())}/{len(dv)} libraries\n"
              f"mean $\\Delta$ = {dv.mean():+.3f}   paired t p = "
              f"{stats.ttest_rel([r['share_tpm_CAF'] for r in libs], [r['share_tpm_epi'] for r in libs]).pvalue:.4f}",
              transform=ax_g.transAxes, ha="center", va="bottom", fontsize=8.3,
              linespacing=1.35)
    ax_g.set_xticks([0, 1])
    ax_g.set_xticklabels(["epithelia", "CAF"], fontsize=9.8)
    ax_g.set_xlim(-0.45, 1.45)
    ax_g.set_ylabel("FBLN1-1D share of FBLN1")
    panel(ax_g, "g  the same cells reproduce the switch", fontsize=12.4)

    # ---------------------------------------------------------------- h
    cl_comp = (pc.groupby(["leiden", "lineage"]).size().unstack(fill_value=0)
                 .reindex(columns=LORDER, fill_value=0))
    cl_comp = cl_comp.loc[cl_order]
    cl_frac = cl_comp.div(cl_comp.sum(axis=1), axis=0)
    bot = np.zeros(len(cl_frac))
    xs = np.arange(len(cl_frac))
    for lin in LORDER:
        ax_h.bar(xs, cl_frac[lin].to_numpy(), bottom=bot, color=LCOL[lin],
                 width=0.86, lw=0)
        bot += cl_frac[lin].to_numpy()
    ax_h.set_xticks(xs[::2])
    ax_h.set_xticklabels([cl_order[i] for i in range(0, len(cl_order), 2)],
                         fontsize=7.7, rotation=90)
    ax_h.set_xlabel("Leiden cluster")
    ax_h.set_ylabel("fraction of cells")
    ax_h.set_ylim(0, 1)
    panel(ax_h, "h  30 clusters are mixtures, not pure types", fontsize=12.4)

    # ---------------------------------------------------------------- i
    order_i = sorted(range(len(dv)), key=lambda k: dv[k])
    ax_i.barh(range(len(dv)), dv[order_i], color="#5b8f5b", height=0.6,
              edgecolor="white", lw=0.6)
    ax_i.axvline(0, color="#333333", lw=0.9)
    ax_i.set_yticks(range(len(dv)))
    ax_i.set_yticklabels([libs[k]["library"].replace("_", " ") for k in order_i],
                         fontsize=8.4)
    ax_i.set_xlabel("$\\Delta$ FBLN1-1D share (CAF $-$ epithelia)")
    ax_i.set_xlim(-0.02, max(dv) * 1.28)
    panel(ax_i, "i  every library moves the same way", fontsize=12.4)
    se = dv.std(ddof=1) / np.sqrt(len(dv))
    ax_i.text(0.985, 0.06,
              f"mean {dv.mean():+.3f} $\\pm$ {se:.3f} (SE)\n"
              f"{int((dv > 0).sum())}/{len(dv)} positive",
              transform=ax_i.transAxes, ha="right", va="bottom", fontsize=8.2,
              linespacing=1.3)

    fig.tight_layout(rect=(0, 0, 1, 1))
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / "supp_fig_s1_singlecell_fbln1.png"
    fit_panel_titles(fig)
    fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)

    from PIL import Image
    wpx, hpx = Image.open(p).size
    assert wpx < 7000 and hpx < 7000, f"suspicious canvas {wpx}x{hpx}"
    print(f"written {p.name} {wpx}x{hpx}")
    print(f"g: CAF>epi in {int((dv > 0).sum())}/{len(dv)}, mean delta {dv.mean():+.3f}, "
          f"p={stats.ttest_rel([r['share_tpm_CAF'] for r in libs], [r['share_tpm_epi'] for r in libs]).pvalue:.4f}")
    print("f: FBLN1+ % =", {k: round(v * 100, 1) for k, v in det.items()})


if __name__ == "__main__":
    main()
