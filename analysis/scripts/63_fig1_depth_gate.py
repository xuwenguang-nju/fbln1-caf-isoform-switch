"""Fig 2 (8 panels): sequencing depth decides which isoform switches are reproducible.

Panels a-d were the original figure. Panels e-h add the depth noise model and the
per-candidate consequence, all from tables that already existed but were unused:

  e  replicate SD as a function of per-gene depth bin      caf_isoform_noise_by_depth.csv
  f  every panel candidate, full vs 10% depth              caf_panel_depth_consistency.csv
  g  how many candidates survive at each dIF/FDR threshold caf_switch_threshold_sweep.csv
  h  reads per gene required to reach a target replicate SD caf_depth_requirement.csv

Why these four: a-d show THAT the gate removes genes; e-h show WHY (noise falls
with depth), WHICH candidates it removes (all but FBLN1), HOW the count depends on
the threshold, and WHAT depth would have been needed to avoid the problem.

Outputs figures/fig_depth_gate.png
"""
from __future__ import annotations
from pathlib import Path

import sys as _sys

_sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figstyle import panel, fit_panel_titles  # noqa: E402

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
T = ROOT / "results/tables"
FIG = ROOT / "figures"
STYLE = Path(__file__).resolve().parent / "_style_runtime.py"
exec(compile(STYLE.read_text(encoding="utf-8"), "runtime.py", "exec"), globals())

CALL_ORDER = [
    ("switch", "switch\n(|dIF|>=0.15, reproducible)", "#b8342b"),
    ("small_shift", "small shift\n(0.05-0.15)", "#c47a2c"),
    ("negligible_but_reproducible", "negligible\nbut reproducible", "#7b4b8a"),
    ("not_significant", "not significant", "#9aa3ad"),
    ("excluded_depth_unstable", "excluded:\ndepth-unstable", "#1f6fb2"),
]


def main() -> None:
    apply_figure_style(sizes=(11, 9.5, 8.5))
    # Fixed-referent inputs: dIF is the same reference transcript's share in both
    # conditions (181_gate_under_fixed_referent.py). The old tables subtracted one
    # transcript's fraction in CAF from a *different* transcript's in NF whenever a
    # gene's per-library dominant isoform changed between conditions.
    full = pd.read_csv(T / "switch_test_full_fixed_referent.tsv", sep="\t")
    sub = pd.read_csv(T / "switch_test_sub10_fixed_referent.tsv", sep="\t")
    gated = pd.read_csv(T / "candidates_gated_fixed_referent.tsv", sep="\t")

    both = full.merge(sub[["gene", "mean_dIF"]], on="gene", suffixes=("_full", "_sub"))
    n_full, n_both = len(full), len(both)
    not_meas = n_full - n_both
    assert len(gated) == n_both, f"gated {len(gated)} != both-depth {n_both}"
    counts = gated["call"].astype(str).value_counts()
    print(f"full={n_full}  both={n_both}  only-full={not_meas}  calls={counts.to_dict()}")

    fig, axes = plt.subplots(3, 3, figsize=(14.8, 13.6))
    ax_b, ax_a, ax_d = axes[0]
    ax_c, ax_e, ax_h = axes[1]
    ax_i, ax_g, ax_f = axes[2]
    # ================================================================ a
    x = both["mean_dIF_full"].abs().to_numpy()
    y = both["mean_dIF_sub"].abs().to_numpy()
    ax_a.scatter(x, y, s=4, alpha=0.16, linewidths=0, color="#8a8f98", rasterized=True)
    lim = 0.40
    ax_a.plot([0, lim], [0, lim], ls="--", lw=1.0, color="#333333")
    for gid, lab, dx, dy, col in (("ENSG00000077942", "FBLN1", 0.045, -0.055, "#1f6fb2"),
                                  ("ENSG00000115884", "SDC1", 0.028, 0.050, "#c47a2c")):
        r = both[both["gene"].astype(str).eq(gid)]
        if not len(r):
            raise KeyError(f"{lab} ({gid}) absent from the both-depth set")
        px, py = abs(float(r.iloc[0]["mean_dIF_full"])), abs(float(r.iloc[0]["mean_dIF_sub"]))
        ax_a.scatter([px], [py], s=70, color=col, edgecolor="white", lw=1.0, zorder=5)
        ax_a.annotate(lab, (px, py), xytext=(px + dx, py + dy), color=col, fontsize=10.0,
                      arrowprops=dict(arrowstyle="-", color=col, lw=0.9))
    ax_a.set_xlim(-0.005, lim)
    ax_a.set_ylim(-0.005, lim)
    ax_a.set_xlabel("|dIF| at full depth")
    ax_a.set_ylabel("|dIF| at 10% depth")
    panel(ax_a, f"B  {n_both:,} genes, both depths", fontsize=11.9, ha="center")

    # ================================================================ b
    items = [(key, lab, col, int(counts.get(key, 0))) for key, lab, col in CALL_ORDER]
    items.append(("__notmeas", "not measurable\nat 10% depth", "#e8e8e8", int(not_meas)))
    items.sort(key=lambda t: t[3], reverse=True)
    row_pos = np.arange(len(items))[::-1]
    # do NOT name the loop variable `y`; it shadows the |dIF| array used by panel c
    for yy, (key, lab, col, v) in zip(row_pos, items):
        edge = "#bbbbbb" if key == "__notmeas" else "white"
        ax_b.barh(yy, v, color=col, height=0.68, edgecolor=edge, lw=0.7)
        ax_b.text(v + 55, yy, f"{v:,}  ({v / n_full:.1%})", va="center", fontsize=8.7,
                  color="#333333")
        ax_b.text(-55, yy, lab, va="center", ha="right", fontsize=8.7, color="#333333")
    ax_b.set_yticks([])
    ax_b.set_xlim(0, n_full * 0.98)
    ax_b.set_xlabel(f"genes (of {n_full:,})")
    ax_b.set_ylim(-0.7, len(items) - 0.3)
    panel(ax_b, "A  where the gate puts every gene", fontsize=11.9)
    for s in ("left", "right", "top"):
        ax_b.spines[s].set_visible(False)

    # ================================================================ c
    floor = 0.01
    xf, ys_ = x[x >= floor], y[y >= floor]
    bins = np.linspace(floor, 0.30, 49)
    ax_c.hist(xf, bins=bins, density=True, histtype="step", lw=1.7, color="#1f6fb2",
              label=f"full  (n={len(xf):,})")
    ax_c.hist(ys_, bins=bins, density=True, histtype="step", lw=1.7, color="#c47a2c",
              label=f"10%  (n={len(ys_):,})")
    ax_c.set_xlim(floor, 0.30)
    ax_c.set_xlabel("|dIF|  (genes with |dIF| $\\geq$ 0.01)")
    ax_c.set_ylabel("density")
    ax_c.legend(frameon=False, fontsize=9.0, loc="upper right")
    panel(ax_c, "D  low depth broadens, not shrinks", fontsize=11.9)
    med_f, med_s = np.median(xf), np.median(ys_)
    ax_c.axvline(med_f, color="#1f6fb2", ls=":", lw=1.1)
    ax_c.axvline(med_s, color="#c47a2c", ls=":", lw=1.1)
    ax_c.text(0.985, 0.56,
              f"median |dIF|\nfull {med_f:.3f}\n10%  {med_s:.3f}",
              transform=ax_c.transAxes, ha="right", va="top", fontsize=8.7,
              color="#333333", linespacing=1.35)

    # ================================================================ d
    g = gated.copy()
    g["abs_full"] = g["mean_dIF_full"].abs()
    edges = np.array([0, 0.05, 0.10, 0.15, 0.20, 0.40])
    g["bin"] = pd.cut(g["abs_full"], edges, include_lowest=True)
    agg = g.groupby("bin", observed=True).apply(
        lambda d: pd.Series({
            "n": len(d),
            "unstable": (d["call"].astype(str) == "excluded_depth_unstable").sum(),
            "med_change": d["abs_change"].median()}),
        include_groups=False)
    # bins can be empty (no gene with |dIF| in that window); drop them explicitly
    agg = agg[agg["n"] > 0]
    agg["rate"] = agg["unstable"] / agg["n"]
    centers = np.array([iv.mid for iv in agg.index])

    ax_d.plot(centers, agg["rate"], "o-", color="#1f6fb2", lw=1.7, ms=6,
              label="fraction excluded")
    for cx, r_, n_ in zip(centers, agg["rate"], agg["n"]):
        ax_d.text(cx, r_ + 0.05, f"{r_:.0%}", ha="center", fontsize=8.7, color="#1f6fb2")
    ax_d.set_xlabel("|dIF| at full depth  (n per bin: "
                    + ", ".join(f"{int(v):,}" for v in agg["n"]) + ")", fontsize=7.9)
    ax_d.set_ylabel("fraction excluded", color="#1f6fb2")
    ax_d.tick_params(axis="y", colors="#1f6fb2")
    ax_d.set_ylim(0, 1.18)
    ax_d2 = ax_d.twinx()
    ax_d2.plot(centers, agg["med_change"], "s--", color="#c47a2c", lw=1.4, ms=5,
               label="median |$\\Delta$dIF|")
    ax_d2.set_ylabel("median |$\\Delta$dIF|", color="#c47a2c")
    ax_d2.tick_params(axis="y", colors="#c47a2c")
    dx_max = max(0.005, float(agg["med_change"].max()) * 1.6)
    ax_d2.set_ylim(0, dx_max)
    h1, l1 = ax_d.get_legend_handles_labels()
    h2, l2 = ax_d2.get_legend_handles_labels()
    ax_d.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8.7, loc="center right")
    panel(ax_d, "C  instability rises with effect size", fontsize=11.9)
    ax_d.text(0.02, 0.86,
              f"even the smallest bin\nloses {agg['rate'].iloc[0]:.0%}",
              transform=ax_d.transAxes, va="top", fontsize=8.4, color="#444444",
              linespacing=1.3)

    # ================================================================ e
    # replicate SD as a function of per-gene depth: the mechanism behind panels a-d
    nd = pd.read_csv(T / "caf_isoform_noise_by_depth.csv")
    order_e = ["[5, 10)", "[10, 20)", "[20, 40)", "[40, 80)", "[80, 160)",
               "[160, 320)", "[320, 1000000)"]
    nd["depth_bin"] = pd.Categorical(nd["depth_bin"], categories=order_e, ordered=True)
    nd = nd.sort_values("depth_bin")
    xe = np.arange(len(nd))
    ax_e.fill_between(xe, nd["median_replicate_sd"], nd["p90_replicate_sd"],
                      color="#1f6fb2", alpha=0.15, lw=0,
                      label="10th-90th pct")
    ax_e.plot(xe, nd["median_replicate_sd"], "o-", color="#1f6fb2", lw=1.8, ms=6,
              label="median replicate SD")
    ax_e.set_xticks(xe)
    ax_e.set_xticklabels([f"{int(r)}" for r in nd["median_reads"]], fontsize=8.4)
    ax_e.set_xlabel("median reads per gene", fontsize=7.8)
    ax_e.set_ylabel("within-group SD of isoform fraction")
    ax_e.legend(frameon=False, fontsize=8.4, loc="upper right")
    panel(ax_e, "E  why depth matters: noise is set by depth", fontsize=11.9)
    ax_e.axhline(0.05, color="#b8342b", ls=":", lw=1.1)
    ax_e.text(0.02, 0.052, "SD = 0.05", color="#b8342b", fontsize=8.2, va="bottom")
    ax_e.text(0.98, 0.50, "group n:\n" + "\n".join(f"{int(v):,}" for v in nd["n_groups"]),
              transform=ax_e.transAxes, ha="right", va="top", fontsize=6.4,
              color="#666666", linespacing=1.15)

    # ================================================================ f
    # the per-candidate consequence: which panel genes survive the gate
    pc = pd.read_csv(T / "caf_panel_depth_consistency_fixed_referent.csv")
    pc = pc.dropna(subset=["mean_dIF", "dIF_sub"]).copy()
    surv = pc["abs_change"] < 0.05            # same cut-off used by the gate
    ax_f.scatter(pc.loc[~surv, "mean_dIF"], pc.loc[~surv, "dIF_sub"], s=44,
                 color="#9aa3ad", edgecolor="white", lw=0.8, label="unstable", zorder=3)
    ax_f.scatter(pc.loc[surv, "mean_dIF"], pc.loc[surv, "dIF_sub"], s=52,
                 color="#b8342b", edgecolor="white", lw=0.8, label="stable", zorder=4)
    mm = max(pc["mean_dIF"].abs().max(), pc["dIF_sub"].abs().max()) * 1.15
    ax_f.plot([-mm, mm], [-mm, mm], ls="--", lw=0.9, color="#333333")
    ax_f.axhline(0, color="#dddddd", lw=0.7)
    ax_f.axvline(0, color="#dddddd", lw=0.7)
    # Place the four labels in genuinely empty space and pull a thin leader
    # line to the point. Fixed offsets put TIMP1 on top of the cluster, so the
    # label anchors are absolute data coordinates in the empty lower-left
    # quadrant (no panel gene sits below x = -0.09 there).
    lab_pos = {
        "FBLN1": None,                      # point sits alone in the upper right
        "MMP2":  None,                      # point sits alone to the right
        "DCN":   (-0.192, -0.118),
        "TIMP1": (-0.192, -0.172),
    }
    for _, r in pc.iterrows():
        sym = r["symbol"]
        if sym not in lab_pos:
            continue
        colr = "#b8342b" if sym == "FBLN1" else "#444444"
        if lab_pos[sym] is None:
            off = {"FBLN1": (-6, 8), "MMP2": (7, 4)}[sym]
            ax_f.annotate(sym, (r["mean_dIF"], r["dIF_sub"]),
                          xytext=off, textcoords="offset points", fontsize=8.4,
                          ha="right" if off[0] < 0 else "left", color=colr)
        else:
            ax_f.annotate(
                sym, xy=(r["mean_dIF"], r["dIF_sub"]),
                xytext=lab_pos[sym], textcoords="data",
                fontsize=8.4, ha="left", va="center", color=colr,
                arrowprops=dict(arrowstyle="-", lw=0.6, color="#999999",
                                shrinkA=1, shrinkB=2),
            )

    ax_f.set_xlim(-mm, mm)
    ax_f.set_ylim(-mm * 1.06, mm * 1.18)   # headroom above for FBLN1, margin below for the leader labels
    ax_f.set_xlabel("dIF at full depth")
    ax_f.set_ylabel("dIF at 10% depth")
    ax_f.legend(frameon=False, fontsize=8.4, loc="upper left")
    panel(ax_f, f"I  {int(surv.sum())}/{len(pc)} panel genes survive", fontsize=11.9)
    # keep this note in the right-hand empty band so it never shares a row
    # with the TIMP1 leader label on the far left.
    # "most depth-stable" was true only of the withdrawn dominant-isoform statistic;
    # under the fixed referent FBLN1 ranks 2nd in relative terms, so the note now
    # reports both ranks without asserting a superlative.
    big = pc[pc["mean_dIF"].abs() > 0.02].copy()
    big["rel"] = big["abs_change"] / big["mean_dIF"].abs()
    if (big.symbol == "FBLN1").any():
        fbl_rel = float(big.loc[big.symbol == "FBLN1", "rel"].iloc[0])
        rank_rel = int(big["rel"].rank().loc[big.symbol == "FBLN1"].iloc[0])
        rank_abs = int(pc["abs_change"].rank().loc[pc.symbol == "FBLN1"].iloc[0])
        note_f = (f"FBLN1: rank {rank_rel}/{len(big)} relative,\n"
                  f"{rank_abs}/{len(pc)} absolute\n"
                  f"($|\\Delta$dIF|/|dIF| = {fbl_rel:.3f})")
    else:
        note_f = "FBLN1 is not measurable at both depths"
    ax_f.text(0.98, 0.10, note_f,
              transform=ax_f.transAxes, fontsize=8.2, color="#b8342b", va="bottom",
              ha="right", linespacing=1.3)

    # ================================================================ g
    # fixed-referent sweep (183_sweep_fixed_referent.py); the replay of the old
    # definition against this table matched the filed sweep in all 15 cells
    sw = pd.read_csv(T / "caf_switch_threshold_sweep_fixed_referent.csv")
    for fdr, col, mk in ((0.10, "#1f6fb2", "o"), (0.05, "#5b8fb9", "s"), (0.01, "#8fb6d4", "^")):
        d = sw[sw["fdr"] == fdr].sort_values("min_abs_dIF")
        if not len(d):
            continue
        ax_g.plot(d["min_abs_dIF"], d["n_candidates"], marker=mk, lw=1.6, ms=5,
                  color=col, label=f"FDR < {fdr}")
    ax_g.set_xlabel("minimum |dIF| required")
    ax_g.set_ylabel("candidate switches")
    ax_g.set_ylim(-0.5, max(sw["n_candidates"].max() + 1, 7))
    ax_g.legend(frameon=False, fontsize=8.4, loc="center right", bbox_to_anchor=(1.0, 0.42))
    panel(ax_g, "H  the count collapses as the threshold tightens", fontsize=11.9)
    ax_g.axvline(0.15, color="#b8342b", ls=":", lw=1.1)
    # anchor from the axis top so the label follows the data instead of assuming
    # the old 6-candidate maximum
    y_ann = ax_g.get_ylim()[1] * 0.74
    ax_g.annotate("0.15\npre-registered", xy=(0.15, y_ann), xytext=(0.108, y_ann),
                  fontsize=7.9, color="#b8342b", ha="right", va="center",
                  linespacing=1.25,
                  arrowprops=dict(arrowstyle="-", color="#b8342b", lw=0.8,
                                  shrinkA=0, shrinkB=2))

    # ================================================================ h
    dr = pd.read_csv(T / "caf_depth_requirement.csv").sort_values("target_replicate_sd")
    ax_h.plot(dr["target_replicate_sd"], dr["reads_per_gene_needed"], "o-",
              color="#1f6fb2", lw=1.8, ms=6)
    for _, r in dr.iterrows():
        ax_h.annotate(f"{int(r['reads_per_gene_needed'])}",
                      (r["target_replicate_sd"], r["reads_per_gene_needed"]),
                      xytext=(5, 5), textcoords="offset points", fontsize=8.4,
                      color="#1f6fb2")
    cur = 465.0     # median reads/gene in the deepest bin of panel e
    # what SD does that depth actually deliver? invert the requirement curve
    sd_at_cur = float(np.interp(cur, dr["reads_per_gene_needed"],
                               dr["target_replicate_sd"]))
    ax_h.axhline(cur, color="#c47a2c", ls="--", lw=1.2)
    # NOTE: x must be in AXES fraction here. Using data coords (x=0.98 while the
    # data range is 0.02-0.05) makes bbox_inches="tight" expand the canvas ~20x,
    # which silently produced a 29,273 px wide PNG.
    ax_h.text(0.97, 0.88, f"deepest bin here: {cur:.0f} reads/gene "
                          f"$\\rightarrow$ SD {sd_at_cur:.2f}",
              transform=ax_h.transAxes, ha="right", va="top",
              fontsize=8.2, color="#c47a2c")
    ax_h.set_xlabel("target replicate SD of isoform fraction")
    ax_h.set_ylabel("reads per gene required")
    panel(ax_h, "F  what depth would have been needed", fontsize=11.9)
    ax_h.invert_xaxis()
    ax_h.set_ylim(0, max(dr["reads_per_gene_needed"].max(), cur) * 1.30)

    # ================================================================ i
    # the other limit on detection: how many multi-isoform genes are even
    # distinguishable by unique k-mers (k=31, >=10% unique within the gene)
    ident = pd.read_csv(T / "isoform_identifiability_genes.tsv", sep="\t")
    bins_n = [1, 2, 3, 4, 6, 10, 10000]
    labs_n = ["2", "3", "4", "5-6", "7-10", ">10"]
    ident["bin"] = pd.cut(ident["n_tx"], bins_n, labels=labs_n)
    rr = ident.groupby("bin", observed=True).agg(n=("identifiable", "size"),
                                                 k=("identifiable", "sum"))
    rr["frac"] = rr["k"] / rr["n"]
    xi = np.arange(len(rr))
    ax_i.bar(xi, rr["frac"], color="#7b4b8a", width=0.66, edgecolor="white", lw=0.6)
    for x, (_, row) in zip(xi, rr.iterrows()):
        ax_i.text(x, row["frac"] + 0.015, f"{row['frac']:.0%}", ha="center", fontsize=8.7,
                  color="#333333")
        ax_i.text(x, 0.02, f"n={int(row['n']):,}".replace(",", "\n"),
                  ha="center", va="bottom", fontsize=6.6, color="white",
                  linespacing=0.95)
    ax_i.axhline(ident["identifiable"].mean(), color="#b8342b", ls=":", lw=1.1)
    # place the mean label clear of the last bar: the bars grow rightwards, so a
    # right-anchored label at the final x collides with it. Anchor above the line
    # at the left end instead, where the bars are still short.
    ax_i.text(-0.42, ident["identifiable"].mean() - 0.055,
              f"overall {ident['identifiable'].mean():.0%}", fontsize=7.9,
              color="#b8342b", ha="left", va="top")
    ax_i.set_xticks(xi)
    ax_i.set_xticklabels(rr.index, fontsize=9.0)
    ax_i.set_xlabel("transcripts in the gene")
    ax_i.set_ylabel("fraction identifiable")
    ax_i.set_ylim(0, 1.16)
    panel(ax_i, "G  the other limit: k-mer identifiability", fontsize=11.9)
    # the note must sit in the empty band above the tallest bar
    ax_i.text(0.015, 1.14,
              f"{len(ident):,} multi-isoform genes; {int(ident['identifiable'].sum()):,} "
              f"with $\\geq$2 distinguishable transcripts",
              transform=ax_i.transAxes, fontsize=7.4, color="#444444", va="top",
              ha="left", linespacing=1.25)

    for ax in (ax_a, ax_b, ax_c, ax_d, ax_e, ax_f, ax_g, ax_h, ax_i):
        ax.tick_params(labelsize=8.5)

    fig.tight_layout(rect=(0, 0, 1, 1), w_pad=1.6, h_pad=1.8)

    OUT = FIG / "fig_depth_gate.png"
    fit_panel_titles(fig)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    plt.close(fig)
    # guard: a mis-placed annotation in data coordinates can make bbox_inches
    # expand the canvas enormously (a 29,273 px wide PNG once slipped through).
    from PIL import Image
    with Image.open(OUT) as im:
        wpx, hpx = im.size
    assert wpx < 7000 and hpx < 7000, f"suspicious canvas {wpx}x{hpx} for {OUT.name}"
    print("written", OUT, f"{wpx}x{hpx}")


if __name__ == "__main__":
    main()
