"""Fig 7 (6 panels): why FBLN1 C:D does not travel between preparations.

# PANEL-SCALE-1.32 (2026-10-22): fonts enlarged to match Fig 2.

a  molar C:D per preparation, against the GSE303690 reference ranges
b  reads ratio vs molar ratio: the length correction that makes cohorts comparable
c  isoform composition per library -- not just C:D, the whole distribution differs
d  library quality metrics that differ between the two cohorts
e  FBLN1 detection level across the two cohorts (is the shift a detection floor?)
f  effective length of C vs D per library (the mechanism behind the correction)

Panels c-f are deliberately facets of ONE difference: the two cohorts differ in
sorting, library preparation and fragment length at the same time. Panels e and f
are not independent evidence about biology; they test two mechanical explanations
(detection floor, length bias) and show that neither is sufficient.

Outputs figures/fig_fbln1_cross_system.png
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
REM = ROOT / "remote/easyai"
FIG = ROOT / "figures"
STYLE = Path(__file__).resolve().parent / "_style_runtime.py"
exec(compile(STYLE.read_text(encoding="utf-8"), "runtime.py", "exec"), globals())

REF_CAF = (0.264, 0.123, 0.598)
REF_NF = (0.727, 0.519, 1.922)
MIN_READS = 10


def main() -> None:
    apply_figure_style(sizes=(11, 9.5, 8.5))
    raw = pd.read_csv(REM / "caf_fbln1_isoform.csv")   # all 11 libraries
    d = raw.dropna(subset=["C_over_D_molar"]).copy()
    usable = d[d["C_plus_D_reads"] >= MIN_READS].copy()
    n, dropped = len(usable), len(d) - len(usable)
    med = float(usable["C_over_D_molar"].median())
    print(f"quantified={len(d)} usable={n} dropped={dropped} median molar={med:.3f}")

    fig, axes = plt.subplots(3, 2, figsize=(11.79, 14.36))
    ax_a, ax_b, ax_f, ax_c, ax_d, ax_e = axes.ravel()
    rng = np.random.default_rng(11)

    # ---------------------------------------------------------------- a
    ax_a.fill_between([-0.45, 0.45], REF_CAF[1], REF_CAF[2], color="#1f6fb2", alpha=0.14, lw=0)
    ax_a.fill_between([0.55, 1.45], REF_NF[1], REF_NF[2], color="#c47a2c", alpha=0.14, lw=0)
    ax_a.hlines([REF_CAF[0]], -0.45, 0.45, colors="#1f6fb2", lw=1.5, ls="--")
    ax_a.hlines([REF_NF[0]], 0.55, 1.45, colors="#c47a2c", lw=1.5, ls="--")
    x = 2.0 + rng.uniform(-0.22, 0.22, size=n)
    ax_a.scatter(x, usable["C_over_D_molar"], s=56, color="#7b4b8a", edgecolor="white",
                 lw=0.8, zorder=4)
    ax_a.set_xticks([0, 1, 2])
    ax_a.set_xticklabels(["GSE303690\nCAF (6 paired)", "GSE303690\nnormal fibro.",
                          f"GSE195832\nsorted CAF ({n})"], fontsize=9.8)
    ax_a.set_ylabel("FBLN1 C:D molar ratio")
    ax_a.set_xlim(-0.6, 2.6)
    ax_a.set_ylim(0, 3.9)
    ax_a.text(0, REF_CAF[2] + 0.09, f"median {REF_CAF[0]:.2f}", ha="center", fontsize=9.2,
              color="#1f6fb2")
    ax_a.text(1, REF_NF[2] + 0.09, f"median {REF_NF[0]:.2f}", ha="center", fontsize=9.2,
              color="#c47a2c")
    below = int((usable["C_over_D_molar"] < REF_CAF[2]).sum())
    ax_a.text(2, 3.62, f"median {med:.2f}\n{below}/{n} below the CAF upper bound",
              ha="center", va="top", fontsize=9.0, color="#5c3669")
    panel(ax_a, "A  sorted CAF do not reproduce the paired-CAF value", fontsize=12.4)

    # ---------------------------------------------------------------- b
    ax_b.scatter(usable["C_over_D_reads"], usable["C_over_D_molar"], s=56,
                 color="#7b4b8a", edgecolor="white", lw=0.8, zorder=4)
    lim = [0, max(usable["C_over_D_reads"].max(), 1.8) * 1.1]
    ax_b.plot(lim, lim, ls="--", lw=1.0, color="#333333")
    for _, r in usable.iterrows():
        ax_b.plot([r["C_over_D_reads"], r["C_over_D_reads"]],
                  [r["C_over_D_reads"], r["C_over_D_molar"]], color="#cccccc", lw=0.8)
    ax_b.set_xlim(lim)
    ax_b.set_ylim(0, max(3.6, usable["C_over_D_molar"].max() * 1.15))
    ax_b.set_xlabel("C:D from raw read counts")
    ax_b.set_ylabel("C:D corrected to molar ratio")
    panel(ax_b, "B  effective length shifts C:D by up to 2$\\times$", fontsize=12.4)
    ratio = (usable["C_over_D_molar"] / usable["C_over_D_reads"])
    ax_b.text(0.03, 0.95,
              f"correction factor {ratio.min():.2f}-{ratio.max():.2f}$\\times$\n"
              f"(C eff.len {usable['C_efflen'].median():.0f} bp vs\n"
              f" D eff.len {usable['D_efflen'].median():.0f} bp)\n"
              f"->  cross-cohort comparison needs the molar scale",
              transform=ax_b.transAxes, va="top", fontsize=8.7, color="#444444",
              linespacing=1.3)

    # ---------------------------------------------------------------- c
    sh = pd.read_csv(T / "fbln1_isoform_shares.tsv", sep="\t")
    iso = [("FBLN1-202", "1D"), ("FBLN1-201", "1C"), ("FBLN1-203", "1A")]
    for k, (cond, col) in enumerate((("NF", "#c47a2c"), ("CAF", "#1f6fb2"))):
        s = sh[sh["condition"].str.upper() == cond]
        for _, r in s.iterrows():
            ax_c.plot(range(len(iso)), [r[c] for c, _ in iso], color=col, alpha=0.35, lw=1.0,
                      marker="o", ms=4)
        ax_c.plot(range(len(iso)), [s[c].mean() for c, _ in iso], color=col, lw=2.4,
                  marker="s", ms=7, label=f"{cond} (GSE303690, mean)")
    for _, r in usable.iterrows():
        prof = [r["D_fraction"], r["C_fraction"], np.nan]
        ax_c.plot([0, 1], [r["D_fraction"], r["C_fraction"]], color="#7b4b8a", alpha=0.45,
                  lw=1.0, marker="^", ms=4)
    ax_c.plot([0, 1], [usable["D_fraction"].mean(), usable["C_fraction"].mean()],
              color="#7b4b8a", lw=2.4, marker="s", ms=7,
              label="sorted CAF (GSE195832, mean)")
    ax_c.set_xticks(range(len(iso)))
    ax_c.set_xticklabels([f"{lab}\n({c.split('-')[1]})" for c, lab in iso], fontsize=9.8)
    ax_c.set_ylabel("share of FBLN1 reads")
    ax_c.legend(frameon=False, fontsize=8.4, loc="upper right",
                bbox_to_anchor=(1.0, 0.99))
    panel(ax_c, "D  the whole isoform profile differs, not just C:D", fontsize=12.4)
    ax_c.text(0.03, 0.14,
              "purple lines are the 7 sorted-CAF libraries;\n"
              "1A is not shown for them (fraction only\n"
              "reported for C and D)",
              transform=ax_c.transAxes, fontsize=8.2, color="#555555", va="bottom",
              linespacing=1.3)

    # ---------------------------------------------------------------- d
    metrics = [
        ("FBLN1 total reads", usable["FBLN1_total_reads"], 1),
        ("C+D reads", usable["C_plus_D_reads"], 1),
        ("FBLN1 TPM", usable["FBLN1_TPM"], 1),
    ]
    ys = np.arange(len(metrics))[::-1]
    for y_, (lab, vals, scale) in zip(ys, metrics):
        v = vals.astype(float)
        ax_d.scatter(v, np.full(len(v), y_), s=46, color="#7b4b8a", alpha=0.8,
                     edgecolor="white", lw=0.7, zorder=4)
        ax_d.text(v.median(), y_ + 0.22, f"median {v.median():.4g}", ha="center",
                  fontsize=8.7, color="#5c3669")
    ax_d.set_yticks(ys)
    ax_d.set_yticklabels([m[0] for m in metrics], fontsize=9.5)
    ax_d.set_xscale("log")
    ax_d.set_xlabel("value per library (log scale)")
    ax_d.set_ylim(-0.6, len(metrics) - 0.2)
    panel(ax_d, f"E  depth is adequate in all {n} usable libraries", fontsize=12.4)
    ax_d.text(0.98, 0.03,
              f"only libraries with C+D $\\geq$ {MIN_READS} reads are shown\n"
              f"({dropped} of {len(d)} quantified libraries excluded)",
              transform=ax_d.transAxes, ha="right", va="bottom", fontsize=8.4,
              color="#555555")

    # ---------------------------------------------------------------- e
    # is the C:D shift simply a detection floor in shallower libraries?
    # NB: uses ALL sorted libraries (raw), not only the ones with usable C:D --
    # the shallow libraries are exactly the evidence about a detection floor.
    paired_tot = sh["total_reads"].astype(float).to_numpy()
    sorted_tot = raw["FBLN1_total_reads"].astype(float).to_numpy()
    groups = [("GSE303690\npaired (12)", paired_tot, "#c47a2c"),
              (f"GSE195832\nsorted ({len(sorted_tot)})", sorted_tot, "#7b4b8a")]
    for k, (lab, vals, col) in enumerate(groups):
        xv = k + rng.uniform(-0.16, 0.16, size=len(vals))
        ax_e.scatter(xv, vals, s=52, color=col, edgecolor="white", lw=0.8, zorder=4)
        ax_e.hlines(np.median(vals), k - 0.34, k + 0.34, colors=col, lw=1.6, ls="--")
        ax_e.text(k + 0.34, np.median(vals) * 1.05, f"median {np.median(vals):,.0f}",
                  ha="left", va="bottom", fontsize=8.6, color=col)
    ax_e.set_yscale("log")
    ax_e.set_xticks([0, 1])
    ax_e.set_xticklabels([lab for lab, _, _ in groups], fontsize=9.8)
    ax_e.set_xlim(-0.6, 1.6)
    ax_e.set_ylim(0.7, max(paired_tot.max(), sorted_tot.max()) * 95)
    ax_e.set_ylabel("FBLN1 reads per library (log)")
    fold = np.median(paired_tot) / np.median(sorted_tot)
    panel(ax_e, "F  the sorted cohort is far shallower at FBLN1", fontsize=12.4)
    # the note has its own band at the top, clear of the two point clouds
    n_below = int((sorted_tot < paired_tot.min()).sum())
    ax_e.text(0.985, 0.995,
              f"median {fold:.1f}$\\times$ lower ({np.median(sorted_tot):,.0f} vs "
              f"{np.median(paired_tot):,.0f} reads).\n"
              f"{n_below} of {len(sorted_tot)} sorted libraries fall below the "
              f"shallowest\npaired library ({paired_tot.min():,.0f} reads), so for those a "
              f"detection\nfloor cannot be excluded; the other "
              f"{len(sorted_tot) - n_below} reach\nnormal depth and still give the same "
              f"high C:D (panel a).",
              transform=ax_e.transAxes, ha="right", va="top", fontsize=7.2,
              color="#444444", linespacing=1.3)

    # ---------------------------------------------------------------- f
    # the mechanism behind panel b: D has roughly twice the effective length of C
    eff = d.dropna(subset=["C_efflen", "D_efflen"]).copy()
    eff["ratio"] = eff["C_efflen"] / eff["D_efflen"]
    for k, r in enumerate(eff.itertuples()):
        ax_f.plot([0, 1], [r.C_efflen, r.D_efflen], color="#c9d2da", lw=1.4, zorder=1)
        ax_f.scatter([0, 1], [r.C_efflen, r.D_efflen], s=44, color="#7b4b8a",
                     edgecolor="white", lw=0.8, zorder=4)
    ax_f.scatter([0, 0], [eff["C_efflen"].median()] * 2, marker="_", s=340,
                 color="#222222", lw=2.0, zorder=5)
    ax_f.scatter([1, 1], [eff["D_efflen"].median()] * 2, marker="_", s=340,
                 color="#222222", lw=2.0, zorder=5)
    ax_f.set_xticks([0, 1])
    ax_f.set_xticklabels(["1C (FBLN1-201)", "1D (FBLN1-202)"], fontsize=9.8)
    ax_f.set_xlim(-0.45, 1.45)
    ax_f.set_ylim(0, eff["D_efflen"].max() * 1.22)
    ax_f.set_ylabel("salmon effective length (bp)")
    panel(ax_f, "C  D is $\\sim$2$\\times$ longer, which biases read ratios",
                   fontsize=12.4)
    ax_f.text(0.5, 0.04,
              f"median {eff['C_efflen'].median():,.0f} vs {eff['D_efflen'].median():,.0f} bp\n"
              f"C/D ratio {eff['ratio'].min():.2f}-{eff['ratio'].max():.2f}\n"
              f"-> reads must be divided by this to get molar C:D",
              transform=ax_f.transAxes, ha="center", va="bottom", fontsize=7.9,
              color="#444444", linespacing=1.3)

    fig.tight_layout(rect=(0, 0, 1, 1), w_pad=1.4, h_pad=1.9)
    fit_panel_titles(fig)
    fig.savefig(FIG / "fig_fbln1_cross_system.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    from PIL import Image
    with Image.open(FIG / "fig_fbln1_cross_system.png") as im:
        wpx, hpx = im.size
    assert wpx < 7000 and hpx < 7000, f"suspicious canvas {wpx}x{hpx}"

    print(f"panel b correction factor: {ratio.min():.2f}-{ratio.max():.2f}x")
    print(f"panel e: median FBLN1 reads {np.median(sorted_tot):,.0f} vs "
          f"{np.median(paired_tot):,.0f} ({fold:.1f}x lower)")
    print(f"panel f: efflen median C={eff['C_efflen'].median():.0f} "
          f"D={eff['D_efflen'].median():.0f}, ratio {eff['ratio'].min():.2f}-"
          f"{eff['ratio'].max():.2f}")
    print(f"written figures/fig_fbln1_cross_system.png {wpx}x{hpx}")


if __name__ == "__main__":
    main()
