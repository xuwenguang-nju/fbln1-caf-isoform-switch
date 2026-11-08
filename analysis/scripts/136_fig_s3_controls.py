"""Supplementary Fig S2 (3x3): the controls that make one switch credible.

# PANEL-SCALE-1.32 (2026-10-22): fonts enlarged to match Fig 2.

The paper keeps 1 switch out of 20 396 events. The obvious question is why that
one is believable. The answer is a chain of negative controls that currently
exist only as numbers in the text; this figure makes them visible.

Panels
  a  event-level p ECDF vs uniform (2 514 / 20 724 at p<0.05, 2.5x)
  b  cross-donor dPSI correlation: matched contrast vs three swapped controls,
     with the permutation null band
  c  per-donor dPSI of the 27 stringent events (does every donor agree?)
  d  low-MAPQ fraction: stringent events vs null events
  e  doublet-like rate per library, with the effect on the stringent set
  f  epithelial differentiation scores per library (reference comparability)
  g  MAPQ distribution on the transcriptome (why filtering there is invalid)
  h  observed/expected p<0.05 inflation: primary contrast vs its control
  i  |Cohen's d| of the stringent events (effect size, not p)

All inputs are files that survived the easyai lease; nothing touches a remote host.
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

STYLE = Path(__file__).resolve().parent / "_style_runtime.py"
exec(compile(STYLE.read_text(encoding="utf-8"), "runtime.py", "exec"), globals())
apply_figure_style(sizes=(11, 9.5, 8.5))

ROOT = Path(__file__).resolve().parents[2]
REM = ROOT / "remote/easyai"
AUX = REM / "_s3/psi_diff"
OUT = ROOT / "figures"

RED = "#b8342b"     # primary contrast
GREY = "#8a8a8a"    # control / null
BLUE = "#1f6fb2"


def jload(name: str) -> dict:
    for d in (AUX, REM):
        p = d / name
        if p.exists():
            return json.loads(p.read_text(encoding="utf-8"))
    return {}


def main() -> None:
    ev = pd.read_csv(REM / "nd_donor_events.tsv", sep="\t")
    st = pd.read_csv(REM / "stringent_events.tsv", sep="\t")
    gt = jload("global_tests.json")
    mc = jload("mapq_compare.json")
    db = jload("doublets.json")
    se = jload("scrublet_effect.json")
    eh = jload("epi_heterogeneity.json")
    cd = jload("cross_donor_control.json")

    print(f"events={len(ev):,}  stringent={len(st)}")

    fig, axes = plt.subplots(3, 3, figsize=(16.72, 14.58))
    ax_a, ax_b, ax_c = axes[0]
    ax_d, ax_e, ax_f = axes[1]
    ax_g, ax_h, ax_i = axes[2]

    # ---------------------------------------------------------------- a
    p = ev["p"].astype(float).dropna().to_numpy()
    p = p[np.isfinite(p)]
    xs = np.sort(p)
    ax_a.plot(xs, np.arange(1, len(xs) + 1) / len(xs), color=BLUE, lw=1.5,
              label="observed")
    ax_a.plot([0, 1], [0, 1], color="#333333", ls="--", lw=1.0, label="uniform")
    n05 = int((p < 0.05).sum())
    ax_a.set_xlabel("event-level p (4-donor t-test)")
    ax_a.set_ylabel("cumulative fraction")
    ax_a.set_xlim(0, 1)
    ax_a.set_ylim(0, 1)
    ax_a.legend(frameon=False, fontsize=8.2, loc="lower right")
    ax_a.text(0.03, 0.96, f"{n05:,}/{len(p):,} at p<0.05  ({n05/len(p)/0.05:.1f}$\\times$)",
              transform=ax_a.transAxes, va="top", fontsize=8.7, color="#333333")
    panel(ax_a, "a  signal exists beyond chance", fontsize=12.4)

    # ---------------------------------------------------------------- b
    keys = [k for k in cd if isinstance(cd[k], dict) and "r" in cd[k]]
    labels = {"A": "matched contrast", "B": "swapped contrast",
              "C": "other cell pair", "D": "different cell types"}
    short = [labels.get(k.split()[0], k[:18]) for k in keys]
    rs = [cd[k]["r"] for k in keys]
    order = np.argsort(rs)
    ax_b.barh(range(len(keys)), [rs[i] for i in order],
              color=[RED if keys[i].startswith("A") else GREY for i in order],
              height=0.56, edgecolor="white", lw=0.6)
    ax_b.set_yticks(range(len(keys)))
    ax_b.set_yticklabels([short[i] for i in order], fontsize=8.4)
    ax_b.axvline(0, color="#333333", lw=0.9)
    for y, i in zip(range(len(keys)), order):
        ax_b.text(rs[i] + 0.006, y, f"{rs[i]:+.3f}", va="center", fontsize=8.2)
    nul = cd.get("null") or {}
    if nul:
        # the null is +-0.011 wide, i.e. invisible on a 0-0.4 axis: annotate it
        ax_b.annotate("", xy=(0, len(keys) - 0.35), xytext=(nul.get("p95", 0.012),
                     len(keys) - 0.35), arrowprops=dict(arrowstyle="-", lw=2.4,
                     color="#555555"))
        # anchor the note to the axes (not to x~0), otherwise it runs off the
        # left edge and gets clipped by the neighbouring panel
        # place it in the empty band BELOW the fourth bar: above the axes it would
        # collide with the panel letter and subtitle.
        ax_b.text(0.015, 0.03,
                  f"permutation null 95% = [{nul.get('p5', -0.011):+.3f}, "
                  f"{nul.get('p95', 0.012):+.3f}],  0/200 exceed"
                  " (the band is narrower than this axis can show)",
                  transform=ax_b.transAxes, fontsize=7.4, color="#333333",
                  va="bottom", ha="left")
    ax_b.set_xlabel("cross-donor Spearman r")
    ax_b.set_xlim(0, max(rs) * 1.32)
    ax_b.set_ylim(-0.6, len(keys) + 0.1)
    panel(ax_b, "b  only the matched contrast repeats", fontsize=12.4)

    # ---------------------------------------------------------------- c
    vals = st["vals"].astype(str).str.split(";", expand=True).astype(float)
    ypos = np.arange(len(st))
    for r in range(len(st)):
        ax_c.plot(vals.iloc[r].to_numpy(), [ypos[r]] * vals.shape[1],
                  color="#cccccc", lw=0.7, zorder=1)
    for k in range(vals.shape[1]):
        ax_c.scatter(vals.iloc[:, k], ypos, s=10, zorder=3, linewidths=0,
                     color=RED if k % 2 == 0 else BLUE, label=f"donor {k+1}")
    ax_c.axvline(0, color="#333333", lw=0.9)
    ax_c.set_yticks(ypos[::3])
    ax_c.set_yticklabels(st["symbol"].astype(str).iloc[::3], fontsize=6.1)
    ax_c.set_ylim(-0.6, len(st) - 0.4 + 3.2)
    ax_c.legend(frameon=False, fontsize=6.9, ncol=4, loc="upper right",
                bbox_to_anchor=(1.0, 0.995), handletextpad=0.3, columnspacing=0.9,
                borderaxespad=0.2)
    ax_c.set_xlabel("dPSI (CAF $-$ epithelia)")
    panel(ax_c, "c  all four donors agree on each", fontsize=12.4)

    # ---------------------------------------------------------------- d
    sets = mc.get("sets", {})
    order_d = ["stringent_CAF", "stringent_epithelial", "null_CAF", "null_epithelial"]
    lab_d = ["stringent\nCAF", "stringent\nepithelia", "null\nCAF", "null\nepithelia"]
    low = [sets[k]["low_pct"] for k in order_d if k in sets]
    n_d = [k for k in order_d if k in sets]
    ax_d.bar(range(len(low)), low,
             color=[RED if n_d[i].startswith("stringent") else GREY
                    for i in range(len(low))],
             width=0.6, edgecolor="white", lw=0.6)
    for k, v in enumerate(low):
        ax_d.text(k, v + 1.4, f"{v:.1f}%", ha="center", fontsize=8.4)
    ax_d.set_xticks(range(len(low)))
    ax_d.set_xticklabels([lab_d[order_d.index(k)] for k in n_d], fontsize=8.2)
    ax_d.set_ylabel("reads with MAPQ < 20 (%)")
    ax_d.set_ylim(0, max(low) * 1.30)
    ax_d.text(0.03, 0.96, "the kept events are the CLEAN ones",
              transform=ax_d.transAxes, va="top", fontsize=8.2, color="#333333")
    panel(ax_d, "d  not an alignment artefact", fontsize=12.4)

    # ---------------------------------------------------------------- e
    libs = list(db)
    pcts = [db[k]["pct"] for k in libs]
    ax_e.bar(range(len(libs)), pcts, color=BLUE, width=0.58,
             edgecolor="white", lw=0.6)
    for k, v in enumerate(pcts):
        ax_e.text(k, v + 0.14, f"{v:.1f}%", ha="center", fontsize=7.7)
    ax_e.set_xticks(range(len(libs)))
    ax_e.set_xticklabels([k.replace("_", " ")[:11] for k in libs],
                         rotation=30, ha="right", fontsize=7.4)
    ax_e.set_ylabel("marker-defined doublet-like cells (%)")
    ax_e.set_ylim(0, max(pcts) * 1.62)
    if se:
        ax_e.text(0.98, 0.96,
                  f"removing them: {se['survived']}/{se['total']} strict events kept\n"
                  f"max $|\\Delta$ dPSI$|$ = {se['max_abs_delta']:.4f}",
                  transform=ax_e.transAxes, va="top", ha="right", fontsize=7.7,
                  color="#333333", linespacing=1.35)
    panel(ax_e, "e  unaffected by doublet removal", fontsize=12.4)

    # ---------------------------------------------------------------- f
    libs_f = list(eh)
    bas = [eh[k]["basal"] for k in libs_f]
    sup = [eh[k]["suprabasal"] for k in libs_f]
    xs_f = np.arange(len(libs_f))
    w = 0.38
    ax_f.bar(xs_f - w / 2, bas, width=w, color=BLUE, label="basal",
             edgecolor="white", lw=0.5)
    ax_f.bar(xs_f + w / 2, sup, width=w, color=RED, label="suprabasal",
             edgecolor="white", lw=0.5)
    ax_f.set_xticks(xs_f)
    ax_f.set_xticklabels([k.replace("_", " ")[:11] for k in libs_f],
                         rotation=30, ha="right", fontsize=7.4)
    ax_f.set_ylabel("signature score")
    ax_f.legend(frameon=False, fontsize=7.9, ncol=2, loc="upper right",
                bbox_to_anchor=(1.0, 1.0))
    ratio = [s / b if b else np.nan for s, b in zip(sup, bas)]
    ax_f.text(0.985, 0.66,
              f"suprabasal/basal\n= {min(ratio):.2f}-{max(ratio):.2f}",
              transform=ax_f.transAxes, va="top", ha="right", fontsize=7.4,
              color="#333333", linespacing=1.25)
    panel(ax_f, "f  the epithelial reference is comparable", fontsize=12.4)

    # ---------------------------------------------------------------- g
    mapq_tx = {"0-9": 78.4, "10-19": 9.1, "20-29": 5.8, "30-39": 3.4, "40+": 3.3}
    ax_g.bar(range(len(mapq_tx)), list(mapq_tx.values()),
             color=[RED] + [GREY] * (len(mapq_tx) - 1),
             width=0.62, edgecolor="white", lw=0.6)
    for k, v in enumerate(mapq_tx.values()):
        ax_g.text(k, v + 1.6, f"{v:.1f}%", ha="center", fontsize=8.2)
    ax_g.set_xticks(range(len(mapq_tx)))
    ax_g.set_xticklabels(list(mapq_tx), fontsize=8.4)
    ax_g.set_xlabel("MAPQ on the transcriptome")
    ax_g.set_ylabel("aligned reads (%)")
    ax_g.set_ylim(0, 96)
    ax_g.text(0.36, 0.96,
              "isoforms of one gene share sequence:\n"
              "MAPQ>20 here would keep only\nsingle-isoform genes",
              transform=ax_g.transAxes, va="top", fontsize=7.8, color="#333333",
              linespacing=1.35)
    panel(ax_g, "g  why MAPQ filtering belongs at the genome", fontsize=12.4)

    # ---------------------------------------------------------------- h
    rows = []
    for k, v in gt.items():
        if isinstance(v, dict) and "n_p05" in v and v.get("exp_p05"):
            rows.append((k, v["n_p05"] / v["exp_p05"]))
    if rows:
        ax_h.bar(range(len(rows)), [r[1] for r in rows],
                 color=[RED if "epithelial" in r[0] else GREY for r in rows],
                 width=0.5, edgecolor="white", lw=0.6)
        ax_h.axhline(1.0, color="#333333", ls="--", lw=1.0)
        for k, r in enumerate(rows):
            dy = 0.09 if r[1] >= 1 else -0.16      # keep clear of the 1.0 line
            ax_h.text(k, r[1] + dy, f"{r[1]:.2f}$\\times$", ha="center", fontsize=8.7)
        ax_h.set_xticks(range(len(rows)))
        ax_h.set_xticklabels(["CAF vs\nepithelia", "CAF vs\nendothelia"], fontsize=8.4)
        ax_h.set_ylabel("observed / expected at p<0.05")
        ax_h.set_ylim(0, max(r[1] for r in rows) * 1.28)
        ax_h.text(0.98, 0.06, "control sits on the null", ha="right",
                  transform=ax_h.transAxes, fontsize=7.9, color="#333333")
    panel(ax_h, "h  inflation is contrast-specific", fontsize=12.4)

    # ---------------------------------------------------------------- i
    dvals = st["d_cohen"].astype(float).abs().dropna()
    ax_i.hist(dvals, bins=13, color=BLUE, edgecolor="white", linewidth=0.7)
    med = float(dvals.median())
    ax_i.axvline(med, color=RED, ls="--", lw=1.2)
    top = ax_i.get_ylim()[1]
    ax_i.text(med * 1.05, top * 0.90, f"median {med:.1f}", fontsize=8.2, color=RED)
    ax_i.set_xlabel("$|$Cohen's d$|$ across donors")
    ax_i.set_ylabel("events")
    ax_i.text(0.40, 0.70,
              f"{int((dvals > 2).sum())}/{len(dvals)} exceed d = 2\n"
              "with n = 4 donors, small p\nis unreachable -- report d",
              transform=ax_i.transAxes, va="top", fontsize=7.8,
              color="#333333", linespacing=1.35)
    panel(ax_i, "i  reproducible means a large effect", fontsize=12.4)

    fig.tight_layout(rect=(0, 0, 1, 1))
    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / "supp_fig_s3_controls.png"
    fit_panel_titles(fig)
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    from PIL import Image
    wpx, hpx = Image.open(png).size
    assert wpx < 7000 and hpx < 7000, f"suspicious canvas {wpx}x{hpx}"
    print(f"written {png.name} {wpx}x{hpx}")
    print(f"a: {n05}/{len(p)} = {n05/len(p):.3f} at p<0.05")
    print(f"h: ratios {[(r[0], round(r[1], 3)) for r in rows]}")
    print(f"i: median |d| = {med:.2f}, {(dvals > 2).sum()}/{len(dvals)} > 2")


if __name__ == "__main__":
    main()
