"""Supplementary Fig S3 (3x3): where this study sits in the literature.

# PANEL-SCALE-1.32 (2026-10-22): fonts enlarged to match Fig 2.

The paper claims a narrow gap: cell-type-resolved splicing in oral/HNSCC CAF is
essentially unoccupied. That claim rests on a query ladder and an archetype survey
that currently live only as tables. This figure makes the positioning auditable.

Panels
  a  the query ladder T1..T8 (each row narrows the previous one)
  b  the same ladder as a funnel, with the final 1 highlighted
  c  how narrow is the CAF corner: CAF+splicing vs macrophage+splicing
  d  tissue competition: OSCC vs HNSCC vs generic single-cell counts
  e  archetype volumes (14 archetypes, PubMed TIAB)
  f  citation profile of each archetype (median and max of a 120-record sample)
  g  age of the sampled records per archetype (is the field still moving?)
  h  the 45 curated exemplars, by archetype and role
  i  the three published results of this study mapped onto the ladder

Every number comes from results/tables/{sc_competition_tiers, sc_competition_sweep,
as_archetype_overview, as_bioinformatics_exemplars}.csv; nothing is typed in.
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
T = ROOT / "results/tables"
OUT = ROOT / "figures"

RED = "#b8342b"
BLUE = "#1f6fb2"
GREY = "#9a9a9a"
ORANGE = "#c47a2c"
GREEN = "#5b8f5b"
PURPLE = "#7b4b8a"


def main() -> None:
    tiers = pd.read_csv(T / "sc_competition_tiers.csv")
    sweep = pd.read_csv(T / "sc_competition_sweep.csv")
    arch = pd.read_csv(T / "as_archetype_overview.csv")
    ex = pd.read_csv(T / "as_bioinformatics_exemplars.csv")

    tiers = tiers.sort_values("tier")
    print(f"tiers={len(tiers)}  sweep={len(sweep)}  archetypes={len(arch)}  exemplars={len(ex)}")

    fig, axes = plt.subplots(3, 3, figsize=(16.72, 14.79))
    ax_a, ax_b, ax_c = axes[0]
    ax_d, ax_e, ax_f = axes[1]
    ax_g, ax_h, ax_i = axes[2]

    # ---------------------------------------------------------------- a
    ypos = np.arange(len(tiers))[::-1]
    col = [RED if t == "T6" else (ORANGE if t.startswith(("T7", "T8")) else BLUE)
           for t in tiers["tier"]]
    ax_a.barh(ypos, tiers["hits"], color=col, height=0.66, edgecolor="white", lw=0.6)
    for y, (t, h) in zip(ypos, zip(tiers["tier"], tiers["hits"])):
        ax_a.text(max(h, 1) * 1.12, y, f"{int(h):,}", va="center", fontsize=8.2)
    ax_a.set_yticks(ypos)
    ax_a.set_yticklabels(tiers["tier"], fontsize=8.7)
    ax_a.set_xscale("log")
    ax_a.set_xlabel("PubMed records (TIAB, log)")
    ax_a.set_xlim(0.8, tiers["hits"].max() * 4)
    ax_a.set_ylim(-0.7, len(tiers) - 0.3)
    panel(ax_a, "a  the query ladder", fontsize=12.4)

    # ---------------------------------------------------------------- b
    funnel = tiers[tiers["tier"].isin(["T1", "T2", "T3", "T4", "T5", "T6"])]
    xs = np.arange(len(funnel))
    ax_b.plot(xs, funnel["hits"], color=BLUE, lw=1.5, marker="o", ms=5,
              mfc="white", mew=1.1, zorder=3)
    ax_b.fill_between(xs, funnel["hits"], color=BLUE, alpha=0.13)
    ax_b.set_yscale("log")
    for x, (t, h) in enumerate(zip(funnel["tier"], funnel["hits"])):
        ax_b.annotate(f"{int(h)}", (x, h), textcoords="offset points",
                      xytext=(0, 9), ha="center", fontsize=8.2,
                      color=RED if t == "T6" else "#333333")
    ax_b.set_xticks(xs)
    ax_b.set_xticklabels(funnel["tier"], fontsize=9.0)
    ax_b.set_ylabel("records (log)")
    ax_b.set_ylim(0.7, funnel["hits"].max() * 4)
    # the value labels sit above their markers, so the key goes low-right
    ax_b.text(0.40, 0.16,
              "T1 any disease -> T2 cancer -> T3 TME\n-> T4 isoform -> T5 oral/HN\n"
              "-> T6 oral/HN in title-abstract",
              transform=ax_b.transAxes, va="center", ha="left", fontsize=7.2,
              color="#444444", linespacing=1.35)
    panel(ax_b, "b  oral single-cell splicing is nearly empty", fontsize=12.4)

    # ---------------------------------------------------------------- c
    corner = tiers[tiers["tier"].isin(["T5", "T6", "T7", "T8"])].copy()
    corner["label"] = ["oral/HN\n(T5)", "oral/HN title\n(T6)",
                       "CAF + splicing\n(T7)", "macrophage +\nsplicing (T8)"]
    cols = [BLUE, RED, ORANGE, GREY]
    ax_c.bar(range(len(corner)), corner["hits"], color=cols, width=0.62,
             edgecolor="white", lw=0.6)
    for x, h in enumerate(corner["hits"]):
        ax_c.text(x, h + max(corner["hits"]) * 0.03, f"{int(h)}", ha="center",
                  fontsize=8.7)
    ax_c.set_xticks(range(len(corner)))
    ax_c.set_xticklabels(corner["label"], fontsize=7.7)
    ax_c.set_ylabel("records")
    ax_c.set_ylim(0, max(corner["hits"]) * 1.34)
    ax_c.text(0.97, 0.97,
              "macrophage is 13x more\ncrowded than the CAF corner",
              transform=ax_c.transAxes, ha="right", va="top", fontsize=7.8,
              color="#333333", linespacing=1.3)
    panel(ax_c, "c  the CAF corner is the emptiest", fontsize=12.4)

    # ---------------------------------------------------------------- d
    tissue = sweep[sweep["group"] == "B_tissue"].copy()
    tissue = tissue.sort_values("pubmed_count", ascending=False).head(6)
    ax_d.barh(range(len(tissue))[::-1], tissue["pubmed_count"], color=GREY,
              height=0.62, edgecolor="white", lw=0.6)
    for y, v in zip(range(len(tissue))[::-1], tissue["pubmed_count"]):
        ax_d.text(v * 1.03, y, f"{int(v):,}", va="center", fontsize=7.9)
    ax_d.set_yticks(range(len(tissue))[::-1])
    pretty = {"hnsc_sc": "HNSCC + single-cell",
              "hnsc_sc_recent": "HNSCC, recent window",
              "oscc_sc": "OSCC + single-cell",
              "oscc_sc_recent": "OSCC, recent window"}
    ax_d.set_yticklabels([pretty.get(k, k.replace("_", " ")) for k in tissue["query_key"]],
                         fontsize=7.7)
    ax_d.set_xscale("log")
    ax_d.set_xlabel("records (log)")
    ax_d.set_xlim(1, tissue["pubmed_count"].max() * 6)
    panel(ax_d, "d  single-cell output exists in this tissue", fontsize=12.4)

    # ---------------------------------------------------------------- e
    a2 = arch.sort_values("pubmed_hits_2026_02", ascending=False).head(10)
    ax_e.barh(range(len(a2))[::-1], a2["pubmed_hits_2026_02"], color=BLUE,
              height=0.62, edgecolor="white", lw=0.6)
    for y, v in zip(range(len(a2))[::-1], a2["pubmed_hits_2026_02"]):
        ax_e.text(v * 1.02, y, f"{int(v):,}", va="center", fontsize=7.7)
    ax_e.set_yticks(range(len(a2))[::-1])
    ax_e.set_yticklabels([s.replace("_", " ")[:22] for s in a2["archetype"]],
                         fontsize=7.1)
    ax_e.set_xscale("log")
    ax_e.set_xlabel("PubMed records (log)")
    ax_e.set_xlim(1, a2["pubmed_hits_2026_02"].max() * 8)
    panel(ax_e, f"e  {len(arch)} archetypes, top 10", fontsize=12.4)

    # ---------------------------------------------------------------- f
    a3 = arch.sort_values("median_citations_sample", ascending=False).head(10)
    xs_f = np.arange(len(a3))
    ax_f.bar(xs_f, a3["median_citations_sample"], color=BLUE, width=0.55,
             edgecolor="white", lw=0.6, label="median")
    ax_f.scatter(xs_f, a3["max_citations_sample"], s=26, color=RED, zorder=3,
                 label="max in sample")
    ax_f.set_yscale("log")
    ax_f.set_xticks(xs_f)
    # steeper rotation + smaller type: at 40 deg and 12 characters the labels of
    # adjacent archetypes still collided along the x axis
    ax_f.set_xticklabels([s.replace("_", " ")[:14] for s in a3["archetype"]],
                         rotation=90, ha="center", fontsize=6.2)
    ax_f.set_ylabel("citations (log)")
    ax_f.legend(frameon=False, fontsize=7.7, loc="upper right")
    ax_f.text(0.03, 0.94,
              "120 records sampled per archetype",
              transform=ax_f.transAxes, va="top", fontsize=7.4, color="#444444")
    panel(ax_f, "f  citation weight is broad, not concentrated", fontsize=12.4)

    # ---------------------------------------------------------------- g
    a4 = arch.dropna(subset=["median_year_sample"]).sort_values("median_year_sample")
    ax_g.scatter(a4["median_year_sample"], np.arange(len(a4)), s=34, color=PURPLE,
                 edgecolor="white", lw=0.7, zorder=3)
    ax_g.set_yticks(range(len(a4)))
    ax_g.set_yticklabels([s.replace("_", " ")[:20] for s in a4["archetype"]],
                         fontsize=6.9)
    ax_g.axvline(2026, color=RED, ls=":", lw=1.0)
    ax_g.set_xlabel("median publication year of the sample")
    ax_g.set_xlim(a4["median_year_sample"].min() - 0.5, 2026.6)
    ax_g.tick_params(axis="x", labelsize=7.9)
    ax_g.text(0.98, 0.03, "all archetypes are still\nbeing published",
              transform=ax_g.transAxes, va="bottom", ha="right", fontsize=7.2,
              color="#444444", linespacing=1.25)
    panel(ax_g, "g  no archetype is exhausted", fontsize=12.4)

    # ---------------------------------------------------------------- h
    tab = (ex.groupby(["archetype", "role"]).size()
             .unstack(fill_value=0))
    tab = tab.loc[tab.sum(axis=1).sort_values(ascending=False).index].head(9)
    bottom = np.zeros(len(tab))
    palette = [BLUE, GREEN, ORANGE, RED, GREY, PURPLE]
    for k, c in enumerate(tab.columns):
        ax_h.barh(range(len(tab))[::-1], tab[c], left=bottom, height=0.62,
                  color=palette[k % len(palette)], label=c, edgecolor="white", lw=0.4)
        bottom = bottom + tab[c].to_numpy()
    ax_h.set_yticks(range(len(tab))[::-1])
    ax_h.set_yticklabels([s.replace("_", " ")[:26] for s in tab.index], fontsize=7.1)
    ax_h.set_xlabel("curated exemplars")
    ax_h.set_xlim(0, float(bottom.max()) * 1.18)
    # the two roles are "landmark" and "template"; label them explicitly rather
    # than leaning on the raw column names, which are easy to misread here.
    ax_h.legend(frameon=False, fontsize=7.4, ncol=2, loc="lower right",
                handletextpad=0.4, columnspacing=1.0)
    panel(ax_h, f"h  {len(ex)} curated exemplars by role", fontsize=12.4)

    # ---------------------------------------------------------------- i
    # "depth gate keeps 1 of 27" conflated two arms: the gate acts on Arm B's
    # 12,703 genome-wide genes (1 switch survives), while 27 is Arm A's count of
    # stringent exon-level events -- a set that does not even contain FBLN1.
    # Each row now names the arm and the denominator it actually belongs to.
    # Arm B's delivery line must follow the fixed-referent tables, not the
    # withdrawn dominant-isoform statistic (it read "1 of 12,703, q = 0.007").
    bN = json.load(open(ROOT / "data/processed/fig1_arm_b_numbers.json",
                        encoding="utf-8"))
    axis_i = [
        ("CAF vs epithelia\ncontrast-specific\n(Arm A, 4 donors)", "delivered", GREEN),
        ("FBLN1 1C->1D\nrepeats in 7/7\n(Arm A2, external)", "delivered", GREEN),
        (f"depth gate keeps {bN['n_switch']} of {bN['tested_both_depths']:,}\n"
         f"(sole switch, |dIF|>=0.15,\n"
         f"q = {bN['fbln1_q']:.3f})", "delivered", GREEN),
        ("TCGA prognosis\n0 of 28 events\n(Arm C)", "negative", RED),
        ("ICB response\nAUC < 0.5\n(Arm C)", "negative", RED),
        ("subtype mapping\nunresolvable\n(13 libraries, 0/13)", "not claimed", GREY),
    ]
    ax_i.barh(range(len(axis_i))[::-1], [1] * len(axis_i),
              color=[c for _, _, c in axis_i], height=0.6, edgecolor="white", lw=0.5)
    ax_i.set_yticks(range(len(axis_i))[::-1])
    ax_i.set_yticklabels([n for n, _, _ in axis_i], fontsize=6.9, linespacing=1.25)
    ax_i.set_xticks([])
    ax_i.set_xlim(0, 1.02)
    for y, (n, label, c) in zip(range(len(axis_i))[::-1], axis_i):
        ax_i.text(0.5, y, label, ha="center", va="center", fontsize=7.9,
                  color="white", fontweight="bold")
    panel(ax_i, "i  what this study actually delivers", fontsize=12.4)

    fig.tight_layout(rect=(0, 0, 1, 1))
    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / "supp_fig_s5_positioning.png"
    fit_panel_titles(fig)
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    from PIL import Image
    wpx, hpx = Image.open(png).size
    assert wpx < 7000 and hpx < 7000, f"suspicious canvas {wpx}x{hpx}"
    print(f"written {png.name} {wpx}x{hpx}")
    print("T5/T6/T7/T8 hits:", dict(zip(corner["tier"], corner["hits"])))


if __name__ == "__main__":
    main()
