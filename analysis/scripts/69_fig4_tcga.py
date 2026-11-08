"""Fig 4 (4 panels): FBLN1 and the 27 CAF events in TCGA-HNSC bulk.

# PANEL-SCALE-1.32 (2026-10-22): fonts enlarged to match Fig 2.

a  FBLN1 terminal-exon dPSI across three cohorts (the reproducible signal)
b  composition correlation: signed rho of each AT event with the NNLS CAF
   fraction (replaces the marker-PSI proxy whose sign was inverted)
c  27-event screen -- observed vs expected counts of p<0.05 (no FDR survives)
d  direction concordance with the single-cell dPSI, by event

Panels c and d are deliberately included: they show what the bulk data do NOT
support. Omitting them would let a reader assume the 27 events validated.

Outputs figures/fig_fbln1_tcga.png
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
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
T = ROOT / "results/tables"
FIG = ROOT / "figures"
STYLE = Path(__file__).resolve().parent / "_style_runtime.py"
exec(compile(STYLE.read_text(encoding="utf-8"), "runtime.py", "exec"), globals())

COH = {"oral_strict": "oral strict", "oral_wide": "oral wide", "all_hnsc": "all HNSC"}
COH_N = {"oral_strict": "16 paired", "oral_wide": "29 paired", "all_hnsc": "42 paired"}
EV = {"62665": "exon 17", "62666": "exon 22"}
N_ICB = "28 patients: 9 R / 19 NR"


def main() -> None:
    apply_figure_style(sizes=(11, 9.5, 8.5))
    fb = pd.read_csv(T / "fbln1_tcga_events.csv")
    fb["as_id"] = fb["as_id"].astype(str)
    tvn = pd.read_csv(T / "tcga_caf_events_tvn.csv")

    fig, axes = plt.subplots(3, 3, figsize=(16.72, 15.44))
    ax_a, ax_b, ax_c = axes[0]
    ax_d, ax_e, ax_h = axes[1]
    ax_f, ax_g, ax_i = axes[2]
    # ---------------------------------------------------------------- a
    cohorts = ["oral_strict", "oral_wide", "all_hnsc"]
    colors = {"62665": "#1f6fb2", "62666": "#c47a2c"}
    xpos, ticks = [], []
    for i, coh in enumerate(cohorts):
        sub = fb[fb["cohort"] == coh].sort_values("as_id")
        for k, (_, r) in enumerate(sub.iterrows()):
            x = i + (k - 0.5) * 0.26
            ax_a.bar(x, r["mean_delta"], width=0.24, color=colors[r["as_id"]],
                     alpha=0.88, edgecolor="white", lw=0.6)
            q = r["q"]
            star = "***" if q < 1e-4 else "**" if q < 1e-2 else "*" if q < 0.05 else "ns"
            ax_a.text(x, r["mean_delta"] + (0.006 if r["mean_delta"] > 0 else -0.006),
                      star, ha="center",
                      va="bottom" if r["mean_delta"] > 0 else "top",
                      fontsize=9.8, color="#333333")
        xpos.append(i)
        ticks.append(f"{COH[coh]}\n({COH_N[coh]})")
    ax_a.axhline(0, color="#333333", lw=0.9)
    ax_a.set_xticks(xpos)
    ax_a.set_xticklabels(ticks, fontsize=9.8)
    ax_a.set_ylabel("paired $\\Delta$PSI (tumour $-$ normal)")
    ax_a.set_ylim(-0.155, 0.178)   # headroom so the ** stars clear the note band
    for aid, lab in EV.items():
        ax_a.bar([np.nan], [np.nan], color=colors[aid], label=f"SpliceSeq {lab} (AT)")
    ax_a.legend(frameon=False, fontsize=8.4, loc="upper center", ncol=2,
                bbox_to_anchor=(0.5, -0.16))
    panel(ax_a, "A  FBLN1 terminal-exon usage is reproducible", fontsize=12.4,
                   pad=14)
    ax_a.text(0.02, 0.985,
              "the two events are competitive ends of one choice\n"
              "(mirror deltas), not two independent findings",
              transform=ax_a.transAxes, fontsize=8.4, color="#555555", va="top",
              linespacing=1.25)
    ax_a.text(0.985, 0.03, "* q<0.05   ** q<0.01   *** q<0.0001",
              transform=ax_a.transAxes, fontsize=8.7, color="#333333", ha="right", va="bottom")

    # ---------------------------------------------------------------- b
    # Composition axis is now the NNLS CAF fraction (148/151_*.py), not a
    # marker-gene PSI proxy. Signed rho, because the two AT events are mirrors
    # of one 3'-end choice -- the sign is the finding, so |rho| would erase it.
    comp = pd.read_csv(T / "fbln1_tcga_composition_by_cohort.csv")
    comp["as_id"] = comp["as_id"].astype(str)
    w = 0.34
    xs = np.arange(len(cohorts))
    for k, (aid, lab) in enumerate(EV.items()):
        vals, ns = [], []
        for coh in cohorts:
            r = comp[(comp["cohort"] == coh) & (comp["as_id"] == aid)].iloc[0]
            vals.append(r["rho_caf"])
            ns.append(int(r["n"]))
        off = (k - 0.5) * w
        ax_b.bar(xs + off, vals, width=w * 0.92, color=colors[aid], alpha=0.88,
                 edgecolor="white", lw=0.6,
                 label=f"SpliceSeq {lab} (AT)")
        for x, v in zip(xs + off, vals):
            ax_b.text(x, v + (0.018 if v > 0 else -0.018), f"{v:+.2f}",
                      ha="center", va="bottom" if v > 0 else "top",
                      fontsize=9.0, color="#333333")
    ax_b.axhline(0, color="#333333", lw=0.9)
    ax_b.set_xticks(xs)
    ax_b.set_xticklabels([COH[c] for c in cohorts], fontsize=9.8)
    ax_b.set_ylabel("Spearman $\\rho$ with CAF fraction")
    # headroom band above the tallest bar (0.52) so the three-line note has its
    # own space instead of sitting on the +0.52 / +0.50 value labels
    ax_b.set_ylim(-0.68, 0.95)
    ax_b.set_yticks([-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6])
    ax_b.legend(frameon=False, fontsize=8.4, loc="lower left", ncol=2,
                handlelength=1.3, columnspacing=1.0)
    panel(ax_b, "B  composition explains only a quarter of it", fontsize=12.4)
    r2lo = comp.groupby("as_id")["r2_caf"].min().min()
    r2hi = comp.groupby("as_id")["r2_caf"].max().max()
    ax_b.text(0.985, 0.995,
              f"$\\rho^2$ = {r2lo:.2f}\u2013{r2hi:.2f} at n = {min(ns)}\u2013{max(ns)} tumours\n"
              "and the paired change runs opposite\nto the stromal gradient (panel a)",
              transform=ax_b.transAxes, ha="right", va="top", fontsize=8.2,
              color="#555555", linespacing=1.28)

    # ---------------------------------------------------------------- c
    strict = tvn[tvn["cohort"] == "oral_strict"]
    n_ev = int(strict["caf_symbol"].nunique())
    n_pairs = len(strict)
    n_p05 = int((strict["p"] < 0.05).sum())
    n_fdr = int(strict["sig_fdr"].sum())
    exp_p05 = n_pairs * 0.05
    bars = [("event pairs tested", n_pairs, "#9aa3ad"),
            ("nominal p<0.05", n_p05, "#c47a2c"),
            ("FDR q<0.05", n_fdr, "#b8342b")]
    yb = np.arange(len(bars))[::-1]
    for y_, (lab, v, col) in zip(yb, bars):
        ax_c.barh(y_, v, color=col, height=0.6, edgecolor="white")
        ax_c.text(v + max(n_pairs * 0.015, 1), y_, f"{v}", va="center", fontsize=9.8)
        ax_c.text(-1.5, y_, lab, va="center", ha="right", fontsize=9.5)
    ax_c.axvline(exp_p05, color="#333333", ls=":", lw=1.1)
    ax_c.text(exp_p05 + max(n_pairs * 0.01, 1), yb[1] + 0.45,
              f"{exp_p05:.0f} expected by chance", fontsize=8.7, color="#333333")
    ax_c.set_yticks([])
    ax_c.set_xlim(0, n_pairs * 1.12)
    ax_c.set_xlabel("count (oral-strict cohort)")
    ax_c.set_ylim(-0.7, len(bars) - 0.3)
    panel(ax_c, f"C  {n_ev} CAF events: nominal signal, no robust hit", fontsize=12.4)
    for s in ("left", "right", "top"):
        ax_c.spines[s].set_visible(False)
    # The bars grow rightwards and are longest on top, so the only clean space is
    # the right side below the shortest bar.
    # The three bars fill the left half; the note goes into the empty right half at
    # the level of the shortest (FDR) bar, so it cannot overlap any bar.
    # This panel's three bars occupy only the left ~40% of the width and the
    # top third of the height is free, so the note goes there, right-aligned.
    ax_c.set_ylim(-0.7, len(bars) - 0.3 + 1.35)   # reserved band for the note
    ax_c.text(0.99, 0.995,
              f"{n_p05} of {n_pairs} pairs reach nominal p<0.05 (3 expected);\n"
              f"the {n_fdr} that survive FDR are MSRB3/COPE pairs whose\n"
              f"direction disagrees with the single-cell dPSI.",
              transform=ax_c.transAxes, ha="right", va="top", fontsize=7.3,
              color="#555555", linespacing=1.25)

    # ---------------------------------------------------------------- d
    # Aggregate to one row per CAF event (median over candidate pairs). Plotting all
    # 60 pairs made the gene axis unreadable and double-counted repeated genes.
    s = strict.dropna(subset=["caf_mean_dpsi", "mean_delta"]).copy()
    agg = s.groupby("caf_symbol").agg(
        caf=("caf_mean_dpsi", "first"),
        bulk=("mean_delta", "median"),
        k=("ss_as_id", "size")).reset_index()
    agg["conc"] = np.sign(agg["caf"]) == np.sign(agg["bulk"])
    agg = agg.sort_values("caf").reset_index(drop=True)
    xpos2 = np.arange(len(agg))
    ax_d.scatter(xpos2, agg["caf"], s=34, color="#7b4b8a", label="single-cell dPSI", zorder=4)
    ax_d.scatter(xpos2, agg["bulk"], s=34, marker="s",
                 c=["#1f6fb2" if c else "#b8342b" for c in agg["conc"]],
                 label="bulk dPSI (blue = same sign)", zorder=4)
    for xi, r in agg.iterrows():
        ax_d.plot([xi, xi], [r["caf"], r["bulk"]], color="#cccccc", lw=0.9, zorder=2)
    ax_d.axhline(0, color="#333333", lw=0.9)
    ax_d.set_xticks(xpos2)
    ax_d.set_xticklabels(agg["caf_symbol"], rotation=90, fontsize=8.7)
    ax_d.set_xlim(-0.8, len(agg) - 0.2)
    ax_d.set_ylabel("dPSI (median over candidate pairs for bulk)")
    n_conc = int(agg["conc"].sum())
    ax_d.legend(frameon=False, fontsize=8.7, loc="upper left")
    panel(ax_d, f"D  bulk direction agrees for {n_conc}/{len(agg)} events", fontsize=12.4)
    from scipy import stats as st
    bt = st.binomtest(n_conc, len(agg), 0.5, alternative="greater")
    ax_d.text(0.985, 0.05,
              f"binomial p = {bt.pvalue:.2f} (n = {len(agg)} events)\n"
              f"-> bulk neither confirms nor refutes",
              transform=ax_d.transAxes, ha="right", va="bottom", fontsize=8.7,
              color="#555555")

    # ================================================================ e
    # prognosis: per-SD hazard ratios for every testable CAF event, none survive.
    # One event (IDH3A) has HR = 1.1e5 with SE = 6.7e4 and p = 0.9999 -- numerical
    # separation, not an effect. It is excluded and the exclusion is stated on the
    # panel, so the axis shows the real spread rather than one non-estimable point.
    cox = pd.read_csv(T / "tcga_caf_events_cox.csv")
    cx_all = cox[cox["cohort"] == "oral_strict"].dropna(subset=["hr_per_sd", "p"])
    bad = cx_all[(cx_all["hr_per_sd"] > 100) | (cx_all["se"] > 10)]
    cx = cx_all.drop(bad.index)
    ax_e.scatter(cx["hr_per_sd"], -np.log10(cx["p"].clip(lower=1e-6)), s=40,
                 color="#5f7a86", edgecolor="white", lw=0.7, zorder=4)
    ax_e.axhline(-np.log10(0.05), color="#b8342b", ls=":", lw=1.1)
    # anchor the label INSIDE the axes just right of the HR = 1 line and just
    # above the dotted p = 0.05 line: the old anchor (min HR x 0.97, ha="right")
    # sat on the left spine and collided with the y tick label.
    ax_e.text(1.02, -np.log10(0.05) + 0.03, "p = 0.05", fontsize=7.6,
              color="#b8342b", va="bottom", ha="left")
    ax_e.axvline(1.0, color="#333333", lw=0.9)
    ax_e.set_xlabel("hazard ratio per SD of PSI")
    ax_e.set_ylabel("$-\\log_{10} p$ (Cox)")
    ax_e.set_xlim(min(0.75, float(cx["hr_per_sd"].min()) * 0.95),
                  max(1.35, float(cx["hr_per_sd"].max()) * 1.08))
    # -log10 p cannot be negative, so the axis keeps a zero floor; the summary
    # note is therefore trimmed to three lines and sits in the empty
    # lower-right corner. The one non-estimable event (IDH3A) is named in the
    # panel title count and in the legend, not in a 5-line in-panel block.
    ax_e.set_ylim(0, max(float(-np.log10(cx["p"].clip(lower=1e-6)).max()) * 1.12, 1.50))
    panel(ax_e, f"E  prognosis: none of {len(cx)} estimable events", fontsize=11.9)
    minq = float(cx["q"].min())
    note = (f"min p = {cx['p'].min():.3f}\nmin q = {minq:.2f}\n"
            f"n = {int(cx['n_used'].median())} cases")
    ax_e.text(0.985, 0.03, note, transform=ax_e.transAxes, ha="right", va="bottom",
              fontsize=8.2, color="#444444", linespacing=1.35)
    if len(bad):
        ax_e.text(0.985, 0.42,
                  f"+{len(bad)} not estimable\n({bad.iloc[0]['caf_symbol']})",
                  transform=ax_e.transAxes, ha="right", va="bottom",
                  fontsize=7.6, color="#8c8c8c", linespacing=1.3)

    # ================================================================ f
    # ICB response in an independent oral cohort: three scoring methods, all <= 0.5
    cv = pd.read_csv(T / "c_arm_scoring_variants.tsv", sep="\t")
    lp = pd.read_csv(T / "c_arm_literature_panels.tsv", sep="\t")
    names = list(cv.iloc[:, 0].astype(str)) + list(lp["panel"].astype(str))
    aucs = list(cv["auc"].astype(float)) + list(lp["auc"].astype(float))
    cols = ["#5f7a86"] * len(cv) + ["#7b4b8a"] * len(lp)
    ypos = np.arange(len(names))[::-1]
    ax_f.barh(ypos, aucs, color=cols, height=0.66, edgecolor="white", lw=0.6)
    for yy, v, nm in zip(ypos, aucs, names):
        # method names go in the left gutter, clear of the bars
        ax_f.text(-0.02, yy, nm if len(nm) <= 44 else nm[:43] + "\u2026",
                  va="center", ha="right", fontsize=7.4, color="#222222")
        ax_f.text(v + 0.012, yy, f"{v:.2f}", va="center", fontsize=7.9, color="#333333")
    ax_f.axvline(0.5, color="#b8342b", ls="--", lw=1.1)
    ax_f.set_yticks([])
    ax_f.set_xlim(0, 1.0)
    ax_f.set_ylim(-0.7, len(names) - 0.3)
    ax_f.set_xlabel("AUC for response (R vs NR)")
    # cohort size belongs in the title: the only free space inside the axes is
    # the top strip, and every bar starts at x = 0, so an in-panel note at the
    # bottom sits on the lowest bar (dark text on a filled patch).
    panel(ax_f, f"G  ICB response: no method exceeds chance (n = {N_ICB})",
          fontsize=11.9)
    ax_f.text(0.52, -0.60, "chance", fontsize=7.9, color="#b8342b", va="bottom",
              transform=ax_f.transData)

    # ================================================================ g
    # splicing factors vs response: 41 tested, 0 significant, direction coherent
    sf = pd.read_csv(T / "c_arm_splicing_factors_response.tsv", sep="\t")
    sf = sf.dropna(subset=["delta_logcpm", "p"]).sort_values("delta_logcpm")
    ax_g.scatter(sf["delta_logcpm"], -np.log10(sf["p"].clip(lower=1e-4)), s=26,
                 color="#4f8a5b", edgecolor="white", lw=0.6, zorder=4)
    ax_g.axhline(-np.log10(0.05), color="#b8342b", ls=":", lw=1.1)
    ax_g.axvline(0, color="#333333", lw=0.9)
    ax_g.set_xlabel("log2 CPM difference (R $-$ NR)")
    ax_g.set_ylabel("$-\\log_{10}$ nominal p")
    panel(ax_g, "H  splicing factors: no signal", fontsize=11.9)
    n_low = int((sf["delta_logcpm"] < 0).sum())
    # the tallest points on the left sit at ~1.06, so the note is pushed right
    # of x = 0.13 axes and given a band of headroom to sit in
    ax_g.set_ylim(0, max(-np.log10(sf["p"].min()) * 1.30, 1.45))
    ax_g.text(0.13, 0.97,
              f"{n_low}/{len(sf)} lower in responders\n"
              f"q < 0.05: {int((sf['q'] < 0.05).sum())}/{len(sf)}\n"
              f"min p = {sf['p'].min():.3f}",
              transform=ax_g.transAxes, va="top", fontsize=8.2, color="#444444",
              linespacing=1.35)

    # ================================================================ h
    # what the cohort could have detected at all: power of the survival analysis
    pw = pd.read_csv(T / "tcga_hnsc_oral_power.csv")
    pw = pw.dropna(subset=["min_detectable_HR"])
    xpos3 = np.arange(len(pw))
    ax_h.bar(xpos3, pw["min_detectable_HR"], color="#5f7a86", width=0.62,
             edgecolor="white", lw=0.6)
    for x, r in zip(xpos3, pw.itertuples()):
        ax_h.text(x, r.min_detectable_HR + 0.02, f"{r.min_detectable_HR:.2f}",
                  ha="center", fontsize=8.4, color="#333333")
    ax_h.axhline(1.0, color="#b8342b", ls="--", lw=1.1)
    ax_h.set_xticks(xpos3)
    ax_h.set_xticklabels([s.replace(", ", ",\n") for s in pw["scenario"]],
                         fontsize=6.4, rotation=90, ha="center",
                         linespacing=1.05)
    ax_h.set_ylabel("minimum detectable HR (80% power)")
    ax_h.set_ylim(1.0, max(pw["min_detectable_HR"].max() * 1.18, 1.9))
    panel(ax_h, "F  what this cohort could detect", fontsize=11.9)
    # the HR = 1 reference is explained in this in-panel note rather than by a
    # red label on the dashed line: every x position along that line carries a
    # bar, so on-line text inevitably sits on a filled patch.
    ax_h.text(0.03, 0.97,
              "only HR > ~1.3 was reachable,\nso a null Cox result is expected\n"
              "(dashed line = HR 1, no effect)",
              transform=ax_h.transAxes, va="top", fontsize=7.6, color="#444444",
              linespacing=1.25)

    # ================================================================ i
    # why the oral-strict cohort is the primary one: events concentrate there
    sb = pd.read_csv(T / "tcga_hnsc_site_breakdown.csv")
    comp = (sb.groupby("compartment", as_index=False)
              .agg(n_cases=("n_cases", "sum"),
                   n_with_os=("n_with_os", "sum"),
                   n_deaths=("n_deaths", "sum")))
    comp = comp.sort_values("n_deaths", ascending=True)
    ypos4 = np.arange(len(comp))
    ax_i.barh(ypos4, comp["n_deaths"], color="#5f7a86", height=0.62,
              edgecolor="white", lw=0.6, label="deaths")
    ax_i.barh(ypos4, comp["n_cases"] - comp["n_deaths"], left=comp["n_deaths"],
              color="#cfd8dd", height=0.62, edgecolor="white", lw=0.6,
              label="cases without a death recorded")
    for y, r in zip(ypos4, comp.itertuples()):
        ax_i.text(r.n_cases + 6, y, f"{r.n_cases} cases / {r.n_deaths} deaths",
                  va="center", fontsize=8.2, color="#333333")
    ax_i.set_yticks(ypos4)
    ax_i.set_yticklabels([c.replace(", ", ",\n") for c in comp["compartment"]],
                         fontsize=8.4)
    ax_i.set_xlim(0, comp["n_cases"].max() * 1.62)
    # headroom so the three-line note has its own band above the bars
    ax_i.set_ylim(-1.55, len(comp) - 0.15 + 0.15)
    ax_i.set_xlabel("TCGA-HNSC cases")
    ax_i.legend(frameon=False, fontsize=7.4, loc="lower right",
                bbox_to_anchor=(1.0, 0.22), ncol=1, handlelength=1.2)
    panel(ax_i, "I  events concentrate in the oral cavity", fontsize=11.9)
    top = comp.sort_values("n_deaths", ascending=False).iloc[0]
    # the note belongs in the empty lower-right quadrant: anchoring it under the
    # bars put its first line on top of the "other" row's value label
    ax_i.text(0.985, 0.03,
              f"{top['compartment']} holds {int(top['n_deaths'])} of "
              f"{int(comp['n_deaths'].sum())} deaths "
              f"({top['n_deaths'] / comp['n_deaths'].sum():.0%}),\n"
              f"which is why it is the primary cohort",
              transform=ax_i.transAxes, va="bottom", ha="right", fontsize=7.9,
              color="#444444", linespacing=1.3)

    fig.tight_layout(rect=(0, 0, 1, 1), w_pad=1.4, h_pad=1.8)
    # guard: a mis-placed annotation in data coordinates can make bbox_inches
    # expand the canvas enormously; see the 29,273 px incident in stage 23.
    fit_panel_titles(fig)

    # --- prove the render, check 1: collision scan ------------------------
    # Must run on the finalised in-memory figure: savefig(bbox_inches="tight")
    # temporarily rescales the canvas and a renderer taken afterwards reports
    # extents that no longer match, which makes every axis label appear to
    # collide with its own tick labels.
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    texts = [t for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible()]
    boxes = {id(t): t.get_window_extent(rend) for t in texts}
    own_ticks = {id(x) for ax in fig.axes
                 for x in ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both")}
    hits = []
    for i, a in enumerate(texts):
        for b in texts[i + 1:]:
            if id(a) in own_ticks and id(b) in own_ticks:
                continue                      # neighbouring tick labels may touch
            if boxes[id(a)].overlaps(boxes[id(b)]):
                hits.append((a.get_text(), b.get_text()))
    for ax in fig.axes:
        for s in ax.spines.values():
            if not s.get_visible():
                continue
            sb = s.get_window_extent(rend)
            own_labels = (ax.xaxis.label, ax.yaxis.label, ax.title)
            for t in texts:
                if id(t) in own_ticks or t in own_labels:
                    continue
                if boxes[id(t)].overlaps(sb):
                    hits.append((t.get_text(), "spine"))
    fig.savefig(FIG / "fig_fbln1_tcga.png", dpi=300, bbox_inches="tight")
    crops = panel_crops(fig)
    save_panel_crops(FIG / "fig_fbln1_tcga.png", crops)
    plt.close(fig)
    if hits:
        raise SystemExit("TEXT COLLISIONS:\n" + "\n".join(
            f"  {a!r} <-> {b!r}" for a, b in hits[:20]))
    print(f"collision scan clean; {len(texts)} text boxes; "
          f"{len(crops)} panel crops in .cache/figure-style/")
    from PIL import Image
    with Image.open(FIG / "fig_fbln1_tcga.png") as im:
        wpx, hpx = im.size
    assert wpx < 7000 and hpx < 7000, f"suspicious canvas {wpx}x{hpx}"

    print(f"panel c: pairs={n_pairs} nominal p<0.05={n_p05} (expected {exp_p05:.0f}) "
          f"FDR={n_fdr}")
    print(f"panel d: concordant {n_conc}/{len(agg)} events  binom p={bt.pvalue:.3f}")
    print(f"panel i: compartments={comp['compartment'].nunique()} "
          f"deaths={int(comp['n_deaths'].sum())}")
    print(f"written figures/fig_fbln1_tcga.png {wpx}x{hpx}")


if __name__ == "__main__":
    main()
