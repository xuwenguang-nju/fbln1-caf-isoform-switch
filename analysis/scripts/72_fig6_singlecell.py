"""Fig 8 (9 panels): the single-cell arm, rebuilt without subtype labels.

# PANEL-SCALE-1.32 (2026-10-22): fonts enlarged to match Fig 2.

a  27 exon-level events x 4 donors, all same sign                (strict panel)
b  observed vs expected counts of p<0.05: main contrast vs control
c  lineage composition of the four libraries
d  effect size vs between-donor variability for the 27 events
e  ECDF of all 20,396 per-event p-values against the uniform null
f  where the signal sits: events tested and events at p<0.05, by splice type
g  distribution of the mean dPSI by splice type
h  pairwise donor-to-donor Spearman correlation over all 20,396 events
i  why no subtype is named: the library-to-subtype map cannot be recovered

Panels c and i are the limits of the arm: the four libraries differ in what they
contain, and the published subtype labels cannot be re-derived from the data.
No subtype vocabulary appears anywhere as a result.

Outputs figures/fig_singlecell_locus.png
"""
from __future__ import annotations
import json
from pathlib import Path

import sys as _sys

_sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figstyle import panel, fit_panel_titles  # noqa: E402

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
REM = ROOT / "remote/easyai"
FIG = ROOT / "figures"
STYLE = Path(__file__).resolve().parent / "_style_runtime.py"
exec(compile(STYLE.read_text(encoding="utf-8"), "runtime.py", "exec"), globals())

TYPE_COLORS = {"SE": "#1f6fb2", "A5": "#c47a2c", "A3": "#7b4b8a",
               "RI": "#b8342b", "MX": "#4f8a5b"}
LIN_COLORS = {"CAF": "#1f6fb2", "epithelial": "#c47a2c", "immune": "#b8342b",
              "endothelial": "#4f8a5b", "unassigned": "#bbbbbb"}
DONORS = ["OCSCC1", "OCSCC13", "OCSCC3", "OCSCC36"]


def main() -> None:
    apply_figure_style(sizes=(11, 9.5, 8.5))
    ev = pd.read_csv(REM / "stringent_events.tsv", sep="\t")
    G = json.loads((REM / "global_tests.json").read_text(encoding="utf-8"))
    nd = json.loads((REM / "nd_donor_summary.json").read_text(encoding="utf-8"))
    allv = pd.read_csv(REM / "nd_donor_events.tsv", sep="\t")
    q0 = pd.read_csv(REM / "caf_q0_assignment.csv")

    assert len(ev) == 27, f"expected 27 events, got {len(ev)}"
    assert (ev["same_sign"] == 4).all()
    assert (ev["q"] > 0.05).all(), "no per-event significance at n=4"
    ev["vals_list"] = ev["vals"].str.split(";").apply(lambda v: [float(x) for x in v])
    ev = ev.sort_values("mean").reset_index(drop=True)

    # 4 donor columns, used by panels e-h
    M = np.array([[float(x) for x in str(v).split(";")] for v in allv["vals"]])
    assert M.shape == (len(allv), 4) and not np.isnan(M).any()

    fig, axes = plt.subplots(3, 3, figsize=(16.72, 15.44))
    ax_e, ax_b, ax_a = axes[0]
    ax_h, ax_i, ax_d = axes[1]
    ax_f, ax_g, ax_c = axes[2]
    # ---------------------------------------------------------------- a
    for i, r in ev.iterrows():
        v = np.array(r["vals_list"])
        col = TYPE_COLORS.get(r["type"], "#666666")
        ax_a.scatter([i] * len(v), v, s=20, color=col, alpha=0.55, linewidths=0, zorder=3)
        ax_a.plot([i - 0.34, i + 0.34], [r["mean"]] * 2, color="#222222", lw=1.5, zorder=4)
    ax_a.axhline(0, color="#333333", lw=0.9)
    ax_a.set_xticks(range(len(ev)))
    ax_a.set_xticklabels(ev["symbol"], rotation=90, fontsize=8.2)
    ax_a.set_xlim(-0.8, len(ev) - 0.2)
    ax_a.set_ylabel("dPSI (CAF $-$ epithelial), per donor")
    panel(ax_a, "C  27 exon-level events, 4 donors, all same sign", fontsize=12.4)
    for t, c in TYPE_COLORS.items():
        ax_a.scatter([], [], color=c, s=20, label=f"{t} ({int((ev['type'] == t).sum())})")
    ax_a.legend(frameon=False, fontsize=8.2, ncol=1, loc="lower right",
                columnspacing=0.8, handletextpad=0.3, labelspacing=0.3)
    ax_a.text(0.015, 0.985,
              "horizontal bar = mean of the four donors\n"
              "no per-event FDR is quoted: at n=4 the BH threshold over ~20k tests\n"
              "needs |Cohen's d| > 35, which no event reaches (all q > 0.23)",
              transform=ax_a.transAxes, fontsize=8.2, color="#444444", va="top")

    # ---------------------------------------------------------------- b
    order = [("CAF_vs_epithelial", "CAF vs epithelial\n(main)", "#1f6fb2"),
             ("CAF_vs_endothelial", "CAF vs endothelial\n(control)", "#9aa3ad")]
    xs = np.arange(len(order))
    obs = [G[k]["n_p05"] for k, _, _ in order]
    exp = [G[k]["exp_p05"] for k, _, _ in order]
    w = 0.34
    ax_b.bar(xs - w / 2, obs, width=w, color=[c for _, _, c in order], alpha=0.9,
             edgecolor="white", label="observed")
    ax_b.bar(xs + w / 2, exp, width=w, color="#dddddd", edgecolor="#999999", hatch="//",
             lw=0.7, label="expected under null")
    for x, o, e in zip(xs, obs, exp):
        ax_b.text(x - w / 2, o + 45, f"{o:,}", ha="center", fontsize=9.2)
        ax_b.text(x + w / 2, e + 45, f"{e:.0f}", ha="center", fontsize=9.2, color="#555555")
        fold = o / e
        ax_b.text(x, max(o, e) + 250, f"{fold:.2f}$\\times$", ha="center", fontsize=10.6,
                  fontweight="bold", color="#b8342b" if fold > 1.5 else "#666666")
    ax_b.set_xticks(xs)
    ax_b.set_xticklabels([lab for _, lab, _ in order], fontsize=9.8)
    ax_b.set_ylabel("events with p < 0.05")
    ax_b.set_ylim(0, max(obs + exp) * 1.42)
    ax_b.legend(frameon=False, fontsize=9.0, loc="upper left")
    panel(ax_b, "B  the signal is contrast-specific", fontsize=12.4)
    ax_b.text(0.985, 0.97,
              f"n = {G['CAF_vs_epithelial']['n']:,} / {G['CAF_vs_endothelial']['n']:,} events\n"
              f"four-donor same-sign:\n"
              f"  main  {G['CAF_vs_epithelial']['all_same_sign']:,} vs "
              f"{G['CAF_vs_epithelial']['exp_allsame']:.0f}\n"
              f"  ctrl  {G['CAF_vs_endothelial']['all_same_sign']:,} vs "
              f"{G['CAF_vs_endothelial']['exp_allsame']:.0f}\n"
              f"same-sign alone is not enough (the control\n"
              f"also exceeds it); the p distribution is the test",
              transform=ax_b.transAxes, ha="right", va="top", fontsize=7.9,
              color="#444444", linespacing=1.3)

    # ---------------------------------------------------------------- c
    donors = nd["donors"]
    lin = {}
    for d in donors:
        p = REM / f"donor_{d}_summary.json"
        if not p.exists():
            raise FileNotFoundError(f"missing {p.name}; copy donor summaries into remote/easyai")
        lin[d] = json.loads(p.read_text(encoding="utf-8"))["lineage"]
    cats = ["CAF", "epithelial", "immune", "endothelial", "unassigned"]
    bottom = np.zeros(len(donors))
    xd = np.arange(len(donors))
    for cat in cats:
        vals = np.array([lin[d].get(cat, 0) for d in donors], dtype=float)
        ax_c.bar(xd, vals, bottom=bottom, color=LIN_COLORS[cat], label=cat,
                 edgecolor="white", lw=0.6)
        bottom += vals
    for x, d in zip(xd, donors):
        tot = sum(lin[d].values())
        caf = lin[d].get("CAF", 0)
        ax_c.text(x, tot + 90, f"CAF {caf / tot:.0%}", ha="center", fontsize=9.0,
                  color="#1f6fb2")
    ax_c.set_xticks(xd)
    ax_c.set_xticklabels([d.replace("_Post", "") for d in donors], fontsize=9.5, rotation=15)
    ax_c.set_ylabel("cells")
    ax_c.set_ylim(0, bottom.max() * 1.62)
    ax_c.legend(frameon=False, fontsize=8.7, ncol=5, loc="upper center",
                bbox_to_anchor=(0.5, 1.0), columnspacing=0.8, handletextpad=0.4)
    panel(ax_c, "I  the four libraries differ in composition", fontsize=12.4)
    # keep the caveat inside the axes: a negative y in axes coords lands in the
    # gap between rows and reads as if it belonged to no panel
    # sit the note in the gap between the OCSCC13 and OCSCC3 stacks: the CAF%
    # labels live just above the bars, so the text must stay below them.
    ax_c.text(0.52, 0.90,
              "all four are labelled CD45neg or Unsorted,\n"
              "yet immune content ranges 0.1-53%\n"
              "-> CAF depth differs ~7$\\times$",
              transform=ax_c.transAxes, fontsize=7.4, color="#444444", va="top",
              ha="left", linespacing=1.3,
              bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.82))

    # ---------------------------------------------------------------- d
    ax_d.scatter(ev["mean"].abs(), ev["sd"], s=46, color="#7b4b8a", alpha=0.8,
                 edgecolor="white", lw=0.7, zorder=4)
    lim = max(ev["mean"].abs().max(), ev["sd"].max()) * 1.12
    ax_d.plot([0, lim], [0, lim], ls="--", lw=1.0, color="#333333")
    n_lab = 0                      # counts annotated labels only, not rows
    for _, r in ev.iterrows():
        if r["d_cohen"] is not None and abs(r["d_cohen"]) > 15:
            # alternate the offset side so neighbouring points do not stack
            # their labels on top of each other
            off = (-11, 10) if n_lab == 0 else (11, -17)
            n_lab += 1
            ax_d.annotate(r["symbol"], (abs(r["mean"]), r["sd"]),
                          xytext=off, textcoords="offset points", fontsize=7.9,
                          ha="right" if n_lab == 1 else "left",
                          color="#5c3669")
    ax_d.set_xlim(0, lim)
    ax_d.set_ylim(0, lim)
    ax_d.set_xlabel("|mean dPSI| across donors")
    ax_d.set_ylabel("SD across donors")
    panel(ax_d, "F  most events are far more consistent than their size", fontsize=12.4)
    ax_d.text(0.985, 0.30,
              f"median SD = {ev['sd'].median():.3f}; {int((ev['sd'] < 0.05).sum())}/27 events "
              f"have SD < 0.05\nCohen's d ranges {ev['d_cohen'].min():.1f} to "
              f"{ev['d_cohen'].max():.1f}",
              transform=ax_d.transAxes, ha="right", va="bottom", fontsize=8.4,
              color="#555555", linespacing=1.3)

    # ---------------------------------------------------------------- e
    # the whole testable universe, not only the p<0.05 count of panel b
    pv = np.sort(allv["p"].to_numpy(dtype=float))
    n_ev = len(pv)
    ecdf = np.arange(1, n_ev + 1) / n_ev
    ax_e.plot([0, 1], [0, 1], ls="--", lw=1.1, color="#333333",
              label="uniform (no signal)")
    ax_e.plot(pv, ecdf, lw=1.8, color="#1f6fb2", label=f"observed (n = {n_ev:,})")
    ax_e.set_xlim(0, 0.5)
    ax_e.set_ylim(0, 1.0)
    ax_e.set_xlabel("per-event p (one-sample t across 4 donors)")
    ax_e.set_ylabel("cumulative fraction of events")
    ax_e.legend(frameon=False, fontsize=8.7, loc="lower right")
    panel(ax_e, "A  the p distribution is shifted below uniform", fontsize=12.4)
    n05 = int((allv["p"] < 0.05).sum())
    n05e = 0.05 * n_ev
    # no arrow: at this scale it crossed both curves. The text sits in the empty
    # wedge below the diagonal, where nothing else is drawn.
    ax_e.text(0.985, 0.30,
              f"p < 0.05: {n05:,} of {n_ev:,}\n({n05 / n05e:.1f}$\\times$ the 5% expected)",
              transform=ax_e.transAxes, ha="right", va="top",
              fontsize=8.4, color="#b8342b", linespacing=1.3)

    # ---------------------------------------------------------------- f
    types = list(TYPE_COLORS)
    cnt = [int((allv["type"] == t).sum()) for t in types]
    sig = [int(((allv["type"] == t) & (allv["p"] < 0.05)).sum()) for t in types]
    xtf = np.arange(len(types))
    # Outline = events tested, filled = events at p<0.05. Using a faint fill for
    # the tested bar made it invisible on a log axis, so it is drawn as an open bar.
    for x, c_, s_, t in zip(xtf, cnt, sig, types):
        ax_f.bar(x, c_, width=0.62, facecolor="none", edgecolor=TYPE_COLORS[t],
                 lw=1.3, label="events tested" if x == 0 else None)
        ax_f.bar(x, s_, width=0.62, facecolor=TYPE_COLORS[t],
                 edgecolor=TYPE_COLORS[t], lw=0.8,
                 label="events at p < 0.05" if x == 0 else None)
        ax_f.text(x, s_ * 1.25, f"{s_:,}", ha="center", fontsize=8.7,
                  color=TYPE_COLORS[t], fontweight="bold")
        ax_f.text(x, c_ * 1.30, f"{c_:,}", ha="center", fontsize=7.9, color="#555555")
    ax_f.set_xticks(xtf)
    ax_f.set_xticklabels(types, fontsize=9.8)
    ax_f.set_yscale("log")
    ax_f.set_ylim(20, max(cnt) * 9.0)
    ax_f.set_ylabel("events (log scale)")
    ax_f.legend(frameon=False, fontsize=8.2, loc="upper center", ncol=2,
                bbox_to_anchor=(0.5, 1.0), columnspacing=1.0, handlelength=1.1)
    panel(ax_f, "G  nominal signal in every splice type", fontsize=12.4)
    ax_f.text(0.985, 0.735,
              "open = tested\nfilled = p<0.05",
              transform=ax_f.transAxes, fontsize=7.2, color="#444444", va="top",
              ha="right", linespacing=1.25)

    # ---------------------------------------------------------------- g
    data_g = [allv.loc[allv["type"] == t, "mean"].to_numpy() for t in types]
    bp = ax_g.boxplot(data_g, orientation="horizontal", widths=0.6, patch_artist=True,
                      showfliers=False, medianprops=dict(color="#222222", lw=1.2),
                      whiskerprops=dict(lw=0.9), capprops=dict(lw=0.9))
    for patch, t in zip(bp["boxes"], types):
        patch.set_facecolor(TYPE_COLORS[t])
        patch.set_alpha(0.55)
        patch.set_edgecolor("white")
    ax_g.axvline(0, color="#333333", lw=0.9)
    ax_g.set_yticks(np.arange(1, len(types) + 1))
    ax_g.set_yticklabels(types, fontsize=9.8)
    ax_g.set_xlabel("mean dPSI (CAF $-$ epithelial)")
    ax_g.set_xlim(-0.45, 0.45)
    panel(ax_g, "H  effect sizes are small in every type", fontsize=12.4)
    med_txt = ", ".join(f"{t} {np.median(v):+.3f}" for t, v in zip(types, data_g))
    ax_g.text(0.5, 0.02, "median dPSI:  " + med_txt,
              transform=ax_g.transAxes, fontsize=7.7, color="#444444", va="bottom",
              ha="center")

    # ---------------------------------------------------------------- h
    R = np.identity(4)
    for i in range(4):
        for j in range(4):
            R[i, j] = stats.spearmanr(M[:, i], M[:, j]).correlation
    im = ax_h.imshow(R, cmap="Blues", vmin=0, vmax=1.0)
    for i in range(4):
        for j in range(4):
            ax_h.text(j, i, f"{R[i, j]:.2f}", ha="center", va="center", fontsize=9.2,
                      color="white" if R[i, j] > 0.7 else "#222222")
    ax_h.set_xticks(range(4))
    ax_h.set_yticks(range(4))
    ax_h.set_xticklabels(DONORS, fontsize=8.7, rotation=30, ha="right")
    ax_h.set_yticklabels(DONORS, fontsize=8.7)
    off = R[~np.eye(4, dtype=bool)]
    panel(ax_h, f"D  donors agree on the same events (rho {off.min():.2f}-{off.max():.2f})",
          fontsize=11.6)

    # ---------------------------------------------------------------- i
    q = q0.sort_values("rho", ascending=True).reset_index(drop=True)
    agree = (q["subtype_inferred"] == q["subtype_hypothesis"]).to_numpy()
    yq = np.arange(len(q))
    ax_i.barh(yq, q["rho"], height=0.66, edgecolor="white", lw=0.6,
              color=["#4f8a5b" if a else "#b8342b" for a in agree])
    ax_i.axvline(0.90, color="#333333", ls="--", lw=1.1)
    ax_i.set_yticks(yq)
    ax_i.set_yticklabels([f"{s} $\\rightarrow$ {t}"
                          for s, t in zip(q["subtype_hypothesis"], q["subtype_inferred"])],
                         fontsize=7.4)
    ax_i.set_xlim(0, 1.02)
    ax_i.set_ylim(-0.9, len(q) - 0.1 + 1.9)      # headroom for the note
    ax_i.set_xlabel("Spearman $\\rho$ to the best-matching published profile")
    panel(ax_i, "E  the subtype map cannot be recovered", fontsize=12.4)
    # the three longest bars occupy the top rows, so the note goes to the lower
    # right, beside the shortest bar, and the criterion label sits on the dashed
    # line itself.
    ax_i.text(0.985, 0.995,
              f"pre-specified criterion $\\rho \\geq$ 0.90 (dashed): "
              f"none of {len(q)} libraries reaches it\n"
              f"reported subtype recovered in {int(agree.sum())}/{len(q)}",
              transform=ax_i.transAxes, ha="right", va="top", fontsize=7.2,
              color="#444444", linespacing=1.25)

    fig.tight_layout(rect=(0, 0, 1, 1), w_pad=1.4, h_pad=2.2)
    fit_panel_titles(fig)
    fig.savefig(FIG / "fig_singlecell_locus.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    from PIL import Image
    with Image.open(FIG / "fig_singlecell_locus.png") as im:
        wpx, hpx = im.size
    assert wpx < 7000 and hpx < 7000, f"suspicious canvas {wpx}x{hpx}"

    print("lineage composition:")
    print(pd.DataFrame(lin).T.fillna(0).astype(int).to_string())
    print(f"\nevent SD median={ev['sd'].median():.3f}  "
          f"n(SD<0.05)={int((ev['sd'] < 0.05).sum())}/27")
    print(f"panel e: p<0.05 {n05:,}/{n_ev:,} = {n05 / n05e:.1f}x expected")
    print(f"panel h: off-diagonal rho {off.min():.3f}-{off.max():.3f}")
    print(f"panel i: rho {q['rho'].min():.3f}-{q['rho'].max():.3f}, "
          f"subtype recovered {int(agree.sum())}/{len(q)}")
    print(f"written figures/fig_singlecell_locus.png {wpx}x{hpx}")


if __name__ == "__main__":
    main()
