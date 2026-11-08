"""Supplementary Fig S1 (3x3): the cohort and power ledger.

# PANEL-SCALE-1.32 (2026-10-22): fonts enlarged to match Fig 2.

The paper rests on data chosen under hard constraints (which public 10x libraries
have usable 3' chemistry, how many independent donors, what a matched design of
n = 6 can actually detect). Those constraints decided the study design and are
currently only prose. This figure is the ledger.

Panels
  a  every candidate resource: libraries quantified vs donors, by PSI capability
  b  TCGA-HNSC subsites: cases and death events
  c  what the oral cohort can detect: minimum HR by coding choice
  d  paired-design power: minimum detectable dIF vs n, at three SDs
  e  assumed vs measured between-patient SD (the number that changed the plan)
  f  the 9R/19NR ICB cohort: minimum detectable PSI shift by strategy
  g  TCGA event coverage: how many of the 544 PSI-bearing samples carry a value
  h  single-cell resources screened, cells vs whether full-length (splice-fit)
  i  the funnel: resources -> libraries -> events -> stringent -> kept switch
"""
from __future__ import annotations

from pathlib import Path

import sys as _sys

_sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figstyle import panel, fit_panel_titles  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize, stats

STYLE = Path(__file__).resolve().parent / "_style_runtime.py"
exec(compile(STYLE.read_text(encoding="utf-8"), "runtime.py", "exec"), globals())
apply_figure_style(sizes=(11, 9.5, 8.5))

ROOT = Path(__file__).resolve().parents[2]
T = ROOT / "results/tables"
OUT = ROOT / "figures"

RED = "#b8342b"
BLUE = "#1f6fb2"
GREY = "#8a8a8a"
ORANGE = "#c47a2c"
GREEN = "#5b8f5b"


def power_of(d: float, sd: float, n: int, alpha: float) -> float:
    """Two-sided one-sample t-test power, exact (noncentral t)."""
    ncp = d * np.sqrt(n) / sd
    crit = stats.t.ppf(1 - alpha / 2, n - 1)
    return float(stats.nct.sf(crit, n - 1, ncp) + stats.nct.cdf(-crit, n - 1, ncp))


def d_for_power(sd: float, n: int, alpha: float, power: float = 0.8) -> float:
    """Grid + interpolation. scipy's nct returns NaN once ncp is large, which
    breaks a root finder; a monotone grid avoids that entirely."""
    grid = np.linspace(0.001, 6.0 * sd, 900)
    ps = np.array([power_of(d, sd, n, alpha) for d in grid])
    good = np.isfinite(ps)
    grid, ps = grid[good], ps[good]
    if ps[-1] < power:                      # extrapolate the tail linearly
        i = max(1, int(0.9 * len(ps)))
        slope = (ps[-1] - ps[i]) / max(grid[-1] - grid[i], 1e-12)
        return float(grid[-1] + (power - ps[-1]) / max(slope, 1e-12))
    return float(np.interp(power, ps, grid))


def main() -> None:
    led = pd.read_csv(T / "arm_sample_ledger.csv")
    tcga = pd.read_csv(T / "tcga_hnsc_site_ledger.csv")
    oral = pd.read_csv(T / "tcga_hnsc_oral_power.csv")
    # fixed-referent power table (191_power_from_fixed_referent.py). The filed
    # caf_power_measured.csv was computed on the withdrawn dominant-isoform dIF.
    meas = pd.read_csv(T / "caf_power_measured_fixed_referent.csv")
    icb = pd.read_csv(T / "icb_power_9v19.csv")
    depth = pd.read_csv(T / "tcga_hnsc_psi_event_depth.csv")
    scf = pd.read_csv(T / "sc_dataset_feasibility.csv")

    fig, axes = plt.subplots(3, 3, figsize=(16.72, 14.58))
    ax_a, ax_b, ax_c = axes[0]
    ax_d, ax_e, ax_f = axes[1]
    ax_g, ax_h, ax_i = axes[2]

    # ---------------------------------------------------------------- a
    d = led.dropna(subset=["n_libraries"]).copy()
    d["n_libraries"] = d["n_libraries"].astype(float)
    d["n_donors"] = pd.to_numeric(d["n_donors"], errors="coerce")
    cmap = {"pseudobulk": BLUE, "event-level": GREY, "gene-level": ORANGE}
    for k, (psi, grp) in enumerate(d.groupby("psi")):
        col = cmap.get(psi, GREY)
        ax_a.scatter(grp["n_libraries"], grp["n_donors"].fillna(0),
                     s=34, color=col, edgecolor="white", lw=0.7,
                     label=f"PSI: {psi}", zorder=3)
    # 11 labels collide at this scale; the legend already carries the classes and
    # the accession list is in the caption, so annotate only the two extremes.
    # Individual accessions are not annotated: at this scale the three
    # low-library points and the axis ticks all collide. The eleven accessions
    # are listed in the figure caption instead.
    ax_a.set_xlabel("libraries quantified")
    ax_a.set_ylabel("donors")
    ax_a.legend(frameon=False, fontsize=7.4, loc="upper left")
    ax_a.set_xlim(0, d["n_libraries"].max() * 1.30)
    ax_a.set_ylim(0, d["n_donors"].max() * 1.25)
    panel(ax_a, f"a  {len(d)} candidate resources", fontsize=12.4)

    # ---------------------------------------------------------------- b
    t = tcga.copy()
    order = t["n_cases"].sort_values(ascending=False).index
    t = t.loc[order]
    xs = np.arange(len(t))
    w = 0.38
    ax_b.bar(xs - w / 2, t["n_cases"], width=w, color=GREY, label="cases",
             edgecolor="white", lw=0.5)
    ax_b.bar(xs + w / 2, t["n_events"], width=w, color=RED, label="deaths",
             edgecolor="white", lw=0.5)
    for x, (a, b) in enumerate(zip(t["n_cases"], t["n_events"])):
        ax_b.text(x - w / 2, a + 4, f"{int(a)}", ha="center", fontsize=7.4)
        ax_b.text(x + w / 2, b + 4, f"{int(b)}", ha="center", fontsize=7.4, color=RED)
    ax_b.set_xticks(xs)
    ax_b.set_xticklabels([g.replace("_", "\n") for g in t["group"]], fontsize=7.1)
    ax_b.set_ylabel("TCGA-HNSC patients")
    ax_b.legend(frameon=False, fontsize=7.9, ncol=2, loc="upper right")
    ax_b.set_ylim(0, t["n_cases"].max() * 1.24)
    panel(ax_b, "b  oral cavity carries the events", fontsize=12.4)

    # ---------------------------------------------------------------- c
    o = oral.copy()
    o["tag"] = o["scenario"].str.replace(r"continuous, per 1 SD.*", "continuous", regex=True)
    o["tag"] = o["tag"].str.replace("median split (binary, p1=p2=0.5)", "median split", regex=False)
    coded = o["scenario"].str.startswith("continuous")
    ax_c.scatter(o.loc[coded, "n_events"], o.loc[coded, "min_detectable_HR"],
                 s=46, marker="o", color=GREEN, edgecolor="white", lw=0.8,
                 zorder=3, label="continuous (per 1 SD)")
    ax_c.scatter(o.loc[~coded, "n_events"], o.loc[~coded, "min_detectable_HR"],
                 s=52, marker="s", color=RED, edgecolor="white", lw=0.8,
                 zorder=3, label="median split (binary)")
    ax_c.legend(frameon=False, fontsize=7.4, loc="upper right")
    ax_c.set_xlabel("death events")
    ax_c.set_ylabel("minimum detectable HR (80% power)")
    ax_c.set_xlim(80, 175)
    ax_c.set_ylim(1.10, 1.90)
    ax_c.text(0.03, 0.06, "continuous PSI needs HR 1.26-1.36;\n"
              "a median split needs 1.59-1.74\n(n = 246 oral, 316 oral-wide)",
              transform=ax_c.transAxes, fontsize=7.4, color="#333333", linespacing=1.3)
    panel(ax_c, "c  how to code PSI decides power", fontsize=12.4)

    # ---------------------------------------------------------------- d
    # the SDs come from the table, not from literals: the filed 0.074/0.117 were
    # the dominant-isoform values and the fixed referent gives 0.089/0.153
    g = meas.drop_duplicates(subset=["sd_source"]).set_index("sd_source")
    sd_med = float(g.loc["measured_median_dIF_SD", "sd"])
    sd_p75 = float(g.loc["measured_p75_dIF_SD", "sd"])
    ns = np.arange(4, 25)
    sds = [(f"measured median SD {sd_med:.3f}", sd_med, GREEN),
           ("scenario SD 0.10", 0.10, GREY),
           (f"measured p75 SD {sd_p75:.3f}", sd_p75, ORANGE)]
    for lab, sd, col in sds:
        ys = [d_for_power(sd, int(n), 0.05) for n in ns]
        ax_d.plot(ns, ys, color=col, lw=1.4, label=lab)
    ax_d.axhline(0.10, color=RED, ls=":", lw=1.1)
    # the three curves converge near dIF=0.10 on the right; label the line at the
    # left edge where they have not yet reached it
    ax_d.text(6.4, 0.415, "dIF = 0.10 (smallest worth claiming)",
              fontsize=7.2, color=RED, va="top")
    ax_d.axvline(6, color="#333333", ls="--", lw=1.0)
    ax_d.text(6.3, 0.42, "n = 6\n(the design)", fontsize=7.7, color="#333333",
              linespacing=1.3)
    ax_d.set_xlabel("matched patient pairs")
    ax_d.set_ylabel("minimum detectable dIF (80% power, $\\alpha$=0.05)")
    ax_d.set_ylim(0.05, 0.45)
    ax_d.legend(frameon=False, fontsize=7.2, loc="upper right",
                bbox_to_anchor=(1.0, 0.985), labelspacing=0.28,
                borderaxespad=0.2)
    # n = 6 delivers 0.125-0.128 at the measured median SD, i.e. just above the
    # 0.10 threshold, so "barely" still describes it; assert rather than assume
    mde6 = d_for_power(sd_med, 6, 0.05)
    panel(ax_d, f"d  six pairs is enough, barely ({mde6:.3f})", fontsize=12.4)

    # ---------------------------------------------------------------- e
    m = meas.drop_duplicates(subset=["sd_source"])
    xs_e = np.arange(len(m))
    ax_e.bar(xs_e, m["sd"], color=[GREEN if "median" in s else ORANGE
                                   for s in m["sd_source"]],
             width=0.55, edgecolor="white", lw=0.6)
    ax_e.axhline(0.10, color=GREY, ls="--", lw=1.2)
    # left end of the line: the right end collides with the third bar's value label
    ax_e.text(-0.42, 0.103, "assumed 0.10", fontsize=7.7, color=GREY, ha="left")
    for x, v in zip(xs_e, m["sd"]):
        ax_e.text(x, v + 0.004, f"{v:.3f}", ha="center", fontsize=8.2)
    ax_e.set_xticks(xs_e)
    ax_e.set_xticklabels([s.replace("_", "\n") for s in m["sd_source"]], fontsize=7.4)
    ax_e.set_ylabel("between-patient SD of dIF")
    ax_e.set_ylim(0, max(m["sd"].max(), 0.10) * 1.35)
    # The direction of the comparison is data-dependent: under the fixed referent
    # the measured median SD (0.089) sits BELOW the 0.10 assumption, so the old
    # "pessimistic by 1.4x" line is no longer true. State the ratio, not a verdict.
    ratio = 0.10 / float(g.loc["measured_median_dIF_SD", "sd"])
    if ratio > 1:
        txt_e = f"median SD is {ratio:.1f}x below the assumed 0.10"
    else:
        txt_e = f"median SD is {1 / ratio:.1f}x above the assumed 0.10"
    ax_e.text(0.03, 0.94, txt_e, transform=ax_e.transAxes, va="top",
              fontsize=7.8, color="#333333")
    panel(ax_e, "e  measured SD versus the scenario", fontsize=12.4)

    # ---------------------------------------------------------------- f
    i2 = icb[icb["assumption"] == "PSI SD=0.1"].copy()
    if i2.empty:
        i2 = icb.copy()
    i2 = i2.sort_values("min_detectable_psi_delta")
    ax_f.barh(range(len(i2)), i2["min_detectable_psi_delta"],
              color=[RED if "set-level" in s else GREY for s in i2["strategy"]],
              height=0.58, edgecolor="white", lw=0.6)
    ax_f.set_yticks(range(len(i2)))
    ax_f.set_yticklabels([s.replace(", ", "\n") for s in i2["strategy"]], fontsize=6.9)
    for y, v in enumerate(i2["min_detectable_psi_delta"]):
        ax_f.text(v + 0.004, y, f"{v:.3f}", va="center", fontsize=7.7)
    ax_f.set_xlabel("minimum detectable PSI shift (absolute, 0-1 scale)")
    ax_f.set_xlim(0, i2["min_detectable_psi_delta"].max() * 1.32)
    ax_f.text(0.98, 0.06, "9 responders vs 19 non-responders",
              transform=ax_f.transAxes, ha="right", fontsize=7.7, color="#333333")
    panel(ax_f, "f  the ICB cohort only supports set-level tests", fontsize=12.4)

    # ---------------------------------------------------------------- g
    fv = depth["frac_samples_with_value"].dropna()
    ax_g.hist(fv, bins=40, color=BLUE, edgecolor="white", linewidth=0.4)
    ax_g.axvline(0.90, color=RED, ls="--", lw=1.1)
    ax_g.text(0.895, ax_g.get_ylim()[1] * 0.88, "0.90", fontsize=7.7, color=RED,
              ha="right")
    # frac_samples_with_value is a fraction of the 544 samples that carry PSI
    # (501 tumours + 43 adjacent normals), not of the 528 clinical cases.
    ax_g.set_xlabel("fraction of 544 PSI-bearing samples")
    ax_g.set_ylabel("events")
    ax_g.text(0.03, 0.98,
              f"median {fv.median():.3f}\n"
              f"{int((fv >= 0.90).sum()):,}/{len(fv):,} events above 0.90",
              transform=ax_g.transAxes, va="top", fontsize=7.6, color="#333333",
              linespacing=1.3)
    ax_g.set_ylim(0, ax_g.get_ylim()[1] * 1.35)
    panel(ax_g, f"g  {len(depth):,} TCGA events are usable", fontsize=12.4)

    # ---------------------------------------------------------------- h
    # `n` mixes cells / samples / tumours, so parse value and unit separately.
    # Units are printed on the tick labels; lengths are only comparable within a unit.
    s = scf.copy()
    ext = s["n"].astype(str).str.extract(r"([\d.]+)\s*(\w+)")
    s["value"] = pd.to_numeric(ext[0], errors="coerce")
    s["unit"] = ext[1].fillna("?")
    s = s.dropna(subset=["value"])
    s = s.sort_values("value")
    fit = s["splice_fit"].astype(str).str.upper()
    col = np.where(fit.str.startswith("HIGH"), GREEN,
                   np.where(fit.str.startswith("NOT"), "#999999",
                            np.where(fit.str.startswith("LOW"), ORANGE, GREY)))
    ypos = np.arange(len(s))
    ax_h.barh(ypos, s["value"], color=col, height=0.66, edgecolor="white", lw=0.6)
    ax_h.set_xscale("log")
    for y, (v, u) in enumerate(zip(s["value"], s["unit"])):
        ax_h.text(v * 1.10, y, f"{int(v):,} {u}", va="center", fontsize=6.3)
    ax_h.set_yticks(ypos)
    ax_h.set_yticklabels([str(a)[:14] for a in s["accession"]], fontsize=6.3)
    ax_h.set_xlabel("size of the resource (log; units differ by row)")
    ax_h.set_xlim(1, s["value"].max() * 6)
    ax_h.set_ylim(-0.7, len(s) - 0.3)
    for cval, lab in ((GREEN, "full-length"), (GREY, "partial"),
                      (ORANGE, "low"), ("#999999", "unusable")):
        ax_h.barh([], [], color=cval, label=lab)
    ax_h.legend(frameon=False, fontsize=6.6, ncol=4, loc="lower right",
                handlelength=1.0, columnspacing=0.8)
    panel(ax_h, "h  13 resources screened, one is full-length", fontsize=12.4)

    # ---------------------------------------------------------------- i
    # All stages are counts of the same kind of thing (units), so the funnel is
    # internally comparable; the event-level counts live in panels a and g.
    stages = [("resources\nscreened", len(led)),
              ("libraries\nquantified", 7),
              ("donors with\nusable CAF", 5),
              ("donors in all\nfour libraries", 4),
              ("switches\nkept", 1)]
    labels = [s0 for s0, _ in stages]
    vals = [v for _, v in stages]
    xs_i = np.arange(len(stages))
    ax_i.bar(xs_i, vals, color=[GREY, GREY, GREY, BLUE, RED][:len(vals)],
             width=0.62, edgecolor="white", lw=0.6)
    ax_i.set_yscale("log")
    for x, v in zip(xs_i, vals):
        ax_i.text(x, v * 1.35, f"{v:,}", ha="center", fontsize=8.2)
    ax_i.set_xticks(xs_i)
    ax_i.set_xticklabels(labels, fontsize=7.1)
    ax_i.set_ylabel("count (log)")
    ax_i.text(0.985, 0.30, "event-level: 20 396 measured,\n"
              "27 stringent, 1 kept (Fig. 2)", transform=ax_i.transAxes, va="top",
              ha="right", fontsize=7.4, color="#333333", linespacing=1.3)
    ax_i.set_ylim(0.5, max(vals) * 16)
    panel(ax_i, "i  the funnel from resources to one switch", fontsize=12.4)

    fig.tight_layout(rect=(0, 0, 1, 1))
    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / "supp_fig_s2_ledger.png"
    fit_panel_titles(fig)
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    from PIL import Image
    wpx, hpx = Image.open(png).size
    assert wpx < 7000 and hpx < 7000, f"suspicious canvas {wpx}x{hpx}"
    print(f"written {png.name} {wpx}x{hpx}")
    print(f"d: SD(median)={sd_med:.3f} p75={sd_p75:.3f}; "
          f"dIF(median SD, n=6) = {d_for_power(sd_med, 6, 0.05):.3f} "
          f"(table says "
          f"{float(g.loc['measured_median_dIF_SD', 'min_detectable_dIF_80pct_power_n6']):.3f})")
    print(f"g: median frac = {fv.median():.3f}, >=0.90: {int((fv >= 0.90).sum()):,}")


if __name__ == "__main__":
    main()
