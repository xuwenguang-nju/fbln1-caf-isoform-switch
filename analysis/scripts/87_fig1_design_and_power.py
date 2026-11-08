"""Fig 1 (2 panels): what the three arms set out to do, and what each delivered.

Revised 2026-10-21. The first version was a *plan* diagram and disagreed with
what was executed in five ways:

  1. Arm A was labelled "junction PSI" but the pipeline is SUPPA2 event-level PSI
     from pseudobulk transcript TPM; the analysis report states these are not
     directly comparable to junction-usage ratios.
  2. Arm B named IsoformSwitchAnalyzeR, which was never run. The delivered test
     is a patient-level paired isoform test with a depth gate.
  3. Arm B promised an "NMD screen"; no NMD analysis was done.
  4. Arm C was drawn as "survival + response" with no outcome. Both were run and
     both are negative (Cox p 0.11-0.32; ICB AUC 0.444-0.468, all < 0.5).
  5. The rule "report when >= 2 of 3 arms agree" does not hold for the headline
     result: FBLN1 is not among Arm A's 27 events.

It also omitted the depth gate, which is the only thing that makes a single-arm
discovery credible (FDR = 0 genome-wide at n = 6; 54% of genes depth-unstable;
exactly 1 event survives).

Panel a therefore now shows, per arm, the executed flow AND the delivered result
with an explicit status. Panel b is unchanged: it is the quantitative reason only
one event survives.

Outputs figures/fig_three_arm_design.png
"""
from __future__ import annotations
from pathlib import Path

import sys as _sys

_sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figstyle import panel, fit_panel_titles  # noqa: E402
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "figures"
PROC = ROOT / "data/processed"
STYLE = Path(__file__).resolve().parent / "_style_runtime.py"
exec(compile(STYLE.read_text(encoding="utf-8"), "runtime.py", "exec"), globals())

AB = json.load(open(ROOT / "data/processed/fig1_arm_b_numbers.json", encoding="utf-8"))
MEAS_SD = AB["measured_sd"]      # paired dIF SD at full depth, fixed-referent statistic
OBS_DIF = AB["fbln1_mean_dIF"]   # FBLN1 mean paired dIF, same statistic
N_PAIRS = 6

COLS = [2.05, 4.25, 6.45, 8.65]
BW = 1.95
REACH = 10.90

ARMS = [
    dict(label="A  single cell\n7 libraries -> 4 donors", colour="#1f6fb2",
         boxes=["CAF annotation", "donor x type\npseudobulk", "SUPPA2\nevent PSI",
                "contrast\nspecificity"],
         gate_idx=None,
         sub="27 events replicate across 4 donors, but none is FBLN1",
         status="delivered", detail="2.47x vs 0.89x control", scolor="#2e7d32"),
    dict(label="B  sorted CAF\n6 paired patients", colour="#c0562b",
         boxes=["Salmon\n(GENCODE v44)", "paired isoform\ntest (dIF)",
                "depth gate\nfull vs 10%", "1 switch\n(FBLN1)"],
         gate_idx=2,
         sub=f"{AB['quantified']:,} quantified -> {AB['depth_unstable']:,} depth-unstable"
             f" -> {AB['n_switch']} survives",
         status="discovery arm",
         detail=f"dIF +{AB['fbln1_mean_dIF']:.3f}, {AB['fbln1_n_pos']}/{AB['fbln1_n']}, "
                f"q {AB['fbln1_q']:.3f}", scolor="#2e7d32"),
    dict(label="C  TCGA-HNSC\n544 samples", colour="#5f7a86",
         boxes=["PSI matrix\n(SpliceSeq)", "tumour vs\nnormal", "survival",
                "ICB response"],
         gate_idx=None,
         sub="42 paired; 501 tumour / 43 normal",
         status="partial", detail="both negative", scolor="#b8342b"),
]


def power_of_d(d, n=N_PAIRS, alpha=0.05):
    df = n - 1
    tc = stats.t.ppf(1 - alpha / 2, df)
    return stats.nct.sf(tc, df, d * np.sqrt(n)) + stats.nct.cdf(-tc, df, d * np.sqrt(n))


def d_for_power(power=0.80, n=N_PAIRS, alpha=0.05, hi=6.0, steps=600):
    ds = np.linspace(1e-4, hi, steps)
    ps = np.array([power_of_d(d, n, alpha) for d in ds])
    ok = np.isfinite(ps)
    ds, ps = ds[ok], ps[ok]
    i = np.where(ps >= power)[0][0]
    if i == 0:
        return float(ds[0])
    d0, d1, p0, p1 = ds[i - 1], ds[i], ps[i - 1], ps[i]
    return float(d0 + (power - p0) * (d1 - d0) / (p1 - p0))


ALPHAS = [(0.05, "#1f6fb2", "a = 0.05"),
          (0.01, "#5b8fb9", "a = 0.01"),
          (0.001, "#8fb6d4", "a = 0.001")]


def panel_a(ax):
    ax.set_xlim(0, 12.55)
    ax.set_ylim(0, 11.1)
    ax.axis("off")
    panel(ax, "a  three arms, and what each actually delivered", fontsize=12.9)

    # column headers
    ax.text(COLS[0], 9.85, "executed pipeline", fontsize=8.7, color="#888888",
            va="center", style="italic")
    ax.text(REACH + 0.30, 9.85, "delivered", fontsize=8.7, color="#888888",
            va="center", style="italic")

    for k, r in enumerate(ARMS):
        y = 7.95 - k * 2.90
        ax.text(0.05, y + 0.52, r["label"], fontsize=10.0, color=r["colour"],
                va="center", ha="left", fontweight="bold", linespacing=1.28)

        for i, (x, txt) in enumerate(zip(COLS, r["boxes"])):
            is_gate = (r["gate_idx"] == i)
            ax.add_patch(FancyBboxPatch(
                (x, y), BW, 1.30,
                boxstyle="round,pad=0.03,rounding_size=0.10",
                linewidth=1.5 if is_gate else 1.05,
                edgecolor=r["colour"], facecolor="#fff3e6" if is_gate else "white",
                zorder=3))
            ax.text(x + BW / 2, y + 0.650, txt, fontsize=8.9, ha="center",
                    va="center", color="#222222", linespacing=1.2, zorder=4,
                    fontweight="bold" if is_gate else "normal")
            if i < 3:
                ax.add_patch(FancyArrowPatch(
                    (x + BW + 0.04, y + 0.650), (COLS[i + 1] - 0.05, y + 0.650),
                    arrowstyle="-|>", mutation_scale=7.5, linewidth=0.85,
                    color=r["colour"], zorder=2))

        ax.text(COLS[0], y - 0.36, r["sub"], fontsize=8.3, color="#666666",
                va="center")

        # status badge
        ax.scatter([REACH + 0.12], [y + 0.80], s=24, color=r["scolor"], zorder=5,
                   edgecolor="white", linewidth=0.7)
        ax.text(REACH + 0.28, y + 0.80, r["status"], fontsize=8.4,
                color=r["scolor"], va="center", fontweight="bold")
        ax.text(REACH + 0.28, y + 0.36, r["detail"], fontsize=8.0,
                color="#555555", va="center")

    # reporting rule, as actually applied
    ax.text(0.05, 0.62,
            "Reporting rule (as applied): a candidate switch must survive the depth gate in the",
            fontsize=9.5, color="#333333", va="center", style="italic")
    ax.text(0.05, 0.24,
            "discovery arm and must not be contradicted by an orthogonal arm; arms that do not "
            "overlap its contrast are context, not confirmation.",
            fontsize=9.5, color="#333333", va="center", style="italic")


def panel_b(ax):
    resolve = d_for_power(0.80, N_PAIRS, ALPHAS[0][0]) * MEAS_SD
    panel(ax, f"b  the paired design can resolve dIF >= {resolve:.3f}", fontsize=12.9)

    x = np.linspace(0.005, 0.42, 400)
    handles = []
    for alpha, colour, lab in ALPHAS:
        df = N_PAIRS - 1
        tc = stats.t.ppf(1 - alpha / 2, df)
        ncp = (x / MEAS_SD) * np.sqrt(N_PAIRS)
        p = stats.nct.sf(tc, df, ncp) + stats.nct.cdf(-tc, df, ncp)
        line, = ax.plot(x, p, color=colour, lw=1.8, label=lab)
        handles.append(line)

    ax.axhline(0.80, color="#999999", lw=0.9, ls=":")
    ax.text(0.008, 0.815, "80% power", fontsize=8.9, color="#777777", va="bottom")

    for alpha, colour, lab in ALPHAS:
        dif = d_for_power(0.80, N_PAIRS, alpha) * MEAS_SD
        ax.plot([dif], [0.80], marker="o", ms=5.2, color=colour,
                mec="white", mew=0.9, zorder=5)
        ax.annotate(f"{dif:.3f}", (dif, 0.80), textcoords="offset points",
                    xytext=(2, -11), fontsize=8.8, color=colour, ha="center")

    ax.axvline(OBS_DIF, color="#b8342b", lw=1.2, ls="--", zorder=1)
    pw_txt = " / ".join(f"{power_of_d(OBS_DIF / MEAS_SD, N_PAIRS, a):.2f}"
                        for a, _, _ in ALPHAS)
    ax.text(OBS_DIF + 0.008, 0.30,
            f"observed FBLN1\ndIF = +{OBS_DIF:.3f}\n(power {pw_txt})",
            fontsize=8.8, color="#b8342b", va="center", linespacing=1.3)

    ax.set_xlim(0, 0.42)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("true paired change in isoform fraction, dIF")
    ax.set_ylabel("power (two-sided paired t, n = 6)")
    ax.legend(handles=handles, frameon=False, fontsize=9.1,
              loc="center left", bbox_to_anchor=(0.012, 0.56), ncol=1,
              title="alpha", title_fontsize=9.1, handlelength=1.6)
    ax.text(0.985, 0.045,
            f"SD of paired dIF = {MEAS_SD:.3f} (measured, B arm)\n"
            "-> why the genome-wide scan finds 0 at n = 6,\n"
            "   yet the pre-registered panel is powered",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8.2, color="#555555", linespacing=1.35)


def main():
    apply_figure_style(sizes=(12.1, 11.5, 10.7))
    # 2026-10-22: canvas was 10.6 x 8.4 and read as too flat -- the roadmap boxes
    # were 4.1:1. Taller canvas gives panel A's boxes and panel B's curves room.
    fig = plt.figure(figsize=(12.2, 12.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.34, 1.0], hspace=0.30,
                          left=0.04, right=0.985, top=0.985, bottom=0.08)
    panel_a(fig.add_subplot(gs[0]))
    panel_b(fig.add_subplot(gs[1]))

    OUT = FIG / "fig_three_arm_design.png"
    fit_panel_titles(fig)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    plt.close(fig)

    det = {f"alpha={a}": round(d_for_power(0.80, N_PAIRS, a) * MEAS_SD, 4)
           for a, _, _ in ALPHAS}
    pw = {f"alpha={a}": round(float(power_of_d(OBS_DIF / MEAS_SD, N_PAIRS, a)), 4)
          for a, _, _ in ALPHAS}
    (PROC / "three_arm_power.json").write_text(json.dumps(dict(
        measured_sd=MEAS_SD, observed_dif=OBS_DIF, n_pairs=N_PAIRS,
        detectable_dif_at_80pct=det, power_at_observed=pw,
        arm_status={
            "A": dict(delivered="contrast specificity, 27 events, 4 donors",
                      contains_fbln1=False),
            "B": dict(delivered="1 switch (FBLN1)", gate="depth full vs 10%"),
            "C": dict(prognosis="negative, Cox p 0.11-0.32",
                      icb_response="negative, AUC 0.444-0.468"),
        }), indent=1), encoding="utf-8")
    print("detectable dIF @80% power:", det)
    print("power at observed dIF   :", pw)
    print("written", OUT)


if __name__ == "__main__":
    main()
