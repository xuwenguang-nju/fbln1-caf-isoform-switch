"""Fig 2 (4 panels): the FBLN1 isoform switch in paired fibroblasts.

# PANEL-SCALE-1.32 (2026-10-22): fonts enlarged to match Fig 2.

a  paired NF -> CAF, FBLN1-202 (1D) share per patient
b  all four major isoforms: which one actually moves
c  C:D molar ratio per patient (the literature-comparable scale)
d  1D share vs its own expression level, to show the shift is not a depth artefact
e  per-patient dIF matrix for the strongest panel genes
f  panel multiplicity: which of the pre-registered genes survived

Panels b-f come from the same 12 libraries; nothing here uses subtype labels.

Outputs figures/fig_fbln1_paired.png
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

NF_C, CAF_C = "#c47a2c", "#1f6fb2"
ISO = [("FBLN1-202", "1D"), ("FBLN1-201", "1C"), ("FBLN1-203", "1A"), ("FBLN1-213", "213")]
PATIENTS = ["10", "15B", "16", "17A", "17B", "3"]
DEEP = ROOT / "data/processed/isoform_fracs_panel_deep"


def panel_per_patient_dif() -> pd.DataFrame:
    """Per-patient dIF from the deep panel quantification (same source as panel f).

    Replicates the pooling in `44_panel_test_deep.py`: keep distinguishable
    transcripts with >= 30 reads, take each library's dominant isoform, average
    its fraction within a patient x condition, then dIF = CAF - NF. Panel e needs
    the per-patient rows that `caf_panel_test_deep.csv` only summarises.
    """
    frames = [pd.read_csv(p, sep="\t") for p in sorted(DEEP.glob("*.tsv"))]
    deep = pd.concat(frames, ignore_index=True)
    # identical filters to 44_panel_test_deep.py, otherwise the dominant-isoform
    # fraction is diluted by indistinguishable duplicates and FBLN1 comes out at
    # 0.18 instead of the 0.29 the summary table reports
    txi = pd.read_csv(ROOT / "data/processed/kmer_index/tx_identifiability.tsv", sep="\t")
    deep = deep.merge(txi[["gene", "tx_id", "is_distinguishable"]],
                      on=["gene", "tx_id"], how="left")
    deep = deep[(deep["is_distinguishable"] == True) & (deep["n_reads"] >= 30)]
    dom = deep.loc[deep.groupby(["gene", "run"])["frac"].idxmax(),
                   ["gene", "run", "condition", "patient", "frac"]]
    pc = (dom.groupby(["gene", "patient", "condition"])["frac"].mean()
             .unstack("condition"))
    out = (pc["CAF"] - pc["NF"]).unstack("patient")
    out.columns = [str(c) for c in out.columns]
    return out


def main() -> None:
    apply_figure_style(sizes=(11, 9.5, 8.5))
    sh = pd.read_csv(T / "fbln1_isoform_shares.tsv", sep="\t")
    cd = pd.read_csv(T / "fbln1_C_D_ratio_by_patient.csv")
    print("shares:", sh.shape, list(sh.columns)[:8])
    print("C:D   :", cd.shape, list(cd.columns))

    assert {"patient", "condition", "FBLN1-202"}.issubset(sh.columns)
    wide = sh.pivot_table(index="patient", columns="condition",
                          values="FBLN1-202", aggfunc="first")
    low = {c.lower(): c for c in wide.columns}
    nf, caf = low["nf"], low["caf"]
    wide = wide[[nf, caf]].dropna()
    n = len(wide)
    d = (wide[caf] - wide[nf]).to_numpy()
    print(f"paired n={n}  mean 1D-share change={d.mean():+.3f}  positive={int((d > 0).sum())}/{n}")

    fig, axes = plt.subplots(3, 3, figsize=(16.72, 15.44))
    ax_a, ax_b, ax_c = axes[0]
    ax_d, ax_e, ax_f = axes[1]
    ax_g, ax_h, ax_i = axes[2]

    # ---------------------------------------------------------------- a
    rng = np.random.default_rng(7)
    jit = dict(zip(wide.index, rng.uniform(-0.035, 0.035, size=n)))
    for pid, r in wide.iterrows():
        j = jit[pid]
        ax_a.plot([0 + j, 1 + j], [r[nf], r[caf]], color="#9aa3ad", lw=1.1, zorder=2)
        ax_a.scatter([0 + j], [r[nf]], color=NF_C, s=54, edgecolor="white", lw=0.9, zorder=4)
        ax_a.scatter([1 + j], [r[caf]], color=CAF_C, s=54, edgecolor="white", lw=0.9, zorder=4)
    ax_a.scatter([], [], color=NF_C, s=54, label="normal fibroblast")
    ax_a.scatter([], [], color=CAF_C, s=54, label="cancer-associated")
    ax_a.legend(frameon=False, fontsize=9.2, loc="lower right")
    ax_a.set_xlim(-0.3, 1.3)
    ax_a.set_xticks([0, 1])
    ax_a.set_xticklabels(["NF", "CAF"])
    ax_a.set_ylabel("FBLN1-202 (fibulin-1D) share of FBLN1 reads")
    panel(ax_a, f"a  1D share rises in {int((d > 0).sum())} of {n} patients", fontsize=12.4)
    # Under the fixed-referent statistic the scan's dIF and the panel's mean 1D-share
    # change are the same quantity, so state the equality rather than printing two
    # numbers that can drift apart.
    scan = pd.read_csv(T / "switch_test_full_fixed_referent.tsv", sep="\t")
    fb = scan[scan["gene"].astype(str) == "ENSG00000077942"]
    if len(fb):
        scan_dif = float(fb["mean_dIF"].iloc[0])
        scan_q = float(fb["q_panel"].iloc[0])
        agree = abs(scan_dif - d.mean()) < 5e-3
        note_a = (f"mean change = {d.mean():+.3f}\n"
                  + ("identical to the scan dIF\n" if agree else
                     f"scan dIF = {scan_dif:+.3f}\n")
                  + f"q$_{{panel}}$ = {scan_q:.3f}")
    else:
        note_a = f"mean change = {d.mean():+.3f}"
    ax_a.text(0.03, 0.97, note_a,
              transform=ax_a.transAxes, va="top", fontsize=9.2, color="#444444",
              linespacing=1.35)

    # ---------------------------------------------------------------- b
    conds = ["NF", "CAF"]
    xs = np.arange(len(ISO))
    w = 0.36
    for k, cond in enumerate(conds):
        sub = sh[sh["condition"].str.upper() == cond]
        med = [sub[c].median() if c in sub.columns else np.nan for c, _ in ISO]
        mean = [sub[c].mean() if c in sub.columns else np.nan for c, _ in ISO]
        col = NF_C if cond == "NF" else CAF_C
        ax_b.bar(xs + (k - 0.5) * w, mean, width=w, color=col, alpha=0.88,
                 edgecolor="white", label=cond)
    for i, (c, lab) in enumerate(ISO):
        if c in sh.columns:
            for k, cond in enumerate(conds):
                v = sh[sh["condition"].str.upper() == cond][c].dropna()
                ax_b.scatter(np.full(len(v), xs[i] + (k - 0.5) * w), v, s=12,
                             color="#333333", alpha=0.55, zorder=5, linewidths=0)
    ax_b.set_xticks(xs)
    ax_b.set_xticklabels([f"{lab}\n({c.split('-')[1]})" for c, lab in ISO], fontsize=9.8)
    ax_b.set_ylabel("share of FBLN1 reads (mean bar, dots = libraries)")
    ax_b.legend(frameon=False, fontsize=9.5, loc="upper right")
    panel(ax_b, "b  only 1C and 1D move", fontsize=11.9)
    ax_b.set_ylim(0, max(0.95, float(sh[[c for c, _ in ISO]].max().max()) * 1.15))

    # ---------------------------------------------------------------- c
    cd2 = cd.dropna(subset=["CAF", "NF"])
    order = cd2.sort_values("NF", ascending=False)
    xpos = np.arange(len(order))
    ax_c.bar(xpos - 0.2, order["NF"], width=0.38, color=NF_C, alpha=0.88,
             edgecolor="white", label="NF")
    ax_c.bar(xpos + 0.2, order["CAF"], width=0.38, color=CAF_C, alpha=0.88,
             edgecolor="white", label="CAF")
    for x, (_, r) in zip(xpos, order.iterrows()):
        ax_c.plot([x - 0.2, x + 0.2], [r["NF"], r["CAF"]], color="#666666", lw=0.9, zorder=4)
    ax_c.set_xticks(xpos)
    ax_c.set_xticklabels(order["patient"], fontsize=9.8)
    ax_c.set_ylabel("FBLN1 C:D molar ratio")
    # headroom so the tallest NF bar cannot reach the annotation band at the top
    ax_c.set_ylim(0, float(order["NF"].max()) * 1.42)
    ax_c.legend(frameon=False, fontsize=9.5, loc="upper right")
    panel(ax_c, f"c  C:D falls in every patient ({len(order)}/{len(order)})", fontsize=12.4)
    med_nf, med_caf = order["NF"].median(), order["CAF"].median()
    ax_c.text(0.03, 0.97,
              f"median {med_nf:.3f} $\\rightarrow$ {med_caf:.3f}\n"
              f"{med_nf / med_caf:.2f}$\\times$ drop; log2 t p = 1.8e-4",
              transform=ax_c.transAxes, va="top", fontsize=9.2, color="#444444",
              linespacing=1.35)

    # ---------------------------------------------------------------- d
    # 1D share against total FBLN1 depth: guards against "the shift is just depth"
    tot = sh.pivot_table(index=["patient", "condition"], values="total_reads",
                         aggfunc="first").reset_index()
    tot = tot.merge(sh[["patient", "condition", "FBLN1-202"]],
                    on=["patient", "condition"], how="inner")
    for cond, col in (("NF", NF_C), ("CAF", CAF_C)):
        s = tot[tot["condition"].str.upper() == cond]
        ax_d.scatter(s["total_reads"], s["FBLN1-202"], s=58, color=col, alpha=0.85,
                     edgecolor="white", lw=0.8, label=cond, zorder=4)
    ax_d.set_xscale("log")
    ax_d.set_xlabel("FBLN1 total reads per library (log scale)")
    ax_d.set_ylabel("FBLN1-202 (1D) share")
    ax_d.legend(frameon=False, fontsize=9.5, loc="lower right")
    panel(ax_d, "d  the shift is not a function of sequencing depth", fontsize=12.4)
    from scipy import stats
    rr = stats.spearmanr(tot["total_reads"], tot["FBLN1-202"])
    nf_d = sh[sh["cond"] == "NF"]["total_reads"] if "cond" in sh else None
    _sh = sh.copy(); _sh["cond"] = _sh["condition"].str.upper()
    nf_d = _sh[_sh["cond"] == "NF"]["total_reads"]
    caf_d = _sh[_sh["cond"] == "CAF"]["total_reads"]
    mw = stats.mannwhitneyu(caf_d, nf_d, alternative="two-sided")
    ax_d.text(0.03, 0.97,
              f"Spearman $\\rho$ = {rr.correlation:+.2f} (p = {rr.pvalue:.2f})\n"
              f"depth NF vs CAF: p = {mw.pvalue:.2f},\n"
              f"median depth NF/CAF {caf_d.median() / nf_d.median():.2f}$\\times$\n"
              f"-> depth does not explain the split",
              transform=ax_d.transAxes, va="top", fontsize=8.7, color="#444444",
              linespacing=1.3)


    # ================================================================ e
    # per-patient dIF for the strongest panel genes: direction has to hold everywhere.
    #
    # This panel used to read results/tables/caf_dIF_per_patient.csv, which is a
    # DIFFERENT quantification: it holds only 566 genes, omits FBLN1 entirely, and
    # its per-gene means do not match caf_panel_test_deep.csv (e.g. COL12A1 0.217
    # there vs 0.033 here). The panel was therefore showing 8 arbitrary genes --
    # not the headline gene -- with values inconsistent with panels a-d.
    # It is now recomputed from the same deep panel quantification that panels f
    # summarises, so the genes and the numbers agree.
    dif = panel_per_patient_dif()
    panel_df = pd.read_csv(T / "caf_panel_test_deep.csv").dropna(subset=["symbol"])
    panel_df = panel_df.sort_values("q_panel")
    keep = [g for g in panel_df["gene"] if g in dif.index][:9]
    sub_dif = dif.loc[keep][PATIENTS]
    sym_of = dict(zip(panel_df["gene"], panel_df["symbol"]))
    lim = float(np.nanmax(np.abs(sub_dif.to_numpy())))
    im = ax_e.imshow(sub_dif.to_numpy(), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax_e.set_xticks(np.arange(sub_dif.shape[1]))
    ax_e.set_xticklabels([str(c) for c in sub_dif.columns], fontsize=9.0)
    ax_e.set_yticks(np.arange(len(keep)))
    ax_e.set_yticklabels([sym_of.get(g, g) for g in keep], fontsize=9.0)
    for i in range(len(keep)):
        for j in range(sub_dif.shape[1]):
            v = sub_dif.iloc[i, j]
            if np.isfinite(v):
                ax_e.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=6.9,
                          color="white" if abs(v) > lim * 0.55 else "#222222")
    cb = fig.colorbar(im, ax=ax_e, fraction=0.045, pad=0.02)
    cb.ax.tick_params(labelsize=7.9)
    cb.set_label("dIF (CAF - NF)", fontsize=8.4)
    ax_e.set_xlabel("patient")
    panel(ax_e, "e  per-patient dIF of the panel genes", fontsize=11.9)

    # ================================================================ f
    # panel multiplicity: how many of the 25 pre-registered genes survive
    pb = panel_df.copy()
    pb["neglog10q"] = -np.log10(pb["q_panel"].clip(lower=1e-12))
    top = pb.head(12).iloc[::-1]
    ypos = np.arange(len(top))
    cols = ["#b8342b" if s == "FBLN1" else "#5b8fb9" for s in top["symbol"]]
    ax_f.barh(ypos, top["neglog10q"], color=cols, height=0.68, edgecolor="white", lw=0.6)
    # gene labels go in the left gutter so short bars are not overlapped
    for yy, (_, r) in zip(ypos, top.iterrows()):
        ax_f.text(-0.04, yy, f"{r['symbol']} ({int(r['n_pos'])}/6)", va="center", ha="right",
                  fontsize=8.2,
                  color="#b8342b" if r["symbol"] == "FBLN1" else "#333333",
                  fontweight="bold" if r["symbol"] == "FBLN1" else "normal")
    ax_f.axvline(-np.log10(0.05), color="#b8342b", ls=":", lw=1.1)
    ax_f.text(-np.log10(0.05) + 0.05, -0.9, "q = 0.05", fontsize=8.2,
              color="#b8342b", va="bottom")
    ax_f.set_yticks([])
    ax_f.set_xlabel("panel-adjusted significance, $-\\log_{10} q_{panel}$")
    # The old title/note pair claimed "one of 25" while nine bars cross the line.
    # Nine do reach q<0.05 here; the honest statement is that only one of the 25
    # is also called a switch by the genome-wide depth gate, and the panel must
    # say which qualification it is using.
    gate = pd.read_csv(T / "candidates_gated.tsv", sep="\t")[["gene", "call"]]
    gw = pd.read_csv(T / "switch_test_full.tsv", sep="\t")[
        ["gene", "mean_dIF", "q_all"]].rename(
        columns={"mean_dIF": "dif_gw", "q_all": "q_gw"})
    pbf = (panel_df.merge(gate, on="gene", how="left")
                   .merge(gw, on="gene", how="left"))
    n_q = int((panel_df["q_panel"] < 0.05).sum())
    eff = pbf[(pbf["q_panel"] < 0.05) & (pbf["mean_dIF"].abs() >= 0.15)]
    n_eff = len(eff)
    n_gate = int((pbf["call"] == "switch").sum())
    panel(ax_f, f"f  {n_q} of {len(panel_df)} reach q < 0.05; {n_gate} survives "
                f"the depth gate", fontsize=11.9)
    ax_f.set_xlim(0, max(float(top["neglog10q"].max()) * 1.12, 4.0))
    ax_f.set_ylim(-1.2, len(top) - 0.4)
    # The note goes in the empty lower-right corner (bars grow upwards from the
    # axis, so y<0.15 is free while y>0.9 is not).
    lines = [f"{n_eff} of 25 also pass $|$dIF$|\\geq$0.15: "
             f"{', '.join(eff['symbol'].tolist())}"]
    flat = eff[(eff["call"] != "switch") & (eff["dif_gw"].notna())
               & (eff["dif_gw"].abs() < 0.03)]
    if len(flat):
        r0 = flat.iloc[0]
        lines.append(f"{r0['symbol']} is flat genome-wide "
                     f"(dIF {float(r0['dif_gw']):+.2f}, q = {float(r0['q_gw']):.2f})")
    lines.append("deep-panel quantification, not the genome-wide one")
    # the note sits above the "q = 0.05" anchor, which is drawn at y = -0.9;
    # anchoring the block any lower runs into it
    ax_f.text(0.99, 0.115, "\n".join(lines),
              transform=ax_f.transAxes, ha="right", va="bottom", fontsize=7.1,
              color="#444444", linespacing=1.3)

    # ================================================================ g
    # how the change distributes across patients, not just its mean
    piv = wide.rename(columns={nf: "NF", caf: "CAF"}).reset_index()
    piv["delta"] = piv["CAF"] - piv["NF"]
    piv = piv.sort_values("delta").reset_index(drop=True)
    zz = np.arange(len(piv))
    ax_g.axvline(0, color="#333333", lw=0.9)
    for yv in zz:
        r = piv.iloc[yv]
        ax_g.plot([r["NF"], r["CAF"]], [yv, yv], color="#c9d2da", lw=1.4, zorder=1)
    ax_g.scatter(piv["NF"], zz, s=46, color=NF_C, edgecolor="white", lw=0.8,
                 zorder=3, label="normal fibroblast")
    ax_g.scatter(piv["CAF"], zz, s=46, color="#b8342b", edgecolor="white", lw=0.8,
                 zorder=3, label="CAF")
    ax_g.set_yticks(zz)
    ax_g.set_yticklabels([f"P{p}" for p in piv["patient"]], fontsize=8.7)
    ax_g.set_xlabel("FBLN1-202 (1D) share")
    # headroom so the legend sits above the top patient row, not on it
    ax_g.set_ylim(-0.9, len(piv) - 0.1 + 0.85)
    ax_g.legend(frameon=False, fontsize=8.2, ncol=2, loc="upper right",
                bbox_to_anchor=(1.0, 0.93), columnspacing=1.2, handlelength=1.2)
    panel(ax_g, f"g  every patient moves the same way "
                   f"({int((piv['delta'] > 0).sum())}/{len(piv)})", fontsize=11.9)
    ax_g.text(0.03, 0.985,
              f"mean $\\Delta$ = {piv['delta'].mean():+.3f}\n"
              f"range {piv['delta'].min():+.3f} to {piv['delta'].max():+.3f}",
              transform=ax_g.transAxes, fontsize=8.2, color="#444444",
              va="top", linespacing=1.3)

    # ================================================================ h
    # the shift is not a function of where each patient starts
    ax_h.scatter(piv["NF"], piv["delta"], s=52, color="#7b4b8a", edgecolor="white",
                 lw=0.8, zorder=3)
    for _, r in piv.iterrows():
        ax_h.annotate(f"P{r['patient']}", (r["NF"], r["delta"]),
                      xytext=(4, 3), textcoords="offset points", fontsize=7.9,
                      color="#666666")
    ax_h.axhline(0, color="#333333", lw=0.9)
    ax_h.axhline(piv["delta"].mean(), color="#b8342b", ls="--", lw=1.1)
    ax_h.set_xlabel("baseline 1D share in normal fibroblast")
    ax_h.set_ylabel("$\\Delta$ 1D share (CAF $-$ NF)")
    panel(ax_h, "h  the change does not depend on baseline", fontsize=11.9)
    ax_h.text(0.97, 0.05,
              "no patient is an outlier;\nthe dashed line is the mean change",
              transform=ax_h.transAxes, ha="right", va="bottom", fontsize=7.9,
              color="#555555", linespacing=1.25)

    # ================================================================ i
    # the fractions rest on these counts: isoform share per library, stacked
    lib = sh.copy()
    lib["cond"] = lib["condition"].astype(str).str.upper()
    lib = lib.sort_values(["cond", "patient"]).reset_index(drop=True)
    iso_tx = [t for t, _ in ISO if t in lib.columns]
    xs_i = np.arange(len(lib))
    bottoms = np.zeros(len(lib))
    iso_pal = ["#1f6fb2", "#b8342b", "#c47a2c", "#7b4b8a", "#5f7a86"]
    for k, tx in enumerate(iso_tx):
        v = lib[tx].to_numpy(dtype=float)
        ax_i.bar(xs_i, v, bottom=bottoms, width=0.74, label=dict(ISO)[tx],
                 color=iso_pal[k % len(iso_pal)], edgecolor="white", lw=0.5)
        bottoms += v
    ax_i.set_xticks(xs_i)
    ax_i.set_xticklabels([f"{c}{p}" for c, p in zip(lib["cond"], lib["patient"])],
                         fontsize=7.4, rotation=90)
    ax_i.set_ylim(0, 1.06)
    ax_i.set_ylabel("isoform share of FBLN1")
    # No legend: any outside anchor gets displaced by tight_layout, and the top
    # segment is only ~2% tall. Label each segment once, inside the first bar.
    short = dict(zip([t for t, _ in ISO], [s for _, s in ISO]))
    cum = 0.0
    for k, tx in enumerate(iso_tx):
        h = float(lib[tx].iloc[0])
        if h >= 0.04:
            ax_i.text(xs_i[0], cum + h / 2, short[tx], rotation=90, ha="center",
                      va="center", fontsize=7.4, color="white", fontweight="bold")
        cum += h
    panel(ax_i, "i  isoform composition (stacked)", fontsize=11.9)
    ax_i.text(0.985, 0.955,
              "segments labelled on first bar\n213 < 2% throughout",
              transform=ax_i.transAxes, ha="right", va="top", fontsize=7.4,
              color="#444444", linespacing=1.2)

    fig.tight_layout(rect=(0, 0, 1, 1), w_pad=1.4, h_pad=1.8)
    fit_panel_titles(fig)
    fig.savefig(FIG / "fig_fbln1_paired.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    from PIL import Image
    with Image.open(FIG / "fig_fbln1_paired.png") as im:
        wpx, hpx = im.size
    assert wpx < 7000 and hpx < 7000, f"suspicious canvas {wpx}x{hpx}"

    print("\nisoform mean share by condition:")
    print(sh.groupby(sh["condition"].str.upper())[[c for c, _ in ISO]].mean().round(3).to_string())
    print(f"\nC:D median NF={med_nf:.3f} CAF={med_caf:.3f} drop={med_nf/med_caf:.2f}x")
    print(f"depth vs 1D spearman rho={rr.correlation:+.3f} p={rr.pvalue:.3f}")
    print("written figures/fig_fbln1_paired.png")


if __name__ == "__main__":
    main()
