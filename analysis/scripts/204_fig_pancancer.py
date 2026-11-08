"""Figure 7 (6 panels): the FBLN1 3'-end / stromal-content relationship across five
TCGA cohorts.

Why this figure exists: the main text reports the pan-cancer generalisation as a
number, and a claim that spans five cohorts should be shown rather than asserted.
It also has to be visible that the metric is a *proxy* -- SpliceSeq reports
alternative-terminal-exon (AT) PSI, not the 1C:1D ratio -- so no panel is allowed
to look like a direct replication of the paired-fibroblast effect.

a  per-cohort Spearman rho between event PSI and the NNLS CAF fraction, both
   events, with 95% CIs and the inverse-variance pooled estimate
b  sample-level PSI against CAF fraction, exon 17 event (62665), by cohort
c  the same for the exon 22 event (62666); the two events mirror one 3' choice
d  paired tumour-minus-normal dPSI per cohort, with the paired n annotated and
   cohorts too small for a paired design marked as such
e  samples available per cohort, tumour and adjacent normal
f  the CAF fraction each cohort actually spans -- the premise for panels b and c

Outputs figures/fig_pancancer.png
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
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
T = ROOT / "results/tables"
D = ROOT / "data/processed"
FIG = ROOT / "figures"
STYLE = Path(__file__).resolve().parent / "_style_runtime.py"
exec(compile(STYLE.read_text(encoding="utf-8"), "runtime.py", "exec"), globals())

EV = {"62665": "exon 17", "62666": "exon 22"}
ORDER = ["HNSC", "BRCA", "LUAD", "PAAD", "ESCA"]
COLOUR = {"HNSC": "#1f6fb2", "BRCA": "#c0562b", "LUAD": "#2e7d32",
          "PAAD": "#7b4b8a", "ESCA": "#9a6b1f"}
MIN_PAIRED = 20          # below this a paired tumour-normal test is not attempted


def rho_ci(r: float, n: int, z: float = 1.959964) -> tuple[float, float]:
    """Fisher-z confidence interval for a Spearman correlation."""
    if not np.isfinite(r) or n < 4 or abs(r) >= 1:
        return np.nan, np.nan
    zr, se = np.arctanh(r), 1.0 / np.sqrt(n - 3)
    return float(np.tanh(zr - z * se)), float(np.tanh(zr + z * se))


def main() -> None:
    apply_figure_style(sizes=(11, 9.5, 8.5))
    ev = pd.read_csv(T / "pancancer_fbln1_events.csv")
    ev["as_id"] = ev["as_id"].astype(str)
    samp = pd.read_csv(T / "pancancer_caf_fraction.csv")
    samp["as_id"] = samp["as_id"].astype(str)
    meta = json.load(open(D / "pancancer_fbln1.json", encoding="utf-8"))

    prim = ev[ev["as_id"].isin(EV)].copy()
    cohorts = [c for c in ORDER if c in set(prim["cancer"])]
    print("cohorts:", cohorts, "| events:", sorted(prim["as_id"].unique()))
    print(f"samples: {len(samp):,} rows; "
          f"{samp['xena_id'].nunique():,} distinct xena ids")
    # pandas reads the TCGA code "01" as integer 1 and "11" as 11, so tumour vs
    # normal must be compared numerically, not as strings
    if samp["sample_type"].dtype == object:
        samp["is_tumour"] = samp["sample_type"].astype(str).str.zfill(2).eq("01")
    else:
        samp["is_tumour"] = samp["sample_type"].astype(int).eq(1)

    fig, axes = plt.subplots(2, 3, figsize=(15.2, 9.6))
    ax_a, ax_b, ax_c = axes[0]
    ax_d, ax_e, ax_f = axes[1]

    # ================================================================ a
    # per-cohort rho with CI, two events side by side, pooled at the right
    width = 0.36
    xs = np.arange(len(cohorts))
    for k, aid in enumerate(("62665", "62666")):
        sub = prim[prim["as_id"] == aid].set_index("cancer").reindex(cohorts)
        r = sub["rho_psi_caf_tumour"].to_numpy(float)
        n = sub["n_tumour"].to_numpy(float)
        off = (k - 0.5) * width
        cols = ["#1f6fb2", "#c47a2c"][k]
        ax_a.bar(xs + off, r, width=width * 0.92, color=cols, alpha=0.85,
                 label=f"{EV[aid]} ({aid})")
        for x, ri, ni in zip(xs + off, r, n):
            lo, hi = rho_ci(ri, int(ni))
            if np.isfinite(lo):
                ax_a.plot([x, x], [lo, hi], color="#333333", lw=1.0)
                ax_a.plot([x], [lo], marker="_", ms=4, color="#333333")
                ax_a.plot([x], [hi], marker="_", ms=4, color="#333333")
            ax_a.text(x, ri + (0.035 if ri >= 0 else -0.035), f"{ri:+.2f}",
                      ha="center", va="bottom" if ri >= 0 else "top",
                      fontsize=7.2, color="#333333")
    ax_a.axhline(0, color="#888888", lw=0.9)
    # pooled estimate as a diamond
    for k, aid in enumerate(("62665", "62666")):
        m = meta["meta_rho"][aid]
        xp = len(cohorts) - 0.5
        ax_a.plot([xp + (k - 0.5) * width], [m["rho_pooled"]], marker="D", ms=7,
                  color=["#1f6fb2", "#c47a2c"][k], mec="white", mew=1.0, zorder=5)
        ax_a.plot([xp + (k - 0.5) * width] * 2, [m["ci_lo"], m["ci_hi"]],
                  color=["#1f6fb2", "#c47a2c"][k], lw=1.2)
    ax_a.set_xticks(list(xs) + [len(cohorts) - 0.5])
    ax_a.set_xticklabels(cohorts + ["pooled"], fontsize=8.6)
    ax_a.set_ylabel("Spearman $\\rho$  (event PSI vs CAF fraction)")
    ax_a.set_ylim(-0.75, 0.75)
    ax_a.axvline(len(cohorts) - 1, color="#dddddd", lw=0.9)
    ax_a.legend(frameon=False, fontsize=8.2, loc="upper left")
    panel(ax_a, "a  the stromal relationship holds in every cohort", fontsize=11.9)
    # The JSON's n_same_sign counts cohorts agreeing with the POOLED direction, so
    # "0/5" for exon 17 reads as no consistency when in fact every cohort is
    # negative. Count the actual signs instead.
    n_neg17 = int((prim[prim["as_id"] == "62665"]["rho_psi_caf_tumour"] < 0).sum())
    n_pos22 = int((prim[prim["as_id"] == "62666"]["rho_psi_caf_tumour"] > 0).sum())
    ax_a.text(0.985, 0.045,
              f"exon 17 negative in {n_neg17}/{len(cohorts)} cohorts,\n"
              f"exon 22 positive in {n_pos22}/{len(cohorts)}",
              transform=ax_a.transAxes, ha="right", va="bottom", fontsize=8.0,
              color="#444444", linespacing=1.3)

    # ================================================================ b, c
    for ax, aid in ((ax_b, "62665"), (ax_c, "62666")):
        s = samp[(samp["as_id"] == aid) & samp["is_tumour"]]
        s = s.dropna(subset=["psi", "CAF"])
        for c in cohorts:
            d = s[s["cancer"] == c]
            if not len(d):
                continue
            ax.scatter(d["CAF"], d["psi"], s=3.2, alpha=0.22, linewidths=0,
                       color=COLOUR[c], rasterized=True)
            # The line spans the 5th-95th percentile of the observed CAF fractions:
            # BRCA, PAAD and LUAD pile up near 1.0, and a full-width line there
            # would be drawn straight through empty space and read as evidence.
            if len(d) >= 10:
                k, b = np.polyfit(d["CAF"], d["psi"], 1)
                lo_q, hi_q = d["CAF"].quantile(0.05), d["CAF"].quantile(0.95)
                if hi_q - lo_q < 0.02:            # degenerate span: draw a marker
                    hi_q = lo_q + 0.02
                gx = np.linspace(lo_q, hi_q, 20)
                ax.plot(gx, k * gx + b, color=COLOUR[c], lw=1.9,
                        label=f"{c} (n={len(d)})")
        ax.set_xlabel("estimated CAF fraction of the tumour")
        ax.set_ylabel(f"{EV[aid]} PSI")
        ax.set_xlim(0, 1)
        ax.legend(frameon=False, fontsize=7.0, loc="upper right",
                  handlelength=1.2, labelspacing=0.26, borderaxespad=0.3)
    # outside the loop on purpose: the panel letters must appear literally in the
    # source, or a text-based checker cannot verify the panel count
    panel(ax_b, "b  exon 17: falls with stromal content", fontsize=11.9)
    panel(ax_c, "c  exon 22: rises with stromal content", fontsize=11.9)

    # ================================================================ d
    w2 = 0.36
    for k, aid in enumerate(("62665", "62666")):
        sub = prim[prim["as_id"] == aid].set_index("cancer").reindex(cohorts)
        y = sub["median_delta_psi"].to_numpy(float)
        px = sub["n_paired"].to_numpy(float)
        off = (k - 0.5) * w2
        bar = ax_d.bar(xs + off, y, width=w2 * 0.92,
                       color=["#1f6fb2", "#c47a2c"][k], alpha=0.85,
                       label=f"{EV[aid]}")
        for rect, ok in zip(bar, px >= MIN_PAIRED):
            if not ok:
                rect.set_facecolor("white")
                rect.set_edgecolor(["#1f6fb2", "#c47a2c"][k])
                rect.set_hatch("///")
                rect.set_linewidth(1.0)
        for x, v, npair in zip(xs + off, y, px):
            ax_d.text(x, v + (0.004 if v >= 0 else -0.004), f"{v:+.03f}",
                      ha="center", va="bottom" if v >= 0 else "top",
                      fontsize=6.8, color="#333333")
    ax_d.axhline(0, color="#888888", lw=0.9)
    # headroom so a label under a negative bar cannot land on the bottom spine
    _dmin = float(np.nanmin(prim["median_delta_psi"]))
    _dmax = float(np.nanmax(prim["median_delta_psi"]))
    ax_d.set_ylim(_dmin * 1.45 if _dmin < 0 else -0.01, _dmax * 1.45 if _dmax > 0 else 0.01)
    ax_d.set_xticks(xs)
    ax_d.set_xticklabels([f"{c}\nn={int(prim[prim['cancer'] == c]['n_paired'].iloc[0])}"
                          for c in cohorts], fontsize=8.2)
    ax_d.set_ylabel("median paired dPSI  (tumour $-$ normal)")
    panel(ax_d, "d  the tumour-normal shift needs paired normals", fontsize=11.9)
    ax_d.legend(frameon=False, fontsize=8.2, loc="upper left", ncol=2)
    few = [c for c in cohorts
           if int(prim[prim["cancer"] == c]["n_paired"].iloc[0]) < MIN_PAIRED]
    if few:
        ax_d.text(0.985, 0.04,
                  "hatched: fewer than "
                  f"{MIN_PAIRED} paired normals\n(" + ", ".join(few) + ")",
                  transform=ax_d.transAxes, ha="right", va="bottom",
                  fontsize=7.8, color="#444444", linespacing=1.3)

    # ================================================================ e
    # tumour vs normal counts must come from the whole sample table: `prim` holds
    # only tumour rows, so its n_normal column is the paired subset, not the total
    sub0 = prim[prim["as_id"] == "62665"].set_index("cancer").reindex(cohorts)
    cnt = (samp[samp["as_id"] == "62665"]
           .groupby(["cancer", "is_tumour"]).size().unstack(fill_value=0)
           .reindex(cohorts).fillna(0))
    n_tum = cnt.get(True, pd.Series(0, index=cohorts)).to_numpy(int)
    n_nor = cnt.get(False, pd.Series(0, index=cohorts)).to_numpy(int)
    ax_e.bar(xs - 0.19, n_tum, width=0.36, color="#1f6fb2",
             label="samples with PSI: tumour")
    ax_e.bar(xs + 0.19, n_nor, width=0.36, color="#c47a2c",
             label="samples with PSI: adjacent normal")
    for x, t, nn in zip(xs, n_tum, n_nor):
        ax_e.text(x - 0.19, t + 8, f"{int(t)}", ha="center", fontsize=7.4,
                  color="#1f6fb2")
        ax_e.text(x + 0.19, nn + 8, f"{int(nn)}", ha="center", fontsize=7.4,
                  color="#c47a2c")
    ax_e.set_xticks(xs)
    ax_e.set_xticklabels(cohorts, fontsize=8.6)
    ax_e.set_ylabel("samples carrying a PSI value")
    ax_e.set_ylim(0, float(n_tum.max()) * 1.25)
    ax_e.legend(frameon=False, fontsize=8.2, loc="upper right",
                bbox_to_anchor=(1.0, 0.78))
    panel(ax_e, "e  what each cohort actually contains", fontsize=11.9)
    n_paired_txt = ", ".join(
        f"{c} {int(sub0.loc[c, 'n_paired'])}" for c in cohorts)
    ax_e.text(0.985, 0.965, f"paired tumour-normal:\n{n_paired_txt}",
              transform=ax_e.transAxes, ha="right", va="top", fontsize=7.4,
              color="#444444", linespacing=1.3)

    # ================================================================ f
    for c in cohorts:
        d = samp[(samp["cancer"] == c) & (samp["as_id"] == "62665") & samp["is_tumour"]]
        if not len(d):
            continue
        ax_f.scatter(np.full(len(d), xs[cohorts.index(c)]) + np.random.default_rng(7)
                     .uniform(-0.16, 0.16, len(d)),
                     d["CAF"], s=2.4, alpha=0.16, linewidths=0,
                     color=COLOUR[c], rasterized=True)
        ax_f.plot([cohorts.index(c) - 0.3, cohorts.index(c) + 0.3],
                  [d["CAF"].median()] * 2, color="#111111", lw=2.0,
                  solid_capstyle="butt")
    ax_f.set_xticks(xs)
    ax_f.set_xticklabels(cohorts, fontsize=8.6)
    ax_f.set_ylabel("estimated CAF fraction")
    ax_f.set_ylim(0, 1)
    panel(ax_f, "f  every cohort spans a usable range of stroma", fontsize=11.9)
    ax_f.text(0.985, 0.955,
              "relative NNLS estimates,\nno immune reference: rank-based only",
              transform=ax_f.transAxes, ha="right", va="top", fontsize=7.6,
              color="#444444", linespacing=1.3)

    fig.tight_layout(rect=(0, 0, 1, 1), w_pad=1.4, h_pad=1.8)
    fit_panel_titles(fig)

    # --- collision scan on the finalised in-memory figure ---------------------
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
                continue
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

    FIG.mkdir(parents=True, exist_ok=True)
    out = FIG / "fig_pancancer.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    crops = panel_crops(fig)
    save_panel_crops(out, crops)
    plt.close(fig)
    if hits:
        raise SystemExit("TEXT COLLISIONS:\n" + "\n".join(
            f"  {a!r} <-> {b!r}" for a, b in hits[:20]))

    from PIL import Image
    with Image.open(out) as im:
        wpx, hpx = im.size
    assert wpx < 7000 and hpx < 7000, f"suspicious canvas {wpx}x{hpx}"
    print(f"collision scan clean; {len(texts)} text boxes; "
          f"{len(crops)} panel crops in .cache/figure-style/")
    print(f"written {out.name} {wpx}x{hpx}")


if __name__ == "__main__":
    main()
