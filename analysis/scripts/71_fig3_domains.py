"""Fig 5 (6 panels): protein-level consequence of the FBLN1 1C/1D switch.

# PANEL-SCALE-1.32 (2026-10-22): fonts enlarged to match Fig 2.

a  domain schematic, 1D vs 1C
b  amino-acid identity along the alignment: where the two proteins agree
c  domain retention matrix: which Pfam domains survive in each isoform
d  alternative-segment composition: length and residue make-up of the divergent
   C-termini
e  NEW: exon structure and the 50-nt NMD rule -- both isoforms qualify, so the
   rule does not distinguish them
f  NEW: physicochemical properties of the two divergent C-termini (GRAVY, pI,
   net charge, sheet propensity) computed from the sequences

Everything is derived from GENCODE v44 sequences + UniProt P23142 + InterPro/Pfam
(`remote/easyai/fbln1_structure.json`, `fbln1_sequences.json`, `fbln1_nmd.json`,
`fbln1_physchem.json`).
No structure, stability, NMD, secretion or binding experiment was performed; the
figure says so on the canvas.

Outputs figures/fig_fbln1_domains.png
"""
from __future__ import annotations
import json
from pathlib import Path

import sys as _sys

_sys.path.insert(0, str(Path(__file__).resolve().parent))
from _figstyle import panel, fit_panel_titles  # noqa: E402

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

ROOT = Path(__file__).resolve().parents[2]
REM = ROOT / "remote/easyai"
FIG = ROOT / "figures"
STYLE = Path(__file__).resolve().parent / "_style_runtime.py"
exec(compile(STYLE.read_text(encoding="utf-8"), "runtime.py", "exec"), globals())

D_COLOR, C_COLOR, ALT_COLOR = "#1f6fb2", "#c47a2c", "#9aa3ad"
PFAM_COLORS = {"PF01821": "#6b8f71", "PF07645": "#7b4b8a",
               "PF12662": "#9a6fb0", "PF22914": "#b8342b"}
AA_GROUPS = {
    "hydrophobic": set("AVILMFWC"),
    "polar": set("STNQ"),
    "positive": set("KRH"),
    "negative": set("DE"),
    "special": set("GP"),
    "aromatic": set("FY"),
}


def main() -> None:
    apply_figure_style(sizes=(11, 9.5, 8.5))
    S = json.loads((REM / "fbln1_structure.json").read_text(encoding="utf-8"))
    seq = json.loads((REM / "fbln1_sequences.json").read_text(encoding="utf-8"))
    nmd = json.loads((REM / "fbln1_nmd.json").read_text(encoding="utf-8"))
    pc = json.loads((REM / "fbln1_physchem.json").read_text(encoding="utf-8"))

    div = S["divergence"]["first_differing_residue_1based"]
    d_len = S["isoforms"]["FBLN1-202"]["aa"]
    c_len = S["isoforms"]["FBLN1-201"]["aa"]
    egf9 = S["divergence"]["egf9"]
    ctd = S["c_terminal_domain"]
    pfam = S["pfam_interpro"]
    pD = seq["FBLN1-202_1D"]["protein"]
    pC = seq["FBLN1-201_1C"]["protein"]
    assert len(pD) == d_len and len(pC) == c_len
    assert S["identity"]["d_matches_uniprot"]

    fig = plt.figure(figsize=(12.22, 15.22))
    gs = fig.add_gridspec(4, 2, height_ratios=[0.85, 0.8, 1.05, 0.95],
                          hspace=0.52, wspace=0.24)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, :])
    ax_c = fig.add_subplot(gs[2, 0])
    ax_d = fig.add_subplot(gs[2, 1])
    ax_e = fig.add_subplot(gs[3, 0])
    ax_f = fig.add_subplot(gs[3, 1])

    # ---------------------------------------------------------------- a
    H, GAP = 0.28, 0.34
    yD, yC = 1.0, 1.0 - (H + GAP)

    def bar(ax, y, x0, x1, color, alpha=1.0, hatch=None, edge=None, h=H):
        ax.add_patch(Rectangle((x0, y), x1 - x0, h, facecolor=color,
                               edgecolor=edge or color, alpha=alpha, lw=0.7, hatch=hatch))

    for y, col in ((yD, D_COLOR), (yC, C_COLOR)):
        bar(ax_a, y, 1, div - 1, col, alpha=0.16, edge=col)
    for d in pfam:
        keep_to = min(d["end"], div - 1)
        if keep_to >= d["start"]:
            for y in (yD, yC):
                bar(ax_a, y, d["start"], keep_to, PFAM_COLORS[d["pfam"]], alpha=0.78)
        if d["start"] < div <= d["end"]:
            for y in (yD, yC):
                bar(ax_a, y, div, d["end"], PFAM_COLORS[d["pfam"]], alpha=0.30,
                    hatch="///", edge=PFAM_COLORS[d["pfam"]])
    bar(ax_a, yD, div, d_len, D_COLOR, alpha=0.85)
    bar(ax_a, yC, div, c_len, ALT_COLOR, alpha=0.85)
    ax_a.add_patch(Rectangle((ctd["start"], yD), ctd["end"] - ctd["start"], H,
                             facecolor=PFAM_COLORS["PF22914"],
                             edgecolor="#7a1f18", lw=0.9, alpha=0.94))
    ax_a.text((ctd["start"] + ctd["end"]) / 2, yD + H / 2, f"{ctd['pfam']}\nC-term Ig-like",
              ha="center", va="center", fontsize=7.9, color="white", linespacing=1.05)
    ax_a.add_patch(Rectangle((ctd["start"], yC), min(ctd["end"], c_len) - ctd["start"], H,
                             facecolor="none", edgecolor=PFAM_COLORS["PF22914"],
                             lw=1.0, ls=(0, (2, 2))))
    ax_a.text((ctd["start"] + min(ctd["end"], c_len)) / 2, yC + H / 2, "absent",
              ha="center", va="center", fontsize=8.2, color=PFAM_COLORS["PF22914"])
    ax_a.axvline(div, color="#222222", lw=1.1, ls=":")
    ax_a.text(div + 4, yD + H + 0.07, f"divergence {div}", fontsize=9.0, va="bottom")
    ax_a.text(-14, yD + H / 2, f"1D  {d_len} aa", ha="right", va="center", fontsize=9.8,
              color=D_COLOR)
    ax_a.text(-14, yC + H / 2, f"1C  {c_len} aa", ha="right", va="center", fontsize=9.8,
              color=C_COLOR)
    ax_a.set_xlim(-150, d_len + 30)
    ax_a.set_ylim(yC - 0.18, yD + H + 0.46)
    ax_a.tick_params(axis="x", pad=2)
    ax_a.set_yticks([])
    ax_a.set_xticks([0, 100, 200, 300, 400, 500, 600, 700])
    ax_a.set_xlabel("residue", labelpad=2)
    for s in ("left", "right", "top"):
        ax_a.spines[s].set_visible(False)
    panel(ax_a, "a  domain architecture", fontsize=12.4)
    ax_a.legend(handles=[Patch(facecolor=PFAM_COLORS[k], label=k) for k in
                         ("PF01821", "PF07645", "PF12662", "PF22914")] +
                        [Patch(facecolor=D_COLOR, label="1D C-terminus"),
                         Patch(facecolor=ALT_COLOR, label="1C C-terminus")],
                frameon=False, fontsize=7.8, ncol=6, loc="upper left",
                bbox_to_anchor=(0.0, 1.02), columnspacing=0.8, handlelength=1.1)

    # ---------------------------------------------------------------- b
    n = min(len(pD), len(pC))
    ident = np.array([pD[i] == pC[i] for i in range(n)])
    # 20-residue windowed identity
    win = 20
    xs, ys_ = [], []
    for i in range(0, n - win, 5):
        xs.append(i + win / 2)
        ys_.append(ident[i:i + win].mean())
    ax_b.plot(xs, ys_, color="#333333", lw=1.3)
    ax_b.axvline(div, color="#222222", lw=1.1, ls=":")
    ax_b.axhline(1.0, color="#cccccc", lw=0.8, ls="--")
    ax_b.fill_between([div, n], 0, 1.05, color="#f2e6e4", alpha=0.7, zorder=0)
    ax_b.text(div + 8, 1.06, "C-termini differ", fontsize=9.0, color="#b8342b",
              va="top")
    ax_b.set_xlim(0, n)
    ax_b.set_ylim(0, 1.12)
    ax_b.set_xlabel("residue")
    ax_b.set_ylabel(f"identity in {win}-aa window")
    panel(ax_b, "b  the two proteins are identical up to the switch point", fontsize=12.4)
    ax_b.text(0.015, 0.30,
              f"EGF-like 9 ({egf9['start']}-{egf9['end']}): only "
              f"{egf9['egf9_residues_kept_by_both']}/{egf9['egf9_total']} residues shared",
              transform=ax_b.transAxes, ha="left", va="top", fontsize=8.4,
              color="#555555")

    # ---------------------------------------------------------------- c
    order = ["PF01821", "PF07645", "PF12662", "PF22914"]
    names = {"PF01821": "Anaphylatoxin", "PF07645": "Ca-binding EGF",
             "PF12662": "Complement Clr-like EGF", "PF22914": "Fibulin C-term Ig-like"}
    mat = np.zeros((len(order), 2))
    for i, acc in enumerate(order):
        for j, iso in enumerate(("FBLN1-202", "FBLN1-201")):
            kept = S["isoforms"][iso]["pfam_kept"]
            mat[i, j] = 1.0 if any(s.startswith(acc) for s in kept) else 0.0
    ax_c.imshow(mat, cmap="Blues", vmin=0, vmax=1.0, aspect="auto")
    for i in range(len(order)):
        for j in range(2):
            ax_c.text(j, i, "present" if mat[i, j] else "absent", ha="center", va="center",
                      fontsize=9.2, color="white" if mat[i, j] else "#b8342b",
                      fontweight="bold" if not mat[i, j] else "normal")
    ax_c.set_xticks([0, 1])
    ax_c.set_xticklabels(["1D (703 aa)", "1C (683 aa)"], fontsize=9.8)
    ax_c.set_yticks(range(len(order)))
    ax_c.set_yticklabels([f"{acc}\n{names[acc]}" for acc in order], fontsize=8.4)
    panel(ax_c, "c  domain retention", fontsize=12.4)
    ax_c.tick_params(length=0)

    # ---------------------------------------------------------------- d
    segD = pD[div - 1:]
    segC = pC[div - 1:]
    groups = list(AA_GROUPS)
    fracD = [sum(aa in AA_GROUPS[g] for aa in segD) / len(segD) for g in groups]
    fracC = [sum(aa in AA_GROUPS[g] for aa in segC) / len(segC) for g in groups]
    xg = np.arange(len(groups))
    w = 0.36
    ax_d.bar(xg - w / 2, fracD, width=w, color=D_COLOR, label=f"1D ({len(segD)} aa)")
    ax_d.bar(xg + w / 2, fracC, width=w, color=ALT_COLOR, label=f"1C ({len(segC)} aa)")
    ax_d.set_xticks(xg)
    ax_d.set_xticklabels(groups, fontsize=9.0, rotation=20)
    ax_d.set_ylabel("fraction of the C-terminal segment")
    ax_d.legend(frameon=False, fontsize=8.6, loc="upper right",
                bbox_to_anchor=(1.0, 0.86))
    panel(ax_d, "d  composition of the divergent segment", fontsize=12.4)
    diff = sum(1 for a, b in zip(segD, segC) if a != b)
    ax_d.set_ylim(0, max(max(fracD), max(fracC)) * 1.42)
    ax_d.text(0.02, 0.98,
              f"{diff} of {min(len(segD), len(segC))} aligned residues differ "
              f"past position {div}",
              transform=ax_d.transAxes, ha="left", va="top", fontsize=8.2,
              color="#444444", linespacing=1.25)

    # ---------------------------------------------------------------- e
    # the 50-nt NMD rule: stop codon versus the last exon-exon junction
    # transcript coordinates are built from the GENCODE exon list
    rows = [("FBLN1-202_1D", "1D", D_COLOR), ("FBLN1-201_1C", "1C", C_COLOR)]
    max_tx = 0
    for k, (lab, short, col) in enumerate(rows):
        d = nmd["transcripts"][lab]
        exons = d["exons_transcript_order"]
        cum = 0
        bounds = []
        for s_, e_ in exons:
            bounds.append((cum, cum + (e_ - s_), s_, e_))
            cum += e_ - s_
        tx_len = cum
        max_tx = max(max_tx, tx_len)
        yv = 1.0 - k * 0.52
        for j, (t0, t1, _, _) in enumerate(bounds):
            term = j == len(bounds) - 1
            ax_e.add_patch(Rectangle((t0, yv), t1 - t0, 0.16,
                                     facecolor=col if term else "#c9d2da",
                                     edgecolor=col if term else "#8a929b",
                                     lw=0.6, alpha=0.95 if term else 0.8))
        ax_e.plot([0, tx_len], [yv + 0.08] * 2, color="#e3e6e9", lw=0.8, zorder=0)
        # the last exon-exon junction is the start of the terminal exon in
        # transcript coordinates; the stop sits inside the terminal exon
        junc = bounds[-1][0]
        stop_tx = junc + (int(d["stop_codon_start"]) - bounds[-1][2])
        dist = int(d["stop_to_last_junction_nt"])
        ax_e.plot([junc, junc], [yv - 0.05, yv + 0.24], color="#222222", lw=1.1, ls=":")
        ax_e.plot([stop_tx, stop_tx], [yv - 0.05, yv + 0.24], color="#b8342b", lw=1.4)
        ax_e.annotate("", xy=(junc, yv - 0.03), xytext=(stop_tx, yv - 0.03),
                      arrowprops=dict(arrowstyle="<->", color="#b8342b", lw=0.9))
        ax_e.text((junc + stop_tx) / 2, yv - 0.09, f"{dist} nt", ha="center",
                  va="top", fontsize=8.2, color="#b8342b")
        ax_e.text(-max_tx * 0.02, yv + 0.08, f"{short}  {d['n_exons']} exons",
                  ha="right", va="center", fontsize=9.0, color=col)
    ax_e.set_xlim(-max_tx * 0.26, max_tx * 1.08)
    ax_e.set_ylim(0.28, 1.72)          # headroom band above the two models
    ax_e.set_yticks([])
    ax_e.set_xlabel("transcript position (nt)")
    for s in ("left", "right", "top"):
        ax_e.spines[s].set_visible(False)
    panel(ax_e, "e  both isoforms pass the 50-nt NMD rule", fontsize=12.4)
    dD = nmd["transcripts"]["FBLN1-202_1D"]
    dC = nmd["transcripts"]["FBLN1-201_1C"]
    ax_e.text(0.015, 0.995,
              f"dotted = last exon-exon junction; red = stop codon\n"
              f"both stops are >50 nt from the junction "
              f"({dD['stop_to_last_junction_nt']} / {dC['stop_to_last_junction_nt']} nt),\n"
              f"so the rule does NOT separate 1C from 1D\n"
              f"3'UTR in the terminal exon: {dD['utr3_nt_in_terminal_exon']} vs "
              f"{dC['utr3_nt_in_terminal_exon']} nt",
              transform=ax_e.transAxes, ha="left", va="top", fontsize=7.9,
              color="#444444", linespacing=1.35)

    # ---------------------------------------------------------------- f
    # physicochemical character of the two divergent C-termini
    segD_k = "1D C-term 567-703"
    segC_k = "1C C-term 567-683"
    metrics = [("GRAVY", "gravy"), ("pI", "pI"), ("net charge\n(pH 7)", "net_charge_pH7"),
               ("sheet\npropensity", "E")]
    xf2 = np.arange(len(metrics))
    vD = [pc["segments"][segD_k][k] for _, k in metrics]
    vC = [pc["segments"][segC_k][k] for _, k in metrics]
    wf = 0.36
    ax_f.bar(xf2 - wf / 2, vD, width=wf, color=D_COLOR,
             label=f"1D {pc['segments'][segD_k]['length']} aa")
    ax_f.bar(xf2 + wf / 2, vC, width=wf, color=ALT_COLOR,
             label=f"1C {pc['segments'][segC_k]['length']} aa")
    for x, a, b in zip(xf2, vD, vC):
        ax_f.text(x - wf / 2, a + (0.12 if a >= 0 else -0.12), f"{a:+.2f}",
                  ha="center", va="bottom" if a >= 0 else "top", fontsize=7.9,
                  color=D_COLOR)
        ax_f.text(x + wf / 2, b + (0.12 if b >= 0 else -0.12), f"{b:+.2f}",
                  ha="center", va="bottom" if b >= 0 else "top", fontsize=7.9,
                  color="#8a5a1e")
    ax_f.axhline(0, color="#333333", lw=0.9)
    ax_f.set_xticks(xf2)
    ax_f.set_xticklabels([lab for lab, _ in metrics], fontsize=8.7)
    ax_f.legend(frameon=False, fontsize=8.4, loc="upper left")
    ax_f.set_ylim(min(min(vC), min(vD)) - 1.6, max(max(vC), max(vD)) + 2.6)
    ax_f.set_ylabel("value")
    panel(ax_f, "f  the 1C tail is more basic", fontsize=12.4)
    ax_f.text(0.985, 0.985,
              f"pI {pc['segments'][segC_k]['pI']:.2f} vs "
              f"{pc['segments'][segD_k]['pI']:.2f}\n"
              f"charge at pH 7: {pc['segments'][segC_k]['net_charge_pH7']:+.2f} vs "
              f"{pc['segments'][segD_k]['net_charge_pH7']:+.2f}\n"
              f"shared 1-566 for reference:\n"
              f"  GRAVY {pc['segments']['shared 1-566']['gravy']:+.3f}, "
              f"q {pc['segments']['shared 1-566']['net_charge_pH7']:+.1f}",
              transform=ax_f.transAxes, ha="right", va="top", fontsize=7.8,
              color="#444444", linespacing=1.35)

    # The disclaimer must sit INSIDE the canvas. A figure-level fig.text at y=0.008
    # is pushed outside when bbox_inches="tight" recomputes the bounding box after
    # the banner was removed, so it is anchored to the last axes instead.
    ax_f.text(0.0, -0.30,
              "Annotation only: sequence from GENCODE v44, domains from UniProt P23142 and "
              "InterPro/Pfam,\nNMD from the 50-nt rule on the GENCODE exon model. No structure, "
              "stability, NMD-reporter, secretion or binding\nassay was performed; any "
              "partner-binding statement is a literature-derived hypothesis.",
              transform=ax_f.transAxes, ha="left", va="top",
              fontsize=7.7, color="#555555", linespacing=1.35)
    fit_panel_titles(fig)
    fig.savefig(FIG / "fig_fbln1_domains.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    from PIL import Image
    with Image.open(FIG / "fig_fbln1_domains.png") as im:
        wpx, hpx = im.size
    assert wpx < 7000 and hpx < 7000, f"suspicious canvas {wpx}x{hpx}"

    print(f"1D={d_len} 1C={c_len} div={div} EGF9={egf9['start']}-{egf9['end']}")
    print(f"identical prefix={div - 1} aa; differing aligned residues past {div}: {diff}")
    print(f"NMD: 1D stop {dD['stop_to_last_junction_nt']} nt from junction, "
          f"1C {dC['stop_to_last_junction_nt']} nt; both >50 -> not discriminating")
    print(f"physchem: 1C pI {pc['segments'][segC_k]['pI']:.2f} vs "
          f"1D {pc['segments'][segD_k]['pI']:.2f}")
    print(f"written figures/fig_fbln1_domains.png {wpx}x{hpx}")


if __name__ == "__main__":
    main()
