"""Does the single-cell CAF>epithelial 1D result survive removal of length normalisation?

The cell-type TPM table was built with bulk-style TPM (count / transcript length).
In 3'-tag 10x data one molecule yields roughly one UMI regardless of transcript
length, so the length-free quantity is the UMI ratio, not the TPM ratio.

TPM and UMI ratios are related exactly:
    TPM_C / TPM_D = (n_C / L_C) / (n_D / L_D)
  => n_C / n_D     = (TPM_C / TPM_D) * (L_C / L_D)

So the UMI-based 1D share can be recovered algebraically from the stored TPMs,
without re-reading the count matrices. If both scales give CAF > epithelial, the
result is not an artefact of the length correction (or of the different 3'UTR
lengths of the two isoforms).
"""
from __future__ import annotations
import glob
import json
import os
from pathlib import Path

EXPR = Path("/root/as_psi/expr")
SEQ = Path("/root/as_psi/ref/pc_transcripts.fa")
OUT = Path("/root/hnscc_caf/results/fbln1_singlecell_validation.json")
TX_C, TX_D = "ENST00000262722", "ENST00000327858"


def tx_lengths():
    want = {TX_C, TX_D}
    out, name, n = {}, None, 0
    for line in open(SEQ):
        if line[0] == ">":
            if name in want:
                out[name] = n
            name = line[1:].split("|")[0].split(".")[0]
            n = 0
        else:
            n += len(line.rstrip("\n"))
    if name in want:
        out[name] = n
    return out


def main() -> None:
    L = tx_lengths()
    assert set(L) == {TX_C, TX_D}, f"lengths not found: {L}"
    lc, ld = L[TX_C], L[TX_D]
    print(f"transcript lengths: 1C={lc} nt, 1D={ld} nt, L_C/L_D={lc / ld:.4f}\n")

    rows = []
    for f in sorted(glob.glob(str(EXPR / "*_celltype_tpm.tsv"))):
        lib = os.path.basename(f).replace("_celltype_tpm.tsv", "")
        hdr = open(f).readline().rstrip("\n").split("\t")
        vals = {}
        for line in open(f):
            p = line.rstrip("\n").split("\t")
            key = p[0][:15]
            if key in (TX_C, TX_D):
                vals[key] = [float(x) for x in p[1:]]
        if len(vals) != 2:
            print(f"{lib}: FBLN1 transcripts missing, skipped")
            continue
        for i, ct in enumerate(hdr):
            tpm_c, tpm_d = vals[TX_C][i], vals[TX_D][i]
            if tpm_c + tpm_d <= 0:
                continue
            share_tpm = tpm_d / (tpm_c + tpm_d)
            # invert to the UMI (length-free) scale
            if tpm_d > 0:
                n_c_over_n_d = (tpm_c / tpm_d) * (lc / ld)
                share_umi = 1.0 / (1.0 + n_c_over_n_d)
            else:
                share_umi = float("nan")
            rows.append(dict(library=lib, cell_type=ct, tpm_c=tpm_c, tpm_d=tpm_d,
                             share_tpm=share_tpm, share_umi=share_umi))

    # per-library CAF vs epithelial comparison on both scales
    by_lib = {}
    for r in rows:
        by_lib.setdefault(r["library"], {})[r["cell_type"]] = r

    print(f"{'library':16s} {'CAF tpm':>8s} {'epi tpm':>8s} {'d_tpm':>7s} "
          f"{'CAF umi':>8s} {'epi umi':>8s} {'d_umi':>7s}")
    res = []
    for lib, d in by_lib.items():
        if "CAF" not in d or "epithelial" not in d:
            continue
        a, b = d["CAF"], d["epithelial"]
        dt = a["share_tpm"] - b["share_tpm"]
        du = a["share_umi"] - b["share_umi"]
        ok = (dt > 0) == (du > 0)
        print(f"{lib:16s} {a['share_tpm']:8.3f} {b['share_tpm']:8.3f} {dt:+7.3f} "
              f"{a['share_umi']:8.3f} {b['share_umi']:8.3f} {du:+7.3f}"
              f"{'' if ok else '   <-- SCALES DISAGREE'}")
        res.append(dict(library=lib, share_tpm_CAF=a["share_tpm"],
                        share_tpm_epi=b["share_tpm"], d_tpm=dt,
                        share_umi_CAF=a["share_umi"], share_umi_epi=b["share_umi"],
                        d_umi=du, scales_agree=ok))

    n_lib = len(res)
    n_tpm = sum(1 for r in res if r["d_tpm"] > 0)
    n_umi = sum(1 for r in res if r["d_umi"] > 0)
    n_both = sum(1 for r in res if r["d_tpm"] > 0 and r["d_umi"] > 0)
    print(f"\nlibraries: {n_lib}")
    print(f"CAF > epithelial on TPM scale: {n_tpm}/{n_lib}")
    print(f"CAF > epithelial on UMI scale: {n_umi}/{n_lib}")
    print(f"agreement on both scales:      {n_both}/{n_lib}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(dict(
        note="1D share = 1D/(1C+1D); UMI scale inverted from TPM via gene lengths",
        transcript_lengths={TX_C: lc, TX_D: ld},
        n_libraries=n_lib, caf_gt_epi_tpm=n_tpm, caf_gt_epi_umi=n_umi,
        agree_both=n_both, per_library=res), indent=1), encoding="utf-8")
    print("written", OUT)


if __name__ == "__main__":
    main()
