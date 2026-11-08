"""Physicochemical properties of the divergent FBLN1 C-termini.

Everything here is computed from the GENCODE v44 protein sequences, not taken
from a database, so the numbers cannot drift from the annotation.

Reported per segment:
  * length, GRAVY (Kyte-Doolittle), isoelectric point, net charge at pH 7.0
  * Chou-Fasman helix/sheet/turn propensities (dimensionless; >1 favours)
  * a 9-residue-window hydrophobicity profile, used by the figure

Segments compared: the full proteins, the shared N-terminal part (1-566), and
the two divergent C-terminal segments (1D 567-703, 1C 567-683).
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SEQ = ROOT / "remote/easyai/fbln1_sequences.json"
STR = ROOT / "remote/easyai/fbln1_structure.json"
OUT = ROOT / "remote/easyai/fbln1_physchem.json"

# Kyte & Doolittle (1982)
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2,
          I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9,
          Y=-1.3, V=4.2)
# Chou & Fasman (1978) conformational parameters
CF = {
    "H": dict(A=1.42, R=0.98, N=0.67, D=1.01, C=0.70, Q=1.11, E=1.51, G=0.57,
              H=1.00, I=1.08, L=1.21, K=1.16, M=1.45, F=1.13, P=0.57, S=0.77,
              T=0.83, W=1.08, Y=0.69, V=1.06),
    "E": dict(A=0.83, R=0.93, N=0.89, D=0.54, C=1.19, Q=1.10, E=0.37, G=0.75,
              H=0.87, I=1.60, L=1.30, K=0.74, M=1.05, F=1.38, P=0.55, S=0.75,
              T=1.19, W=1.37, Y=1.47, V=1.70),
    "T": dict(A=0.66, R=0.95, N=1.56, D=1.46, C=1.19, Q=0.98, E=0.74, G=1.56,
              H=0.95, I=0.47, L=0.59, K=1.01, M=0.60, F=0.60, P=1.52, S=1.43,
              T=0.96, W=0.96, Y=1.14, V=0.50),
}
# Emboss/ExPASy pKa set
PKA = dict(Cterm=3.55, Nterm=7.50, D=4.05, E=4.45, C=9.00, Y=10.00,
           H=5.98, K=10.00, R=12.00)


def gravy(seq: str) -> float:
    return float(np.mean([KD[a] for a in seq]))


def net_charge(seq: str, ph: float = 7.0) -> float:
    q = 0.0
    # terminal charges
    q += -1.0 / (1.0 + 10 ** (PKA["Cterm"] - ph))
    q += 1.0 / (1.0 + 10 ** (ph - PKA["Nterm"]))
    for a in seq:
        if a in ("D", "E", "C", "Y"):
            q -= 1.0 / (1.0 + 10 ** (PKA[a] - ph))
        elif a in ("H", "K", "R"):
            q += 1.0 / (1.0 + 10 ** (ph - PKA[a]))
    return float(q)


def isoelectric(seq: str) -> float:
    lo, hi = 0.0, 14.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if net_charge(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return float((lo + hi) / 2)


def cf_propensity(seq: str) -> dict:
    return {k: float(np.mean([tab[a] for a in seq])) for k, tab in CF.items()}


def kd_profile(seq: str, w: int = 9) -> list:
    """Sliding-window mean hydrophobicity, centred."""
    if len(seq) < w:
        return []
    arr = np.array([KD[a] for a in seq])
    kern = np.ones(w) / w
    conv = np.convolve(arr, kern, mode="valid")
    return conv.tolist()


def describe(seq: str) -> dict:
    return dict(length=len(seq), gravy=gravy(seq), pI=isoelectric(seq),
                net_charge_pH7=net_charge(seq), **cf_propensity(seq))


def main() -> None:
    seqs = json.loads(SEQ.read_text(encoding="utf-8"))
    st = json.loads(STR.read_text(encoding="utf-8"))
    p1d = seqs["FBLN1-202_1D"]["protein"]
    p1c = seqs["FBLN1-201_1C"]["protein"]
    assert p1d[:566] == p1c[:566], "shared N-terminus assumption broken"
    div = st["divergence"]
    split = int(div["shared_prefix"])          # 566 shared residues -> slice index

    segs = {
        "1D full": p1d,
        "1C full": p1c,
        "shared 1-566": p1d[:split],
        "1D C-term 567-703": p1d[split:],
        "1C C-term 567-683": p1c[split:],
    }
    out = {"source": str(SEQ), "split_residue": split + 1, "segments": {}}
    for name, s in segs.items():
        d = describe(s)
        out["segments"][name] = d
        print(f"{name:20s} n={d['length']:4d}  GRAVY={d['gravy']:+.3f}  "
              f"pI={d['pI']:5.2f}  q(pH7)={d['net_charge_pH7']:+7.2f}  "
              f"<H>={d['H']:.3f} <E>={d['E']:.3f} <T>={d['T']:.3f}")

    # hydrophobicity profiles for the figure
    out["kd_profile"] = {
        "1D": {"window": 9, "values": kd_profile(p1d)},
        "1C": {"window": 9, "values": kd_profile(p1c)},
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("\nwritten", OUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
