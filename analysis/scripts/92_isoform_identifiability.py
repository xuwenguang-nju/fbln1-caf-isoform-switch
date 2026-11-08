"""Measure short-read isoform identifiability in this project's own annotation.

Why
---
The Methods text currently cites a sibling-session figure: only 57.5% of
multi-isoform genes (2 301 / 4 000 sampled) have >= 2 transcripts distinguishable
by unique in-gene k-mers, and only 8 103 / 27 856 transcripts qualify. That number
bounds what "isoform switch" can mean on short reads, so it should be measured
here rather than quoted.

Criterion (matching the sibling run so the numbers are comparable)
-----------------------------------------------------------------
  k = 31 (salmon default used for the B arm)
  a transcript is distinguishable if >= 10% of its k-mers are unique within its gene
  a gene is identifiable if >= 2 of its transcripts are distinguishable

Implementation
--------------
Rolling polynomial hash mod 2^59-1 over base-4 encoded sequence, vectorised with
numpy. Chosen over Python string sets because 89 067 transcripts x ~2 900 bp is
~2.6e8 k-mers; the expected collision count in a 2^59 space is ~0.06, so
collisions cannot meaningfully bias the fractions.

Extrapolating a k-mer hash is exact for set membership up to collisions, so no
approximation beyond that.

Outputs data/processed/isoform_identifiability.json and
results/tables/isoform_identifiability_genes.tsv (genes with >=2 transcripts).
"""
from __future__ import annotations
from pathlib import Path
import gzip
import json
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PROC = ROOT / "data/processed"
T = ROOT / "results/tables"

K = 31
MOD = (1 << 59) - 1
UNIQUE_FRAC = 0.10
SAMPLE_GENES = 4000
SEED = 42

# Optional argv[1] = FASTA path. The default is GENCODE's official
# *.pc_transcripts.fa.gz, which holds ALL transcripts of protein-coding genes
# (~111k). The salmon index, however, was built from a biotype-filtered file
# (~89k, biotype == protein_coding). The identifiability bound must be computed
# on the same set that was quantified, so pass that file explicitly.
FA = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    ROOT / "data/raw/reference/gencode.v44.pc_transcripts.fa.gz"

_BASE = np.full(256, -1, dtype=np.int64)
for c, v in zip(b"ACGT", range(4)):
    _BASE[c] = v
for c, v in zip(b"acgt", range(4)):
    _BASE[c] = v


def load_transcripts(path, biotype=None):
    """yield (tx_id, gene_id, seq_bytes) from a GENCODE FASTA

    biotype: if given, keep only records whose 8th pipe field equals it
    (GENCODE pc_transcripts files annotate it as `protein_coding`).
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    tx, gene, buf = None, None, []
    with opener(path, "rb") as fh:
        for line in fh:
            if line[:1] == b">":
                if tx is not None:
                    yield tx, gene, b"".join(buf)
                head = line[1:].decode().strip()
                parts = head.split("|")
                tx = parts[0].split(".")[0]
                gene = None
                for p in parts:
                    if p.startswith("ENSG"):
                        gene = p.split(".")[0]
                        break
                buf = []
                if biotype is not None and len(parts) > 7 and parts[7] != biotype:
                    tx = None            # skip this record
            elif tx is not None:
                buf.append(line.strip())
    if tx is not None:
        yield tx, gene, b"".join(buf)


def kmer_hashes(seq: bytes, k: int = K):
    """Sorted unique rolling-hash k-mers. Returns None if the sequence is too short."""
    if len(seq) < k:
        return None
    codes = _BASE[np.frombuffer(seq, dtype=np.uint8)]
    if (codes < 0).any():
        return None                      # ambiguous base
    n = len(codes) - k + 1
    h = np.zeros(n, dtype=np.int64)
    for i in range(k):
        h = (h * 4 + codes[i:i + n]) % MOD
    return np.unique(h)


def main():
    biotype = sys.argv[2] if len(sys.argv) > 2 else None
    print(f"FASTA    : {FA}")
    print(f"biotype  : {biotype or '(no filter)'}")

    # ---- group transcripts by gene ----------------------------------------
    genes = {}
    n_tx_total = 0
    for tx, g, seq in load_transcripts(FA, biotype):
        if g is None:
            continue
        genes.setdefault(g, []).append((tx, seq))
        n_tx_total += 1
    multi = {g: v for g, v in genes.items() if len(v) >= 2}
    print(f"genes: {len(genes)}  transcripts: {n_tx_total}  multi-isoform genes: {len(multi)}")

    # ---- per gene: which transcripts are distinguishable -------------------
    rows = []
    for gi, (g, txs) in enumerate(sorted(multi.items()), 1):
        hashes = {}
        for tx, seq in txs:
            h = kmer_hashes(seq)
            if h is not None and len(h):
                hashes[tx] = h
        if len(hashes) < 2:
            continue
        # union of all k-mers in the gene, and per-transcript membership counts
        allh = np.unique(np.concatenate(list(hashes.values())))
        distinguishable = []
        for tx, h in hashes.items():
            others = [v for k2, v in hashes.items() if k2 != tx]
            other = np.unique(np.concatenate(others)) if others else np.empty(0, np.int64)
            # k-mers present in this tx but in no other tx of the gene
            uniq = h[~np.isin(h, other, assume_unique=True)]
            frac = len(uniq) / len(h)
            distinguishable.append((tx, frac, len(h)))
        n_ok = sum(1 for _, f, _ in distinguishable if f >= UNIQUE_FRAC)
        rows.append(dict(gene=g, n_tx=len(txs), n_scored=len(hashes),
                         n_distinguishable=n_ok, identifiable=int(n_ok >= 2),
                         fracs=";".join(f"{t}:{f:.4f}" for t, f, _ in distinguishable)))
        if gi % 2000 == 0:
            print(f"  {gi}/{len(multi)} genes ...")

    import pandas as pd
    df = pd.DataFrame(rows)
    print(f"\nscored genes: {len(df)}  transcripts scored: {int(df['n_scored'].sum())}")

    # ---- headline numbers, full population -------------------------------
    n_gene = len(df)
    n_ident = int(df["identifiable"].sum())
    n_tx = int(df["n_scored"].sum())
    # transcripts that themselves pass the 10% bar
    tx_ok = 0
    for s in df["fracs"]:
        for item in s.split(";"):
            if item and float(item.split(":")[1]) >= UNIQUE_FRAC:
                tx_ok += 1

    print("\n=== full population ===")
    print(f"  identifiable genes : {n_ident:,}/{n_gene:,} = {n_ident/n_gene:.1%}")
    print(f"  distinguishable tx : {tx_ok:,}/{n_tx:,} = {tx_ok/n_tx:.1%}")

    # ---- the sibling's 4 000-gene sample, same criterion ------------------
    rng = np.random.default_rng(SEED)
    if len(df) > SAMPLE_GENES:
        idx = rng.choice(len(df), SAMPLE_GENES, replace=False)
        s = df.iloc[idx]
    else:
        s = df
    s_ident = int(s["identifiable"].sum())
    s_tx = int(s["n_scored"].sum())
    s_ok = 0
    for v in s["fracs"]:
        for item in v.split(";"):
            if item and float(item.split(":")[1]) >= UNIQUE_FRAC:
                s_ok += 1
    print(f"\n=== random {SAMPLE_GENES}-gene sample (seed {SEED}) ===")
    print(f"  identifiable genes : {s_ident:,}/{len(s):,} = {s_ident/len(s):.1%}")
    print(f"  distinguishable tx : {s_ok:,}/{s_tx:,} = {s_ok/s_tx:.1%}")

    out = dict(
        fasta=str(FA), biotype_filter=biotype,
        k=K, unique_kmer_frac_threshold=UNIQUE_FRAC,
        n_genes_total=len(genes), n_transcripts_total=n_tx_total,
        n_multi_isoform_genes=len(multi), n_genes_scored=n_gene,
        n_transcripts_scored=n_tx,
        n_identifiable_genes=n_ident,
        frac_identifiable_genes=round(n_ident / n_gene, 4),
        n_distinguishable_transcripts=tx_ok,
        frac_distinguishable_transcripts=round(tx_ok / n_tx, 4),
        sample=dict(n=len(s), seed=SEED,
                    n_identifiable=int(s_ident),
                    frac_identifiable=round(s_ident / len(s), 4),
                    n_distinguishable=s_ok,
                    frac_distinguishable=round(s_ok / s_tx, 4)),
        note="hash-based k-mer sets, mod 2^59-1; expected collisions ~0.06",
    )
    PROC.mkdir(parents=True, exist_ok=True)
    (PROC / "isoform_identifiability.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8")
    df.to_csv(T / "isoform_identifiability_genes.tsv", sep="\t", index=False)
    print(f"\nwrote {PROC/'isoform_identifiability.json'} and "
          f"{T/'isoform_identifiability_genes.tsv'}")


if __name__ == "__main__":
    main()
