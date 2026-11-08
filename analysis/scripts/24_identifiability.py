"""Validation step 4b: which isoform fractions are even identifiable?

A read can only distinguish two isoforms where their sequences differ. For nested isoforms
(same exons, different start/end) the fraction of the shorter one is not identifiable from
short reads, and any estimate of it is prior-driven rather than data-driven. Isoform-switch
claims must therefore be restricted to genes with distinguishable isoforms.

This script computes, per transcript, the fraction of its 25-mers that no other isoform of
the same gene carries ("unique k-mer fraction"). Genes where fewer than two transcripts have
a substantial unique fraction are flagged non-identifiable and excluded downstream.

Inputs : data/raw/reference/gencode.v44.pc_transcripts.fa.gz
         data/processed/kmer_index/gene_tx_table.tsv
Outputs: data/processed/kmer_index/tx_identifiability.tsv
         data/processed/kmer_index/identifiability_summary.json
"""
import gzip, json
import numpy as np
import pandas as pd

FASTA = "data/raw/reference/gencode.v44.pc_transcripts.fa.gz"
IDX = "data/processed/kmer_index"
K = 25
UNIQUE_FRAC_MIN = 0.10      # a transcript needs this share of gene-unique k-mers to count as distinguishable

_BASE = np.full(256, 255, dtype=np.uint8)
for ch, v in zip("ACGT", range(4)):
    _BASE[ord(ch)] = v
    _BASE[ord(ch.lower())] = v


def kmer_set(seq, k=K):
    arr = _BASE[np.frombuffer(seq.encode(), dtype=np.uint8)]
    n = len(arr) - k + 1
    if n <= 0:
        return np.empty(0, np.uint64)
    val = np.zeros(n, dtype=np.uint64)
    for j in range(k):
        val <<= np.uint64(2)
        val += arr[j:j + n].astype(np.uint64)
    cum = np.concatenate(([0], np.cumsum((arr == 255).astype(np.int32))))
    ok = (cum[k:] - cum[:-k]) == 0
    return np.unique(val[ok])


def main():
    tab = pd.read_csv(f"{IDX}/gene_tx_table.tsv", sep="\t")
    want = set(tab["tx_id"])
    seqs, cur = {}, None
    with gzip.open(FASTA, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                tid = line[1:].split("|")[0].split(".")[0]
                cur = tid if tid in want else None
                if cur:
                    seqs[cur] = []
            elif cur:
                seqs[cur].append(line.strip())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    print("sequences:", len(seqs), flush=True)

    rows = []
    for gene, sub in tab.groupby("gene"):
        txs = [t for t in sub["tx_id"] if t in seqs]
        if len(txs) < 2:
            for t in txs:
                rows.append({"gene": gene, "tx_id": t, "n_tx_in_gene": len(txs),
                             "unique_kmer_frac": 1.0, "is_distinguishable": True})
            continue
        sets = {t: set(kmer_set(seqs[t]).tolist()) for t in txs}
        union_counts = {}
        for t in txs:
            for km in sets[t]:
                union_counts[km] = union_counts.get(km, 0) + 1
        for t in txs:
            own = sets[t]
            uniq = sum(1 for km in own if union_counts[km] == 1)
            frac = uniq / max(len(own), 1)
            rows.append({"gene": gene, "tx_id": t, "n_tx_in_gene": len(txs),
                         "unique_kmer_frac": frac,
                         "is_distinguishable": bool(frac >= UNIQUE_FRAC_MIN)})
    out = pd.DataFrame(rows)
    out.to_csv(f"{IDX}/tx_identifiability.tsv", sep="\t", index=False)

    gene_level = out.groupby("gene").agg(n_tx=("tx_id", "size"),
                                         n_distinguishable=("is_distinguishable", "sum")) \
        .reset_index()
    gene_level["gene_is_identifiable"] = gene_level["n_distinguishable"] >= 2
    gene_level.to_csv(f"{IDX}/gene_identifiability.tsv", sep="\t", index=False)
    summary = {
        "k": K, "unique_frac_min": UNIQUE_FRAC_MIN,
        "genes": int(len(gene_level)),
        "genes_identifiable": int(gene_level["gene_is_identifiable"].sum()),
        "frac_genes_identifiable": float(gene_level["gene_is_identifiable"].mean()),
        "transcripts": int(len(out)),
        "transcripts_distinguishable": int(out["is_distinguishable"].sum()),
    }
    json.dump(summary, open(f"{IDX}/identifiability_summary.json", "w", encoding="utf-8"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
