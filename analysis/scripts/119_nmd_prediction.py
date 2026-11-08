"""NMD prediction for the two FBLN1 isoforms, from the GENCODE v44 GTF.

The 50-nt rule (Longman et al.; reviewed in Lykke-Andersen & Jensen 2015):
a stop codon more than 50 nt upstream of the LAST exon-exon junction marks the
transcript as an NMD substrate. Position of the stop within the final exon is
what decides it, so this needs exon-level annotation, not just sequences.

Output: JSON with, per transcript, the exon list, the terminal exon length,
the number of coding nucleotides in the terminal exon, and the distance from the
stop codon to the last junction.
"""
from __future__ import annotations
import gzip
import json
import re
from pathlib import Path

GTF = Path("/root/as_psi/ref/gencode.v44.annotation.gtf")
OUT = Path("/root/hnscc_caf/results/fbln1_nmd.json")
TX = {"FBLN1-202_1D": "ENST00000327858", "FBLN1-201_1C": "ENST00000262722"}


def parse_attrs(s: str) -> dict:
    out = {}
    for m in re.finditer(r'(\w+) "([^"]*)"', s):
        out[m.group(1)] = m.group(2)
    return out


def main() -> None:
    assert GTF.exists(), f"missing {GTF}"
    want = set(TX.values())
    tx = {t: {"exons": [], "cds": [], "strand": None, "chrom": None,
              "tx_start": None, "tx_end": None} for t in want}

    op = gzip.open if GTF.suffix == ".gz" else open
    with op(GTF, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "CDS", "transcript"):
                continue
            a = parse_attrs(f[8])
            tid = a.get("transcript_id", "").split(".")[0]
            if tid not in want:
                continue
            s, e = int(f[3]), int(f[4])
            rec = tx[tid]
            rec["chrom"] = f[0]
            rec["strand"] = f[6]
            if f[2] == "exon":
                rec["exons"].append((s, e))
            elif f[2] == "CDS":
                rec["cds"].append((s, e))
            else:
                rec["tx_start"], rec["tx_end"] = s, e

    out = {"rule": "stop codon >50 nt upstream of the last exon-exon junction "
                   "marks an NMD substrate (50-nt rule)",
           "source": str(GTF), "transcripts": {}}

    for label, tid in TX.items():
        r = tx[tid]
        if not r["exons"]:
            raise KeyError(f"{tid} ({label}): no exons parsed — check the transcript id")
        strand = r["strand"]
        exons = sorted(r["exons"], key=lambda x: x[0], reverse=(strand == "-"))
        cds = sorted(r["cds"], key=lambda x: x[0], reverse=(strand == "-"))
        # transcript order: 5'->3'
        cds_span = (min(c[0] for c in cds), max(c[1] for c in cds))
        stop = cds_span[0] - 1 if strand == "-" else cds_span[1] + 1

        # last exon-exon junction in transcript order = boundary between the
        # penultimate and the terminal exon
        term = exons[-1]
        penult = exons[-2] if len(exons) > 1 else None
        if strand == "-":
            # terminal exon is the leftmost; junction is at its 3' end (higher coord)
            junction = term[1]
            dist = abs(junction - stop)
        else:
            junction = term[0]
            dist = abs(stop - junction)

        term_len = term[1] - term[0] + 1
        # coding nt of the terminal exon
        if strand == "-":
            cds_in_term = max(0, min(term[1], cds_span[1]) - max(term[0], cds_span[0]) + 1)
        else:
            cds_in_term = max(0, min(term[1], cds_span[1]) - max(term[0], cds_span[0]) + 1)

        coding_nt = sum(c[1] - c[0] + 1 for c in cds)
        out["transcripts"][label] = {
            "transcript_id": tid,
            "chrom": r["chrom"], "strand": strand,
            "n_exons": len(exons),
            "n_coding_exons": len(cds),
            "coding_nt": coding_nt,
            # GENCODE CDS spans exclude the stop codon, so coding_nt/3 is the
            # protein length exactly (verified: 2109/3 = 703 aa for 1D)
            "protein_aa_from_cds": coding_nt // 3,
            "terminal_exon": [term[0], term[1]],
            "terminal_exon_len": term_len,
            "coding_nt_in_terminal_exon": cds_in_term,
            "stop_codon_start": stop,
            "last_junction": junction,
            "stop_to_last_junction_nt": dist,
            "utr3_nt_in_terminal_exon": term_len - cds_in_term,
            "nmd_substrate_50nt_rule": bool(dist > 50),
            "exons_transcript_order": exons,
        }
        print(f"{label}: {len(exons)} exons, protein {coding_nt // 3} aa, "
              f"3'UTR in terminal exon {term_len - cds_in_term} nt, "
              f"stop {dist} nt from the last junction -> "
              f"{'NMD substrate' if dist > 50 else 'not an NMD substrate'}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("written", OUT)


if __name__ == "__main__":
    main()
