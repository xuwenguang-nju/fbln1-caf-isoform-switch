#!/usr/bin/env bash
# Step 3: Salmon quantification of all 36 libraries, plus a 10%-subsample pass that serves as
# the depth-sensitivity gate.
#
# Why the subsample pass is part of the pipeline rather than an optional check: on the local
# machine three panel genes (DCN, IGFBP7, MYL12B) flipped the sign of their CAF-NF isoform
# difference between 0.6 M reads and the deeper slice pass. A candidate that cannot survive two
# depths must never reach a figure, so both depths are quantified here and compared in step 5.
#
# Usage:  bash 03_quant.sh <workdir> <fastq_dir> <manifest.tsv> [threads]
# Output: <workdir>/quant_full/<run>/quant.sf
#         <workdir>/quant_sub10/<run>/quant.sf
#         <workdir>/quant_run_log.tsv
set -euo pipefail
WORK=${1:?usage: 03_quant.sh <workdir> <fastq_dir> <manifest.tsv> [threads]}
FASTQ_DIR=${2:?usage: 03_quant.sh <workdir> <fastq_dir> <manifest.tsv> [threads]}
MANIFEST=${3:?usage: 03_quant.sh <workdir> <fastq_dir> <manifest.tsv> [threads]}
THREADS=${4:-16}
IDX="$WORK/ref/salmon_index"
LOG="$WORK/quant_run_log.tsv"
SUBSAMPLE_FRAC=0.1
SEED=42

[ -s "$IDX/info.json" ] || { echo "missing index $IDX (run 02_salmon_index.sh)"; exit 1; }
mkdir -p "$WORK/sub10"
[ -f "$LOG" ] || printf 'run\tpatient\tcondition\treplicate\tmode\tseconds\tstatus\n' > "$LOG"

meta_for() {  # meta_for <run> -> patient<TAB>condition<TAB>replicate
  awk -F'\t' -v r="$1" '$1==r{printf "%s\t%s\t%s", $2, $3, $4}' "$MANIFEST"
}

run_salmon() {  # run_salmon <r1> <r2> <outdir> <logfile>
  salmon quant \
    --index "$IDX" \
    --libType A \
    -1 "$1" -2 "$2" \
    --threads "$THREADS" \
    --gcBias --seqBias --posBias \
    --minAssignedFrags 10 \
    --output "$3" > "$4" 2>&1
}

subsample() {  # subsample <in.fastq.gz> <out.fastq.gz>
  if command -v seqtk >/dev/null 2>&1; then
    seqtk sample -s "$SEED" "$1" "$SUBSAMPLE_FRAC" | gzip -1 > "$2"
  else
    # deterministic fallback: keep read records whose index is divisible by 10
    zcat "$1" | awk 'BEGIN{n=0; buf=""} {buf = buf $0 "\n"; if (NR % 4 == 0) {n++; if (n % 10 == 0) printf "%s", buf; buf=""}}' | gzip -1 > "$2"
  fi
}

for r1 in "$FASTQ_DIR"/*_1.fastq.gz; do
  run=$(basename "$r1" _1.fastq.gz)
  r2="$FASTQ_DIR/${run}_2.fastq.gz"
  [ -s "$r2" ] || { echo "missing mate2 for $run"; continue; }
  IFS=$'\t' read -r patient condition replicate <<< "$(meta_for "$run")"

  outf="$WORK/quant_full/$run"
  if [ -s "$outf/quant.sf" ]; then
    echo "[$run] full: cached"
  else
    t0=$(date +%s)
    run_salmon "$r1" "$r2" "$outf" "$outf.log"
    secs=$(( $(date +%s) - t0 ))
    printf '%s\t%s\t%s\t%s\tfull\t%s\tok\n' "$run" "$patient" "$condition" "$replicate" "$secs" >> "$LOG"
    echo "[$run] full: ${secs}s  ($condition HT$patient rep$replicate)"
  fi

  outs="$WORK/quant_sub10/$run"
  if [ -s "$outs/quant.sf" ]; then
    echo "[$run] sub10: cached"
  else
    s1="$WORK/sub10/${run}_1.sub.fastq.gz"; s2="$WORK/sub10/${run}_2.sub.fastq.gz"
    [ -s "$s1" ] || subsample "$r1" "$s1"
    [ -s "$s2" ] || subsample "$r2" "$s2"
    t0=$(date +%s)
    run_salmon "$s1" "$s2" "$outs" "$outs.log"
    secs=$(( $(date +%s) - t0 ))
    printf '%s\t%s\t%s\t%s\tsub10\t%s\tok\n' "$run" "$patient" "$condition" "$replicate" "$secs" >> "$LOG"
    echo "[$run] sub10: ${secs}s"
    rm -f "$s1" "$s2"          # 15 GB of subsampled fastq is not worth keeping
  fi
done

echo
echo "=== summary ==="
echo "full : $(find "$WORK/quant_full" -name quant.sf | wc -l) / 36"
echo "sub10: $(find "$WORK/quant_sub10" -name quant.sf | wc -l) / 36"
find "$WORK/quant_full"  -name quant.sf | wc -l | grep -qx 36 \
  && echo "QUANT COMPLETE" || { echo "INCOMPLETE - rerun this script"; exit 1; }
