#!/usr/bin/env bash
# Driver for step 1 on the analysis node (runs headless; re-runnable).
#
#   phase 1  Salmon index over GENCODE v44 transcripts (k=31, --gencode, no decoys:
#            these are sorted fibroblasts, not tumour tissue, so decoys buy little)
#   phase 2  download 36 libraries x 2 mates with aria2c, md5-verified per file
#   phase 3  salmon quant as soon as both mates of a library verify (full depth),
#            then a 10% subsample pass that serves as the depth-consistency gate
#
# Download and quantification run as two concurrent loops because they are I/O- and
# CPU-bound respectively, so the first quant results appear hours before the last file lands.
#
# Concurrency numbers were measured on this node (see README): one aria2c file at 16
# connections sustains ~1.3 MB/s and three files reach ~2.7 MB/s, so the pool matters.
#
# Usage:  bash run_step1.sh          (env: DOWN_CONNS, MAX_DL_JOBS, QUANT_THREADS, MAX_QUANT_JOBS)
set -uo pipefail
WORK=/root/hnscc_caf
SCRIPTS="$WORK/remote_step1"
FASTQ="$WORK/fastq"
LOG="$WORK/run_step1.log"
MANIFEST="$SCRIPTS/gse303690_ena_manifest.tsv"
GENCODE_REL=44
FA="gencode.v$GENCODE_REL.transcripts.fa.gz"
IDX="$WORK/ref/salmon_index"
DOWN_CONNS=${DOWN_CONNS:-10}
MAX_DL_JOBS=${MAX_DL_JOBS:-4}
QUANT_THREADS=${QUANT_THREADS:-14}
MAX_QUANT_JOBS=${MAX_QUANT_JOBS:-2}

mkdir -p "$FASTQ" "$WORK/ref" "$WORK/results" "$WORK/logs"
exec > >(tee -a "$LOG") 2>&1
ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { printf '[%s] %s\n' "$(ts)" "$*"; }

log "=== step 1 start (workdir $WORK) ==="
log "disk: $(df -h "$WORK" | awk 'NR==2{print $4" free"}')"

# ---------------------------------------------------------------- phase 1: index
if [ ! -s "$IDX/info.json" ]; then
  if [ ! -s "$WORK/ref/$FA" ]; then
    log "downloading GENCODE v$GENCODE_REL transcripts (EBI throttles this host to ~70 KB/s)"
    aria2c -x 8 -s 8 --file-allocation=none --console-log-level=warn -d "$WORK/ref" -o "$FA" \
      "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_REL}/${FA}" \
      || { log "FATAL: gencode download failed"; exit 1; }
  fi
  log "building salmon index"
  salmon index --transcripts "$WORK/ref/$FA" --index "$IDX" --kmerLen 31 --threads 30 --gencode \
    || { log "FATAL: salmon index failed"; exit 1; }
  { echo "annotation: GENCODE v$GENCODE_REL transcripts"
    echo "salmon: $(salmon --version | head -1)"
    echo "index: salmon index --kmerLen 31 --gencode (no decoys)"
    echo "built_utc: $(ts)"; } > "$IDX/index_meta.txt"
  log "index ready"
else
  log "index already present"
fi

# ---------------------------------------------------------------- helpers
verified() {  # verified <run> <mate> <bytes> <md5>
  local f="$FASTQ/$1_$2.fastq.gz"
  [ -s "$f" ] || return 1
  [ "$(stat -c%s "$f")" = "$3" ] || return 1
  [ "$(md5sum "$f" | cut -d' ' -f1)" = "$4" ] || return 1
}

fetch_file() {
  aria2c -x "$DOWN_CONNS" -s "$DOWN_CONNS" -k 1M --file-allocation=none \
         --max-tries=5 --retry-wait=5 --console-log-level=warn --summary-interval=0 \
         --auto-file-renaming=false --allow-overwrite=false \
         -d "$(dirname "$2")" -o "$(basename "$2")" "$1"
}

download_one() {  # <run> <url1> <md5_1> <bytes_1> <url2> <md5_2> <bytes_2>
  local run=$1 u1=$2 m1=$3 b1=$4 u2=$5 m2=$6 b2=$7
  verified "$run" 1 "$b1" "$m1" || { log "  [$run] mate1 downloading"; fetch_file "$u1" "$FASTQ/${run}_1.fastq.gz" || return 1; }
  verified "$run" 1 "$b1" "$m1" || { log "  [$run] mate1 FAILED verification"; return 1; }
  verified "$run" 2 "$b2" "$m2" || { log "  [$run] mate2 downloading"; fetch_file "$u2" "$FASTQ/${run}_2.fastq.gz" || return 1; }
  verified "$run" 2 "$b2" "$m2" || { log "  [$run] mate2 FAILED verification"; return 1; }
  return 0
}

# ---------------------------------------------------------------- quant loop
quant_worker() {
  while true; do
    local pending=0
    while IFS=$'\t' read -r run patient condition replicate _ _ _ _ _ b1 _ _ b2; do
      [ "$run" = "run" ] && continue
      [ -s "$WORK/quant_full/$run/quant.sf" ] && continue
      [ -s "$FASTQ/${run}_1.fastq.gz" ] && [ -s "$FASTQ/${run}_2.fastq.gz" ] || continue
      [ "$(stat -c%s "$FASTQ/${run}_1.fastq.gz")" = "$b1" ] || continue
      [ "$(stat -c%s "$FASTQ/${run}_2.fastq.gz")" = "$b2" ] || continue
      while [ "$(jobs -rp | wc -l)" -ge "$MAX_QUANT_JOBS" ]; do sleep 20; done
      pending=1
      (
        t0=$(date +%s)
        salmon quant --index "$IDX" --libType A \
          -1 "$FASTQ/${run}_1.fastq.gz" -2 "$FASTQ/${run}_2.fastq.gz" \
          --threads "$QUANT_THREADS" --gcBias --seqBias --posBias --minAssignedFrags 10 \
          --output "$WORK/quant_full/$run" > "$WORK/logs/${run}_full.log" 2>&1
        rc=$?
        secs=$(( $(date +%s) - t0 ))
        printf '%s\t%s\t%s\t%s\tfull\t%s\t%s\n' "$run" "$patient" "$condition" "$replicate" "$secs" "$rc" >> "$WORK/quant_run_log.tsv"
        log "  [$run] salmon full rc=$rc ${secs}s ($condition HT$patient rep$replicate)"
        if [ "$rc" = "0" ] && [ ! -s "$WORK/quant_sub10/$run/quant.sf" ]; then
          s1="$WORK/sub10/${run}_1.sub.fastq.gz"; s2="$WORK/sub10/${run}_2.sub.fastq.gz"
          mkdir -p "$WORK/sub10"
          [ -s "$s1" ] || seqtk sample -s 42 "$FASTQ/${run}_1.fastq.gz" 0.1 | pigz -p 4 > "$s1"
          [ -s "$s2" ] || seqtk sample -s 42 "$FASTQ/${run}_2.fastq.gz" 0.1 | pigz -p 4 > "$s2"
          salmon quant --index "$IDX" --libType A -1 "$s1" -2 "$s2" \
            --threads "$QUANT_THREADS" --gcBias --seqBias --posBias --minAssignedFrags 10 \
            --output "$WORK/quant_sub10/$run" > "$WORK/logs/${run}_sub10.log" 2>&1
          rc2=$?
          printf '%s\t%s\t%s\t%s\tsub10\t-\t%s\n' "$run" "$patient" "$condition" "$replicate" "$rc2" >> "$WORK/quant_run_log.tsv"
          log "  [$run] salmon sub10 rc=$rc2"
          rm -f "$s1" "$s2"
        fi
      ) &
      sleep 5
    done < "$MANIFEST"
    n_done=$(find "$WORK/quant_full" -name quant.sf 2>/dev/null | wc -l)
    [ "$n_done" -ge 36 ] && break
    [ -f "$WORK/.download_finished" ] && [ "$pending" = "0" ] && break
    sleep 60
  done
  wait
}

# ---------------------------------------------------------------- download loop
download_worker() {
  local pass
  for pass in 1 2 3; do
    while IFS=$'\t' read -r run patient condition replicate _ _ _ u1 m1 b1 u2 m2 b2; do
      [ "$run" = "run" ] && continue
      [ "$replicate" = "$pass" ] || continue
      if [ -s "$WORK/quant_full/$run/quant.sf" ]; then
        log "  [$run] already quantified, skipping download"; continue
      fi
      while [ "$(jobs -rp | wc -l)" -ge "$MAX_DL_JOBS" ]; do sleep 5; done
      (
        if download_one "$run" "$u1" "$m1" "$b1" "$u2" "$m2" "$b2"; then
          log "  [$run] DONE $condition HT$patient rep$replicate | disk $(df -h "$WORK" | awk 'NR==2{print $4}')"
        else
          log "  [$run] FAILED $condition HT$patient rep$replicate (retried next pass)"
        fi
      ) &
      sleep 2
    done < "$MANIFEST"
    wait
  done
  touch "$WORK/.download_finished"
  log "download loop finished"
}

quant_worker &
Q_PID=$!
download_worker
wait "$Q_PID" 2>/dev/null
touch "$WORK/.download_finished"

log "=== step 1 finished ==="
log "fastq files: $(ls "$FASTQ"/*.fastq.gz 2>/dev/null | wc -l) / 72"
log "quant full:  $(find "$WORK/quant_full" -name quant.sf 2>/dev/null | wc -l) / 36"
log "quant sub10: $(find "$WORK/quant_sub10" -name quant.sf 2>/dev/null | wc -l) / 36"
