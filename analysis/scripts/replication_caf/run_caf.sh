#!/usr/bin/env bash
# Replication of the FBLN1 C/D result in an independent, pure-CAF cohort.
#
# GSE195832 ships 16 "Sorted CAF" bulk libraries. GEO titles only say "Rep 1..16",
# but the authors' Mendeley matrix labels them by subtype S1-S4 (S1=5, S2=3, S3=3,
# S4=5). The Rep->subtype mapping is a HYPOTHESIS (matrix column order == GSM
# order) and is verified empirically after quantification by correlating our TPMs
# against that matrix; it is never assumed downstream.
#
# No paired normal fibroblasts exist in this cohort, so this cannot reproduce the
# paired fold change. It tests two things instead:
#   (1) is the C:D ratio in pure sorted CAF in the GSE303690 CAF range (~0.19)
#       rather than the NF range (~0.50)?
#   (2) does the C:D ratio differ between CAF subtypes S1-S4?
#
# Reuses the GSE303690 Salmon index; no rebuild needed.
set -uo pipefail
WORK=/root/hnscc_caf
CAF=$WORK/gse195832_caf
FASTQ=$CAF/fastq
MAN=$CAF/sorted_caf_manifest.tsv
IDX=$WORK/ref/salmon_index
LOG=$CAF/run_caf.log
LOCKS=$CAF/.locks
THREADS=${THREADS:-14}
DOWN_CONNS=${DOWN_CONNS:-8}
MAX_DL_JOBS=${MAX_DL_JOBS:-3}
# Which replicates to fetch. Pilot default spans all four hypothesised subtypes.
REPS=${REPS:-1 2 6 8}

mkdir -p "$FASTQ" "$CAF/quant" "$CAF/logs" "$CAF/sub10" "$LOCKS"
rm -rf "$LOCKS"/* 2>/dev/null || true
exec > >(tee -a "$LOG") 2>&1
ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }
log() { printf '[%s] %s\n' "$(ts)" "$*"; }
echo "$REPS" > "$CAF/.reps"

[ -s "$IDX/info.json" ] || { log "FATAL: salmon index missing at $IDX"; exit 1; }
log "=== replication start | reps: $REPS | disk $(df -h "$CAF" | awk 'NR==2{print $4}') ==="

# completion test learned the hard way: apparent size is meaningless mid-download
verified() {  # <run> <mate> <bytes> <md5>
  local run=$1 mate=$2 bytes=$3 md5=$4 f
  f="$FASTQ/${run}_${mate}.fastq.gz"
  [ -f "$f" ] || return 1
  [ -f "$f.aria2" ] && return 1
  [ "$(stat -c%s "$f")" = "$bytes" ] || return 1
  [ "$(md5sum "$f" | cut -d' ' -f1)" = "$md5" ] || return 1
  return 0
}

download_one() {  # <run> <u1> <m1> <b1> <u2> <m2> <b2>
  local run=$1 u1=$2 m1=$3 b1=$4 u2=$5 m2=$6 b2=$7 attempt mate url bytes md5
  for mate in 1 2; do
    case $mate in
      1) url=$u1; bytes=$b1; md5=$m1 ;;
      2) url=$u2; bytes=$b2; md5=$m2 ;;
    esac
    verified "$run" "$mate" "$bytes" "$md5" && continue
    for attempt in 1 2 3 4; do
      aria2c -x "$DOWN_CONNS" -s "$DOWN_CONNS" -k 1M \
        --file-allocation=none --allow-overwrite=true --auto-file-renaming=false \
        --console-log-level=warn --summary-interval=0 \
        -d "$FASTQ" -o "${run}_${mate}.fastq.gz" "$url" \
        > "$CAF/logs/dl_${run}_${mate}.log" 2>&1
      verified "$run" "$mate" "$bytes" "$md5" && break
      log "  [$run] mate$mate attempt $attempt failed, backing off 30s"
      sleep 30
    done
    verified "$run" "$mate" "$bytes" "$md5" || { log "  [$run] mate$mate FAILED"; return 1; }
  done
  return 0
}

quant_one() {  # <run> <rep> <sub>
  local run=$1 rep=$2 sub=$3 t0 rc secs s1 s2 rc2
  t0=$(date +%s)
  salmon quant --index "$IDX" --libType A \
    -1 "$FASTQ/${run}_1.fastq.gz" -2 "$FASTQ/${run}_2.fastq.gz" \
    --threads "$THREADS" --gcBias --seqBias --posBias --minAssignedFrags 10 \
    --output "$CAF/quant/$run" > "$CAF/logs/${run}_full.log" 2>&1
  rc=$?
  secs=$(( $(date +%s) - t0 ))
  printf '%s\t%s\t%s\tfull\t%s\t%s\n' "$run" "$rep" "$sub" "$secs" "$rc" >> "$CAF/quant_log.tsv"
  log "  [$run] salmon full rc=$rc ${secs}s (rep$rep $sub)"
  if [ "$rc" = "0" ] && [ ! -s "$CAF/sub10/$run/quant.sf" ]; then
    s1="$CAF/sub10/${run}_1.sub.fastq.gz"; s2="$CAF/sub10/${run}_2.sub.fastq.gz"
    [ -s "$s1" ] || seqtk sample -s 42 "$FASTQ/${run}_1.fastq.gz" 0.1 | pigz -p 4 > "$s1"
    [ -s "$s2" ] || seqtk sample -s 42 "$FASTQ/${run}_2.fastq.gz" 0.1 | pigz -p 4 > "$s2"
    salmon quant --index "$IDX" --libType A -1 "$s1" -2 "$s2" \
      --threads "$THREADS" --gcBias --seqBias --posBias --minAssignedFrags 10 \
      --output "$CAF/sub10/$run" > "$CAF/logs/${run}_sub10.log" 2>&1
    rc2=$?
    printf '%s\t%s\t%s\tsub10\t-\t%s\n' "$run" "$rep" "$sub" "$rc2" >> "$CAF/quant_log.tsv"
    log "  [$run] salmon sub10 rc=$rc2"
    rm -f "$s1" "$s2"
  fi
  rm -rf "$LOCKS/$run.lock"
}

quant_worker() {
  local pending run rep sub _rest
  local want
  want=" $(echo $REPS) "
  while true; do
    pending=0
    while IFS=$'\t' read -r run rep sub _rest; do
      [ "$run" = "run" ] && continue
      [ "${want/ $rep /}" = "$want" ] && continue     # not selected
      [ -s "$CAF/quant/$run/quant.sf" ] && continue
      [ -s "$FASTQ/${run}_1.fastq.gz" ] && [ -s "$FASTQ/${run}_2.fastq.gz" ] || continue
      mkdir "$LOCKS/$run.lock" 2>/dev/null || continue
      while [ "$(jobs -rp | wc -l)" -ge 2 ]; do sleep 15; done
      pending=1
      quant_one "$run" "$rep" "$sub" &
      sleep 5
    done < "$MAN"
    [ "$(find "$CAF/quant" -name quant.sf 2>/dev/null | wc -l)" -ge "$(echo $REPS | wc -w)" ] && break
    [ -f "$CAF/.dl_done" ] && [ "$pending" = "0" ] && break
    sleep 60
  done
  wait
  log "quant loop finished"
}

dl_worker() {
  local run rep sub u1 m1 b1 u2 m2 b2 _t _l
  local want; want=" $(echo $REPS) "
  while IFS=$'\t' read -r run rep sub _reads _layout u1 m1 b1 u2 m2 b2 _title; do
    [ "$run" = "run" ] && continue
    [ "${want/ $rep /}" = "$want" ] && continue
    [ -s "$CAF/quant/$run/quant.sf" ] && { log "  [$run] already quantified, skip"; continue; }
    while [ "$(jobs -rp | wc -l)" -ge "$MAX_DL_JOBS" ]; do sleep 5; done
    (
      if download_one "$run" "$u1" "$m1" "$b1" "$u2" "$m2" "$b2"; then
        log "  [$run] DONE rep$rep $sub | disk $(df -h "$CAF" | awk 'NR==2{print $4}')"
      else
        log "  [$run] FAILED rep$rep $sub"
      fi
    ) &
    sleep 2
  done < "$MAN"
  wait
  touch "$CAF/.dl_done"
  log "download loop finished"
}

quant_worker &
QP=$!
dl_worker
wait "$QP" 2>/dev/null
touch "$CAF/.dl_done"
N=$(echo $REPS | wc -w)
log "=== replication finished ==="
log "fastq files: $(ls "$FASTQ"/*.fastq.gz 2>/dev/null | wc -l) / $((N*2))"
log "quant:       $(find "$CAF/quant" -name quant.sf 2>/dev/null | wc -l) / $N"
log "quant sub10: $(find "$CAF/sub10" -name quant.sf 2>/dev/null | wc -l) / $N"
