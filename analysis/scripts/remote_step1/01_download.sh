#!/usr/bin/env bash
# Step 1: fetch the 36 GSE303690 paired-end libraries (152.3 GB) and verify checksums.
#
# Why a manifest instead of an API call: the ENA run table was frozen locally
# (gse303690_ena_manifest.tsv, generated 2026-10-18 from the ENA Portal API), so this step
# needs no metadata call and the file list stays reproducible.
#
# Files are downloaded one pair at a time and verified by exact size then md5, so an
# interrupted run can simply be restarted (finished files are re-verified and skipped).
#
# Usage:  bash 01_download.sh <fastq_dir> [manifest.tsv]
# Output: <fastq_dir>/<run>_1.fastq.gz, <run>_2.fastq.gz, <fastq_dir>/download_report.tsv
set -uo pipefail
FASTQ_DIR=${1:?usage: 01_download.sh <fastq_dir> [manifest]}
MANIFEST=${2:-gse303690_ena_manifest.tsv}
mkdir -p "$FASTQ_DIR"
REPORT="$FASTQ_DIR/download_report.tsv"
[ -f "$REPORT" ] || printf 'run\tpatient\tcondition\treplicate\tfile\tbytes_expected\tbytes_actual\tstatus\n' > "$REPORT"

fetch_one() {  # fetch_one <url> <md5> <bytes> <dest> -> prints a status word
  local url=$1 md5=$2 bytes=$3 dest=$4 have got
  if [ -f "$dest" ]; then
    have=$(stat -c%s "$dest" 2>/dev/null || echo 0)
    if [ "$have" -eq "$bytes" ]; then
      got=$(md5sum "$dest" | awk '{print $1}')
      if [ "$got" = "$md5" ]; then echo "cached-ok"; return 0; fi
      rm -f "$dest"            # corrupt; download again
    fi
  fi
  local attempt
  for attempt in 1 2 3; do
    curl -L --fail --silent --show-error --retry 2 --retry-delay 5 \
         --continue-at - --output "$dest" "$url" && break
    sleep 5
  done
  have=$(stat -c%s "$dest" 2>/dev/null || echo 0)
  if [ "$have" -ne "$bytes" ]; then echo "size-mismatch($have/$bytes)"; return 1; fi
  got=$(md5sum "$dest" | awk '{print $1}')
  if [ "$got" != "$md5" ]; then echo "md5-mismatch"; return 1; fi
  echo "downloaded-ok"
}

n_ok=0
# manifest columns: run patient condition replicate layout instrument reads url_1 md5_1 bytes_1 url_2 md5_2 bytes_2
while IFS=$'\t' read -r run patient condition replicate layout instrument reads url1 md5_1 bytes1 url2 md5_2 bytes2; do
  [ "$run" = "run" ] && continue
  echo "== $run  $condition  HT$patient  rep$replicate  ($(( bytes1 / 1000000 )) MB + $(( bytes2 / 1000000 )) MB)"
  for mate in 1 2; do
    if [ "$mate" = "1" ]; then url=$url1; md5=$md5_1; bytes=$bytes1; else url=$url2; md5=$md5_2; bytes=$bytes2; fi
    [ -z "${url:-}" ] && continue
    dest="$FASTQ_DIR/${run}_${mate}.fastq.gz"
    status=$(fetch_one "$url" "$md5" "$bytes" "$dest")
    actual=$(stat -c%s "$dest" 2>/dev/null || echo 0)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$run" "$patient" "$condition" "$replicate" "${run}_${mate}.fastq.gz" \
      "$bytes" "$actual" "$status" >> "$REPORT"
    echo "   mate$mate  $status"
    case "$status" in downloaded-ok|cached-ok) n_ok=$((n_ok + 1));; esac
  done
done < "$MANIFEST"

echo
echo "=== summary ==="
awk -F'\t' 'NR>1{c[$8]++} END{for(k in c) printf "%-28s %d\n", k, c[k]}' "$REPORT"
echo "verified files: $n_ok / 72"
if [ "$n_ok" -eq 72 ]; then
  echo "DOWNLOAD COMPLETE"
  du -sh "$FASTQ_DIR"
else
  echo "INCOMPLETE - rerun this script (finished files are skipped), inspect $REPORT"
  exit 1
fi
