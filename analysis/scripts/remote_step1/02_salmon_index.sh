#!/usr/bin/env bash
# Step 2: Salmon transcript index (GENCODE v44).
#
# Version pinning matters here: every downstream isoform fraction is relative to this
# annotation, so the release and the index parameters are recorded in index_meta.txt next to
# the index and must be quoted in Methods.
#
# Usage:  bash 02_salmon_index.sh <workdir> [threads]
# Output: <workdir>/ref/gencode.v44.transcripts.fa.gz
#         <workdir>/ref/salmon_index/  (+ index_meta.txt)
set -euo pipefail
WORK=${1:?usage: 02_salmon_index.sh <workdir> [threads]}
THREADS=${2:-16}
REF="$WORK/ref"
IDX="$REF/salmon_index"
GENCODE_REL=44
FA_GZ="gencode.v44.transcripts.fa.gz"
mkdir -p "$REF"

if [ ! -s "$REF/$FA_GZ" ]; then
  echo "== downloading GENCODE v$GENCODE_REL transcripts"
  curl -L --fail --retry 3 -o "$REF/$FA_GZ" \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_REL}/${FA_GZ}"
fi
ls -l "$REF/$FA_GZ"

if [ -s "$IDX/info.json" ]; then
  echo "== index already present: $IDX"
else
  echo "== building index (k=31, no decoys: these are sorted fibroblasts, not a tumour genome)"
  salmon index \
    --transcripts "$REF/$FA_GZ" \
    --index "$IDX" \
    --kmerLen 31 \
    --threads "$THREADS" \
    --gencode
fi

{
  echo "annotation: GENCODE v${GENCODE_REL} (transcripts, all biotypes)"
  echo "source: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_REL}/${FA_GZ}"
  echo "salmon: $(salmon --version 2>&1 | head -1)"
  echo "index_build: salmon index --kmerLen 31 --gencode"
  echo "built_utc: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host: $(hostname)"
} | tee "$IDX/index_meta.txt"

echo "INDEX READY: $IDX"
