#!/usr/bin/env bash
# Step 0: environment gate. Exits non-zero if anything the pipeline needs is missing, so a
# 152 GB download never starts on a node that cannot finish the job.
#
# Usage:  bash 00_env_check.sh [workdir]
# Exit:   0 ready, 1 not ready
set -uo pipefail
WORK=${1:-$PWD}
FAIL=0
WARN=0

say()  { printf '  %-22s %s\n' "$1" "$2"; }
head1() { printf '\n=== %s ===\n' "$1"; }

head1 "host"
say hostname "$(hostname)"
say kernel "$(uname -sr)"
say user "$(whoami)"
say workdir "$WORK"

head1 "resources"
CPUS=$(nproc 2>/dev/null || echo 0); say cpu_cores "$CPUS"
MEMG=$(free -g 2>/dev/null | awk '/^Mem:|^内存/{print $2}'); say mem_total_GB "${MEMG:-?}"
DISKG=$(df -BG "$WORK" 2>/dev/null | awk 'NR==2{gsub("G","",$4); print $4}')
say disk_avail_GB "${DISKG:-?}"
DISKREQ=220
[ "${DISKG:-0}" -lt "$DISKREQ" ] && { say disk "FAIL: need >= ${DISKREQ} GB (152 GB FASTQ + intermediates)"; FAIL=1; }
[ "${CPUS:-0}" -lt 8 ] && { say cpu "WARN: < 8 cores, salmon will be slow"; WARN=1; }
[ "${MEMG:-0}" -lt 30 ] && { say mem "WARN: < 30 GB, salmon index may not fit"; WARN=1; }

head1 "required tools"
ver() {  # ver <binary> <version-args...>
  local b=$1; shift
  if ! command -v "$b" >/dev/null 2>&1; then say "$b" "FAIL: missing"; FAIL=1; return 1; fi
  local v; v=$("$@" 2>&1 | head -1)
  [ -z "$v" ] && v="(present, version unavailable)"
  say "$b" "$v"
}
ver salmon  salmon --version
ver seqtk   sh -c 'seqtk 2>&1 | head -2 | tail -1'
ver curl    curl --version
ver gzip    gzip --version
ver md5sum  md5sum --version
ver awk     awk -W version
ver python3 python3 --version

head1 "python packages"
python3 - <<'PY' 2>&1 | sed 's/^/  /'
import importlib
for m in ("numpy", "pandas", "scipy"):
    try:
        mod = importlib.import_module(m)
        print("%-22s %s" % (m, getattr(mod, "__version__", "?")))
    except Exception as e:
        print("%-22s FAIL: %s" % (m, e))
PY
python3 -c "import numpy, pandas, scipy" 2>/dev/null || FAIL=1

head1 "salmon option support"
# capture help output first: piping into `grep -q` sends SIGPIPE to salmon, and with
# pipefail that turns a successful match into a non-zero pipeline (bug found on first run)
SALMON_QUANT_HELP=$(salmon quant --help-reads 2>&1 || true)
SALMON_INDEX_HELP=$(salmon index --help 2>&1 || true)
for opt in --gcBias --seqBias --posBias --minAssignedFrags; do
  case "$SALMON_QUANT_HELP" in
    *"$opt"*) say "$opt" ok ;;
    *)        say "$opt" "FAIL: unsupported by this salmon"; FAIL=1 ;;
  esac
done
case "$SALMON_INDEX_HELP" in
  *--gencode*) say --gencode ok ;;
  *)           say --gencode "FAIL"; FAIL=1 ;;
esac

head1 "optional (speed)"
for b in pigz aria2c seqtk fastqc multiqc Rscript; do
  if command -v "$b" >/dev/null 2>&1; then say "$b" present; else say "$b" "absent (ok)"; fi
done

head1 "network reachability"
for host in ftp.sra.ebi.ac.uk ftp.ebi.ac.uk; do
  if curl -sI --max-time 20 "https://$host/" >/dev/null 2>&1; then say "$host" reachable
  else say "$host" "WARN: not reachable"; WARN=1; fi
done

head1 "verdict"
if [ "$FAIL" -eq 0 ]; then
  [ "$WARN" -eq 1 ] && echo "READY (with warnings above)"
  [ "$WARN" -eq 0 ] && echo "READY"
  exit 0
fi
echo "NOT READY (see FAIL lines)"
exit 1
