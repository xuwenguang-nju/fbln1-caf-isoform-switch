#!/usr/bin/env bash
# One-screen status of step 1. Read-only, safe to run at any time.
#
# Two traps this script avoids, both found the hard way on 2026-10-19:
#   * `du -sb` and `find -size` use the APPARENT size, and aria2 writes at offsets across its
#     connections, so a half-finished file already reports its full size. Real progress needs
#     block counts (`du -s --block-size=1`) plus the absence of the `.aria2` control file.
#   * `x.fastq.gz` and `x.fastq.gz.aria2` are one file, not two: counting both doubles the total.
W=/root/hnscc_caf
MAN=$W/remote_step1/gse303690_ena_manifest.tsv
echo "== $(date -u +%FT%TZ) =="
echo "in flight        : $(pgrep -c aria2c 2>/dev/null || echo 0) aria2, $(pgrep -c salmon 2>/dev/null || echo 0) salmon"

python3 - "$W" "$MAN" <<'PY'
import os, subprocess, sys, time
W, man = sys.argv[1], sys.argv[2]
exp = {}
with open(man) as fh:
    next(fh)
    for line in fh:
        p = line.rstrip("\n").split("\t")
        for mate, url, b in ((1, p[7], int(p[9])), (2, p[10], int(p[12]))):
            if url:
                exp[f"{p[0]}_{mate}.fastq.gz"] = (b, p[2], p[1], p[3])
have = done = 0
parts = []
for name, (b, cond, pat, rep) in exp.items():
    f = os.path.join(W, "fastq", name)
    if not os.path.exists(f):
        continue
    blocks = int(subprocess.run(["du", "-s", "--block-size=1", f],
                                capture_output=True, text=True).stdout.split()[0])
    have += blocks
    if not os.path.exists(f + ".aria2") and blocks == b:
        done += 1
    else:
        parts.append((100.0 * blocks / b, blocks / 1e9, name, cond, pat, rep))
exp_total = sum(v[0] for v in exp.values())
rx1 = int(open("/sys/class/net/eth0/statistics/rx_bytes").read())
time.sleep(5)
rx2 = int(open("/sys/class/net/eth0/statistics/rx_bytes").read())
rate = (rx2 - rx1) / 5
print(f"progress         : {have/1e9:6.1f} / {exp_total/1e9:.1f} GB  ({100*have/exp_total:5.1f}%)")
print(f"files complete   : {done} / {len(exp)}")
print(f"throughput       : {rate/1e6:.2f} MB/s")
if rate > 100_000:
    print(f"ETA at this rate : {(exp_total-have)/rate/3600:.1f} h")
for pct, gb, name, cond, pat, rep in sorted(parts, reverse=True)[:8]:
    print(f"   {pct:5.1f}%  {gb:5.2f} GB  {name}  {cond} HT{pat} rep{rep}")
PY

echo "quant full       : $(find "$W/quant_full"  -name quant.sf 2>/dev/null | wc -l) / 36"
echo "quant sub10      : $(find "$W/quant_sub10" -name quant.sf 2>/dev/null | wc -l) / 36"
echo "disk free        : $(df -h "$W" | awk 'NR==2{print $4}')"
echo "-- recent events --"
grep -E "DONE|FAILED|salmon full rc|index ready" "$W/run_step1.log" 2>/dev/null | tail -4 || true
