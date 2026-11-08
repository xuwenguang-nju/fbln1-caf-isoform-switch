"""Smoke test: run every figure script from this release and report which reproduce.

Usage, from the release root (the folder containing this README):

    python analysis/scripts/482_smoke_test_release.py            # all 12 figures
    python analysis/scripts/482_smoke_test_release.py 3 9 S2     # a subset

Each script is run in a subprocess with this repository root as the working
directory. A pass means: exit status 0, and the expected PNG exists in figures/
with a plausible size. Stale PNGs are removed first, so a pass cannot be masked
by a file left over from an earlier run. Exit status is 1 if anything failed.

No third-party imports: this must run before the analysis dependencies are
confirmed to work.
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# (label as used in the manuscript, script, expected output, pass/fail above this size)
FIGURES = [
    ("1", "87_fig1_design_and_power.py", "fig_three_arm_design.png", 100_000),
    ("2", "63_fig1_depth_gate.py", "fig_depth_gate.png", 100_000),
    ("3", "67_fig2_fbln1_isoforms.py", "fig_fbln1_paired.png", 100_000),
    ("4", "126_fig_s1_singlecell.py", "supp_fig_s1_singlecell_fbln1.png", 100_000),
    ("5", "71_fig3_domains.py", "fig_fbln1_domains.png", 100_000),
    ("6", "69_fig4_tcga.py", "fig_fbln1_tcga.png", 100_000),
    ("7", "204_fig_pancancer.py", "fig_pancancer.png", 100_000),
    ("8", "70_fig5_cross_system.py", "fig_fbln1_cross_system.png", 100_000),
    ("9", "72_fig6_singlecell.py", "fig_singlecell_locus.png", 100_000),
    ("S1", "138_fig_s2_ledger.py", "supp_fig_s2_ledger.png", 100_000),
    ("S2", "136_fig_s3_controls.py", "supp_fig_s3_controls.png", 100_000),
    ("S3", "140_fig_s5_positioning.py", "supp_fig_s5_positioning.png", 100_000),
]


def main(argv: list[str]) -> int:
    wanted = {a.lstrip("S").lstrip("s").upper().replace("FIGURE", "") for a in argv[1:]}
    jobs = [f for f in FIGURES if not wanted or f[0].replace("S", "").upper() in wanted]
    if not jobs:
        print(f"no figure matched {argv[1:]}; labels are {[f[0] for f in FIGURES]}")
        return 2

    (ROOT / "figures").mkdir(exist_ok=True)
    failures = []
    for label, script, png, min_bytes in jobs:
        out = ROOT / "figures" / png
        if out.exists():
            out.unlink()
        start = time.time()
        proc = subprocess.run([sys.executable, f"analysis/scripts/{script}"],
                              cwd=ROOT, capture_output=True, text=True)
        elapsed = time.time() - start
        wrote = out.exists() and out.stat().st_size >= min_bytes
        ok = proc.returncode == 0 and wrote
        detail = (f"{out.stat().st_size / 1e6:.1f} MB" if out.exists()
                  else "no output")
        print(f"{'PASS' if ok else 'FAIL'}  Figure {label:<2s} {script:<31s} "
              f"{elapsed:5.1f}s  {detail}")
        if not ok:
            failures.append(label)
            tail = (proc.stderr or proc.stdout).strip().splitlines()[-6:]
            for line in tail:
                print(f"        {line}")

    print(f"\n{len(jobs) - len(failures)} of {len(jobs)} figure scripts reproduced.")
    if failures:
        print("failed: Figure " + ", Figure ".join(failures))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
