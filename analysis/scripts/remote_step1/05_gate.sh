#!/usr/bin/env bash
# Step 5: the depth-consistency gate and the final candidate table.
#
# Rule fixed before looking at the numbers (it comes from the local validation, where DCN
# flipped from +0.001 to -0.122 between depths):
#   * a candidate is only reported if the two depths agree in sign AND |delta dIF| <= 0.04
#   * effect sizes are reported in three tiers: switch (>=0.15), small shift (0.05-0.15),
#     negligible (<0.05); only the first may be called an isoform switch
#
# Usage:  bash 05_gate.sh <outdir> <panel.tsv>
set -euo pipefail
OUT=${1:?usage: 05_gate.sh <outdir> <panel.tsv>}
PANEL=${2:?usage: 05_gate.sh <outdir> <panel.tsv>}

F="$OUT/switch_test_full.tsv"
S="$OUT/switch_test_sub10.tsv"
[ -s "$F" ] && [ -s "$S" ] || { echo "missing switch_test_full.tsv / switch_test_sub10.tsv"; exit 1; }

python3 - "$F" "$S" "$PANEL" "$OUT" <<'PY'
import sys
import pandas as pd
import numpy as np

f, s, panel_path, out = sys.argv[1:5]
full = pd.read_csv(f, sep="\t")
sub = pd.read_csv(s, sep="\t", dtype={"gene": str})
panel = pd.read_csv(panel_path, sep="\t")

keep = ["gene", "n_patients", "mean_dIF", "sd_dIF", "p", "q_panel", "q_all", "n_pos"]
full = full[[c for c in keep if c in full.columns]].add_suffix("_full").rename(columns={"gene_full": "gene"})
sub = sub[[c for c in keep if c in sub.columns]].add_suffix("_sub").rename(columns={"gene_sub": "gene"})
m = full.merge(sub, on="gene", how="inner")
m = m.merge(panel, on="gene", how="left")

m["abs_change"] = (m["mean_dIF_full"] - m["mean_dIF_sub"]).abs()
m["sign_same"] = np.sign(m["mean_dIF_full"]) == np.sign(m["mean_dIF_sub"])
m["reproducible"] = m["sign_same"] & (m["abs_change"] <= 0.04)

def tier(r):
    if not r["reproducible"]:
        return "excluded_depth_unstable"
    a = abs(r["mean_dIF_full"])
    q = r.get("q_panel_full", np.nan)
    if a >= 0.15 and q < 0.05:
        return "switch"
    if a >= 0.05 and q < 0.05:
        return "small_shift"
    if a < 0.05 and q < 0.05:
        return "negligible_but_reproducible"
    return "not_significant"
m["call"] = m.apply(tier, axis=1)

m = m.sort_values(["call", "mean_dIF_full"], ascending=[True, False])
m.to_csv(f"{out}/candidates_gated.tsv", sep="\t", index=False)

print("== calls ==")
print(m["call"].value_counts().to_string())
print("\n== panel candidates (reproducible, q_panel<0.05) ==")
c = m[m["call"].isin(["switch", "small_shift"])]
if len(c):
    print(c[["symbol", "mean_dIF_full", "sd_dIF_full", "q_panel_full", "mean_dIF_sub",
             "q_panel_sub", "abs_change", "n_pos_full"]].round(4).to_string(index=False))
else:
    print("(none)")
print("\n== excluded for depth instability ==")
e = m[m["call"] == "excluded_depth_unstable"]
if len(e):
    print(e[["symbol", "mean_dIF_full", "mean_dIF_sub", "abs_change"]].round(4).to_string(index=False))
else:
    print("(none)")
PY

echo
echo "WROTE $OUT/candidates_gated.tsv"
