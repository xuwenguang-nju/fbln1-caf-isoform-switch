"""Rebuild `caf_panel_depth_consistency.csv` under the fixed-referent statistic.

The filed version was produced from the old dominant-isoform dIF, so its numbers
(FBLN1 |delta dIF| = 0.0034, rank 1 of 14 in relative terms, 8 of 26 in absolute
terms) belong to the withdrawn definition. Figure 2 panel f sits beside panels
a-d, which now use the fixed referent, so this table must use it too.

Same schema as before so the plotting code is unchanged.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
T = ROOT / "results/tables"
PANEL = ROOT / "analysis/scripts/remote_step1/panel_genes.tsv"
FBLN1 = "ENSG00000077942"


def main() -> None:
    full = pd.read_csv(T / "switch_test_full_fixed_referent.tsv", sep="\t")
    gated = pd.read_csv(T / "candidates_gated_fixed_referent.tsv", sep="\t").rename(
        columns={"mean_dIF_full": "g_dIF_full", "mean_dIF_sub": "dIF_sub",
                 "q_panel_full": "q_panel_g"})
    panel = pd.read_csv(PANEL, sep="\t")
    panel["gene"] = panel["gene"].astype(str)
    print(f"panel genes in file: {len(panel)}")

    m = (panel.merge(full[["gene", "mean_dIF", "sd_dIF", "p", "n_pos", "n_patients",
                           "q_panel", "q_all"]], on="gene", how="left")
              .merge(gated[["gene", "dIF_sub", "abs_change", "reproducible", "call"]],
                     on="gene", how="left"))
    m["dIF_full"] = m["mean_dIF"]

    meas = m[m["mean_dIF"].notna() & m["dIF_sub"].notna()]
    print(f"measurable at both depths: {len(meas)}")
    survive = int((meas["abs_change"] <= 0.05).sum())
    print(f"within the 0.05 gate tolerance: {survive}")
    print(f"beyond it: {len(meas) - survive}")

    f = meas.copy()
    f["rel"] = f["abs_change"] / f["mean_dIF"].abs()
    r = f[f["gene"] == FBLN1]
    if len(r):
        x = r.iloc[0]
        print(f"\nFBLN1: dIF {x['mean_dIF']:+.4f}  dIF_sub {x['dIF_sub']:+.4f}  "
              f"|delta| {x['abs_change']:.4f}  rel {x['rel']:.4f}")
        sub = f[f["mean_dIF"].abs() > 0.02].reset_index(drop=True)
        print(f"  relative rank: {int(sub['rel'].rank().loc[sub['gene'].eq(FBLN1)].iloc[0])}"
              f" of {len(sub)} (genes with |dIF| > 0.02)")
        print(f"  absolute rank: {int(f['abs_change'].rank().loc[f['gene'].eq(FBLN1)].iloc[0])}"
              f" of {len(f)}")

    out = m[["gene", "symbol", "mean_dIF", "sd_dIF", "p", "n_pos", "n_patients",
             "q_panel", "q_all", "dIF_sub", "abs_change", "reproducible", "call",
             "dIF_full"]]
    out.to_csv(T / "caf_panel_depth_consistency_fixed_referent.csv", index=False,
               encoding="utf-8-sig")
    print(f"\nwritten: {T / 'caf_panel_depth_consistency_fixed_referent.csv'}")


if __name__ == "__main__":
    main()
