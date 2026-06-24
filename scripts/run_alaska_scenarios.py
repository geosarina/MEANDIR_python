"""Run all five published Alaska scenarios (Kemeny et al., 2023) end-to-end and
write a summary table of basin-median results.

This is the validation harness: it exercises the whole pipeline across every
published scenario and reports the headline quantities (per-end-member
fractional contribution to the cation budget, and the R/Z/C/W/Y weathering
metrics). Pair the emitted CSV with the paper's reported values to validate.

Usage:
    python -m scripts.run_alaska_scenarios [--success N] [--out FILE.csv]
"""

from __future__ import annotations

import argparse
import csv
import statistics

import numpy as np

from meandir.run import run

SCENARIOS = [
    "AK_scenario1_carbonate_slctindi",
    "AK_scenario2_carbonate_slctindi_degas_2p5",
    "AK_scenario3_carbonate_slctindi_degas_25",
    "AK_scenario4_carbonate_slctindi_degas_2p5_highfrac",
    "AK_scenario5_carbonate_slctindi_degas_2p5_lowfrac",
]

UE = "reference/data/MEANDIR_UserEntries.xlsx"
RD = "reference/data/RiverDataSpreadsheet_Kemeny_etal_2023.xlsx"


def _basin_median(per_sample_median):
    vals = [v for v in per_sample_median if v == v]  # drop NaN
    return statistics.median(vals) if vals else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--success", type=int, default=50,
                    help="successful simulations per sample (default 50)")
    ap.add_argument("--max-iter", type=int, default=20000,
                    help="iteration ceiling per sample (default 20000)")
    ap.add_argument("--zerohits", type=int, default=5000,
                    help="give up on a sample after this many zero-success iters")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="reference/alaska_results_python.csv")
    args = ap.parse_args()

    rows = []
    for name in SCENARIOS:
        res, ctx = run(name, UE, RD, seed=args.seed, max_success=args.success,
                       max_iter=args.max_iter, max_zerohits=args.zerohits)
        S = ctx["summary"]
        p = ctx["params"]
        n_ok = sum(1 for v in S.successes.values() if v > 0)
        print(f"\n{name}\n  samples inverted: {n_ok}")

        # basin-median fractional contribution to the normalization per EM
        for em in p.EMList0:
            bm = _basin_median(S.fraction["norm"][em]["median"])
            rows.append({"scenario": name, "quantity": f"f_norm[{em}]",
                         "basin_median": bm})
            print(f"    f_norm[{em:8s}] = {bm:+.3f}")

        # basin-median RZCWY (gross, unscaled)
        for V in "RZCWY":
            arr = S.rzcwy.get(V, {}).get("gross", {}).get("unscaled", {}).get("median")
            bm = _basin_median(arr) if arr is not None else float("nan")
            rows.append({"scenario": name, "quantity": f"{V}_gross_unscaled",
                         "basin_median": bm})
            print(f"    {V}_gross_unscaled = {bm:+.3f}")

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["scenario", "quantity", "basin_median"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {args.out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
