"""Validate the Python port against Table S2 of Kemeny et al. (2023) supplement.

Table S2 reports, per scenario, the *mean over rivers of the median fractional
contribution* of each end-member to each dissolved ion (the "All samples"
column). That is exactly ``mean_over_samples( fraction[ion][em]["median"] )``
in this port, so we run scenarios 1 and 2 over every functional sample and diff
against the published percentages.

Usage:  python -m scripts.validate_against_tableS2 [--success N]
"""

from __future__ import annotations

import argparse

import numpy as np

from meandir.run import run

UE = "reference/data/MEANDIR_UserEntries.xlsx"
RD = "reference/data/RiverDataSpreadsheet_Kemeny_etal_2023.xlsx"

# Published "All samples" mean-of-median contributions (%), Table S2.
# (ion, paper_end_member, port_em) -> {scenario: published_pct}
# Silicate maps to the ion-specific resolved silicate end-member.
PUBLISHED = {
    "DIC": [("Carbonate", "carb", {1: 68.5, 2: 72.7}),
            ("Corg oxidation", "corg", {1: 31.5, 2: 38.5}),
            ("Degassing", "degas", {1: None, 2: -10.7})],
    "Ca": [("Carbonate", "carb", {1: 69.0, 2: 75.2}),
           ("Evaporite", "evap", {1: 15.8, 2: 15.7}),
           ("Silicate", "slct_Ca", {1: 10.3, 2: 4.6}),
           ("Precipitation", "prec", {1: 0.5, 2: 0.5})],
    "Mg": [("Carbonate", "carb", {1: 67.4, 2: 70.4}),
           ("Silicate", "slct_Mg", {1: 32.4, 2: 29.4}),
           ("Precipitation", "prec", {1: 0.2, 2: 0.2})],
    "Na": [("Silicate", "slct_Na", {1: 94.2, 2: 94.2}),
           ("Precipitation", "prec", {1: 5.8, 2: 5.8})],
    "K": [("Silicate", "slct_K", {1: 90.8, 2: 90.6}),
          ("Precipitation", "prec", {1: 9.2, 2: 9.4})],
    "Cl": [("Precipitation", "prec", {1: 100.0, 2: 100.0})],
    "SO4": [("H2SO4 production", "pyri", {1: 79.4, 2: 79.6}),
            ("Evaporite", "evap", {1: 20.0, 2: 19.8}),
            ("Precipitation", "prec", {1: 0.6, 2: 0.6})],
}

SCENARIOS = {1: "AK_scenario1_carbonate_slctindi",
             2: "AK_scenario2_carbonate_slctindi_degas_2p5"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--success", type=int, default=30)
    ap.add_argument("--max-iter", type=int, default=12000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    for sc, name in SCENARIOS.items():
        res, ctx = run(name, UE, RD, seed=args.seed, max_success=args.success,
                       max_iter=args.max_iter, max_zerohits=4000)
        S = ctx["summary"]
        print(f"\n================ Scenario {sc}: {name} ================")
        print(f"{'ion':4s} {'end-member':16s} {'published':>10s} {'python':>9s} {'Δ':>8s}")
        worst = 0.0
        for ion, rows in PUBLISHED.items():
            for paper_em, em, pub in rows:
                target = pub.get(sc)
                med = S.fraction[ion][em]["median"]
                val = float(np.nanmean(med)) * 100 if np.any(~np.isnan(med)) else float("nan")
                if target is None:
                    print(f"{ion:4s} {paper_em:16s} {'—':>10s} {val:>8.1f}%   (n/a sc{sc})")
                    continue
                diff = val - target
                worst = max(worst, abs(diff))
                flag = "" if abs(diff) <= 3 else "  <-- check"
                print(f"{ion:4s} {paper_em:16s} {target:>9.1f}% {val:>8.1f}% {diff:>+7.1f}{flag}")
        print(f"  worst absolute deviation: {worst:.1f} percentage points")


if __name__ == "__main__":
    main()
