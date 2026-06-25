"""Focused, fast A/B: SLSQP vs trust-constr on scenario 1 only (the square
system where the Ca/Mg residual is purely the optimizer's feasible-point
choice). Prints the Mg and Ca carbonate/silicate split for each optimizer
against the published Table S2 values.
"""
import argparse
import time
import numpy as np

import meandir.engine.inversion as inv
from meandir.run import run

UE = "reference/data/MEANDIR_UserEntries.xlsx"
RD = "reference/data/RiverDataSpreadsheet_Kemeny_etal_2023.xlsx"
NAME = "AK_scenario1_carbonate_slctindi"

# published Table S2 (scenario 1, All samples) for the cells of interest
PUB = {("Mg", "carb"): 67.4, ("Mg", "slct_Mg"): 32.4,
       ("Ca", "carb"): 69.0, ("Ca", "slct_Ca"): 10.3, ("Ca", "evap"): 15.8,
       ("DIC", "carb"): 68.5, ("DIC", "corg"): 32.4,
       ("SO4", "pyri"): 79.4, ("Na", "slct_Na"): 94.2}


def score(success, max_iter, seed):
    res, ctx = run(NAME, UE, RD, seed=seed, max_success=success,
                   max_iter=max_iter, max_zerohits=4000)
    S = ctx["summary"]
    out = {}
    for (ion, em), target in PUB.items():
        med = S.fraction[ion][em]["median"]
        out[(ion, em)] = float(np.nanmean(med)) * 100 if np.any(~np.isnan(med)) else float("nan")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--success", type=int, default=20)
    ap.add_argument("--max-iter", type=int, default=8000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    results = {}
    for method in ("SLSQP", "trust-constr"):
        inv.OPTIMIZER_METHOD = method
        t0 = time.time()
        results[method] = score(args.success, args.max_iter, args.seed)
        print(f"[{method}] done in {time.time()-t0:.0f}s", flush=True)

    print(f"\n{'ion':4s} {'em':9s} {'pub':>6s} {'SLSQP':>8s} {'Δsl':>6s} {'trustc':>8s} {'Δtc':>6s}")
    for (ion, em), target in PUB.items():
        vs = results["SLSQP"][(ion, em)]
        vt = results["trust-constr"][(ion, em)]
        ds, dt = vs - target, vt - target
        mark = "  *" if abs(dt) + 0.2 < abs(ds) else ("  x" if abs(dt) > abs(ds) + 0.2 else "")
        print(f"{ion:4s} {em:9s} {target:6.1f} {vs:7.1f}% {ds:+6.1f} {vt:7.1f}% {dt:+6.1f}{mark}")


if __name__ == "__main__":
    main()
