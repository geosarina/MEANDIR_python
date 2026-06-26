"""Scenario-1 Table S2 comparison for three solver modes:
  SLSQP      - current nonlinear optimizer
  lsq_linear - exact convex bounded least-squares (true global minimum)
  x0clip     - X0 clipped to bounds, no optimization (maximal under-convergence)

Tests whether the Ca/Mg residual is the optimizer (it isn't, if SLSQP==exact)
and whether MATLAB-style under-convergence (x0clip) moves toward published.
"""
import numpy as np

import meandir.engine.inversion as inv
from meandir.run import run

UE = "reference/data/MEANDIR_UserEntries.xlsx"
RD = "reference/data/RiverDataSpreadsheet_Kemeny_etal_2023.xlsx"
NAME = "AK_scenario1_carbonate_slctindi"
PUB = {("Mg", "carb"): 67.4, ("Mg", "slct_Mg"): 32.4,
       ("Ca", "carb"): 69.0, ("Ca", "slct_Ca"): 10.3, ("Ca", "evap"): 15.8,
       ("DIC", "carb"): 68.5, ("SO4", "pyri"): 79.4, ("Na", "slct_Na"): 94.2}


def score(mode, success, max_iter):
    inv.SOLVER_MODE = mode
    _r, ctx = run(NAME, UE, RD, seed=1, max_success=success, max_iter=max_iter,
                  max_zerohits=4000)
    S = ctx["summary"]
    out = {}
    for (ion, em), _t in PUB.items():
        med = S.fraction[ion][em]["median"]
        out[(ion, em)] = float(np.nanmean(med)) * 100 if np.any(~np.isnan(med)) else float("nan")
    inv.SOLVER_MODE = "default"
    return out


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--success", type=int, default=40)
    ap.add_argument("--max-iter", type=int, default=8000)
    args = ap.parse_args()

    res = {m: score(m, args.success, args.max_iter)
           for m in ("default", "lsq_linear", "x0clip")}
    print(f"{'ion':4s} {'em':9s} {'pub':>6s} {'SLSQP':>8s} {'exact':>8s} {'x0clip':>8s}"
          f"   {'Δslsqp':>7s} {'Δexact':>7s} {'Δx0clip':>8s}")
    for (ion, em), t in PUB.items():
        vs, ve, vx = res["default"][(ion, em)], res["lsq_linear"][(ion, em)], res["x0clip"][(ion, em)]
        print(f"{ion:4s} {em:9s} {t:6.1f} {vs:7.1f}% {ve:7.1f}% {vx:7.1f}%"
              f"   {vs-t:+7.1f} {ve-t:+7.1f} {vx-t:+8.1f}")


if __name__ == "__main__":
    main()
