"""Sweep SLSQP ftol on scenario 1, scored against Table S2.

The x0clip experiment showed that *less* optimization lands closer to the
published values, implying MATLAB fmincon (TolFun=1e-6) stops earlier than our
ftol=1e-10. This sweeps ftol to find the value that best matches fmincon's
effective convergence.
"""
import numpy as np

import meandir.engine.inversion as inv
from meandir.run import run

UE = "reference/data/MEANDIR_UserEntries.xlsx"
RD = "reference/data/RiverDataSpreadsheet_Kemeny_etal_2023.xlsx"
NAME = "AK_scenario1_carbonate_slctindi"
PUB = {("Mg", "carb"): 67.4, ("Mg", "slct_Mg"): 32.4,
       ("Ca", "carb"): 69.0, ("Ca", "slct_Ca"): 10.3, ("Ca", "evap"): 15.8,
       ("DIC", "carb"): 68.5, ("DIC", "corg"): 31.5,
       ("SO4", "pyri"): 79.4, ("SO4", "evap"): 20.0, ("Na", "slct_Na"): 94.2}
FTOLS = [1e-10, 1e-6, 1e-4, 1e-3]


def score(ftol, success, max_iter, maxidx):
    inv.SLSQP_FTOL = ftol
    _r, ctx = run(NAME, UE, RD, seed=1, max_success=success, max_iter=max_iter,
                  max_zerohits=4000, samples=range(0, maxidx))
    S = ctx["summary"]
    out = {}
    for (ion, em), _t in PUB.items():
        med = S.fraction[ion][em]["median"]
        out[(ion, em)] = float(np.nanmean(med)) * 100 if np.any(~np.isnan(med)) else float("nan")
    return out


def main():
    import argparse
    import time
    ap = argparse.ArgumentParser()
    ap.add_argument("--success", type=int, default=30)
    ap.add_argument("--max-iter", type=int, default=6000)
    ap.add_argument("--maxidx", type=int, default=199)
    args = ap.parse_args()

    res = {}
    for ft in FTOLS:
        t0 = time.time()
        res[ft] = score(ft, args.success, args.max_iter, args.maxidx)
        worst = max(abs(res[ft][k] - t) for k, t in PUB.items())
        rms = np.sqrt(np.mean([(res[ft][k] - t) ** 2 for k, t in PUB.items()]))
        print(f"ftol={ft:.0e}  worstΔ={worst:5.1f}  rmsΔ={rms:5.2f}  ({time.time()-t0:.0f}s)", flush=True)

    print(f"\n{'ion':4s} {'em':9s} {'pub':>6s}   " + "  ".join(f"{ft:.0e}" for ft in FTOLS))
    for (ion, em), t in PUB.items():
        cells = "  ".join(f"{res[ft][(ion, em)]:5.1f}" for ft in FTOLS)
        print(f"{ion:4s} {em:9s} {t:6.1f}   {cells}")


if __name__ == "__main__":
    main()
