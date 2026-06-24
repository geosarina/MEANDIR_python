"""A/B test: SLSQP vs trust-constr optimizer, scored against Table S2.

The optimizer-activity diagnostic showed the Ca/Mg residual is the bounded
optimizer's choice of feasible point. MATLAB fmincon defaults to interior-point;
this checks whether SciPy's trust-constr (interior-point-style) matches the
published contributions better than SLSQP.
"""
import argparse
import time
import numpy as np

import meandir.engine.inversion as inv
from meandir.run import run
from scripts.validate_against_tableS2 import PUBLISHED, SCENARIOS, UE, RD


def score(sc, name, success, max_iter, seed):
    res, ctx = run(name, UE, RD, seed=seed, max_success=success,
                   max_iter=max_iter, max_zerohits=4000)
    S = ctx["summary"]
    rows, worst = [], 0.0
    for ion, items in PUBLISHED.items():
        for paper_em, em, pub in items:
            target = pub.get(sc)
            if target is None or em not in S.fraction[ion]:
                continue
            med = S.fraction[ion][em]["median"]
            val = float(np.nanmean(med)) * 100 if np.any(~np.isnan(med)) else float("nan")
            rows.append((ion, paper_em, target, val, val - target))
            worst = max(worst, abs(val - target))
    return rows, worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--success", type=int, default=40)
    ap.add_argument("--max-iter", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    for sc, name in SCENARIOS.items():
        print(f"\n############ Scenario {sc}: {name} ############")
        out = {}
        for method in ("SLSQP", "trust-constr"):
            inv.OPTIMIZER_METHOD = method
            t0 = time.time()
            out[method] = (score(sc, name, args.success, args.max_iter, args.seed),
                           time.time() - t0)
        (rows_s, worst_s), ts = out["SLSQP"]
        (rows_t, worst_t), tt = out["trust-constr"]
        print(f"{'ion':4s} {'end-member':16s} {'pub':>7s} {'SLSQP':>9s} {'Δsl':>6s} "
              f"{'trustc':>9s} {'Δtc':>6s}")
        for (ion, pe, t, vs, ds), (_, _, _, vt, dt) in zip(rows_s, rows_t):
            mark = "  *" if abs(dt) + 0.2 < abs(ds) else ("  x" if abs(dt) > abs(ds) + 0.2 else "")
            print(f"{ion:4s} {pe:16s} {t:7.1f} {vs:8.1f}% {ds:+6.1f} {vt:8.1f}% {dt:+6.1f}{mark}")
        print(f"  worst:  SLSQP={worst_s:.1f} ({ts:.0f}s)   trust-constr={worst_t:.1f} ({tt:.0f}s)")


if __name__ == "__main__":
    main()
