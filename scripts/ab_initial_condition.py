"""A/B test: minimum-norm vs MATLAB-mldivide basic initial condition, scored
against Table S2. Runs scenarios 1 and 2 under each mode and prints the per-cell
deviation so we can see whether the basic solution tightens the Ca/Mg split.
"""
import argparse
import numpy as np

import meandir.engine.inversion as inv
from meandir.run import run
from scripts.validate_against_tableS2 import PUBLISHED, SCENARIOS, UE, RD


def score(sc, name, success, max_iter, seed):
    res, ctx = run(name, UE, RD, seed=seed, max_success=success,
                   max_iter=max_iter, max_zerohits=4000)
    S = ctx["summary"]
    rows = []
    worst = 0.0
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
        for mode in ("minnorm", "basic"):
            inv.MLDIVIDE_INITIAL_CONDITION = mode
            out[mode] = score(sc, name, args.success, args.max_iter, args.seed)
        rows_m, worst_m = out["minnorm"]
        rows_b, worst_b = out["basic"]
        print(f"{'ion':4s} {'end-member':16s} {'pub':>7s} {'minnorm':>9s} {'Δmn':>6s} "
              f"{'basic':>9s} {'Δba':>6s}")
        for (ion, pe, t, vm, dm), (_, _, _, vb, db) in zip(rows_m, rows_b):
            mark = "  *" if abs(db) + 0.2 < abs(dm) else ("  x" if abs(db) > abs(dm) + 0.2 else "")
            print(f"{ion:4s} {pe:16s} {t:7.1f} {vm:8.1f}% {dm:+6.1f} {vb:8.1f}% {db:+6.1f}{mark}")
        print(f"  worst deviation:  minnorm={worst_m:.1f}   basic={worst_b:.1f}")


if __name__ == "__main__":
    main()
