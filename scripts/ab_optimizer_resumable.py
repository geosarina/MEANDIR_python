"""Resumable, full-scale optimizer A/B: SLSQP vs trust-constr, scored against
Table S2. Checkpoints every (method, scenario, sample) result to disk, so a
worker restart loses at most one sample and the next invocation resumes.

Run repeatedly (foreground in ~10-min chunks, or detached); each invocation
makes progress and, once a (method, scenario) is fully covered, prints its
aggregated comparison against the published values.

  python -m scripts.ab_optimizer_resumable --success 100
"""
import argparse
import json
import os
import time

import numpy as np

import meandir.engine.inversion as inv
from meandir.run import run
from meandir.scenarios import find_scenario_parameters
from meandir.io.user_entries import load_conc2equi
from meandir.io.river_data import read_river_observations, build_model_variables
from scripts.validate_against_tableS2 import PUBLISHED, SCENARIOS, UE, RD

CACHE = "/tmp/claude-0/-home-user-MEANDIR-python/feec57dc-e16b-557c-afd4-47f9242d36db/scratchpad/ab_opt_cache.json"
METHODS = ("SLSQP", "trust-constr")


def functional_indices():
    p = find_scenario_parameters(SCENARIOS[1], UE)
    river = read_river_observations(RD, load_conc2equi(UE))
    mask = build_model_variables(river, p.ObsInNormalization, p.EMUnits, p.ObsList)
    return [int(i) for i in np.where(mask)[0]]


def load_cache():
    if os.path.exists(CACHE):
        with open(CACHE) as f:
            return json.load(f)
    return {}


def save_cache(c):
    tmp = CACHE + ".tmp"
    with open(tmp, "w") as f:
        json.dump(c, f)
    os.replace(tmp, CACHE)


def cells_for(sc):
    """list of (ion, em, target%) with a published value for this scenario."""
    out = []
    for ion, items in PUBLISHED.items():
        for _paper_em, em, pub in items:
            t = pub.get(sc)
            if t is not None:
                out.append((ion, em, t))
    return out


def run_one(method, sc, name, sample, success, max_iter):
    inv.OPTIMIZER_METHOD = method
    _res, ctx = run(name, UE, RD, seed=1, max_success=success, max_iter=max_iter,
                    max_zerohits=4000, samples=[sample])
    S = ctx["summary"]
    out = {}
    for ion, em, _t in cells_for(sc):
        med = S.fraction[ion][em]["median"][sample]
        out[f"{ion}|{em}"] = None if np.isnan(med) else float(med) * 100
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--success", type=int, default=100)
    ap.add_argument("--max-iter", type=int, default=15000)
    ap.add_argument("--budget", type=int, default=10**9, help="seconds before stopping")
    args = ap.parse_args()

    idx = functional_indices()
    cache = load_cache()
    t_start = time.time()
    done_this_run = 0

    for method in METHODS:
        for sc, name in SCENARIOS.items():
            for sample in idx:
                if time.time() - t_start > args.budget:
                    print(f"time budget reached; {done_this_run} samples this run", flush=True)
                    _report(cache, idx)
                    return
                key = f"{method}|{sc}|{sample}"
                if key in cache:
                    continue
                t0 = time.time()
                cache[key] = run_one(method, sc, name, sample, args.success, args.max_iter)
                save_cache(cache)
                done_this_run += 1
                print(f"{key}  ({time.time()-t0:.0f}s)", flush=True)

    print("\nALL DONE.", flush=True)
    _report(cache, idx)


def _report(cache, idx):
    for sc in SCENARIOS:
        cells = cells_for(sc)
        print(f"\n===== Scenario {sc} =====")
        print(f"{'ion':4s} {'em':9s} {'pub':>6s} {'SLSQP':>8s} {'Δsl':>6s} {'trustc':>8s} {'Δtc':>6s}")
        worst = {"SLSQP": 0.0, "trust-constr": 0.0}
        complete = {m: all(f"{m}|{sc}|{s}" in cache for s in idx) for m in METHODS}
        for ion, em, target in cells:
            vals = {}
            for m in METHODS:
                per = [cache.get(f"{m}|{sc}|{s}", {}).get(f"{ion}|{em}") for s in idx]
                per = [v for v in per if v is not None]
                vals[m] = float(np.mean(per)) if per else float("nan")
                if per and not np.isnan(vals[m]):
                    worst[m] = max(worst[m], abs(vals[m] - target))
            ds = vals["SLSQP"] - target
            dt = vals["trust-constr"] - target
            mark = "  *" if abs(dt) + 0.2 < abs(ds) else ("  x" if abs(dt) > abs(ds) + 0.2 else "")
            print(f"{ion:4s} {em:9s} {target:6.1f} {vals['SLSQP']:7.1f}% {ds:+6.1f} "
                  f"{vals['trust-constr']:7.1f}% {dt:+6.1f}{mark}")
        for m in METHODS:
            tag = "complete" if complete[m] else "PARTIAL"
            print(f"  worst {m}: {worst[m]:.1f} ({tag})")


if __name__ == "__main__":
    main()
