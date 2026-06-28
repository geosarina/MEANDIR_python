"""Is the scenario-1 Mg-carbonate residual Monte Carlo noise or systematic?

Run scenario 1 under several RNG seeds (same config) and report the spread of the
All-samples and Mainstem mean-of-median Mg<-carbonate contribution. If the
published 67.4% lies within the seed-to-seed spread, the ~3-point gap is
realization noise; if every seed clusters well below it, the gap is systematic.

Resumable: caches the per-(seed, sample) Mg-carbonate median.
"""
import json
import os

import numpy as np

import meandir.run as R
from scripts.validate_groups import classify, SCEN

UE = "reference/data/MEANDIR_UserEntries.xlsx"
RD = "reference/data/RiverDataSpreadsheet_Kemeny_etal_2023.xlsx"
NAME = SCEN[1]
SEEDS = [1, 2, 3, 4, 5, 6]
SUCCESS = 120
CACHE = "/tmp/claude-0/-home-user-MEANDIR-python/feec57dc-e16b-557c-afd4-47f9242d36db/scratchpad/seed_spread.json"


def load():
    return json.load(open(CACHE)) if os.path.exists(CACHE) else {}


def save(c):
    json.dump(c, open(CACHE + ".tmp", "w")); os.replace(CACHE + ".tmp", CACHE)


def main():
    main_idx, slough_idx = classify()[:2]
    samples = main_idx          # mainstem only (clean test; slough is slow + noisy)
    cache = load()
    for seed in SEEDS:
        for s in samples:
            key = f"{seed}|{s}"
            if key in cache:
                continue
            r, ctx = R.run(NAME, UE, RD, seed=seed, max_success=SUCCESS,
                           max_iter=12000, max_zerohits=4000, samples=[s])
            ci = ctx["params"].EMList0.index("carb")
            med = ctx["summary"].fraction["Mg"]["carb"]["median"][s]
            cache[key] = None if np.isnan(med) else float(med) * 100
            save(cache)
        print(f"seed {seed} done", flush=True)

    PUB_MAIN = 64.3   # Table S2 Mainstem Mg<-carbonate, scenario 1
    print(f"\nMainstem (n=26) Mg<-carbonate mean-of-median (%), published = {PUB_MAIN}")
    mainvals = []
    for seed in SEEDS:
        m = [cache[f"{seed}|{s}"] for s in main_idx if cache.get(f"{seed}|{s}") is not None]
        if len(m) < len(main_idx):
            print(f"  seed {seed}: incomplete ({len(m)}/{len(main_idx)})")
            continue
        mv = np.mean(m)
        mainvals.append(mv)
        print(f"  seed {seed}: {mv:.1f}")
    if len(mainvals) >= 2:
        sd = np.std(mainvals, ddof=1)
        print(f"\nmean={np.mean(mainvals):.1f}  sd={sd:.2f}  "
              f"range=[{min(mainvals):.1f},{max(mainvals):.1f}]")
        print(f"published {PUB_MAIN} is {(PUB_MAIN-np.mean(mainvals))/sd:+.1f} seed-sd from the python mean")


if __name__ == "__main__":
    main()
