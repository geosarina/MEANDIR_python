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
    samples = main_idx + slough_idx
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

    print(f"\nMg<-carbonate mean-of-median (%), published = 67.4")
    print(f"{'seed':>4}  {'All(n=30)':>10}  {'Mainstem(n=26)':>14}")
    allvals, mainvals = [], []
    for seed in SEEDS:
        a = [cache[f"{seed}|{s}"] for s in samples if cache.get(f"{seed}|{s}") is not None]
        m = [cache[f"{seed}|{s}"] for s in main_idx if cache.get(f"{seed}|{s}") is not None]
        av, mv = np.mean(a), np.mean(m)
        allvals.append(av); mainvals.append(mv)
        print(f"{seed:>4}  {av:>10.1f}  {mv:>14.1f}")
    print(f"\nAll      : mean={np.mean(allvals):.1f}  sd={np.std(allvals, ddof=1):.2f}  "
          f"range=[{min(allvals):.1f},{max(allvals):.1f}]  -> published 67.4 is "
          f"{(67.4-np.mean(allvals))/np.std(allvals, ddof=1):+.1f} sd away")
    print(f"Mainstem : mean={np.mean(mainvals):.1f}  sd={np.std(mainvals, ddof=1):.2f}  "
          f"range=[{min(mainvals):.1f},{max(mainvals):.1f}]")


if __name__ == "__main__":
    main()
