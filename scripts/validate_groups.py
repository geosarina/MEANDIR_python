"""Full Table S2 comparison: scenarios 1/2/3 x groups All/Mainstem/Slough.

Mainstem = ExtraField1 'river', Slough = 'slough'; the paper's "All samples" is
those two together. For each scenario we run once over the mainstem+slough
samples (per-sample checkpointed to disk for restart-resilience), then aggregate
the mean-of-median contributions for each group and diff against the published
Table S2 values. Emits markdown tables to scratchpad/validation_groups.md.

  python -m scripts.validate_groups --success 60
"""
import argparse
import json
import os

import numpy as np
import openpyxl

import meandir.engine.inversion as inv  # noqa: F401 (kept for parity)
from meandir.run import run
from meandir.scenarios import find_scenario_parameters
from meandir.io.user_entries import load_conc2equi
from meandir.io.river_data import read_river_observations, build_model_variables

UE = "reference/data/MEANDIR_UserEntries.xlsx"
RD = "reference/data/RiverDataSpreadsheet_Kemeny_etal_2023.xlsx"
SCEN = {1: "AK_scenario1_carbonate_slctindi",
        2: "AK_scenario2_carbonate_slctindi_degas_2p5",
        3: "AK_scenario3_carbonate_slctindi_degas_25"}
CACHE = "/tmp/claude-0/-home-user-MEANDIR-python/feec57dc-e16b-557c-afd4-47f9242d36db/scratchpad/validate_groups_cache.json"
OUT = "/tmp/claude-0/-home-user-MEANDIR-python/feec57dc-e16b-557c-afd4-47f9242d36db/scratchpad/validation_groups.md"

# Published Table S2 (mean-of-median %, from coordinate extraction of the PDF).
# (ion, paper label, port end-member) -> {group: [sc1, sc2, sc3]}
PUB = [
    ("DIC", "Carbonate", "carb", {"All": [68.5, 72.7, 72.2], "Main": [70.9, 75.6, 75.0], "Slough": [52.8, 53.5, 53.9]}),
    ("DIC", "Corg oxidation", "corg", {"All": [31.5, 38.5, 37.5], "Main": [29.1, 36.6, 35.6], "Slough": [47.3, 50.7, 49.7]}),
    ("DIC", "Degassing", "degas", {"All": [None, -10.7, -9.1], "Main": [None, -11.7, -10.1], "Slough": [None, -3.9, -2.6]}),
    ("Ca", "Carbonate", "carb", {"All": [69.0, 75.2, 75.5], "Main": [66.5, 73.4, 73.7], "Slough": [85.2, 86.7, 87.0]}),
    ("Ca", "Evaporite", "evap", {"All": [15.8, 15.7, 15.5], "Main": [17.0, 16.9, 16.6], "Slough": [8.1, 8.1, 8.3]}),
    ("Ca", "Silicate", "slct_Ca", {"All": [10.3, 4.6, 4.1], "Main": [11.2, 5.0, 4.5], "Slough": [1.9, 1.9, 1.9]}),
    ("Ca", "Precipitation", "prec", {"All": [0.5, 0.5, 0.5], "Main": [0.5, 0.5, 0.5], "Slough": [1.0, 1.0, 1.0]}),
    ("Mg", "Carbonate", "carb", {"All": [67.4, 70.4, 69.5], "Main": [64.3, 67.6, 66.5], "Slough": [87.6, 88.6, 89.2]}),
    ("Mg", "Silicate", "slct_Mg", {"All": [32.4, 29.4, 30.3], "Main": [35.6, 32.2, 33.4], "Slough": [11.9, 10.9, 10.3]}),
    ("Mg", "Precipitation", "prec", {"All": [0.2, 0.2, 0.2], "Main": [0.2, 0.2, 0.2], "Slough": [0.5, 0.5, 0.5]}),
    ("Na", "Silicate", "slct_Na", {"All": [94.2, 94.2, 94.2], "Main": [95.2, 95.2, 95.1], "Slough": [88.0, 88.0, 88.2]}),
    ("Na", "Precipitation", "prec", {"All": [5.8, 5.8, 5.8], "Main": [4.8, 4.8, 4.9], "Slough": [12.0, 12.0, 11.8]}),
    ("K", "Silicate", "slct_K", {"All": [90.8, 90.6, 90.8], "Main": [90.8, 90.5, 90.7], "Slough": [91.5, 91.6, 91.6]}),
    ("K", "Precipitation", "prec", {"All": [9.2, 9.4, 9.2], "Main": [9.2, 9.5, 9.3], "Slough": [8.5, 8.4, 8.4]}),
    ("Cl", "Precipitation", "prec", {"All": [100.0, 100.0, 100.0], "Main": [100.0, 100.0, 100.0], "Slough": [100.0, 100.0, 100.0]}),
    ("SO4", "H2SO4 production", "pyri", {"All": [79.4, 79.6, 79.7], "Main": [80.7, 80.7, 81.0], "Slough": [71.2, 71.9, 70.7]}),
    ("SO4", "Evaporite", "evap", {"All": [20.0, 19.8, 19.8], "Main": [19.0, 18.9, 18.6], "Slough": [26.5, 25.7, 27.2]}),
    ("SO4", "Precipitation", "prec", {"All": [0.6, 0.6, 0.6], "Main": [0.3, 0.3, 0.3], "Slough": [2.2, 2.2, 2.1]}),
]


def classify():
    p = find_scenario_parameters(SCEN[1], UE)
    river = read_river_observations(RD, load_conc2equi(UE))
    mask = build_model_variables(river, p.ObsInNormalization, p.EMUnits, p.ObsList)
    ws = openpyxl.load_workbook(RD, read_only=True)["MyRiverData"]
    rows = list(ws.iter_rows(values_only=True))
    ti = rows[0].index("ExtraField1")
    types = [str(r[ti]).strip() if r[ti] is not None else "" for r in rows[1:]]
    func = np.where(mask)[0]
    main = [int(i) for i in func if types[i] == "river"]
    slough = [int(i) for i in func if types[i] == "slough"]
    return main, slough


def load_cache():
    return json.load(open(CACHE)) if os.path.exists(CACHE) else {}


def save_cache(c):
    json.dump(c, open(CACHE + ".tmp", "w"))
    os.replace(CACHE + ".tmp", CACHE)


def run_sample(sc, sample, success, max_iter):
    _r, ctx = run(SCEN[sc], UE, RD, seed=1, max_success=success, max_iter=max_iter,
                  max_zerohits=4000, samples=[sample])
    S = ctx["summary"]
    out = {}
    for ion, _lab, em, _pub in PUB:
        if em not in S.fraction[ion]:
            out[f"{ion}|{em}"] = None
            continue
        med = S.fraction[ion][em]["median"][sample]
        out[f"{ion}|{em}"] = None if np.isnan(med) else float(med) * 100
    return out


def aggregate(cache, sc, samples):
    res = {}
    for ion, _lab, em, _pub in PUB:
        vals = [cache.get(f"{sc}|{s}", {}).get(f"{ion}|{em}") for s in samples]
        vals = [v for v in vals if v is not None]
        res[(ion, em)] = float(np.mean(vals)) if vals else float("nan")
    return res


def emit(cache, main, slough):
    groups = {"All samples": main + slough, "Mainstem": main, "Slough": slough}
    lines = ["## Table S2 comparison — all groups and scenarios\n",
             f"Mainstem = {len(main)} `river` samples, Slough = {len(slough)} "
             "`slough` samples; All = the two together (the paper's grouping). "
             "Values are mean-of-median % contributions; **py** is this port, "
             "**Δ** = py − published.\n"]
    for gname, gidx in groups.items():
        agg = {sc: aggregate(cache, sc, gidx) for sc in SCEN}
        lines.append(f"\n### {gname} (n={len(gidx)})\n")
        lines.append("| ion | end-member | S1 pub | S1 py | S1 Δ | S2 pub | S2 py | S2 Δ | S3 pub | S3 py | S3 Δ |")
        lines.append("|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
        for ion, lab, em, pub in PUB:
            cells = [ion, lab]
            for sc in (1, 2, 3):
                t = pub[{"All samples": "All", "Mainstem": "Main", "Slough": "Slough"}[gname]][sc - 1]
                v = agg[sc][(ion, em)]
                if t is None:
                    cells += ["—", f"{v:.1f}" if not np.isnan(v) else "—", "—"]
                else:
                    cells += [f"{t:.1f}", f"{v:.1f}" if not np.isnan(v) else "—",
                              f"{v - t:+.1f}" if not np.isnan(v) else "—"]
            lines.append("| " + " | ".join(cells) + " |")
    open(OUT, "w").write("\n".join(lines) + "\n")
    print(f"wrote {OUT}")
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--success", type=int, default=60)
    ap.add_argument("--max-iter", type=int, default=8000)
    args = ap.parse_args()
    main_idx, slough_idx = classify()
    samples = main_idx + slough_idx
    cache = load_cache()
    for sc in SCEN:
        for s in samples:
            key = f"{sc}|{s}"
            if key in cache:
                continue
            cache[key] = run_sample(sc, s, args.success, args.max_iter)
            save_cache(cache)
            print(f"{key} done", flush=True)
    emit(cache, main_idx, slough_idx)


if __name__ == "__main__":
    main()
