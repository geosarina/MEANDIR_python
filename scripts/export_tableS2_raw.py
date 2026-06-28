"""Capture and export the Table S2 validation at three granularities:

  validation_outputs/tableS2_summary_200successes.csv       (mean-of-median per group/scenario/cell)
  validation_outputs/tableS2_persample_200successes.csv     (per-sample percentiles)
  validation_outputs/tableS2_persimulation_200successes.csv (every successful simulation)

It re-runs the inversion capturing every successful simulation's fractional
contribution for each of the 18 Table S2 (ion, end-member) cells, caching the
raw per-simulation values per (scenario, sample) so a restart resumes.

  python -m scripts.export_tableS2_raw --success 200       # run + cache
  python -m scripts.export_tableS2_raw --write-only        # CSVs from cache
"""
import argparse
import csv
import json
import os

import numpy as np
import openpyxl

from meandir.run import run
from meandir.scenarios import find_scenario_parameters
from meandir.io.user_entries import load_conc2equi
from meandir.io.river_data import read_river_observations, build_model_variables
from scripts.validate_groups import PUB, SCEN

UE = "reference/data/MEANDIR_UserEntries.xlsx"
RD = "reference/data/RiverDataSpreadsheet_Kemeny_etal_2023.xlsx"
RAW = "/tmp/claude-0/-home-user-MEANDIR-python/feec57dc-e16b-557c-afd4-47f9242d36db/scratchpad/tableS2_raw_s200.json"
OUTDIR = "validation_outputs"
DEGAS = {1: "none", 2: "<2.5x DIC", 3: "<25x DIC"}
CELLS = [(ion, lab, em) for ion, lab, em, _pub in PUB]          # 18 cells, in order
CELLKEY = [f"{ion}|{lab}" for ion, lab, _em in CELLS]


def river_names():
    ws = openpyxl.load_workbook(RD, read_only=True)["MyRiverData"]
    rows = list(ws.iter_rows(values_only=True))
    ni = rows[0].index("Name"); ti = rows[0].index("ExtraField1")
    names = [r[ni] for r in rows[1:]]
    types = [str(r[ti]).strip() if r[ti] is not None else "" for r in rows[1:]]
    return names, types


def classify():
    p = find_scenario_parameters(SCEN[1], UE)
    river = read_river_observations(RD, load_conc2equi(UE))
    mask = build_model_variables(river, p.ObsInNormalization, p.EMUnits, p.ObsList)
    names, types = river_names()
    func = np.where(mask)[0]
    main = [int(i) for i in func if types[i] == "river"]
    slough = [int(i) for i in func if types[i] == "slough"]
    return main, slough, names


def load_raw():
    return json.load(open(RAW)) if os.path.exists(RAW) else {}


def save_raw(c):
    json.dump(c, open(RAW + ".tmp", "w"))
    os.replace(RAW + ".tmp", RAW)


def run_sample(sc, sample, success, max_iter):
    r, ctx = run(SCEN[sc], UE, RD, seed=1, max_success=success, max_iter=max_iter,
                 max_zerohits=4000, samples=[sample])
    ems = ctx["params"].EMList0
    sims = []
    for s in r.sample_results.get(sample, []):
        row = []
        for ion, _lab, em in CELLS:
            if em in ems:   # degas is absent in scenario 1 -> None
                row.append(round(float(s["fractions"][ion][ems.index(em)]) * 100, 4))
            else:
                row.append(None)
        sims.append(row)
    return sims


def run_and_cache(success, max_iter):
    main, slough, _ = classify()
    samples = main + slough
    cache = load_raw()
    for sc in SCEN:
        for s in samples:
            key = f"{sc}|{s}"
            if key in cache:
                continue
            cache[key] = run_sample(sc, s, success, max_iter)
            save_raw(cache)
            print(f"{key}: {len(cache[key])} sims", flush=True)
    return cache


def write_csvs(success):
    cache = load_raw()
    main, slough, names = classify()
    group_of = {**{s: "Mainstem" for s in main}, **{s: "Slough" for s in slough}}
    os.makedirs(OUTDIR, exist_ok=True)

    # (1) per-simulation (long)
    p_sim = os.path.join(OUTDIR, f"tableS2_persimulation_s{success}.csv")
    with open(p_sim, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["scenario", "degassing_constraint", "group", "sample_index",
                    "river_name", "sim_index", "ion", "end_member", "contribution_pct"])
        for sc in SCEN:
            for s in (main + slough):
                sims = cache.get(f"{sc}|{s}", [])
                for j, row in enumerate(sims):
                    for ci, (ion, lab, _em) in enumerate(CELLS):
                        if row[ci] is None:   # cell absent (degas in S1)
                            continue
                        w.writerow([f"S{sc}", DEGAS[sc], group_of[s], s, names[s],
                                    j, ion, lab, row[ci]])

    # (2) per-sample (percentiles)
    p_smp = os.path.join(OUTDIR, f"tableS2_persample_s{success}.csv")
    with open(p_smp, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["scenario", "degassing_constraint", "group", "sample_index",
                    "river_name", "ion", "end_member", "n_successes",
                    "median_pct", "p05_pct", "p25_pct", "p75_pct", "p95_pct"])
        for sc in SCEN:
            for s in (main + slough):
                sims = cache.get(f"{sc}|{s}", [])
                n = len(sims)
                for ci, (ion, lab, _em) in enumerate(CELLS):
                    col = np.array([r[ci] for r in sims if r[ci] is not None], float)
                    if not col.size:   # cell absent (degas in S1)
                        continue
                    def q(p):
                        return round(float(np.percentile(col, p)), 4)
                    w.writerow([f"S{sc}", DEGAS[sc], group_of[s], s, names[s],
                                ion, lab, n, round(float(np.median(col)), 4),
                                q(5), q(25), q(75), q(95)])

    # (3) summary (mean-of-median per group)
    groups = {"All samples": main + slough, "Mainstem": main, "Slough": slough}
    gk = {"All samples": "All", "Mainstem": "Main", "Slough": "Slough"}
    p_sum = os.path.join(OUTDIR, f"tableS2_summary_s{success}.csv")
    with open(p_sum, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["group", "n_samples", "scenario", "degassing_constraint",
                    "ion", "end_member", "published_pct", "python_pct", "delta_pct"])
        for gname, gidx in groups.items():
            for sc in SCEN:
                for ci, (ion, lab, _em) in enumerate(CELLS):
                    permed = []
                    for s in gidx:
                        col = [r[ci] for r in cache.get(f"{sc}|{s}", []) if r[ci] is not None]
                        if col:
                            permed.append(np.median(col))
                    py = float(np.mean(permed)) if permed else float("nan")
                    t = PUB[ci][3][gk[gname]][sc - 1]
                    w.writerow([gname, len(gidx), f"S{sc}", DEGAS[sc], ion, lab,
                                "" if t is None else t,
                                "" if py != py else round(py, 2),
                                "" if (t is None or py != py) else round(py - t, 2)])
    for p in (p_sum, p_smp, p_sim):
        print(f"wrote {p} ({sum(1 for _ in open(p)) - 1} rows)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--success", type=int, default=200)
    ap.add_argument("--max-iter", type=int, default=15000)
    ap.add_argument("--write-only", action="store_true")
    args = ap.parse_args()
    if not args.write_only:
        run_and_cache(args.success, args.max_iter)
    write_csvs(args.success)


if __name__ == "__main__":
    main()
