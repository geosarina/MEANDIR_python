"""Capture and export the Table S2 validation at multiple granularities.

Captures every successful simulation's FULL end-member -> observation
contribution matrix (resumable per-sample cache), then writes:

  validation_outputs/tableS2_summary_s200.csv        mean-of-median vs published (18 Table S2 cells)
  validation_outputs/tableS2_persample_s200.csv      per-sample percentiles      (18 Table S2 cells)
  validation_outputs/tableS2_persimulation_s200.csv  raw per simulation, long    (18 Table S2 cells)
  validation_outputs/full_persimulation_s200.csv     raw per simulation, wide    (full 9x10 matrix)

  python -m scripts.export_tableS2_raw --success 200    # run + cache, then write
  python -m scripts.export_tableS2_raw --write-only      # write CSVs from cache
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
RAW = "/tmp/claude-0/-home-user-MEANDIR-python/feec57dc-e16b-557c-afd4-47f9242d36db/scratchpad/tableS2_rawfull_s200.json"
OUTDIR = "validation_outputs"
DEGAS = {1: "none", 2: "<2.5x DIC", 3: "<25x DIC"}

# fixed observation order (ObsList) and the union of end-members (EMList0 of S2/S3)
OBS = ["Ca", "Mg", "Na", "K", "Cl", "SO4", "DIC", "d34S", "d13C"]
EMALL = ["prec", "carb", "slct_Ca", "slct_Mg", "slct_Na", "slct_K", "pyri",
         "evap", "corg", "degas"]
CELLS = [(ion, lab, em) for ion, lab, em, _pub in PUB]   # 18 Table S2 cells


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
        flat = []                                  # obs-major, em-minor
        for o in OBS:
            frac = s["fractions"][o]
            for k in range(len(ems)):
                flat.append(round(float(frac[k]) * 100, 4))
        sims.append(flat)
    return {"ems": ems, "sims": sims}


def _cell(sim, ems, ion, em):
    """contribution of `em` to `ion` for one flattened simulation, or None."""
    if em not in ems:
        return None
    return sim[OBS.index(ion) * len(ems) + ems.index(em)]


def run_and_cache(success, max_iter):
    main, slough, _ = classify()
    cache = load_raw()
    for sc in SCEN:
        for s in (main + slough):
            key = f"{sc}|{s}"
            if key in cache:
                continue
            cache[key] = run_sample(sc, s, success, max_iter)
            save_raw(cache)
            print(f"{key}: {len(cache[key]['sims'])} sims", flush=True)
    return cache


def write_csvs(success):
    cache = load_raw()
    main, slough, names = classify()
    group_of = {**{s: "Mainstem" for s in main}, **{s: "Slough" for s in slough}}
    allsamples = main + slough
    os.makedirs(OUTDIR, exist_ok=True)

    def entry(sc, s):
        e = cache.get(f"{sc}|{s}")
        return (e["ems"], e["sims"]) if e else (None, [])

    # (1) Table S2 18-cell: per-simulation (long)
    with open(os.path.join(OUTDIR, f"tableS2_persimulation_s{success}.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["scenario", "degassing_constraint", "group", "sample_index",
                    "river_name", "sim_index", "ion", "end_member", "contribution_pct"])
        for sc in SCEN:
            for s in allsamples:
                ems, sims = entry(sc, s)
                for j, sim in enumerate(sims):
                    for ion, lab, em in CELLS:
                        v = _cell(sim, ems, ion, em)
                        if v is None:
                            continue
                        w.writerow([f"S{sc}", DEGAS[sc], group_of[s], s, names[s], j, ion, lab, v])

    # (2) Table S2 18-cell: per-sample percentiles
    with open(os.path.join(OUTDIR, f"tableS2_persample_s{success}.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["scenario", "degassing_constraint", "group", "sample_index",
                    "river_name", "ion", "end_member", "n_successes",
                    "median_pct", "p05_pct", "p25_pct", "p75_pct", "p95_pct"])
        for sc in SCEN:
            for s in allsamples:
                ems, sims = entry(sc, s)
                for ion, lab, em in CELLS:
                    col = np.array([_cell(sim, ems, ion, em) for sim in sims
                                    if _cell(sim, ems, ion, em) is not None], float)
                    if not col.size:
                        continue
                    w.writerow([f"S{sc}", DEGAS[sc], group_of[s], s, names[s], ion, lab,
                                len(sims), round(float(np.median(col)), 4),
                                *[round(float(np.percentile(col, p)), 4) for p in (5, 25, 75, 95)]])

    # (3) Table S2 18-cell: summary (mean-of-median vs published)
    groups = {"All samples": allsamples, "Mainstem": main, "Slough": slough}
    gk = {"All samples": "All", "Mainstem": "Main", "Slough": "Slough"}
    with open(os.path.join(OUTDIR, f"tableS2_summary_s{success}.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["group", "n_samples", "scenario", "degassing_constraint",
                    "ion", "end_member", "published_pct", "python_pct", "delta_pct"])
        for gname, gidx in groups.items():
            for sc in SCEN:
                for ci, (ion, lab, em) in enumerate(CELLS):
                    permed = []
                    for s in gidx:
                        ems, sims = entry(sc, s)
                        col = [_cell(sim, ems, ion, em) for sim in sims]
                        col = [v for v in col if v is not None]
                        if col:
                            permed.append(np.median(col))
                    py = float(np.mean(permed)) if permed else float("nan")
                    t = PUB[ci][3][gk[gname]][sc - 1]
                    w.writerow([gname, len(gidx), f"S{sc}", DEGAS[sc], ion, lab,
                                "" if t is None else t,
                                "" if py != py else round(py, 2),
                                "" if (t is None or py != py) else round(py - t, 2)])

    # (4) FULL 9x10 matrix: per-simulation (wide)
    cols = [f"{o}.{e}" for o in OBS for e in EMALL]
    with open(os.path.join(OUTDIR, f"full_persimulation_s{success}.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["scenario", "degassing_constraint", "group", "sample_index",
                    "river_name", "sim_index"] + cols)
        for sc in SCEN:
            for s in allsamples:
                ems, sims = entry(sc, s)
                for j, sim in enumerate(sims):
                    vals = []
                    for o in OBS:
                        for e in EMALL:
                            v = _cell(sim, ems, o, e)
                            vals.append("" if v is None else v)
                    w.writerow([f"S{sc}", DEGAS[sc], group_of[s], s, names[s], j] + vals)

    for p in sorted(os.listdir(OUTDIR)):
        if p.endswith(f"s{success}.csv"):
            full = os.path.join(OUTDIR, p)
            print(f"wrote {full} ({sum(1 for _ in open(full)) - 1} rows)")


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
