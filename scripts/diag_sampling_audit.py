"""Direct sampling audit for the Ca/Mg split (scenario 1).

For the carbonate and silicate Mg ratios (and evaporite Ca/Mg), compare:
  (a) the spreadsheet specification (Min/Max/type),
  (b) the RAW per-cell draw (before closure/rejection),
  (c) the ACCEPTED draw (after mass-balance closure + reject/retry),
and report the reject/retry rate. This tests whether the sampler reproduces the
intended distributions and how much the acceptance conditioning shifts them.
"""
import numpy as np

from meandir.scenarios import find_scenario_parameters
from meandir.io.user_entries import load_conc2equi, load_delta2r, load_endmember_group
from meandir.io.river_data import read_river_observations, build_model_variables
from meandir.engine.endmembers import read_endmembers
from meandir.engine.distributions import make_em_distributions
from meandir.engine.runtime_data import adjust_river_data
from meandir.engine.sampling import pull_end_member_ratios

UE = "reference/data/MEANDIR_UserEntries.xlsx"
RD = "reference/data/RiverDataSpreadsheet_Kemeny_etal_2023.xlsx"
NAME = "AK_scenario1_carbonate_slctindi"

delta2r = load_delta2r(UE)
params = find_scenario_parameters(NAME, UE)
river = read_river_observations(RD, load_conc2equi(UE))
build_model_variables(river, params.ObsInNormalization, params.EMUnits, params.ObsList)
em_data = read_endmembers(load_endmember_group(UE, params.EMdatasource))
params.EMHaveSO4 = em_data.all_have(params.EMList0, "SO4")
dists = make_em_distributions(params, em_data, delta2r)

obs, ems = params.ObsList, params.EMList0
Mg, Ca = obs.index("Mg"), obs.index("Ca")
carb, slctMg, evap = ems.index("carb"), ems.index("slct_Mg"), ems.index("evap")


def spec(em, num):
    out = {}
    for stat in ("Min", "Max", "Men", "Sig", "Mod"):
        if em_data.has(em, num, stat):
            out[stat] = round(float(np.real(em_data.get(em, num, stat))), 5)
    return out


def stats(a):
    a = np.asarray(a, float)
    return (f"median={np.median(a):.4f}  mean={np.mean(a):.4f}  "
            f"min={np.min(a):.4f}  max={np.max(a):.4f}")


print("=== spreadsheet specification (normalized ratios) ===")
for label, (o, e) in [("carb -> Mg", (Mg, carb)), ("slct_Mg -> Mg", (Mg, slctMg)),
                      ("carb -> Ca (closure ion)", (Ca, carb)),
                      ("evap -> Ca", (Ca, evap)), ("evap -> Mg", (Mg, evap))]:
    print(f"  {label:26s} code={dists.distcode[o][e]:7s} spec={spec(ems[e], obs[o])}")

print("\n=== RAW per-cell draws (before closure/rejection), 20000 draws ===")
rng = np.random.default_rng(0)
for label, (o, e) in [("carb -> Mg", (Mg, carb)), ("slct_Mg -> Mg", (Mg, slctMg)),
                      ("evap -> Ca", (Ca, evap)), ("evap -> Mg", (Mg, evap))]:
    raw = np.array([dists.dist[o][e].rvs(rng) for _ in range(20000)])
    print(f"  {label:14s} {stats(raw)}")

print("\n=== ACCEPTED draws (after closure + reject/retry), 8000 matrices ===")
rng = np.random.default_rng(1)
inv, rc0 = adjust_river_data(river, params, 3, rng, delta2r)   # a representative river column
acc = {k: [] for k in ("carb_Mg", "slctMg_Mg", "carb_Ca", "evap_Ca", "evap_Mg")}
tries, nan = [], 0
for _ in range(8000):
    M, fail, tc = pull_end_member_ratios(params, em_data, dists, rc0, rng,
                                         delta2r=delta2r, sample_index=3, return_trycount=True)
    tries.append(tc)
    if np.any(np.isnan(M)):
        nan += 1
        continue
    acc["carb_Mg"].append(M[Mg, carb]);   acc["slctMg_Mg"].append(M[Mg, slctMg])
    acc["carb_Ca"].append(M[Ca, carb]);   acc["evap_Ca"].append(M[Ca, evap])
    acc["evap_Mg"].append(M[Mg, evap])
for label, key in [("carb -> Mg", "carb_Mg"), ("slct_Mg -> Mg", "slctMg_Mg"),
                   ("carb -> Ca", "carb_Ca"), ("evap -> Ca", "evap_Ca"),
                   ("evap -> Mg", "evap_Mg")]:
    print(f"  {label:14s} {stats(acc[key])}")

tries = np.array(tries)
print(f"\n=== reject/retry ===")
print(f"  mean attempts/accept = {tries.mean():.2f}   "
      f"frac needing >1 try = {100*np.mean(tries > 1):.0f}%   "
      f"frac >5 tries = {100*np.mean(tries > 5):.0f}%   NaN(maxtry) = {nan}/8000")
