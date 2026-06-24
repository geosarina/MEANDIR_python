"""Diagnostic: for the square scenario-1 system, is the optimizer actually
moving the solution, or is X0 (the unique exact solve) already feasible?

If X0 is almost always feasible and the optimizer barely moves it, the Ca/Mg
residual lives in the end-member sampling, not the solver. If X0 is often
infeasible, the residual is an SLSQP-vs-fmincon boundary effect.
"""
import numpy as np

from meandir.scenarios import find_scenario_parameters
from meandir.io.user_entries import load_conc2equi, load_delta2r, load_endmember_group
from meandir.io.river_data import read_river_observations, build_model_variables
from meandir.engine.endmembers import read_endmembers
from meandir.engine.distributions import make_em_distributions
from meandir.engine.runtime_data import adjust_river_data, prepare_updated_data
from meandir.engine.sampling import pull_end_member_ratios
from meandir.engine.inversion import invert_active_simulation, _initial_condition
from meandir.engine.fractionation import find_fractionation_pairs

UE = "reference/data/MEANDIR_UserEntries.xlsx"
RD = "reference/data/RiverDataSpreadsheet_Kemeny_etal_2023.xlsx"
name = "AK_scenario1_carbonate_slctindi"

delta2r = load_delta2r(UE)
params = find_scenario_parameters(name, UE)
river = read_river_observations(RD, load_conc2equi(UE))
funct = build_model_variables(river, params.ObsInNormalization, params.EMUnits, params.ObsList)
em_data = read_endmembers(load_endmember_group(UE, params.EMdatasource))
params.EMHaveSO4 = em_data.all_have(params.EMList0, "SO4")
dists = make_em_distributions(params, em_data, delta2r)

ems = params.EMList0
obs = params.ObsList
relpos = np.array([c == "rel" for c in params.CostFunType])
abspos = np.array([c == "abs" for c in params.CostFunType])
weighting = np.asarray(params.WeightingList, float)
sources = np.array([e in params.EMsources for e in ems])
minf0 = np.array(params.MinFractionalContribution0, float)
maxf0 = np.array(params.MaxFractionalContribution0, float)
rng = np.random.default_rng(7)

ci = ems.index("carb")
si = ems.index("slct_Mg")
n_inf = 0
n_tot = 0
dx_list = []
split_shift = []   # change in (carb - slct_Mg) X from X0 to final

for i in np.where(funct)[0]:
    got = 0
    for _ in range(4000):
        if got >= 25:
            break
        inv, rc0 = adjust_river_data(river, params, i, rng, delta2r)
        em_inst0, fail = pull_end_member_ratios(params, em_data, dists, rc0, rng,
                                                delta2r=delta2r, sample_index=i)
        if fail or np.any(np.isnan(em_inst0)):
            continue
        xd = np.full(len(ems), np.nan)
        rc_r, em_r, ems_r, obs_r, xd, dcode_r = prepare_updated_data(
            xd, rc0, obs, ems, ems, relpos, em_inst0, params.Solver, dists.distcode)
        if np.any(np.isnan(rc_r)) or em_r.size == 0:
            continue
        frac = find_fractionation_pairs(dists.distcode, dcode_r, obs, obs_r, ems, ems_r,
                                        params.carbonisotopematch)
        nan = np.isnan(xd)
        minf_r, maxf_r = minf0[nan], maxf0[nan]
        X0 = _initial_condition(params.Solver, em_r.copy(), rc_r)
        viol = np.maximum(np.maximum(minf_r - X0, X0 - maxf_r), 0)
        infeasible = np.any(viol > 1e-9)
        in_r = np.array([o in obs_r for o in obs])
        X, _, _ = invert_active_simulation(
            params.Solver, em_inst0, em_r, rc0, rc_r, ems, xd, minf_r, maxf_r,
            np.array([o not in params.nCFList for o in obs_r]),
            abspos[in_r], relpos[in_r], weighting[in_r], frac, sources)
        # map reduced X0 back to full index for carb/slct comparison
        full0 = np.full(len(ems), np.nan); full0[nan] = X0
        n_tot += 1
        got += 1
        if infeasible:
            n_inf += 1
        dx_list.append(np.nanmax(np.abs(X - full0)))
        split_shift.append((X[ci] - X[si]) - (full0[ci] - full0[si]))

dx = np.array(dx_list); ss = np.array(split_shift)
print(f"instances: {n_tot}")
print(f"X0 infeasible (optimizer must act): {n_inf}/{n_tot} = {100*n_inf/n_tot:.0f}%")
print(f"max |X_final - X0| : median={np.nanmedian(dx):.4f}  90th={np.nanpercentile(dx,90):.4f}")
print(f"optimizer shift of (carb - slct_Mg): median={np.nanmedian(ss):+.4f}  mean={np.nanmean(ss):+.4f}")
