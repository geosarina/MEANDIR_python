"""Is SLSQP under-converging on the convex scenario-1 problem?

Scenario 1 has no fractionation, so the cost is a bound-constrained linear
least-squares problem: strictly convex, unique global minimum. We compare the
cost achieved by SLSQP vs scipy.optimize.lsq_linear (an exact bounded-LSQ
solver). If lsq_linear reaches a *lower* cost, SLSQP is stopping short of the
true optimum -- and lsq_linear is a cheap, exact fix.
"""
import numpy as np
from scipy.optimize import lsq_linear

from meandir.scenarios import find_scenario_parameters
from meandir.io.user_entries import load_conc2equi, load_delta2r, load_endmember_group
from meandir.io.river_data import read_river_observations, build_model_variables
from meandir.engine.endmembers import read_endmembers
from meandir.engine.distributions import make_em_distributions
from meandir.engine.runtime_data import adjust_river_data, prepare_updated_data
from meandir.engine.sampling import pull_end_member_ratios
from meandir.engine.inversion import invert_active_simulation, cost_function
from meandir.engine.fractionation import find_fractionation_pairs

UE = "reference/data/MEANDIR_UserEntries.xlsx"
RD = "reference/data/RiverDataSpreadsheet_Kemeny_etal_2023.xlsx"
NAME = "AK_scenario1_carbonate_slctindi"

delta2r = load_delta2r(UE)
params = find_scenario_parameters(NAME, UE)
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
ci, si = ems.index("carb"), ems.index("slct_Mg")
rng = np.random.default_rng(11)


def lsq_linear_solve(em_r, rc_r, solvecf_r, relpos_r, abspos_r, w_r, minf_r, maxf_r):
    # rows used in the cost; Alaska uses all-relative
    rel = solvecf_r & relpos_r
    C = (np.sqrt(w_r[rel])[:, None] / rc_r[rel][:, None]) * em_r[rel, :]
    d = np.sqrt(w_r[rel])  # = (sqrt(w)/b)*b
    res = lsq_linear(C, d, bounds=(minf_r, maxf_r), method="bvls")
    return res.x


n = 0
slsqp_cost, lsq_cost = [], []
split_slsqp, split_lsq = [], []
lsq_better = 0
for i in np.where(funct)[0]:
    got = 0
    for _ in range(3000):
        if got >= 30:
            break
        inv, rc0 = adjust_river_data(river, params, i, rng, delta2r)
        em0, fail = pull_end_member_ratios(params, em_data, dists, rc0, rng,
                                           delta2r=delta2r, sample_index=i)
        if fail or np.any(np.isnan(em0)):
            continue
        xd = np.full(len(ems), np.nan)
        rc_r, em_r, ems_r, obs_r, xd, dcode_r = prepare_updated_data(
            xd, rc0, obs, ems, ems, relpos, em0, params.Solver, dists.distcode)
        if np.any(np.isnan(rc_r)) or em_r.size == 0:
            continue
        frac = find_fractionation_pairs(dists.distcode, dcode_r, obs, obs_r, ems, ems_r,
                                        params.carbonisotopematch)
        nan = np.isnan(xd)
        minf_r, maxf_r = minf0[nan], maxf0[nan]
        in_r = np.array([o in obs_r for o in obs])
        scf_r = np.array([o not in params.nCFList for o in obs_r])
        rp_r, ap_r, w_r = relpos[in_r], abspos[in_r], weighting[in_r]

        Xs, _, _ = invert_active_simulation(params.Solver, em0, em_r, rc0, rc_r, ems,
                                            xd, minf_r, maxf_r, scf_r, ap_r, rp_r, w_r,
                                            frac, sources)
        Xl = lsq_linear_solve(em_r, rc_r, scf_r, rp_r, ap_r, w_r, minf_r, maxf_r)

        args = (em0, em_r, rc0, rc_r, scf_r, ap_r, rp_r, w_r, frac, xd, sources)
        cs = cost_function(Xs[nan], *args)
        cl = cost_function(Xl, *args)
        slsqp_cost.append(cs); lsq_cost.append(cl)
        if cl < cs - 1e-12:
            lsq_better += 1
        # full split (carb - slct_Mg)
        Xl_full = np.full(len(ems), np.nan); Xl_full[nan] = Xl
        split_slsqp.append(Xs[ci] - Xs[si])
        split_lsq.append(Xl_full[ci] - Xl_full[si])
        n += 1
        got += 1

sc = np.array(slsqp_cost); lc = np.array(lsq_cost)
print(f"instances: {n}")
print(f"mean cost  SLSQP={sc.mean():.6e}   lsq_linear={lc.mean():.6e}")
print(f"lsq_linear achieves STRICTLY lower cost in {lsq_better}/{n} = {100*lsq_better/n:.0f}% of instances")
print(f"median cost ratio (slsqp/lsq): {np.median(sc/np.maximum(lc,1e-30)):.4f}")
print(f"mean (carb - slct_Mg):  SLSQP={np.mean(split_slsqp):+.4f}   lsq_linear={np.mean(split_lsq):+.4f}")
print(f"  -> lsq_linear shifts split by {np.mean(split_lsq)-np.mean(split_slsqp):+.4f} vs SLSQP")
