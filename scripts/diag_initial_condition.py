"""Diagnostic: inspect the inversion linear system for one instance and compare
how np.linalg.lstsq (minimum-norm) vs a MATLAB-mldivide basic solution behave as
the optimizer's initial condition.
"""
import numpy as np
import scipy.linalg as sla
from scipy.optimize import nnls

from meandir.run import run  # noqa
from meandir.scenarios import find_scenario_parameters
from meandir.io.user_entries import load_conc2equi, load_delta2r, load_endmember_group
from meandir.io.river_data import read_river_observations, build_model_variables
from meandir.engine.endmembers import read_endmembers
from meandir.engine.distributions import make_em_distributions
from meandir.engine.runtime_data import adjust_river_data, prepare_updated_data
from meandir.engine.sampling import pull_end_member_ratios

UE = "reference/data/MEANDIR_UserEntries.xlsx"
RD = "reference/data/RiverDataSpreadsheet_Kemeny_etal_2023.xlsx"
name = "AK_scenario2_carbonate_slctindi_degas_2p5"

delta2r = load_delta2r(UE)
params = find_scenario_parameters(name, UE)
river = read_river_observations(RD, load_conc2equi(UE))
build_model_variables(river, params.ObsInNormalization, params.EMUnits, params.ObsList)
em_data = read_endmembers(load_endmember_group(UE, params.EMdatasource))
params.EMHaveSO4 = em_data.all_have(params.EMList0, "SO4")
dists = make_em_distributions(params, em_data, delta2r)

ems = params.EMList0
obs = params.ObsList
relpos = np.array([c == "rel" for c in params.CostFunType])
rng = np.random.default_rng(1)

# pull one valid instance for sample 3
i = 3
for _ in range(50):
    inv, rc0 = adjust_river_data(river, params, i, rng, delta2r)
    em_inst0, fail = pull_end_member_ratios(params, em_data, dists, rc0, rng,
                                            delta2r=delta2r, sample_index=i)
    if not fail and not np.any(np.isnan(em_inst0)):
        break
xdirect = np.full(len(ems), np.nan)
rc_r, em_r, ems_r, obs_r, xdirect, distcode_r = prepare_updated_data(
    xdirect, rc0, obs, ems, ems, relpos, em_inst0, params.Solver, dists.distcode)

A = em_r
b = rc_r
print(f"system shape A = {A.shape}  (obs_r={len(obs_r)}, em_r={len(ems_r)})")
print("rank(A) =", np.linalg.matrix_rank(A))
print("ems_r =", ems_r)

# (1) numpy lstsq -> minimum-norm
x_ls, *_ = np.linalg.lstsq(A, b, rcond=None)
# (2) MATLAB-style basic solution via column-pivoted QR
Q, R, piv = sla.qr(A, mode="economic", pivoting=True)
rank = np.linalg.matrix_rank(A)
x_basic = np.zeros(A.shape[1])
# solve using the first `rank` pivoted columns
Rb = R[:rank, :rank]
Qb = Q[:, :rank]
x_basic[piv[:rank]] = sla.solve_triangular(Rb, Qb.T @ b)

print("\n               " + "  ".join(f"{e:>8s}" for e in ems_r))
print("lstsq (min-norm):", "  ".join(f"{v:8.3f}" for v in x_ls))
print("mldivide (basic):", "  ".join(f"{v:8.3f}" for v in x_basic))
print("\nresidual norms:  lstsq=%.2e  basic=%.2e" % (
    np.linalg.norm(A @ x_ls - b), np.linalg.norm(A @ x_basic - b)))
print("nonzeros: lstsq=%d  basic=%d" % (np.sum(np.abs(x_ls) > 1e-9),
                                        np.sum(np.abs(x_basic) > 1e-9)))
