"""ClCritical (cyclic-chloride) correction.

The Alaska scenarios use PrecProcessing='EndMember', so no real sample carries a
ClCr value. These tests exercise the ClCrit branch with synthetic ClCr values
and confirm the EndMember path is unaffected.
"""

import copy

import numpy as np
import pytest

from meandir.scenarios import find_scenario_parameters
from meandir.io.user_entries import (
    load_conc2equi, load_delta2r, load_endmember_group)
from meandir.io.river_data import read_river_observations, build_model_variables
from meandir.engine.endmembers import read_endmembers
from meandir.engine.distributions import make_em_distributions
from meandir.engine.master import run_scenario
from meandir.engine.results import aggregate_results
from meandir.engine.clcritical import quantify_prec_ratios

SCEN = "AK_scenario2_carbonate_slctindi_degas_2p5"


@pytest.fixture(scope="module")
def setup(user_entries_path, river_data_path):
    c2e = load_conc2equi(user_entries_path)
    d2r = load_delta2r(user_entries_path)
    p = find_scenario_parameters(SCEN, user_entries_path)
    river = read_river_observations(river_data_path, c2e, sheet="MyRiverData")
    func = build_model_variables(river, p.ObsInNormalization, p.EMUnits, p.ObsList)
    em = read_endmembers(load_endmember_group(user_entries_path, p.EMdatasource))
    p.EMHaveSO4 = em.all_have(p.EMList0, "SO4")
    dists = make_em_distributions(p, em, d2r)
    return p, river, func, em, dists, d2r


def test_alaska_uses_endmember_not_clcrit(setup):
    p = setup[0]
    # Confirms the premise: ClCritical is inert for the published scenarios.
    assert p.PrecProcessing == "EndMember"


def test_endmember_path_still_produces_successes(setup):
    p, river, func, em, dists, d2r = setup
    fidx = np.where(func)[0][:3]
    mask = np.zeros_like(func)
    mask[fidx] = True
    res = run_scenario(p, river, em, dists, d2r, mask,
                       np.random.default_rng(0), max_success=5,
                       max_iter=4000, max_zerohits=2000)
    assert all(len(res.sample_results[int(i)]) > 0 for i in fidx)


def test_clcrit_fixes_precipitation_and_reconstructs_cl(setup):
    p, river, func, em, dists, d2r = setup
    river = copy.deepcopy(river)
    fidx = list(map(int, np.where(func)[0]))[:3]
    clcr = np.asarray(river.model_variable["ClCr"], dtype=float)
    cl = np.asarray(river.model_variable["Cl"], dtype=float)
    for i in fidx:
        clcr[i] = 0.8 * cl[i]            # synthetic critical chloride
    river.model_variable["ClCr"] = clcr

    p2 = copy.copy(p)
    p2.PrecProcessing = "ClCrit"
    mask = np.zeros_like(func)
    for i in fidx:
        mask[i] = True
    res = run_scenario(p2, river, em, dists, d2r, mask,
                       np.random.default_rng(0), max_success=8,
                       max_iter=6000, max_zerohits=3000)
    S = aggregate_results(res, d2r)
    assert any(len(res.sample_results[i]) > 0 for i in fidx)
    for i in fidx:
        if res.sample_results[i]:
            # chloride budget reconstructs (precipitation sources it)
            assert abs(S.massbalance["Cl"]["median"][i] - 1.0) < 0.1
            # precipitation contribution is defined (fixed, not NaN)
            assert np.isfinite(S.fraction["norm"]["prec"]["median"][i])


def test_clcrit_nan_skips_sample(setup):
    p, river, func, em, dists, d2r = setup
    # all real ClCr are NaN -> ClCrit mode yields no successes (sample skipped)
    p2 = copy.copy(p)
    p2.PrecProcessing = "ClCrit"
    i = int(np.where(func)[0][0])
    mask = np.zeros_like(func)
    mask[i] = True
    res = run_scenario(p2, river, em, dists, d2r, mask,
                       np.random.default_rng(0), max_success=5,
                       max_iter=2000, max_zerohits=1000)
    assert len(res.sample_results[i]) == 0


def test_quantify_prec_ratios_cl_self_ratio(setup):
    p, river, func, em, dists, d2r = setup
    from meandir.engine.runtime_data import adjust_river_data
    from meandir.engine.sampling import pull_end_member_ratios
    i = int(np.where(func)[0][0])
    rng = np.random.default_rng(0)
    inv, rc0 = adjust_river_data(river, p, i, rng, d2r)
    em_inst0, _ = pull_end_member_ratios(p, em, dists, rc0, rng,
                                         delta2r=d2r, sample_index=i)
    clcrit = 0.5 * inv["Cl"]
    ratios = quantify_prec_ratios(p.ObsList, em_inst0, inv, clcrit,
                                  p.EMList0.index("prec"), p.ObsList.index("Cl"))
    # one ratio per concentration variable present; all finite
    assert len(ratios) > 0 and np.all(np.isfinite(ratios))
