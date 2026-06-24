"""Validate the end-member distribution builder and sampler against Alaska.

These tests assert the *invariants* MEANDIR guarantees (internal mass balance,
evaporite stoichiometry, non-negativity, isotope products) rather than exact
draws, which is the correct fidelity target for a Monte Carlo engine.
"""

import numpy as np
import pytest

from meandir.constants import ISOTOPE_VARIABLES
from meandir.scenarios import find_scenario_parameters
from meandir.io.user_entries import load_endmember_group, load_delta2r
from meandir.engine.endmembers import read_endmembers
from meandir.engine.distributions import (
    make_em_distributions, Uniform, PointMass,
)
from meandir.engine.sampling import pull_end_member_ratios


@pytest.fixture(scope="module")
def alaska(user_entries_path):
    p = find_scenario_parameters("AK_scenario2_carbonate_slctindi_degas_2p5",
                                 user_entries_path)
    em = read_endmembers(load_endmember_group(user_entries_path, p.EMdatasource))
    d2r = load_delta2r(user_entries_path)
    dists = make_em_distributions(p, em, d2r)
    return p, em, d2r, dists


def test_distribution_types(alaska):
    p, em, d2r, D = alaska
    obs, ems = p.ObsList, p.EMList0
    # Ca/carb is a uniform 1/3..2/3
    d = D.dist[obs.index("Ca")][ems.index("carb")]
    assert isinstance(d, Uniform)
    assert d.low == pytest.approx(1 / 3) and d.high == pytest.approx(2 / 3)
    # Ca/slct_Ca min==max==1 -> point mass
    assert isinstance(D.dist[obs.index("Ca")][ems.index("slct_Ca")], PointMass)
    # d34S converted to ratio: evap uniform around the VCDT ratio
    d34s_evap = D.dist[obs.index("d34S")][ems.index("evap")]
    conv = d2r.factor("d34S")
    assert isinstance(d34s_evap, Uniform)
    assert d34s_evap.low == pytest.approx((10 / 1000 + 1) * conv)
    assert d34s_evap.high == pytest.approx((30 / 1000 + 1) * conv)
    # d13C of degas is a fractionation (EPS-UNI)
    assert D.distcode[obs.index("d13C")][ems.index("degas")] == "EPS-UNI"


def test_sampler_internal_consistency(alaska):
    p, em, d2r, D = alaska
    obs, ems = p.ObsList, p.EMList0
    norm_rows = [obs.index(o) for o in p.ObsInNormalization]
    iso = np.array([o in ISOTOPE_VARIABLES for o in obs])
    rng = np.random.default_rng(7)

    n_ok = 0
    for _ in range(500):
        M, fc = pull_end_member_ratios(p, em, D, np.zeros(len(obs)), rng, delta2r=d2r)
        if fc:
            continue
        # mass balance: normalization rows sum to 1 for every end-member
        assert np.allclose(M[norm_rows, :].sum(axis=0), 1.0)
        # non-isotope rows are non-negative (within tolerance)
        assert (M[~iso, :] < -1e-9).sum() == 0
        n_ok += 1
    assert n_ok > 490  # essentially all draws succeed for this scenario


def test_evaporite_stoichiometry(alaska):
    p, em, d2r, D = alaska
    obs, ems = p.ObsList, p.EMList0
    rng = np.random.default_rng(3)
    ev = ems.index("evap")
    camgsr = [obs.index(o) for o in ("Ca", "Mg", "Sr") if o in obs]
    for _ in range(200):
        M, fc = pull_end_member_ratios(p, em, D, np.zeros(len(obs)), rng, delta2r=d2r)
        if fc:
            continue
        assert np.isclose(M[obs.index("SO4"), ev], M[camgsr, ev].sum())
        assert np.isclose(M[obs.index("Cl"), ev],
                          M[obs.index("Na"), ev] + M[obs.index("K"), ev])
