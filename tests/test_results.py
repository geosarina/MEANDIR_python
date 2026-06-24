"""Results aggregation: percentile summaries, mass balance, reconstruction."""

import numpy as np
import pytest

from meandir.run import run


@pytest.fixture(scope="module")
def summary(user_entries_path, river_data_path):
    res, ctx = run(
        "AK_scenario2_carbonate_slctindi_degas_2p5",
        user_entries_path, river_data_path,
        seed=3, max_success=25, max_iter=20000, max_zerohits=5000,
        samples=range(0, 8),
    )
    return ctx["summary"], ctx


def test_fraction_percentiles_ordered(summary):
    S, ctx = summary
    i = next(k for k, v in S.successes.items() if v > 0)
    for em in ctx["params"].EMList0:
        f = S.fraction["norm"][em]
        assert f["pct05"][i] <= f["pct25"][i] <= f["median"][i] <= f["pct75"][i] <= f["pct95"][i] + 1e-9


def test_mass_balance_near_one(summary):
    S, ctx = summary
    # every inverted non-isotope observation reconstructs to ~100%
    for i in S.sample_indices:
        if S.successes[i] == 0:
            continue
        for o in ("Ca", "Mg", "Na", "K", "SO4", "DIC"):
            assert abs(S.massbalance[o]["median"][i] - 1.0) < 0.05


def test_reconstructed_isotopes_track_observations(summary):
    S, ctx = summary
    mv = ctx["river"].model_variable
    for i in S.sample_indices:
        if S.successes[i] == 0:
            continue
        # reconstructed delta values land within the inversion window of the obs
        assert abs(S.reconstructed["d13C"]["median"][i] - mv["d13C"][i]) < 2.0
        assert abs(S.reconstructed["d34S"]["median"][i] - mv["d34S"][i]) < 3.0


def test_rzcwy_sane_and_net_equals_gross(summary):
    S, ctx = summary
    Z = S.rzcwy
    assert set("RZCWY").issubset(Z)
    i = next(k for k, v in S.successes.items() if v > 0)
    # carbonate-dominated catchment: most cation weathering is carbonate (R high)
    assert 0.0 <= Z["R"]["gross"]["unscaled"]["median"][i] <= 1.5
    # most river SO4 is pyrite-derived in this scenario
    assert Z["Y"]["gross"]["unscaled"]["median"][i] > 0.3
    # no sinks in the Alaska scenarios -> net recycling factor is 1
    for V in "RZCWY":
        g = Z[V]["gross"]["unscaled"]["median"][i]
        nt = Z[V]["net"]["unscaled"]["median"][i]
        if not (np.isnan(g) or np.isnan(nt)):
            assert abs(g - nt) < 1e-9
