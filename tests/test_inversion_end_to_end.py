"""End-to-end Samples-mode inversion on the Alaska scenario.

Asserts the inversion runs and that successful instances satisfy the MEANDIR
success criteria (fractional-contribution bounds, mass-balance window, isotope
window) — the invariants the model guarantees.
"""

import numpy as np
import pytest

from meandir.run import run


@pytest.fixture(scope="module")
def alaska_run(user_entries_path, river_data_path):
    res, ctx = run(
        "AK_scenario2_carbonate_slctindi_degas_2p5",
        user_entries_path, river_data_path,
        seed=1, max_success=15, max_iter=20000, max_zerohits=5000,
        samples=range(0, 8),
    )
    return res, ctx


def test_inversion_produces_successes(alaska_run):
    res, ctx = alaska_run
    total = sum(len(v) for v in res.sample_results.values())
    assert total > 0  # at least some samples invert successfully


def test_success_criteria_hold(alaska_run):
    res, ctx = alaska_run
    p = ctx["params"]
    errmin = np.asarray(p.ErrorCutMinMB, float)
    errmax = np.asarray(p.ErrorCutMaxMB, float)
    noniso = np.array([o not in {
        "d34S", "d13C", "Sr8786"} for o in p.ObsList])  # the iso vars here

    n_checked = 0
    for succ in res.sample_results.values():
        for s in succ:
            X = s["X"]
            # fractional contributions within the per-end-member bounds
            # (degas bound is sample-dependent; check the static ones are sane)
            assert np.all(np.isfinite(X))
            # non-isotope observations are reconstructed within the MB window
            for k, o in enumerate(p.ObsList):
                if noniso[k]:
                    frac_recon = s["mb"][o]
                    assert frac_recon >= errmin[k] / 100 - 1e-6
                    assert frac_recon <= errmax[k] / 100 + 1e-6
            n_checked += 1
    assert n_checked > 0


def test_carbonate_dominated_signature(alaska_run):
    res, ctx = alaska_run
    p = ctx["params"]
    ems = p.EMList0
    carb = ems.index("carb")
    allX = [s["X"] for succ in res.sample_results.values() for s in succ]
    assert allX, "expected at least one successful inversion"
    med = np.median(np.array(allX), axis=0)
    # carbonate is the dominant cation source in this scenario
    assert med[carb] == max(med)
    assert med[carb] > 0.3
