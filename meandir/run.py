"""Top-level entry point — the Python equivalent of running ``MEANDIR_Master``.

Loads the workbooks, resolves the scenario, builds the river matrix and
end-member distributions, and runs the Monte Carlo inversion.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .scenarios import find_scenario_parameters
from .io.user_entries import (
    load_conc2equi, load_delta2r, load_endmember_group,
)
from .io.river_data import read_river_observations, build_model_variables
from .engine.endmembers import read_endmembers
from .engine.distributions import make_em_distributions
from .engine.master import run_scenario
from .engine.results import aggregate_results
from .engine.rzcwy import calculate_rzcwy


def run(scenario_name: str, user_entries_path, river_data_path, *,
        seed: int = 0, max_iter=None, max_success=None, max_zerohits=None,
        samples=None, river_sheet="MyRiverData", progress=False):
    """Run a scenario end-to-end and return (ScenarioResults, context dict).

    ``samples`` optionally restricts the run to a subset of sample indices
    (useful for quick checks); otherwise every sample with complete data runs.
    """
    user_entries_path = str(user_entries_path)
    river_data_path = str(river_data_path)

    conc2equi = load_conc2equi(user_entries_path)
    delta2r = load_delta2r(user_entries_path)
    params = find_scenario_parameters(scenario_name, user_entries_path)

    river = read_river_observations(river_data_path, conc2equi, sheet=river_sheet)
    functional = build_model_variables(
        river, params.ObsInNormalization, params.EMUnits, params.ObsList)

    em_table = load_endmember_group(user_entries_path, params.EMdatasource)
    em_data = read_endmembers(em_table)
    params.EMHaveSO4 = em_data.all_have(params.EMList0, "SO4")
    dists = make_em_distributions(params, em_data, delta2r)

    if samples is not None:
        mask = np.zeros_like(functional)
        mask[list(samples)] = functional[list(samples)]
        functional = mask

    rng = np.random.default_rng(seed)
    results = run_scenario(
        params, river, em_data, dists, delta2r, functional, rng,
        max_iter=max_iter, max_success=max_success, max_zerohits=max_zerohits,
        progress=progress)

    summary = aggregate_results(results, delta2r)
    summary.rzcwy = calculate_rzcwy(results, conc2equi)

    context = {
        "params": params, "river": river, "em_data": em_data,
        "dists": dists, "delta2r": delta2r, "functional": functional,
        "summary": summary,
    }
    return results, context
