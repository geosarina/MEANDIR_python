"""Monte Carlo driver — port of the Samples-mode loop of ``MEANDIR_Master.m``.

For each sample with complete data, repeatedly: draw the river observations
within analytical error, draw an internally-consistent end-member matrix, reduce
for zero entries, locate fractionation pairs, optionally reset the degassing
bound, invert, and evaluate. Successful instances are collected until
``max_success`` is reached or the iteration budget is exhausted.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .runtime_data import adjust_river_data, prepare_updated_data
from .sampling import pull_end_member_ratios
from .fractionation import find_fractionation_pairs
from .inversion import invert_active_simulation
from .evaluate import evaluate_inversion_instance


@dataclass
class ScenarioResults:
    params: object
    river: object
    sample_results: dict = field(default_factory=dict)   # sample index -> [SimResult]
    iterations: dict = field(default_factory=dict)        # sample index -> attempts


def _reset_degas(params, rc):
    minf = np.array(params.MinFractionalContribution0, dtype=float)
    maxf = np.array(params.MaxFractionalContribution0, dtype=float)
    if params.ResetDegasDICContribution == 1 and "degas" in params.EMList0:
        dic = rc[params.ObsList.index("DIC")]
        di = params.EMList0.index("degas")
        minf[di] = params.DegasDICContributionMin * dic
        maxf[di] = params.DegasDICContributionMax * dic
    return minf, maxf


def run_scenario(params, river, em_data, dists, delta2r, functional_mask, rng,
                 max_iter=None, max_success=None, max_zerohits=None,
                 progress=False):
    """Run the Samples-mode inversion. Returns :class:`ScenarioResults`."""
    if params.IterateOver != "Samples":
        raise NotImplementedError(
            "run_scenario currently implements IterateOver='Samples' "
            "(the Alaska scenarios). End-members mode is a later phase.")

    obs = params.ObsList
    ems = params.EMList0
    solver = params.Solver
    relpos = np.array([c == "rel" for c in params.CostFunType])
    weighting = np.asarray(params.WeightingList, dtype=float)
    abspos = np.array([c == "abs" for c in params.CostFunType])
    sources = np.array([e in params.EMsources for e in ems])

    max_iter = int(params.maxiterations if max_iter is None else max_iter)
    max_success = int(params.maxsuccess if max_success is None else max_success)
    max_zerohits = int(params.maxzerohits if max_zerohits is None else max_zerohits)

    results = ScenarioResults(params=params, river=river)
    sample_indices = np.where(functional_mask)[0]

    for i in sample_indices:
        success = []
        it = 0
        while it < max_iter and len(success) < max_success:
            it += 1
            inv, rc0 = adjust_river_data(river, params, i, rng, delta2r)
            em_inst0, failcase = pull_end_member_ratios(
                params, em_data, dists, rc0, rng, delta2r=delta2r, sample_index=i)
            if failcase or np.any(np.isnan(em_inst0)):
                if len(success) == 0 and it >= max_zerohits:
                    break
                continue

            xdirect = np.full(len(ems), np.nan)
            rc_r, em_r, ems_r, obs_r, xdirect, distcode_r = prepare_updated_data(
                xdirect, rc0, obs, ems, ems, relpos, em_inst0, solver, dists.distcode)
            frac = find_fractionation_pairs(
                dists.distcode, distcode_r, obs, obs_r, ems, ems_r,
                params.carbonisotopematch)

            minf, maxf = _reset_degas(params, rc0)
            nanmask = np.isnan(xdirect)
            minf_r = minf[nanmask]
            maxf_r = maxf[nanmask]
            solvecf_r = np.array([o not in params.nCFList for o in obs_r])
            in_r = np.array([o in obs_r for o in obs])
            abspos_r = abspos[in_r]
            relpos_r = relpos[in_r]
            weighting_r = weighting[in_r]

            X, em_updated, fcost = invert_active_simulation(
                solver, em_inst0, em_r, rc0, rc_r, ems, xdirect,
                minf_r, maxf_r, solvecf_r, abspos_r, relpos_r, weighting_r,
                frac, sources)

            res = evaluate_inversion_instance(
                params, X, em_inst0, em_updated, rc0, inv, fcost,
                minf, maxf, delta2r)
            if res is not None:
                success.append(res)

            if len(success) == 0 and it >= max_zerohits:
                break

        results.sample_results[int(i)] = success
        results.iterations[int(i)] = it
        if progress:
            print(f"sample {int(i)}: {len(success)} successes in {it} iters")

    return results
