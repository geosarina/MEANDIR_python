"""Aggregate inversion results — ports the summarizing scripts of MEANDIR_Master.

Covers ``UnpackInversionResults`` / ``SaveFractionalContributions`` /
``CalculateMassBalance`` / ``CalculateReconstructedObservations`` /
``CalculateExcessSO4``: turns the per-success records collected by the engine
into per-sample percentile summaries.

MATLAB's ``prctile`` uses Hazen plotting positions ``(k-0.5)/N``; NumPy's
``method="hazen"`` reproduces it, so percentiles match the reference.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..constants import ISOTOPE_VARIABLES, isotope_ion

_STATS = ("median", "mean", "pct05", "pct25", "pct75", "pct95")


def _summ(mat, axis):
    """median/mean/percentiles over ``axis`` (NaN-aware), matching MATLAB."""
    with np.errstate(all="ignore"):
        out = {
            "median": np.nanmedian(mat, axis=axis),
            "mean": np.nanmean(mat, axis=axis),
        }
        for name, p in (("pct05", 5), ("pct25", 25), ("pct75", 75), ("pct95", 95)):
            out[name] = np.nanpercentile(mat, p, axis=axis, method="hazen")
    return out


@dataclass
class RiverResults:
    params: object
    n_samples: int
    sample_indices: list
    successes: dict                       # sample -> n successful instances
    fraction: dict = field(default_factory=dict)        # target -> em -> stat -> array(n_samples)
    massbalance: dict = field(default_factory=dict)     # obs -> stat -> array(n_samples)
    reconstructed: dict = field(default_factory=dict)   # obs -> stat -> array(n_samples)
    excess_so4: dict = field(default_factory=dict)
    misfit_model: dict = field(default_factory=dict)    # sample -> array(n instances)
    rzcwy: dict = field(default_factory=dict)           # V -> gross/net -> scaled/unscaled -> stat -> array


def aggregate_results(scenario_results, delta2r):
    p = scenario_results.params
    river = scenario_results.river
    obs = p.ObsList
    ems = p.EMList0
    nEM = len(ems)
    n_samples = len(river.model_variable["norm"])
    sample_idx = sorted(scenario_results.sample_results)

    rr = RiverResults(params=p, n_samples=n_samples, sample_indices=sample_idx,
                      successes={})

    # storage: target -> em -> stat -> array(n_samples) prefilled with NaN
    targets = list(obs) + ["norm"]
    for t in targets:
        rr.fraction[t] = {e: {st: np.full(n_samples, np.nan) for st in _STATS}
                          for e in ems}
    for o in obs:
        rr.massbalance[o] = {st: np.full(n_samples, np.nan) for st in _STATS}
        rr.reconstructed[o] = {st: np.full(n_samples, np.nan) for st in _STATS}

    for i in sample_idx:
        succ = scenario_results.sample_results[i]
        rr.successes[i] = len(succ)
        if not succ:
            continue

        # stack instances: (nEM, ninst) per target; (nOL, ninst) reconstructed
        norm_mat = np.stack([s["X"] for s in succ], axis=1)
        frac_mat = {o: np.stack([s["fractions"][o] for s in succ], axis=1)
                    for o in obs}
        recon_mat = np.stack([s["reconstructed"] for s in succ], axis=1)
        rr.misfit_model[i] = np.array([s["misfit_model"] for s in succ])

        # (a) fractional contributions per end-member
        for j, e in enumerate(ems):
            st = _summ(norm_mat[j, :], axis=0)
            for k in _STATS:
                rr.fraction["norm"][e][k][i] = st[k]
            for o in obs:
                st = _summ(frac_mat[o][j, :], axis=0)
                for k in _STATS:
                    rr.fraction[o][e][k][i] = st[k]

        # (b) mass balance: sum of contributions across end-members
        for o in obs:
            st = _summ(np.sum(frac_mat[o], axis=0), axis=0)
            for k in _STATS:
                rr.massbalance[o][k][i] = st[k]

        # (c) reconstructed observations
        recon_obs = {o: recon_mat[k, :].copy() for k, o in enumerate(obs)}
        for o in obs:
            if o in ISOTOPE_VARIABLES:
                ion = isotope_ion(o, p.carbonisotopematch)
                ratio = recon_obs[o] / recon_obs[ion]
                if o in p.ConvertDelta2RList:
                    ratio = (ratio / delta2r.factor(o) - 1) * 1000
                recon_obs[o] = ratio
            st = _summ(recon_obs[o], axis=0)
            for k in _STATS:
                rr.reconstructed[o][k][i] = st[k]

    _excess_so4(rr, scenario_results)
    return rr


def _excess_so4(rr, scenario_results):
    """CalculateExcessSO4 — only when SO4 is not itself inverted."""
    p = scenario_results.params
    if ("SO4" in p.ObsList or p.NormalizationType == "SO4"
            or not getattr(p, "EMHaveSO4", False)):
        return
    river = scenario_results.river
    so4 = river.model_variable["SO4"]
    for stat in ("absolute", "fraction"):
        rr.excess_so4[stat] = {st: np.full(rr.n_samples, np.nan) for st in _STATS}
    for i in rr.sample_indices:
        succ = scenario_results.sample_results[i]
        if not succ:
            continue
        frac_so4_sum = np.array([np.sum(s["fractions"]["SO4"]) for s in succ])
        f_excess = 1 - frac_so4_sum
        abs_excess = so4[i] * f_excess
        for st, val in (("absolute", abs_excess), ("fraction", f_excess)):
            summ = _summ(val, axis=0)
            for k in _STATS:
                rr.excess_so4[st][k][i] = summ[k]
