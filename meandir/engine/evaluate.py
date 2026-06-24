"""Evaluate one inversion instance — port of ``MEANDIR_EvaluateInversionInstance.m``.

Checks that the fractional contributions lie in the user range, reconstructs the
fraction of each observation explained by the end-members, enforces the
mass-balance (and optional normalization) window and the isotope-precision
window, and — on success — returns the per-observation fractional contributions
and misfits. Returns ``None`` for a failed instance.
"""

from __future__ import annotations

import math

import numpy as np

from ..constants import ISOTOPE_VARIABLES, isotope_ion


def evaluate_inversion_instance(params, X, em_inst0, em_updated, river_col0,
                                inv, functioncost, minfrac, maxfrac, delta2r):
    obs = params.ObsList
    nOL = len(obs)
    nEM = len(params.EMList0)
    X = np.asarray(X, dtype=float)
    if np.any(np.isnan(X)):
        return None

    # fractional-contribution range check
    if not np.all((X >= np.asarray(minfrac)) & (X <= np.asarray(maxfrac))):
        return None

    isopos = np.array([o in ISOTOPE_VARIABLES for o in obs])
    errmin = np.asarray(params.ErrorCutMinMB, dtype=float)
    errmax = np.asarray(params.ErrorCutMaxMB, dtype=float)

    # fractional contribution of each end-member to each observation
    fractions = {}
    mb = np.full(nOL, np.nan)
    for k in range(nOL):
        working = X * em_updated[k, :] / river_col0[k]
        fractions[obs[k]] = working
        mb[k] = np.sum(working)

    # optional normalization-sum check row
    mb_check = list(mb)
    iso_check = list(isopos)
    emin_check = list(errmin)
    emax_check = list(errmax)
    if params.NormalizationType == "SumObs" and params.ImposeNormalizationCheck == 1:
        mb_check.append(np.sum(X))
        iso_check.append(False)
        noniso = ~isopos
        emin_check.append(np.min(errmin[noniso]))
        emax_check.append(np.max(errmax[noniso]))
    mb_check = np.array(mb_check)
    iso_check = np.array(iso_check)
    emin_check = np.array(emin_check)
    emax_check = np.array(emax_check)

    noniso = ~iso_check
    within = ((mb_check[noniso] >= emin_check[noniso] / 100)
              & (mb_check[noniso] <= emax_check[noniso] / 100))
    if not np.all(within):
        return None

    # isotope-precision window
    reconstructed = em_updated @ X
    if isopos.any():
        for ll in np.where(isopos)[0]:
            name = obs[ll]
            ion = isotope_ion(name, params.carbonisotopematch)
            ionpos = obs.index(ion)
            recon_iso = reconstructed[ll] / reconstructed[ionpos]
            river_iso = river_col0[ll] / river_col0[ionpos]
            lo = errmin[ll]
            hi = errmax[ll]
            if name in params.ConvertDelta2RList:
                conv = delta2r.factor(name)
                lo = lo / 1000 * conv
                hi = hi / 1000 * conv
            if not ((recon_iso >= river_iso + lo) and (recon_iso <= river_iso + hi)):
                return None

    # success: compute the model-observation misfit
    misfit = _misfit(params, river_col0, reconstructed)
    return {
        "X": X,
        "fractions": fractions,       # obs -> array(nEM)
        "mb": {obs[k]: mb[k] for k in range(nOL)},
        "reconstructed": reconstructed,
        "em_updated": em_updated,
        "em_inst0": em_inst0,
        "river_col0": river_col0,
        "inv": inv,
        "misfit_model": misfit,
        "functioncost": functioncost,
    }


def _misfit(params, river_col0, reconstructed):
    obs = params.ObsList
    solvecf = np.array([o not in params.nCFList for o in obs])
    relpos = np.array([c == "rel" for c in params.CostFunType])
    abspos = np.array([c == "abs" for c in params.CostFunType])
    weighting = np.asarray(params.WeightingList, dtype=float)
    if params.Solver in ("mldivide", "lsqnonneg"):
        misfitvec = (river_col0 - reconstructed) ** 2
        return math.sqrt(np.sum(misfitvec[solvecf]))
    rel = ((river_col0 - reconstructed) / river_col0) ** 2
    ab = (river_col0 - reconstructed) ** 2
    return math.sqrt(np.sum(weighting[solvecf & relpos] * rel[solvecf & relpos])
                     + np.sum(weighting[solvecf & abspos] * ab[solvecf & abspos]))
