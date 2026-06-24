"""ClCritical (cyclic-chloride) correction.

Ports ``MEANDIR_ClCriticalCorrection.m`` and ``MEANDIR_QuantifyPrecRatios.m``.
Used only when a scenario sets ``PrecProcessing == 'ClCrit'`` (the published
Alaska scenarios use ``'EndMember'``, so this path is inert for them). It fixes
the precipitation contribution from the chloride budget, subtracts precipitation
from the river column, and removes ``prec`` from the free inversion.
"""

from __future__ import annotations

import numpy as np

from .sampling import pull_end_member_ratios

# Dissolved concentration variables tested for negativity (QuantifyPrecRatios);
# isotopes and the normalization are excluded, in MATLAB order.
_CONC_VARS = ["ALK", "DIC", "Ca", "Mg", "Na", "K", "Sr", "Fe", "Cl", "SO4",
              "NO3", "PO4", "Si", "Ge", "Li", "F", "B", "Re", "Mo", "Os", "HCO3"]


def quantify_prec_ratios(obs, em_inst0, inv, clcrit, prec_idx, cl_idx):
    """Cl removed from precipitation relative to what each ion can support.

    Any entry > 1 means subtracting the ClCritical chloride would drive that
    ion negative.
    """
    cl_river = inv["Cl"]
    cl_test = cl_river if clcrit >= cl_river else clcrit
    ratios = []
    cl_prec = em_inst0[cl_idx, prec_idx]
    for o in _CONC_VARS:
        if o not in obs:
            continue
        prec_cl_over_o = cl_prec / em_inst0[obs.index(o), prec_idx]
        ratios.append(cl_test / (inv[o] * prec_cl_over_o))
    return np.array(ratios)


# Sentinel returned when ClCrit is NaN for a sample -> skip the whole sample.
SKIP_SAMPLE = object()


def cl_critical_correction(params, river, inv, em_inst0, rc0, i, rng, em_data,
                           dists, delta2r, minfrac):
    """Apply the ClCritical correction for sample ``i``.

    Returns a tuple
    ``(river_column, em_noprec, em_inst0, xdirect, ems_noprec, distcode_noprec)``
    or ``SKIP_SAMPLE`` if ClCrit is NaN, or ``None`` if the instance failed
    (river forced to NaN / end-member pull failed) and should be retried.
    """
    obs = params.ObsList
    ems = params.EMList0
    nOL = len(obs)
    prec_idx = ems.index("prec")
    cl_idx = obs.index("Cl")

    rc = rc0.copy()
    cl_over_norm_river = rc0[cl_idx]
    cl_river_conc = inv["Cl"]

    if params.ClCriticalValuesGiven == 1:
        clcrit = river.model_variable["ClCr"][i]
    else:
        clcrit = cl_river_conc

    if np.isnan(clcrit):
        return SKIP_SAMPLE

    # (1) re-pull the precipitation chemistry until removing ClCrit no longer
    # forces a negative ion -- only for Samples mode with no negative min bound.
    minfrac = np.asarray(minfrac, dtype=float)
    precratios = quantify_prec_ratios(obs, em_inst0, inv, clcrit, prec_idx, cl_idx)
    if (np.sum(precratios > 1) > 0 and params.IterateOver == "Samples"
            and np.sum(minfrac < 0) == 0):
        counter = 0
        while True:
            counter += 1
            em_inst0, failcase = pull_end_member_ratios(
                params, em_data, dists, rc0, rng, delta2r=delta2r, sample_index=i)
            if failcase or np.any(np.isnan(em_inst0)):
                return None
            precratios = quantify_prec_ratios(obs, em_inst0, inv, clcrit, prec_idx, cl_idx)
            if np.sum(precratios > 1) == 0:
                break
            if counter > 100:
                return None  # could not realize ClCrit physically; retry instance

    # (2) fractional contribution of precipitation to the normalization ion
    cl_over_norm_prec = em_inst0[cl_idx, prec_idx]
    if cl_river_conc <= clcrit:
        f_norm_prec = cl_over_norm_river / cl_over_norm_prec
    else:
        f_norm_prec = (clcrit / cl_river_conc) * cl_over_norm_river / cl_over_norm_prec

    # (3) subtract precipitation's contribution from the river column
    for j in range(nOL):
        rc[j] = rc0[j] - em_inst0[j, prec_idx] * f_norm_prec
    rc[np.abs(rc) < 5e-16] = 0.0

    # (4) drop precipitation from the inverted matrices and fix it in Xdirect
    keep = [k for k in range(len(ems)) if k != prec_idx]
    em_noprec = em_inst0[:, keep]
    distcode_noprec = [[row[k] for k in keep] for row in dists.distcode]
    ems_noprec = [ems[k] for k in keep]
    xdirect = np.full(len(ems), np.nan)
    xdirect[prec_idx] = f_norm_prec

    return rc, em_noprec, em_inst0, xdirect, ems_noprec, distcode_noprec
