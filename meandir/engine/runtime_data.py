"""Per-simulation river-data handling.

Ports of ``MEANDIR_AdjustRiverDataForAnalyticalError.m`` (draw the sample's
observations within analytical error and rebuild the normalized river column)
and ``MEANDIR_PrepareUpdatedDataForInversion.m`` (drop zero-valued, relatively-
weighted entries and the end-members that source them).
"""

from __future__ import annotations

import numpy as np

from ..constants import ISOTOPE_VARIABLES, isotope_ion

# Isotopes whose delta values get converted to ratios in AdjustRiverData (the
# intrinsic-ratio variables Sr8786/Os8788/Fmod are never converted here).
_DELTA_CONVERTIBLE = {"d7Li", "d13C", "d18O", "d26Mg", "d30Si", "d34S",
                      "d42Ca", "d44Ca", "d56Fe", "d98Mo"}


def adjust_river_data(river, params, i, rng, delta2r):
    """Return (inv_riv_dat dict, river_column0 array) for sample ``i``."""
    obs = params.ObsList
    mv = river.model_variable
    inv = {}

    # observations needed: ObsList plus SO4 (always carried forward)
    needed = list(obs)
    if "SO4" not in needed:
        needed.append("SO4")

    if params.AdjustRiverObs == 1:
        for o in needed:
            inv[o] = mv[o][i] + mv[f"{o}_asd"][i] * rng.standard_normal()
        norm = 0.0
        for o in params.ObsInNormalization:
            norm += inv[o]
    else:
        for o in needed:
            inv[o] = mv[o][i]
        norm = mv["norm"][i]
    inv["norm"] = norm

    # delta -> ratio conversion for the convertible isotopes
    for iso in params.ConvertDelta2RList:
        if iso in _DELTA_CONVERTIBLE and iso in inv:
            conv = delta2r.factor(iso)
            inv[iso] = (inv[iso] / 1000 + 1) * conv

    # rebuild normalized river column in ObsList order
    rc = np.full(len(obs), np.nan)
    for k, o in enumerate(obs):
        if o not in ISOTOPE_VARIABLES:
            rc[k] = inv[o] / norm
        else:
            ion = isotope_ion(o, params.carbonisotopematch)
            rc[k] = inv[o] * (inv[ion] / norm)
    return inv, rc


def prepare_updated_data(Xdirect, river_column, obs, ems, ems0, relpos,
                         em_inst, solver, distcode):
    """Reduce the inversion for zero, relatively-weighted river entries."""
    if solver in ("mldivide", "lsqnonneg"):
        return (river_column.copy(), em_inst.copy(), list(ems), list(obs),
                Xdirect, [row[:] for row in distcode])

    # optimize-family solvers: drop zero river entries with relative weighting
    zero_remove = (river_column == 0) & relpos
    keep_rows = ~zero_remove

    em_remove = np.zeros(len(ems), dtype=bool)
    for pos in np.where(zero_remove)[0]:
        nonzero_cols = np.where(em_inst[pos, :] != 0)[0]
        for col in nonzero_cols:
            Xdirect[ems0.index(ems[col])] = 0.0
            em_remove[col] = True
    keep_cols = ~em_remove

    rc_r = river_column[keep_rows]
    em_r = em_inst[np.ix_(keep_rows, keep_cols)]
    ems_r = [e for e, k in zip(ems, keep_cols) if k]
    obs_r = [o for o, k in zip(obs, keep_rows) if k]
    distcode_r = [[distcode[ri][ci] for ci in range(len(ems)) if keep_cols[ci]]
                  for ri in range(len(obs)) if keep_rows[ri]]
    return rc_r, em_r, ems_r, obs_r, Xdirect, distcode_r
