"""Draw an internally-consistent end-member matrix — port of
``MEANDIR_PullEndMemberRatios.m``.

For each (observation, end-member) a value is drawn from the distribution built
by :mod:`meandir.engine.distributions`, then the matrix is made internally
consistent: mass-balance closure to the normalization, optional evaporite
stoichiometry, optional charge-balance closure, sign checks, and (when
normalizing to a sum) a reject/retry loop up to ``maxtrycount`` attempts.
Finally isotope rows are converted to the product of isotope ratio and carrier
concentration (the quantity that mixes ~linearly).
"""

from __future__ import annotations

import math

import numpy as np

from ..constants import ISOTOPE_VARIABLES, isotope_ion

MAXTRYCOUNT = 500


def pull_end_member_ratios(params, em_data, dists, river_column0, rng,
                           delta2r=None, sample_index=None, return_trycount=False):
    """Return (EMIterInst [nOL x nEM], fail_case). If ``return_trycount`` is set,
    also return the number of draw attempts (for sampling diagnostics)."""
    obs = params.ObsList
    ems = params.EMList0
    nIR, nEM = dists.nOL, dists.nEM
    charges = params.IonCharges
    norm_type = params.NormalizationType
    obs_in_norm = params.ObsInNormalization
    cats_norm = params.CationsInNormalization
    ans_norm = params.AnionsInNormalization
    neu_norm = params.NeutralInNormalization
    carbonmatch = params.carbonisotopematch
    neg_ems = params.EndMembersWithNegativeRatios

    iso_pos = np.array([o in ISOTOPE_VARIABLES for o in obs])
    indx_iso = np.where(iso_pos)[0]
    river_column0 = np.asarray(river_column0, dtype=float).reshape(-1)
    if river_column0.size != nIR:  # End-members mode passes scalar 0
        river_column0 = np.zeros(nIR)

    EMIterInst = np.full((nIR, nEM), np.nan)
    fail = 0
    trycount = 0
    accepted = False

    while not accepted:
        trycount += 1

        # (1) draw a value for each (observation, end-member)
        for ii in range(nIR):
            for jj in range(nEM):
                code = dists.distcode[ii][jj]
                if code in ("UNI", "NOR", "EPS-UNI", "EPS-NOR"):
                    EMIterInst[ii, jj] = dists.dist[ii][jj].rvs(rng)
                elif code in ("LGU", "EPS-LGU"):
                    EMIterInst[ii, jj] = math.exp(dists.dist[ii][jj].rvs(rng)) * dists.lgu_sign[ii, jj]
                elif code in ("SMP-NOR", "SMP-UNI", "SMP-LGU"):
                    EMIterInst[ii, jj] = _sample_relative(
                        code, ii, jj, obs, river_column0, dists, carbonmatch, rng)

        # (2) evaporite stoichiometry (pre-closure)
        if params.BalanceEvaporite == 1 and "evap" in ems:
            _balance_evaporite_pre(EMIterInst, obs, ems, charges)

        # (3) charge-balance closure when all ions explicitly resolved
        if params.AllIonsExplicitlyResolved == 1:
            _charge_closure(EMIterInst, obs, ems, charges, obs_in_norm,
                            cats_norm, ans_norm, neu_norm, params.ListChargeClosure)

        # (4) single-observation normalization is already consistent
        if norm_type != "SumObs":
            accepted = True

        # (5) sum normalization: enforce internal consistency
        if norm_type == "SumObs":
            # (5.1) one ratio per end-member by mass balance
            for i, emc in enumerate(ems):
                offion = obs.index(params.ListNormClosure[i])
                onions = [obs.index(o) for o in obs_in_norm if obs.index(o) != offion]
                EMIterInst[offion, i] = 1 - np.sum(EMIterInst[onions, i])

            # (5.2) evaporite closure
            if (params.BalanceEvaporite == 1 and "evap" in ems
                    and "SO4" in obs and "Cl" in obs):
                _balance_evaporite_post(EMIterInst, obs, ems, charges,
                                        obs_in_norm, params.ListNormClosure)

            # (5.3) charge closure
            if params.AllIonsExplicitlyResolved == 1:
                _charge_closure(EMIterInst, obs, ems, charges, obs_in_norm,
                                cats_norm, ans_norm, neu_norm, params.ListChargeClosure)

            # (5.4) clamp tiny negatives to zero
            mask = (EMIterInst > -1e-12) & (EMIterInst < 0)
            EMIterInst[mask] = 0

            # (5.5.x) consistency checks
            if not _check_mass_balance(EMIterInst, obs, ems, obs_in_norm):
                fail = 1
            if params.AllIonsExplicitlyResolved == 1:
                if not _check_charge_balance(EMIterInst, charges):
                    fail = 1

            # (5.6) accept if nothing that must be positive is negative
            nonnegrows = ~iso_pos
            negem = np.array([e in neg_ems for e in ems])
            sub = EMIterInst[np.ix_(nonnegrows, ~negem)]
            if np.sum(sub < 0) == 0:
                accepted = True

        # (6) give up after maxtrycount attempts
        if trycount == MAXTRYCOUNT:
            accepted = True
            EMIterInst = np.full((nIR, nEM), np.nan)

    # (7) couple FeS2 d34S across end-members
    if ("SO4" in obs and "d34S" in obs and params.CoupleFeS2d34SintoEM):
        _couple_fes2_d34s(EMIterInst, obs, ems, em_data, params, delta2r, rng)

    # (8) isotope rows -> isotope ratio * carrier concentration
    for ii in indx_iso:
        ion = isotope_ion(obs[ii], carbonmatch)
        ionpos = obs.index(ion)
        EMIterInst[ii, :] = EMIterInst[ii, :] * EMIterInst[ionpos, :]

    if return_trycount:
        return EMIterInst, fail, trycount
    return EMIterInst, fail


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _river_ratio(ii, obs, river_column0, carbonmatch):
    """Isolate the river chemical parameter (ratio for isotopes)."""
    name = obs[ii]
    if name not in ISOTOPE_VARIABLES:
        return river_column0[ii]
    ion = isotope_ion(name, carbonmatch)
    return river_column0[ii] / river_column0[obs.index(ion)]


def _sample_relative(code, ii, jj, obs, river_column0, dists, carbonmatch, rng):
    riverentry = _river_ratio(ii, obs, river_column0, carbonmatch)
    if code == "SMP-NOR":
        return riverentry + dists.smp_mean[ii, jj] + rng.standard_normal() * dists.smp_std[ii, jj]
    if code == "SMP-UNI":
        minval = (riverentry + dists.smp_min[ii, jj]) if dists.distmin[ii][jj] == "SMP" else dists.smp_min[ii, jj]
        maxval = (riverentry + dists.smp_max[ii, jj]) if dists.distmax[ii][jj] == "SMP" else dists.smp_max[ii, jj]
        if minval <= maxval:
            return minval + rng.random() * (maxval - minval)
        return np.nan
    if code == "SMP-LGU":
        minval = (riverentry + dists.smp_min[ii, jj]) if dists.distmin[ii][jj] == "SMP" else dists.smp_min[ii, jj]
        maxval = (riverentry + dists.smp_max[ii, jj]) if dists.distmax[ii][jj] == "SMP" else dists.smp_max[ii, jj]
        if minval <= maxval:
            minval = 1e-10 if minval <= 0 else minval
            maxval = 1e-10 if maxval <= 0 else maxval
            v1, v2 = math.log(minval), math.log(maxval)
            low, high = (v1, v2) if v1 < v2 else (v2, v1)
            return math.exp(low + rng.random() * (high - low))
        return np.nan
    return np.nan


def _cation_mask(obs, charges):
    return np.array([c == "+" for c in charges])


def _balance_evaporite_pre(EMIterInst, obs, ems, charges):
    ev = ems.index("evap")
    camgsr = np.array([o in ("Ca", "Mg", "Sr") for o in obs])
    if "SO4" in obs and camgsr.any():
        EMIterInst[obs.index("SO4"), ev] = np.sum(EMIterInst[camgsr, ev])
    if "Cl" in obs:
        catpos = _cation_mask(obs, charges) & ~camgsr
        EMIterInst[obs.index("Cl"), ev] = np.sum(EMIterInst[catpos, ev])


def _balance_evaporite_post(EMIterInst, obs, ems, charges, obs_in_norm, list_norm_closure):
    ev = ems.index("evap")
    closureion = list_norm_closure[ev]
    camgsr_names = ("Ca", "Mg", "Sr")
    if closureion in camgsr_names:
        if "SO4" not in obs_in_norm:
            camgsr = np.array([o in camgsr_names for o in obs])
            EMIterInst[obs.index("SO4"), ev] = np.sum(EMIterInst[camgsr, ev])
        else:
            others = [o for o in obs_in_norm if o not in ("SO4", closureion)]
            massremaining = 1 - np.sum([EMIterInst[obs.index(o), ev] for o in others])
            chargeaccounted = 0.0
            for cat in camgsr_names:
                if cat != closureion and cat in obs:
                    chargeaccounted += EMIterInst[obs.index(cat), ev]
            closure = 0.5 * (massremaining - chargeaccounted)
            EMIterInst[obs.index(closureion), ev] = closure
            EMIterInst[obs.index("SO4"), ev] = chargeaccounted + closure
    else:
        catpos = _cation_mask(obs, charges) & ~np.array([o in camgsr_names for o in obs])
        if "Cl" not in obs_in_norm:
            EMIterInst[obs.index("Cl"), ev] = np.sum(EMIterInst[catpos, ev])
        else:
            others = [o for o in obs_in_norm if o not in ("Cl", closureion)]
            massremaining = 1 - np.sum([EMIterInst[obs.index(o), ev] for o in others])
            catpos2 = _cation_mask(obs, charges) & ~np.array(
                [o in camgsr_names or o == closureion for o in obs])
            chargeaccounted = np.sum(EMIterInst[catpos2, ev])
            closure = 0.5 * (massremaining - chargeaccounted)
            EMIterInst[obs.index(closureion), ev] = closure
            EMIterInst[obs.index("Cl"), ev] = chargeaccounted + closure


def _charge_closure(EMIterInst, obs, ems, charges, obs_in_norm,
                    cats_norm, ans_norm, neu_norm, list_charge_closure):
    for i in range(len(ems)):
        ratio = np.sum([EMIterInst[obs.index(o), i] for o in obs_in_norm])
        offpos = obs.index(list_charge_closure[i])
        for j in range(len(obs)):
            if j == offpos:
                continue
            cj, oj = charges[j], obs[j]
            if cj == "+":
                if oj not in cats_norm:
                    ratio = ratio + EMIterInst[j, i]
            elif cj == "-":
                if oj in ans_norm:
                    ratio = ratio - 2 * EMIterInst[j, i]
                else:
                    ratio = ratio - EMIterInst[j, i]
            elif cj == "0":
                if oj in neu_norm:
                    ratio = ratio - EMIterInst[j, i]
        sign = 1 if charges[offpos] == "-" else -1
        EMIterInst[offpos, i] = ratio * sign


def _check_mass_balance(EMIterInst, obs, ems, obs_in_norm):
    rows = [obs.index(o) for o in obs_in_norm]
    sums = np.sum(EMIterInst[rows, :], axis=0)
    return np.all(np.isclose(sums, 1.0))


def _check_charge_balance(EMIterInst, charges):
    catpos = np.array([c == "+" for c in charges])
    anipos = np.array([c == "-" for c in charges])
    cb = np.sum(EMIterInst[catpos, :], axis=0) - np.sum(EMIterInst[anipos, :], axis=0)
    return np.all(np.isclose(cb, 0.0))


def _couple_fes2_d34s(EMIterInst, obs, ems, em_data, params, delta2r, rng):
    """Overwrite the d34S of the CoupleFeS2d34SintoEM end-members with the
    pyrite FeS2 d34S distribution (MEANDIR_PullEndMemberRatios step 7)."""
    pos_d34s = obs.index("d34S")
    targets = [k for k, e in enumerate(ems) if e in params.CoupleFeS2d34SintoEM]
    dt = em_data.disttype.get("FeS2d34S")
    do_convert = "d34S" in params.ConvertDelta2RList
    conv = delta2r.factor("d34S") if do_convert else None
    if dt == "NOR":
        mean = float(em_data.get("pyri", "FeS2d34S", "Men").real)
        std = float(em_data.get("pyri", "FeS2d34S", "Sig").real)
        if do_convert:
            rmean = (mean / 1000 + 1) * conv
            rstd = (std / 1000) * conv
        else:
            rmean, rstd = mean, std
        for k in targets:
            EMIterInst[pos_d34s, k] = rmean + rstd * rng.standard_normal()
    elif dt == "UNI":
        lo = float(em_data.get("pyri", "FeS2d34S", "Min").real)
        hi = float(em_data.get("pyri", "FeS2d34S", "Max").real)
        if do_convert:
            lo = (lo / 1000 + 1) * conv
            hi = (hi / 1000 + 1) * conv
        for k in targets:
            EMIterInst[pos_d34s, k] = lo + rng.random() * (hi - lo)
