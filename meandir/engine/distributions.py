"""End-member chemical distributions — port of ``MEANDIR_makeEMdistributions.m``.

For every (observation, end-member) pair this builds the probability
distribution that the Monte Carlo sampler draws from, plus the bookkeeping
arrays that flag which entries are defined relative to a sample
(``SMP-*``) or as a fractionation (``EPS-*``), and the sign array for
log-uniform distributions.

Distribution primitives expose ``.rvs(rng)`` returning a single draw. Truncated
normals use ``scipy.stats.truncnorm`` (inverse-CDF, like MATLAB's
``truncate``+``random``); a zero-variance normal collapses to a point mass.
Exact draws differ from MATLAB (different RNG), but the distributions are
identical — the basis of the statistical-fidelity guarantee.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy import stats

from ..constants import ISOTOPE_VARIABLES
from .endmembers import EMData


# ---------------------------------------------------------------------------
# Distribution primitives
# ---------------------------------------------------------------------------
class Dist:
    def rvs(self, rng) -> float:  # pragma: no cover - interface
        raise NotImplementedError


@dataclass
class PointMass(Dist):
    mu: float

    def rvs(self, rng) -> float:
        return self.mu


@dataclass
class Uniform(Dist):
    low: float
    high: float

    def rvs(self, rng) -> float:
        return self.low + rng.random() * (self.high - self.low)


@dataclass
class TruncNormal(Dist):
    mu: float
    sigma: float
    low: float = -math.inf
    high: float = math.inf

    def rvs(self, rng) -> float:
        if self.sigma == 0:
            return self.mu
        a = (self.low - self.mu) / self.sigma
        b = (self.high - self.mu) / self.sigma
        return float(stats.truncnorm.rvs(a, b, loc=self.mu, scale=self.sigma,
                                         random_state=rng))


@dataclass
class EMDistributions:
    """Holds the per-(obs, em) distribution machinery from makeEMdistributions."""

    nOL: int
    nEM: int
    distcode: list = field(default_factory=list)     # [ii][jj] -> str (distEMdist0)
    dist: list = field(default_factory=list)         # [ii][jj] -> Dist (TrEMdist)
    distmin: list = field(default_factory=list)      # [ii][jj] -> str|None
    distmax: list = field(default_factory=list)      # [ii][jj] -> str|None
    smp_mean: np.ndarray = None
    smp_std: np.ndarray = None
    smp_min: np.ndarray = None
    smp_max: np.ndarray = None
    lgu_sign: np.ndarray = None


# Per-mil delta variables get a lower truncation bound of -1000 (a delta cannot
# be below -1000 permil); converted ratios and intrinsic ratios are bounded at 0.
_DELTA_PERMIL = {"d7Li", "d13C", "d18O", "d26Mg", "d30Si", "d34S",
                 "d42Ca", "d44Ca", "d56Fe", "d98Mo"}
_INTRINSIC_RATIO = {"Fmod", "Sr8786", "Os8788"}


def make_em_distributions(params, em: EMData, delta2r) -> EMDistributions:
    """Build :class:`EMDistributions` for the scenario (makeEMdistributions)."""
    obs_list = params.ObsList
    em_list = params.EMList0
    nOL, nEM = len(obs_list), len(em_list)
    norm_type = params.NormalizationType
    convert = set(params.ConvertDelta2RList)
    obs_in_norm = set(params.ObsInNormalization)
    neg_ems = set(params.EndMembersWithNegativeRatios)
    couple_so4 = set(params.CoupleFeS2SO4intoEM)
    indx_iso = {ii for ii, o in enumerate(obs_list) if o in ISOTOPE_VARIABLES}

    out = EMDistributions(nOL=nOL, nEM=nEM)
    out.distcode = [[None] * nEM for _ in range(nOL)]
    out.dist = [[None] * nEM for _ in range(nOL)]
    out.distmin = [[None] * nEM for _ in range(nOL)]
    out.distmax = [[None] * nEM for _ in range(nOL)]
    out.smp_mean = np.full((nOL, nEM), np.nan)
    out.smp_std = np.full((nOL, nEM), np.nan)
    out.smp_min = np.full((nOL, nEM), np.nan)
    out.smp_max = np.full((nOL, nEM), np.nan)
    out.lgu_sign = np.full((nOL, nEM), np.nan)

    for ii in range(nOL):
        obs = obs_list[ii]
        dist_type = em.disttype[obs]
        do_convert = obs in convert
        conv = delta2r.factor(obs) if do_convert else None
        for jj in range(nEM):
            emc = em_list[jj]
            pd = None
            stdentry = None
            out.distcode[ii][jj] = dist_type

            if dist_type == "NOR":
                pd, stdentry = _build_normal(out, ii, jj, em, emc, obs,
                                             norm_type, do_convert, conv,
                                             params.IterateOver)
            elif dist_type == "UNI":
                pd, stdentry = _build_uniform(out, ii, jj, em, emc, obs,
                                              norm_type, do_convert, conv)
            elif dist_type == "LGU":
                pd, stdentry = _build_loguniform(out, ii, jj, em, emc, obs,
                                                 norm_type, do_convert, conv)

            # Truncate normal distributions (only for NOR-declared numerators).
            if dist_type == "NOR" and stdentry not in (0, None):
                low = _lower_bound(obs, ii, indx_iso, do_convert, emc, neg_ems)
                high = 1.0 if (norm_type == "SumObs" and obs in obs_in_norm) else math.inf
                td = _truncate(pd, low, high)
            else:
                td = pd

            out.dist[ii][jj] = td

            # FeS2-coupled SO4: override the SO4 distribution with the pyrite one.
            if obs == "SO4" and emc in couple_so4:
                out.dist[ii][jj] = _pyrite_so4_dist(em, norm_type)

    return out


# ---------------------------------------------------------------------------
# distribution constructors per type
# ---------------------------------------------------------------------------
def _build_normal(out, ii, jj, em, emc, obs, norm_type, do_convert, conv, iterate_over):
    mean = em.get(emc, obs, "Men")
    std = em.get(emc, obs, "Sig")
    SAMPLE_MEN = 123456789j
    FRAC_MEN = 314159265j
    if mean != SAMPLE_MEN and mean != FRAC_MEN:
        if do_convert:
            rmean = (mean.real / 1000 + 1) * conv
            rstd = (std.real / 1000) * conv
        else:
            rmean, rstd = mean.real, std.real
        return TruncNormal(rmean, rstd), rstd
    if mean == SAMPLE_MEN and iterate_over == "Samples":
        out.distcode[ii][jj] = "SMP-NOR"
        pull_mean = em.get_addvalue(emc, obs, "Men")
        pull_std = std.real
        if do_convert:
            out.smp_mean[ii, jj] = (pull_mean / 1000) * conv
            out.smp_std[ii, jj] = (pull_std / 1000) * conv
        else:
            out.smp_mean[ii, jj] = pull_mean
            out.smp_std[ii, jj] = pull_std
        return PointMass(0.0), 0.0
    if mean == FRAC_MEN and iterate_over == "Samples":
        pull_mean = em.get_addvalue(emc, obs, "Men")
        pull_std = std.real
        if do_convert:
            rmean = (pull_mean / 1000) * conv
            rstd = (pull_std / 1000) * conv
        else:
            rmean, rstd = pull_mean, pull_std
        out.distcode[ii][jj] = "EPS-NOR"
        return TruncNormal(rmean, rstd), rstd
    return PointMass(0.0), 0.0


def _build_uniform(out, ii, jj, em, emc, obs, norm_type, do_convert, conv):
    SMIN, SMAX, FMIN, FMAX = -123456789.0, 123456789.0, -314159265.0, 314159265.0
    minentry = em.get(emc, obs, "Min")
    maxentry = em.get(emc, obs, "Max")
    minentry = minentry.real if isinstance(minentry, complex) else minentry
    maxentry = maxentry.real if isinstance(maxentry, complex) else maxentry

    if minentry not in (SMIN, FMIN) and maxentry not in (SMAX, FMAX):
        if minentry == maxentry:
            if do_convert:
                rmin = (minentry / 1000 + 1) * conv
                rmax = (maxentry / 1000 + 1) * conv
            else:
                rmin, rmax = minentry, maxentry
            out.distcode[ii][jj] = "NOR"
            return PointMass((rmin + rmax) / 2), 0.0
        if minentry < maxentry:
            if do_convert:
                rmin = (minentry / 1000 + 1) * conv
                rmax = (maxentry / 1000 + 1) * conv
            else:
                rmin, rmax = minentry, maxentry
            return Uniform(rmin, rmax), None
        return None, None  # min > max: invalid (MATLAB warns)

    # sample-relative / fractionation cases
    if minentry == SMIN and maxentry != SMAX:
        out.distmin[ii][jj], out.distmax[ii][jj] = "SMP", "UNI"
        pull_min = em.get_addvalue(emc, obs, "Min")
        pull_max = maxentry
        if do_convert:
            out.smp_min[ii, jj] = (pull_min / 1000) * conv
            out.smp_max[ii, jj] = (pull_max / 1000 + 1) * conv
        else:
            out.smp_min[ii, jj], out.smp_max[ii, jj] = pull_min, pull_max
    elif minentry != SMIN and maxentry == SMAX:
        out.distmin[ii][jj], out.distmax[ii][jj] = "UNI", "SMP"
        pull_min = minentry
        pull_max = em.get_addvalue(emc, obs, "Max")
        if do_convert:
            out.smp_min[ii, jj] = (pull_min / 1000 + 1) * conv
            out.smp_max[ii, jj] = (pull_max / 1000) * conv
        else:
            out.smp_min[ii, jj], out.smp_max[ii, jj] = pull_min, pull_max
    elif minentry == SMIN and maxentry == SMAX:
        out.distmin[ii][jj], out.distmax[ii][jj] = "SMP", "SMP"
        pull_min = em.get_addvalue(emc, obs, "Min")
        pull_max = em.get_addvalue(emc, obs, "Max")
        if do_convert:
            out.smp_min[ii, jj] = (pull_min / 1000) * conv
            out.smp_max[ii, jj] = (pull_max / 1000) * conv
        else:
            out.smp_min[ii, jj], out.smp_max[ii, jj] = pull_min, pull_max
    elif minentry == FMIN and maxentry == FMAX:
        pull_min = em.get_addvalue(emc, obs, "Min")
        pull_max = em.get_addvalue(emc, obs, "Max")
        if do_convert:
            rmin = (pull_min / 1000) * conv
            rmax = (pull_max / 1000) * conv
        else:
            rmin, rmax = pull_min, pull_max
        if pull_min == pull_max:
            out.distcode[ii][jj] = "EPS-NOR"
            return PointMass((rmin + rmax) / 2), 0.0
        if pull_min < pull_max:
            out.distcode[ii][jj] = "EPS-UNI"
            return Uniform(rmin, rmax), None
        return None, None

    # if either endpoint is sample-floating, set placeholder
    if minentry == SMIN or maxentry == SMAX:
        out.distcode[ii][jj] = "SMP-UNI"
        return PointMass(0.0), None
    return None, None


def _build_loguniform(out, ii, jj, em, emc, obs, norm_type, do_convert, conv):
    SMIN, SMAX, FMIN, FMAX = -123456789.0, 123456789.0, -314159265.0, 314159265.0
    minentry = em.get(emc, obs, "Min")
    maxentry = em.get(emc, obs, "Max")
    minentry = minentry.real if isinstance(minentry, complex) else minentry
    maxentry = maxentry.real if isinstance(maxentry, complex) else maxentry

    if minentry not in (SMIN, FMIN) and maxentry not in (SMAX, FMAX):
        if minentry == maxentry:
            if do_convert:
                rmin = (minentry / 1000 + 1) * conv
                rmax = (maxentry / 1000 + 1) * conv
            else:
                rmin, rmax = minentry, maxentry
            out.distcode[ii][jj] = "NOR"
            return PointMass((rmin + rmax) / 2), 0.0
        if do_convert:
            rmin = (minentry / 1000 + 1) * conv
            rmax = (maxentry / 1000 + 1) * conv
        else:
            rmin, rmax = minentry, maxentry
        if rmin == 0 or rmax == 0:
            return None, None
        if rmin < rmax:
            isneg = rmin < 0 or rmax < 0
            out.lgu_sign[ii, jj] = -1 if isneg else 1
            if isneg:
                v1, v2 = math.log(-rmin), math.log(-rmax)
            else:
                v1, v2 = math.log(rmin), math.log(rmax)
            low, high = (v1, v2) if v1 < v2 else (v2, v1)
            return Uniform(low, high), None
        return None, None

    if minentry == SMIN and maxentry != SMAX:
        out.distmin[ii][jj], out.distmax[ii][jj] = "SMP", "LGU"
        pull_min = em.get_addvalue(emc, obs, "Min")
        pull_max = maxentry
        if do_convert:
            out.smp_min[ii, jj] = (pull_min / 1000) * conv
            out.smp_max[ii, jj] = (pull_max / 1000 + 1) * conv
        else:
            out.smp_min[ii, jj], out.smp_max[ii, jj] = pull_min, pull_max
    elif minentry != SMIN and maxentry == SMAX:
        out.distmin[ii][jj], out.distmax[ii][jj] = "LGU", "SMP"
        pull_min = minentry
        pull_max = em.get_addvalue(emc, obs, "Max")
        if do_convert:
            out.smp_min[ii, jj] = (pull_min / 1000 + 1) * conv
            out.smp_max[ii, jj] = (pull_max / 1000) * conv
        else:
            out.smp_min[ii, jj], out.smp_max[ii, jj] = pull_min, pull_max
    elif minentry == SMIN and maxentry == SMAX:
        out.distmin[ii][jj], out.distmax[ii][jj] = "SMP", "SMP"
        pull_min = em.get_addvalue(emc, obs, "Min")
        pull_max = em.get_addvalue(emc, obs, "Max")
        if do_convert:
            out.smp_min[ii, jj] = (pull_min / 1000) * conv
            out.smp_max[ii, jj] = (pull_max / 1000) * conv
        else:
            out.smp_min[ii, jj], out.smp_max[ii, jj] = pull_min, pull_max
    elif minentry == FMIN and maxentry == FMAX:
        pull_min = em.get_addvalue(emc, obs, "Min")
        pull_max = em.get_addvalue(emc, obs, "Max")
        if do_convert:
            rmin = (pull_min / 1000) * conv
            rmax = (pull_max / 1000) * conv
        else:
            rmin, rmax = pull_min, pull_max
        if pull_min == pull_max:
            out.distcode[ii][jj] = "EPS-NOR"
            return PointMass((rmin + rmax) / 2), 0.0
        if rmin == 0 or rmax == 0:
            return None, None
        if rmin < rmax:
            isneg = rmin < 0 or rmax < 0
            out.lgu_sign[ii, jj] = -1 if isneg else 1
            if isneg:
                v1, v2 = math.log(-rmin), math.log(-rmax)
            else:
                v1, v2 = math.log(rmin), math.log(rmax)
            low, high = (v1, v2) if v1 < v2 else (v2, v1)
            out.distcode[ii][jj] = "EPS-LGU"
            return Uniform(low, high), None
        return None, None

    if minentry == SMIN or maxentry == SMAX:
        out.distcode[ii][jj] = "SMP-LGU"
        return PointMass(0.0), None
    return None, None


def _lower_bound(obs, ii, indx_iso, do_convert, emc, neg_ems):
    if ii in indx_iso:
        if obs in _INTRINSIC_RATIO:
            return 0.0
        if obs in _DELTA_PERMIL:
            return 0.0 if do_convert else -1000.0
        return 0.0
    if emc in neg_ems:
        return -math.inf
    return 0.0


def _truncate(pd, low, high):
    if isinstance(pd, PointMass):
        return pd
    if isinstance(pd, TruncNormal):
        return TruncNormal(pd.mu, pd.sigma, low, high)
    return pd


def _pyrite_so4_dist(em: EMData, norm_type):
    dt = em.disttype.get("FeS2SO4")
    if dt == "NOR":
        mu = em.get("pyri", "FeS2SO4", "Men")
        sg = em.get("pyri", "FeS2SO4", "Sig")
        return TruncNormal(float(mu.real), float(sg.real))
    if dt == "UNI":
        lo = float(em.get("pyri", "FeS2SO4", "Min").real)
        hi = float(em.get("pyri", "FeS2SO4", "Max").real)
        if lo == hi:
            return PointMass((lo + hi) / 2)
        return Uniform(lo, hi)
    return PointMass(0.0)
