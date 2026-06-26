"""The per-simulation inversion — ports of ``MEANDIR_InvertActiveSimulation.m``
and ``MEANDIR_CostFunction.m``.

Solver mapping (statistical-fidelity port):
  * ``mldivide``  -> ``numpy.linalg.lstsq``
  * ``lsqnonneg`` -> ``scipy.optimize.nnls``  (same Lawson-Hanson algorithm)
  * ``optimize``  family -> ``scipy.optimize.minimize(method="SLSQP")`` with
    bound constraints, in place of MATLAB's ``fmincon``.
"""

from __future__ import annotations

import math

import numpy as np
import scipy.linalg as sla
from scipy.optimize import nnls, minimize, Bounds, lsq_linear

# Bounded optimizer for the *_optimize solvers, in place of MATLAB fmincon.
#   "SLSQP"        -> active-set SQP (fast; the default here)
#   "trust-constr" -> interior-point-style method, closer to fmincon's default
# Toggleable so we can A/B which better matches the published results.
OPTIMIZER_METHOD = "SLSQP"

# Experimental solver modes for the convex (no-fractionation) case:
#   "default"    -> the OPTIMIZER_METHOD nonlinear optimizer
#   "lsq_linear" -> exact bounded least-squares (true convex minimum, fast)
#   "x0clip"     -> X0 clipped to bounds, no optimization (maximal under-convergence)
# Only affects instances without fractionation; others fall back to "default".
SOLVER_MODE = "default"


def _bvls(A, b, solvecf_r, relpos_r, abspos_r, w_r, lo, hi):
    """Exact bound-constrained least squares for the MEANDIR relative/absolute
    cost (no fractionation). Builds the row-scaled design matrix and calls
    scipy.optimize.lsq_linear (BVLS)."""
    rel = solvecf_r & relpos_r
    rows = [(np.sqrt(w_r[rel]) / b[rel])[:, None] * A[rel, :]]
    d = [np.sqrt(w_r[rel])]
    ab = solvecf_r & abspos_r
    if ab.any():
        rows.append(np.sqrt(w_r[ab])[:, None] * A[ab, :])
        d.append(np.sqrt(w_r[ab]) * b[ab])
    C = np.vstack(rows)
    dvec = np.concatenate(d)
    res = lsq_linear(C, dvec, bounds=(lo, hi), method="bvls")
    return res.x

# Initial-condition strategy for the mldivide-family solvers:
#   "minnorm"-> numpy.linalg.lstsq (minimum-2-norm solution)
#   "basic"  -> literal MATLAB backslash semantics (column-pivoted QR; for an
#               underdetermined system a basic solution with <=rank nonzeros)
#
# These differ only for the underdetermined scenarios (more end-members than
# observations, e.g. scenario 2 at 9x10). Empirically (scripts/ab_initial_condition.py)
# "minnorm" reproduces the published Table S2 contributions far better: the
# literal "basic" solution zeros a structural end-member (e.g. carbonate) at the
# start and drives the optimizer to a wrong minimum (worst Table S2 deviation
# 3.4 -> 42 points on scenario 2). MATLAB's effective behaviour matches the
# minimum-norm solution here, so "minnorm" is the default; "basic" is retained
# for reproducibility of that finding.
MLDIVIDE_INITIAL_CONDITION = "minnorm"


def _mldivide(A, b):
    """Reproduce MATLAB ``A\\b`` for the inversion's possibly-underdetermined
    linear system. Square/overdetermined full-rank -> ordinary least squares;
    underdetermined or rank-deficient -> a basic solution via column-pivoted QR
    (at most rank(A) nonzero entries), as MATLAB's backslash returns."""
    m, n = A.shape
    rank = np.linalg.matrix_rank(A)
    if rank == n and m >= n:
        x, *_ = np.linalg.lstsq(A, b, rcond=None)
        return x
    Q, R, piv = sla.qr(A, mode="economic", pivoting=True)
    x = np.zeros(n)
    x[piv[:rank]] = sla.solve_triangular(R[:rank, :rank], Q[:, :rank].T @ b)
    return x



def cost_function(X0, em_inst0, em_inst_r_base, river_col0, river_col_r,
                  solvecf_r, abspos_r, relpos_r, weighting_r, frac, xdirect,
                  sources):
    """Optimization objective (MEANDIR_CostFunction)."""
    em_inst_r = em_inst_r_base.copy()
    nanmask = np.isnan(xdirect)
    # End-members that do not source the carrier ion give 0/0 entries in em_iso;
    # they are excluded from `src` so the NaN/Inf never propagate (MATLAB
    # produces the same values silently).
    with np.errstate(divide="ignore", invalid="ignore"):
        for ii in range(frac.red.n):
            Xfull = xdirect.copy()
            Xfull[nanmask] = X0
            ionpos0 = frac.red.ionpos0[ii]
            isopos0 = frac.red.isopos0[ii]
            empos0 = frac.red.empos0[ii]
            Xactive = em_inst0[ionpos0, :] * Xfull / river_col0[ionpos0]
            em_iso = em_inst0[isopos0, :] / em_inst0[ionpos0, :]
            src = sources & (Xactive != 0)
            denom = np.sum(Xactive[src])
            source_iso = np.sum(Xactive[src] * em_iso[src]) / denom
            active_frac = em_inst0[isopos0, empos0] / em_inst0[ionpos0, empos0]
            new_iso = (source_iso + active_frac) + active_frac * Xactive[empos0] / denom
            ir, jr = frac.red.isopos_r[ii], frac.red.empos_r[ii]
            kr = frac.red.ionpos_r[ii]
            em_inst_r[ir, jr] = em_inst_r[kr, jr] * new_iso

    model = em_inst_r @ X0
    cost_rel = ((river_col_r - model) / river_col_r) ** 2
    cost_abs = (river_col_r - model) ** 2
    mask_rel = solvecf_r & relpos_r
    mask_abs = solvecf_r & abspos_r
    total = (np.sum(weighting_r[mask_rel] * cost_rel[mask_rel])
             + np.sum(weighting_r[mask_abs] * cost_abs[mask_abs]))
    return math.sqrt(total)


def _initial_condition(solver, initial, river_col_r):
    if solver in ("mldivide", "mldivide_optimize"):
        if MLDIVIDE_INITIAL_CONDITION == "minnorm":
            X0, *_ = np.linalg.lstsq(initial, river_col_r, rcond=None)
            return X0
        return _mldivide(initial, river_col_r)
    if solver in ("lsqnonneg", "lsqnonneg_optimize"):
        X0, _ = nnls(initial, river_col_r)
        return X0
    # 'optimize': equal fractional contributions
    n = initial.shape[1]
    return np.ones(n) / n


def _bounds(minfrac_r, maxfrac_r):
    out = []
    for lo, hi in zip(minfrac_r, maxfrac_r):
        out.append((None if lo == -math.inf else lo,
                    None if hi == math.inf else hi))
    return out


def invert_active_simulation(solver, em_inst0, em_inst_r, river_col0,
                             river_col_r, ems0, xdirect, minfrac_r, maxfrac_r,
                             solvecf_r, abspos_r, relpos_r, weighting_r, frac,
                             sources):
    nEM = len(ems0)
    if np.any(np.isnan(river_col_r)) or em_inst_r.size == 0:
        return np.full(nEM, np.nan), em_inst0, np.nan

    # (2) initial condition; encode half the fractionation offset
    initial = em_inst_r.copy()
    for ii in range(frac.red.n):
        ir, jr, kr = frac.red.isopos_r[ii], frac.red.empos_r[ii], frac.red.ionpos_r[ii]
        ionvalue = initial[kr, jr]
        fracvalue = initial[ir, jr] / ionvalue
        riveriso = river_col_r[ir] / river_col_r[kr]
        initial[ir, jr] = ionvalue * (riveriso + 0.5 * fracvalue)

    X0 = _initial_condition(solver, initial, river_col_r)
    functioncost = math.sqrt(np.sum((river_col_r - initial @ X0) ** 2))
    if math.isnan(functioncost):
        functioncost = 999999.0

    optimize = solver in ("optimize", "mldivide_optimize", "lsqnonneg_optimize")
    if optimize and not np.any(np.isnan(X0)):
        args = (em_inst0, em_inst_r, river_col0, river_col_r, solvecf_r,
                abspos_r, relpos_r, weighting_r, frac, xdirect, sources)
        bounds = _bounds(minfrac_r, maxfrac_r)
        if SOLVER_MODE != "default" and frac.red.n == 0:
            # convex linear case: exact bounded LSQ or unoptimized X0-clip
            if SOLVER_MODE == "lsq_linear":
                Xtemp = _bvls(em_inst_r, river_col_r, solvecf_r, relpos_r,
                              abspos_r, weighting_r, minfrac_r, maxfrac_r)
            else:  # x0clip
                Xtemp = np.clip(X0, minfrac_r, maxfrac_r)
            functioncost = cost_function(Xtemp, *args)
            nanmask = np.isnan(xdirect)
            X = np.full(nEM, np.nan)
            if not np.all(nanmask):
                X[~nanmask] = xdirect[~nanmask]
                X[nanmask] = Xtemp
            else:
                X = Xtemp.copy()
            return X, _update_fractionation(em_inst0, X, river_col0, frac,
                                            sources), functioncost
        if OPTIMIZER_METHOD == "trust-constr":
            lb = np.array([b[0] if b[0] is not None else -np.inf for b in bounds])
            ub = np.array([b[1] if b[1] is not None else np.inf for b in bounds])
            res = minimize(
                cost_function, X0, args=args, method="trust-constr",
                bounds=Bounds(lb, ub),
                options={"maxiter": 1000 * nEM, "gtol": 1e-10, "xtol": 1e-12},
            )
        else:
            res = minimize(
                cost_function, X0, args=args,
                method="SLSQP", bounds=bounds,
                options={"maxiter": 1000 * nEM, "ftol": 1e-10},
            )
        Xtemp = res.x
        functioncost = res.fun
    else:
        Xtemp = X0

    # (4) reassemble the full solution vector
    nanmask = np.isnan(xdirect)
    X = np.full(nEM, np.nan)
    if not np.all(nanmask):
        X[~nanmask] = xdirect[~nanmask]
        X[nanmask] = Xtemp
    else:
        X = Xtemp.copy()

    # (5) update end-member matrix for fractionation-derived isotope values
    em_updated = _update_fractionation(em_inst0, X, river_col0, frac, sources)
    return X, em_updated, functioncost


def _update_fractionation(em_inst0, X, river_col0, frac, sources):
    """Post-inversion substitution of fractionation-derived isotope values
    into the end-member matrix (MEANDIR_InvertActiveSimulation step 5)."""
    em_updated = em_inst0.copy()
    with np.errstate(divide="ignore", invalid="ignore"):
        for i in range(frac.all.n):
            isopos0 = frac.all.isopos0[i]
            empos0 = frac.all.empos0[i]
            ionpos0 = frac.all.ionpos0[i]
            if frac.all.isopos_r[i] is None or frac.all.empos_r[i] is None:
                em_updated[isopos0, empos0] = em_updated[ionpos0, empos0] * (-1001)
            else:
                Xactive = em_inst0[ionpos0, :] * X / river_col0[ionpos0]
                em_iso = em_inst0[isopos0, :] / em_inst0[ionpos0, :]
                src = sources & (Xactive != 0)
                denom = np.sum(Xactive[src])
                source_iso = np.sum(Xactive[src] * em_iso[src]) / denom
                active_frac = em_inst0[isopos0, empos0] / em_inst0[ionpos0, empos0]
                new_iso = (source_iso + active_frac) + active_frac * Xactive[empos0] / denom
                em_updated[isopos0, empos0] = em_updated[ionpos0, empos0] * new_iso
    return em_updated
