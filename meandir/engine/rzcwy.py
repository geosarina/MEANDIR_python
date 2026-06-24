"""R, Z, C, W, Y weathering metrics — port of ``MEANDIR_CalculateRZCWY.m``.

For each sample/instance these combine the inversion-constrained river chemistry
with the fractional contributions to express weathering ratios:

  * R = (carbonate-sourced cations) / (all silicate+carbonate cations)
  * Z = pyrite-derived SO4 / cation denominator (sulfuric-acid weathering)
  * C = organic-carbon-derived DIC / cation denominator
  * W = pyrite SO4 / normalization;  Y = pyrite SO4 / total river SO4

Each is produced gross/net (net adds sink/source recycling) and
scaled/unscaled (scaled divides each contribution by the total reconstructed
fraction). Units follow ``EMUnits`` ('equi' -> multiplier 1).
"""

from __future__ import annotations

import numpy as np

from .results import _summ, _STATS

_VARIANTS = (("gross", "unscaled"), ("gross", "scaled"),
             ("net", "unscaled"), ("net", "scaled"))


def calculate_rzcwy(scenario_results, conc2equi):
    p = scenario_results.params
    if not p.CalculateRZCWY:
        return {}

    ems = p.EMList0
    river = scenario_results.river
    n = len(river.model_variable["norm"])
    equi = p.EMUnits == "equi"
    src = np.array([e in p.EMsources for e in ems])
    snk = np.array([e in p.EMsinks for e in ems])

    def mult(ion):
        return 1.0 if equi else conc2equi[ion]

    so4mult = 1.0 if equi else conc2equi["SO4"]

    # output skeleton: V -> gross/net -> scaled/unscaled -> stat -> array(n)
    out = {V: {g: {s: {k: np.full(n, np.nan) for k in _STATS}
                   for s in ("unscaled", "scaled")}
               for g in ("gross", "net")} for V in "RZCWY"}

    for i in sorted(scenario_results.sample_results):
        succ = scenario_results.sample_results[i]
        if not succ:
            continue
        ninst = len(succ)

        def fr(ion):
            return np.stack([s["fractions"][ion] for s in succ], axis=1)  # (nEM,ninst)

        def inv(ion):
            return np.array([s["inv"][ion] for s in succ])

        def accumulate(ionlist, emlist, mfn):
            acc = {v: np.zeros(ninst) for v in _VARIANTS}
            emidx = [ems.index(e) for e in emlist]
            for ion in ionlist:
                f = fr(ion)
                rv = inv(ion)
                m = mfn(ion)
                allnet = f.sum(0)
                allsrc = f[src].sum(0) if src.any() else np.zeros(ninst)
                allsnk = f[snk].sum(0) if snk.any() else np.zeros(ninst)
                with np.errstate(divide="ignore", invalid="ignore"):
                    netfac = 1 + (allsnk / allsrc if src.any() else 0.0)
                    for k in emidx:
                        fk = f[k]
                        acc[("gross", "unscaled")] += m * (rv * fk)
                        acc[("gross", "scaled")] += m * (rv * fk / allnet)
                        acc[("net", "unscaled")] += m * (rv * fk * netfac)
                        acc[("net", "scaled")] += m * (rv * fk / allnet * netfac)
            return acc

        den = accumulate(p.RZC_Denominator_IonList, p.RZC_Denominator_EMList, mult)
        rnum = accumulate(p.R_Numerator_IonList, p.R_Numerator_EMList, mult)

        # --- Z numerator: SO4 from pyrite ---
        znum = {v: np.full(ninst, np.nan) for v in _VARIANTS}
        inv_so4 = np.full(ninst, np.nan)
        if p.ZfromEM:
            znum = accumulate(["SO4"], p.Z_Numerator_EMList, lambda _ion: 1.0)
            inv_so4 = inv("SO4")
        elif p.ZfromriverSO4:
            inv_so4 = inv("SO4")
            znum[("gross", "unscaled")] = inv_so4.copy()
            znum[("net", "unscaled")] = inv_so4.copy()
        elif p.ZfromSO4excess:
            so4_obs = river.model_variable["SO4"][i]
            frac_src_so4 = fr("SO4")[src].sum(0) if "SO4" in succ[0]["fractions"] else 0.0
            excess = so4_obs * (1 - frac_src_so4)
            inv_so4 = np.full(ninst, so4_obs)
            znum[("gross", "unscaled")] = excess
            znum[("net", "unscaled")] = excess
        # Znotcalculated -> all NaN

        # --- C numerator: DIC (or HCO3) from organic-carbon oxidation ---
        cnum = None
        if p.C_Numerator_EMList:
            carbon = "DIC" if p.carbonisotopematch == "DIC" else "HCO3"
            if carbon in succ[0]["fractions"]:
                cnum = accumulate([carbon], p.C_Numerator_EMList, lambda _ion: 1.0)

        # normalization (denominator of W)
        if equi:
            inv_norm = inv("norm")
        else:
            inv_norm = np.zeros(ninst)
            for ion in p.ObsInNormalization:
                inv_norm = inv_norm + conc2equi[ion] * inv(ion)

        # assemble the five metrics, per variant
        with np.errstate(divide="ignore", invalid="ignore"):
            for v in _VARIANTS:
                g, s = v
                metrics = {
                    "R": rnum[v] / den[v],
                    "Z": so4mult * znum[v] / den[v],
                    "Y": znum[v] / inv_so4,
                    "W": so4mult * znum[v] / inv_norm,
                }
                if cnum is not None:
                    metrics["C"] = cnum[v] / den[v]
                for V, arr in metrics.items():
                    summ = _summ(arr, axis=0)
                    for k in _STATS:
                        out[V][g][s][k][i] = summ[k]
    return out
