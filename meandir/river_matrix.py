"""Port of ``MEANDIR_GenerateRiverMatrix.m``.

Builds ``RiverMatrix0``: rows are the normalized elemental/isotopic ratios in
``ObsList`` order, columns are samples. Isotopic rows are the (optionally
delta->ratio converted) isotope value times its carrier-ion concentration,
divided by the normalization.
"""

from __future__ import annotations

import numpy as np

from .constants import ISOTOPE_VARIABLES, isotope_ion


def generate_river_matrix(river, obs_list, carbonisotopematch,
                          convert_delta2r_list, delta2r) -> np.ndarray:
    n = river.n_samples
    rivernorm = river.model_variable["norm"]
    matrix = np.full((len(obs_list), n), np.nan)

    for i, name in enumerate(obs_list):
        if name not in ISOTOPE_VARIABLES:
            matrix[i, :] = river.model_variable[name] / rivernorm
            continue

        deltavalue = river.model_variable[name]
        ion = isotope_ion(name, carbonisotopematch)
        if name in convert_delta2r_list:
            conversionfactor = delta2r.factor(name)
            scaled = (deltavalue / 1000.0 + 1.0) * conversionfactor
        else:
            # User-supplied value used directly (delta value or raw ratio).
            scaled = deltavalue
        # The carrier ion is itself in ObsList (MATLAB find(ismember(...))),
        # so it is present in model_variable.
        numerator = river.model_variable[ion]
        matrix[i, :] = scaled * numerator / rivernorm

    return matrix
