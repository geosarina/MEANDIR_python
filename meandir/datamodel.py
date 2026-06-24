"""The ``river`` data container.

MEANDIR's MATLAB code accumulates everything into one big ``river`` struct with
nested fields (``river.observations.Ca_conc``, ``river.model_variable.norm`` …).
We mirror that with a light class whose sections are plain dicts of NumPy
arrays, so dynamic name-based access (``model_variable[obs]``) stays simple.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class River:
    """Container mirroring the MATLAB ``river`` struct."""

    n_samples: int = 0
    info: dict = field(default_factory=dict)
    observations: dict = field(default_factory=dict)
    model_variable: dict = field(default_factory=dict)
    settings: dict = field(default_factory=dict)
    # Result sections (populated by the calculation layer) live here too:
    results: dict = field(default_factory=dict)

    def obs(self, key: str) -> np.ndarray:
        """Return an observations array, or an all-NaN array if absent."""
        if key in self.observations:
            return self.observations[key]
        return np.full(self.n_samples, np.nan)
