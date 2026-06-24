"""MEANDIR — Mixing Elements ANd Dissolved Isotopes in Rivers (Python port).

A faithful Python port of the MEANDIR MATLAB model of Kemeny & Torres for
Monte Carlo inversion of dissolved river chemistry.

Original MATLAB model: Preston Cosslett Kemeny & Mark Albert Torres,
"Presentation and applications of Mixing Elements ANd Dissolved Isotopes in
Rivers (MEANDIR), a customizable MATLAB model for Monte Carlo inversion of
dissolved river chemistry", American Journal of Science (2021); v1.3 (2023).

This port targets *algorithmic and statistical* fidelity: identical equations
and model structure, with results validated to match the MATLAB output within
Monte Carlo sampling noise (bit-for-bit reproduction is not possible across the
MATLAB/Python RNG and optimizer boundary).
"""

from .constants import ALL_VARIABLES, ISOTOPE_VARIABLES
from .go import go

__all__ = ["ALL_VARIABLES", "ISOTOPE_VARIABLES", "go"]
__version__ = "0.1.0"
