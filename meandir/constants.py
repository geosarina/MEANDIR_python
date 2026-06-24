"""Canonical variable definitions for MEANDIR.

This module centralizes the lists of dissolved observations that MEANDIR can
handle, the mapping between spreadsheet column headers and internal variable
names, and the isotope bookkeeping. It is the Python equivalent of the
information that is spread across ``MEANDIR_Go.m``, ``MEANDIR_StandardInput.m``,
``MEANDIR_DefineBasicFields.m`` and ``MEANDIR_GenerateRiverMatrix.m``.

Fidelity note: the *order* of ``ALL_VARIABLES`` matches the order in which the
MATLAB code declares fields, and the spreadsheet column headers / unit scales
reproduce ``MEANDIR_StandardInput.m`` exactly (including the ``/1e6`` rescaling
of the picomolar trace metals Fe, B, Re, Mo, Os).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class VariableSpec:
    """Describes one dissolved observation MEANDIR understands."""

    name: str          # internal name, e.g. "Ca", "d34S", "Sr8786"
    column: str        # spreadsheet header for the value, e.g. "Ca_uM"
    scale: float       # multiply the raw spreadsheet value by this (unit fix)
    kind: str          # "conc" | "isotope_delta" | "isotope_ratio"

    @property
    def unc_column(self) -> str:
        """Header of the 1-sigma absolute uncertainty column."""
        return f"Unc_{self.column}_1sigma"


# ---------------------------------------------------------------------------
# The full, ordered list of variables (matches MEANDIR_DefineBasicFields.m /
# MEANDIR_Go.m). Concentrations are carried internally in micromolar (uM).
# ---------------------------------------------------------------------------
_CONC_UM = [
    "ALK", "DIC", "Ca", "Mg", "Na", "K", "Sr", "Cl", "SO4",
    "NO3", "PO4", "Si", "Ge", "Li", "F", "HCO3",
]
_CONC_PM = ["Fe", "B", "Re", "Mo", "Os"]  # reported picomolar, stored as uM (/1e6)
_ISO_DELTA = [
    "d7Li", "d13C", "d18O", "d26Mg", "d30Si", "d34S",
    "d42Ca", "d44Ca", "d56Fe", "d98Mo",
]
_ISO_RATIO = {  # name -> spreadsheet column header
    "Sr8786": "Sr8786_ratio",
    "Os8788": "Os187188_ratio",
}

_SPECS: dict[str, VariableSpec] = {}
for _n in _CONC_UM:
    _SPECS[_n] = VariableSpec(_n, f"{_n}_uM", 1.0, "conc")
for _n in _CONC_PM:
    _SPECS[_n] = VariableSpec(_n, f"{_n}_pM", 1.0e-6, "conc")
for _n in _ISO_DELTA:
    _SPECS[_n] = VariableSpec(_n, f"{_n}_permil", 1.0, "isotope_delta")
for _n, _col in _ISO_RATIO.items():
    _SPECS[_n] = VariableSpec(_n, _col, 1.0, "isotope_ratio")
# Fmod is a special "modern fraction" isotope-like variable with no unit suffix.
_SPECS["Fmod"] = VariableSpec("Fmod", "Fmod", 1.0, "isotope_ratio")

# Canonical declaration order (used for deterministic iteration / output).
ALL_VARIABLES: list[str] = (
    ["ALK", "DIC", "Ca", "Mg", "Na", "K", "Sr", "Fe", "Cl", "SO4",
     "NO3", "PO4", "Si", "Ge", "Li", "F", "B", "Re", "Mo", "Os", "HCO3"]
    + ["d7Li", "d13C", "d18O", "d26Mg", "d30Si", "d34S", "d42Ca", "d44Ca",
       "d56Fe", "Sr8786", "d98Mo", "Os8788", "Fmod"]
)

# Subset that MEANDIR treats as isotopic (matches the ``isotopeposition``
# logic at the top of MEANDIR_Master.m).
ISOTOPE_VARIABLES: frozenset[str] = frozenset({
    "d34S", "d18O", "Sr8786", "d7Li", "d26Mg", "d30Si", "d42Ca", "d44Ca",
    "d56Fe", "d98Mo", "Os8788", "d13C", "Fmod",
})

# For each isotopic variable, which elemental concentration it rides on when
# building the river matrix (MEANDIR_GenerateRiverMatrix.m). d13C and Fmod are
# resolved at runtime against ``carbonisotopematch`` ("HCO3" or "DIC").
ISOTOPE_ION: dict[str, str] = {
    "d7Li": "Li",
    "d18O": "SO4",
    "d26Mg": "Mg",
    "d30Si": "Si",
    "d34S": "SO4",
    "d42Ca": "Ca",
    "d44Ca": "Ca",
    "d56Fe": "Fe",
    "Sr8786": "Sr",
    "Os8788": "Os",
    "d98Mo": "Mo",
    # "d13C" and "Fmod" depend on carbonisotopematch -> handled in code.
}


def spec(name: str) -> VariableSpec:
    """Return the :class:`VariableSpec` for ``name`` (KeyError if unknown)."""
    return _SPECS[name]


def is_known(name: str) -> bool:
    return name in _SPECS


def isotope_ion(name: str, carbonisotopematch: str | None = None) -> str:
    """Return the concentration variable an isotopic variable normalizes to."""
    if name in ("d13C", "Fmod"):
        if carbonisotopematch not in ("HCO3", "DIC"):
            raise ValueError(
                f"{name} requires carbonisotopematch in ('HCO3','DIC'), "
                f"got {carbonisotopematch!r}"
            )
        return carbonisotopematch
    return ISOTOPE_ION[name]
