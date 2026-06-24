"""Readers for the ``MEANDIR_UserEntries`` workbook.

Three sheets:

* ``MEANDIR_conc2equi``        -> :func:`load_conc2equi`  (``MEANDIR_DefineConc2Equi.m``)
* ``MEANDIR_DeltaNotationToR`` -> :func:`load_delta2r`    (``MEANDIR_ReadConvertDelta2R.m``)
* ``MEANDIR_Endmembers``       -> :func:`load_endmember_group` (raw table; the
  full distribution structure is built in :mod:`meandir.engine`)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .excel import read_grid, header_index

USER_ENTRIES_DEFAULT = "MEANDIR_UserEntries.xlsx"


def _is_number(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


# ---------------------------------------------------------------------------
# conc2equi
# ---------------------------------------------------------------------------
def load_conc2equi(path: str | Path) -> dict[str, float]:
    """Return ``{variable: conversion_factor}`` from the conc2equi sheet."""
    grid = read_grid(path, "MEANDIR_conc2equi")
    header = grid[0]
    ci_var = header_index(header, "Variable")
    ci_fac = header_index(header, "ConversionFactor")
    out: dict[str, float] = {}
    for row in grid[1:]:
        var = row[ci_var] if ci_var < len(row) else None
        fac = row[ci_fac] if ci_fac < len(row) else None
        if var is None:
            continue
        out[str(var)] = float(fac) if _is_number(fac) else float("nan")
    return out


# Map the spreadsheet "Charge" words to MEANDIR's single-character codes.
_CHARGE_CODE = {"positive": "+", "negative": "-", "neutral": "0"}


def load_charges(path: str | Path) -> dict[str, str]:
    """Return ``{variable: '+'/'-'/'0'}`` from the conc2equi sheet's Charge column."""
    grid = read_grid(path, "MEANDIR_conc2equi")
    header = grid[0]
    ci_var = header_index(header, "Variable")
    ci_chg = header_index(header, "Charge")
    out: dict[str, str] = {}
    for row in grid[1:]:
        var = row[ci_var] if ci_var < len(row) else None
        chg = row[ci_chg] if ci_chg < len(row) else None
        if var is None:
            continue
        out[str(var)] = _CHARGE_CODE.get(str(chg), "NaN")
    return out


# ---------------------------------------------------------------------------
# delta -> ratio conversion
# ---------------------------------------------------------------------------
@dataclass
class Delta2R:
    """Conversion factors from delta notation to absolute isotope ratios."""

    isotope: list[str] = field(default_factory=list)
    ratio: list[float] = field(default_factory=list)

    def factor(self, isotope: str) -> float:
        """Return the absolute ratio of the standard for ``isotope``."""
        return self.ratio[self.isotope.index(isotope)]

    def has(self, isotope: str) -> bool:
        return isotope in self.isotope


def load_delta2r(path: str | Path) -> Delta2R:
    grid = read_grid(path, "MEANDIR_DeltaNotationToR")
    header = grid[0]
    ci_iso = header_index(header, "IsotopeSystem")
    ci_val = header_index(header, "Values")
    out = Delta2R()
    for row in grid[1:]:
        iso = row[ci_iso] if ci_iso < len(row) else None
        val = row[ci_val] if ci_val < len(row) else None
        if iso is None:
            continue
        out.isotope.append(str(iso))
        out.ratio.append(float(val) if _is_number(val) else float("nan"))
    return out


# ---------------------------------------------------------------------------
# End-member raw table for one group
# ---------------------------------------------------------------------------
@dataclass
class EndMemberTable:
    """Raw end-member entries for a single end-member group.

    This is a faithful, lightly-structured view of the ``MEANDIR_Endmembers``
    sheet rows whose group (column A) matches ``group``. The MEANDIR
    distribution machinery (``MEANDIR_makeEMdistributions.m`` /
    ``MEANDIR_ReadEndMembersIntoMatlab.m``) consumes this in
    :mod:`meandir.engine`.

    Attributes
    ----------
    group:
        The end-member group name (matches ``EMdatasource``).
    em_codes:
        Short end-member codes from sheet row 2 (e.g. ``['carb', 'slct_Ca', ...]``).
    rows:
        One :class:`EndMemberRow` per data row for this group.
    """

    group: str
    em_codes: list[str] = field(default_factory=list)
    rows: list["EndMemberRow"] = field(default_factory=list)


@dataclass
class EndMemberRow:
    numerator: str          # column C, e.g. "Ca", "d34S"
    normalization: str      # column D, e.g. "SumObs", "Na"
    nature: str             # column E, e.g. "Min", "Max", "Mean", "Std"
    matlab_name: str        # column F, the constructed variable name
    values: dict[str, object]  # em_code -> numeric or string entry (frac#/sample#)


# Fixed sheet layout (0-based): the first 6 columns are metadata, EM columns
# start at index 6. Row 0 = headers, Row 1 = short codes, Row 2+ = data.
_META_COLS = 6
_COL_GROUP = 0
_COL_DIST = 1
_COL_NUM = 2
_COL_NORM = 3
_COL_NATURE = 4
_COL_NAME = 5


def list_endmember_groups(path: str | Path) -> list[str]:
    """Return the distinct end-member group names present in the workbook."""
    grid = read_grid(path, "MEANDIR_Endmembers")
    seen: list[str] = []
    for row in grid[2:]:
        g = row[_COL_GROUP] if _COL_GROUP < len(row) else None
        if g is not None and str(g) not in seen:
            seen.append(str(g))
    return seen


def load_endmember_group(path: str | Path, group: str) -> EndMemberTable:
    """Load the raw end-member table for ``group`` (== ``EMdatasource``)."""
    grid = read_grid(path, "MEANDIR_Endmembers")
    code_row = grid[1]
    em_codes = [str(code_row[i]) for i in range(_META_COLS, len(code_row))
                if code_row[i] is not None]
    n_em = len(em_codes)

    table = EndMemberTable(group=group, em_codes=em_codes)
    for row in grid[2:]:
        if (_COL_GROUP >= len(row)) or (str(row[_COL_GROUP]) != group):
            continue
        values: dict[str, object] = {}
        for j, code in enumerate(em_codes):
            col = _META_COLS + j
            values[code] = row[col] if col < len(row) else None
        table.rows.append(EndMemberRow(
            numerator=_cell_str(row, _COL_NUM),
            normalization=_cell_str(row, _COL_NORM),
            nature=_cell_str(row, _COL_NATURE),
            matlab_name=_cell_str(row, _COL_NAME),
            values=values,
        ))
    if not table.rows:
        groups = list_endmember_groups(path)
        raise KeyError(
            f"End-member group {group!r} not found. Available groups: {groups}"
        )
    return table


def _cell_str(row: list, idx: int) -> str:
    return str(row[idx]) if idx < len(row) and row[idx] is not None else ""
