"""Reader for the river-observation spreadsheet.

Port of ``MEANDIR_StandardInput.m`` (column parsing + unit conversion) plus the
normalization / ``model_variable`` construction at the bottom of
``MEANDIR_ReadRiverDataIntoMatlab.m``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..constants import ALL_VARIABLES, spec
from ..datamodel import River
from .excel import read_grid, header_index

RIVER_SHEET_DEFAULT = "MyRiverData"

# Info columns carried alongside the observations.
_INFO_TEXT = ["Name", "SampleNumber", "ExtraField1", "ExtraField2",
              "ExtraField3", "ExtraField4", "ExtraField5"]
_INFO_NUMERIC = {"ClCritical": "ClCr", "Latitude_oN": "latd", "Longitude_oE": "long"}


def _is_number(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def _column(data_rows, idx) -> list:
    return [(row[idx] if idx < len(row) else None) for row in data_rows]


def _numeric_array(data_rows, idx, scale=1.0) -> np.ndarray:
    """Build a float array from a column, NaN where the cell is non-numeric."""
    vals = _column(data_rows, idx)
    out = np.full(len(vals), np.nan)
    for i, v in enumerate(vals):
        if _is_number(v):
            out[i] = float(v) * scale
    return out


def read_river_observations(path: str | Path, conc2equi: dict,
                            sheet: str = RIVER_SHEET_DEFAULT) -> River:
    """Read river observations into a :class:`River` (StandardInput port)."""
    grid = read_grid(path, sheet)
    header = grid[0]
    data_rows = grid[1:]
    n = len(data_rows)

    river = River(n_samples=n)

    # (1) info / text fields
    for col in _INFO_TEXT:
        idx = header_index(header, col)
        river.info[col] = _column(data_rows, idx) if idx != -1 else [None] * n
    for col, name in _INFO_NUMERIC.items():
        idx = header_index(header, col)
        arr = _numeric_array(data_rows, idx) if idx != -1 else np.full(n, np.nan)
        river.info[name] = arr
    river.model_variable["ClCr"] = river.info["ClCr"]

    # (2) observations + (3) uncertainties, with conc/equi pairs
    for name in ALL_VARIABLES:
        s = spec(name)
        factor = conc2equi.get(name, np.nan)
        vi = header_index(header, s.column)
        if vi != -1:
            conc = _numeric_array(data_rows, vi, s.scale)
            river.observations[f"{name}_conc"] = conc
            river.observations[f"{name}_equi"] = conc * factor
        ui = header_index(header, s.unc_column)
        if ui != -1:
            conc_asd = _numeric_array(data_rows, ui, s.scale)
            river.observations[f"{name}_conc_asd"] = conc_asd
            river.observations[f"{name}_equi_asd"] = conc_asd * factor

    # (4) convenience sums used for ternary plots (only when both present)
    def _sum(a, b, fa, fb, key):
        ka, kb = f"{a}_conc", f"{b}_conc"
        if ka in river.observations and kb in river.observations:
            river.observations[f"{key}_conc"] = (
                river.observations[ka] + river.observations[kb])
            river.observations[f"{key}_equi"] = (
                river.observations[ka] * fa + river.observations[kb] * fb)

    cq = conc2equi
    _sum("Na", "K", cq.get("Na", np.nan), cq.get("K", np.nan), "NaK")
    _sum("Ca", "Mg", cq.get("Ca", np.nan), cq.get("Mg", np.nan), "CaMg")
    _sum("Mg", "Na", cq.get("Mg", np.nan), cq.get("Na", np.nan), "MgNa")
    _sum("Ca", "Na", cq.get("Ca", np.nan), cq.get("Na", np.nan), "CaNa")

    return river


def build_model_variables(river: River, obs_in_normalization, em_units: str,
                          obs_list) -> np.ndarray:
    """Construct ``norm`` and ``model_variable`` and the functional-sample mask.

    Mirrors the bottom of ``MEANDIR_ReadRiverDataIntoMatlab.m``. Returns the
    boolean ``functionalsamplelist`` (samples having every ObsList variable).
    """
    n = river.n_samples
    norm_conc = np.zeros(n)
    norm_equi = np.zeros(n)
    for v in obs_in_normalization:
        norm_conc = norm_conc + river.obs(f"{v}_conc")
        norm_equi = norm_equi + river.obs(f"{v}_equi")
    river.observations["norm_conc"] = norm_conc
    river.observations["norm_equi"] = norm_equi
    river.model_variable["norm"] = norm_conc if em_units == "conc" else norm_equi

    for v in obs_list:
        river.model_variable[v] = river.obs(f"{v}_{em_units}")
        river.model_variable[f"{v}_asd"] = river.obs(f"{v}_{em_units}_asd")

    # Always carry SO4 forward even if not inverted.
    if "SO4" not in obs_list:
        river.model_variable["SO4"] = river.obs(f"SO4_{em_units}")
        river.model_variable["SO4_asd"] = river.obs(f"SO4_{em_units}_asd")

    present = np.ones(n, dtype=bool)
    for v in obs_list:
        present &= ~np.isnan(river.model_variable[v])
    river.info["numbersampleswithdata"] = int(present.sum())
    return present
