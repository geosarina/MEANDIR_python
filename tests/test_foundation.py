"""Validate the foundational I/O + matrix layer against the Alaska workbooks."""

import numpy as np
import pytest

from meandir.go import go
from meandir.constants import ALL_VARIABLES
from meandir.io.user_entries import (
    load_conc2equi, load_delta2r, load_endmember_group, list_endmember_groups,
)
from meandir.io.river_data import read_river_observations, build_model_variables
from meandir.river_matrix import generate_river_matrix


# ---- conc2equi -------------------------------------------------------------
def test_conc2equi_known_factors(user_entries_path):
    c2e = load_conc2equi(user_entries_path)
    assert c2e["Ca"] == 2
    assert c2e["Mg"] == 2
    assert c2e["Na"] == 1
    assert c2e["SO4"] == 2
    assert c2e["Cl"] == 1
    assert c2e["PO4"] == 3
    assert c2e["d34S"] == 1  # isotopes are neutral, factor 1


# ---- delta -> ratio --------------------------------------------------------
def test_delta2r_factors(user_entries_path):
    d2r = load_delta2r(user_entries_path)
    assert d2r.has("d34S")
    assert d2r.factor("d34S") == pytest.approx(0.04416258898761681)
    assert d2r.factor("d13C") == pytest.approx(0.01118)
    assert d2r.factor("d7Li") == pytest.approx(12.175596926415688)


# ---- end-member table ------------------------------------------------------
def test_endmember_group(user_entries_path):
    groups = list_endmember_groups(user_entries_path)
    assert "PCK23_Alaska_SumCatSO4DIC_uniform" in groups

    tbl = load_endmember_group(user_entries_path, "PCK23_Alaska_SumCatSO4DIC_uniform")
    assert tbl.em_codes == [
        "carb", "slct_Ca", "slct_Mg", "slct_Na", "slct_K",
        "evap", "prec", "pyri", "degas", "corg",
    ]
    # A known Ca/SumObs Min entry: carbonate min Ca/SumObs = 1/3.
    ca_min = [r for r in tbl.rows if r.numerator == "Ca" and r.nature == "Min"][0]
    assert ca_min.values["carb"] == pytest.approx(1 / 3)
    assert ca_min.values["slct_Ca"] == 1
    # d13C of corg carries a fractionation string entry on degas.
    d13c_min = [r for r in tbl.rows if r.numerator == "d13C" and r.nature == "Min"][0]
    assert str(d13c_min.values["degas"]).startswith("frac#")


# ---- river data ------------------------------------------------------------
def test_read_river_data(river_data_path, user_entries_path):
    c2e = load_conc2equi(user_entries_path)
    river = read_river_observations(river_data_path, c2e)
    assert river.n_samples == 198  # 199 rows - header
    # equi = conc * factor for a charged species
    ca_conc = river.observations["Ca_conc"]
    ca_equi = river.observations["Ca_equi"]
    finite = ~np.isnan(ca_conc)
    assert np.allclose(ca_equi[finite], ca_conc[finite] * 2)


def test_go_flags():
    g = go(["Na", "Ca", "Mg", "SO4", "d34S"])
    assert g["Na"] and g["Ca"] and g["SO4"] and g["d34S"]
    assert not g["Sr"]
    # every canonical variable has a flag
    assert set(g) == set(ALL_VARIABLES)


# ---- river matrix end-to-end ----------------------------------------------
def test_river_matrix_shapes(river_data_path, user_entries_path):
    c2e = load_conc2equi(user_entries_path)
    d2r = load_delta2r(user_entries_path)
    river = read_river_observations(river_data_path, c2e)

    obs_list = ["Ca", "Mg", "Na", "K", "SO4", "DIC"]
    build_model_variables(river, obs_in_normalization=obs_list,
                          em_units="equi", obs_list=obs_list)
    M = generate_river_matrix(river, obs_list, carbonisotopematch="DIC",
                              convert_delta2r_list=[], delta2r=d2r)
    assert M.shape == (len(obs_list), river.n_samples)

    # Each non-isotopic row is variable/norm; columns with full data sum sensibly.
    norm = river.model_variable["norm"]
    row_ca = M[0]
    finite = ~np.isnan(row_ca)
    assert np.allclose(row_ca[finite],
                       (river.model_variable["Ca"] / norm)[finite])
