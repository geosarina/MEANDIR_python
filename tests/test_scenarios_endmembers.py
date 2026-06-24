"""Tests for the scenario-parameter layer and the end-member parser."""

import math

import pytest

from meandir.scenarios import find_scenario_parameters, available_scenarios
from meandir.io.user_entries import load_endmember_group
from meandir.engine.endmembers import read_endmembers, FRAC_MIN, FRAC_MAX


def test_alaska_scenarios_registered():
    names = available_scenarios()
    assert "AK_scenario1_carbonate_slctindi" in names
    assert len([n for n in names if n.startswith("AK_scenario")]) == 5


def test_scenario1_derived(user_entries_path):
    p = find_scenario_parameters("AK_scenario1_carbonate_slctindi", user_entries_path)
    assert p.ObsList == ["Ca", "Mg", "Na", "K", "Cl", "SO4", "DIC", "d34S", "d13C"]
    assert p.NormalizationType == "SumObs"
    assert p.carbonisotopematch == "DIC"
    # IonCharges in ObsList order
    assert p.IonCharges == ["+", "+", "+", "+", "-", "-", "0", "0", "0"]
    assert p.CationsInNormalization == ["Ca", "Mg", "Na", "K"]
    assert p.AnionsInNormalization == ["SO4"]
    assert p.NeutralInNormalization == ["DIC"]
    assert p.nEM == 9
    assert p.ResetDegasDICContribution == 0


def test_scenario2_degas(user_entries_path):
    p = find_scenario_parameters("AK_scenario2_carbonate_slctindi_degas_2p5", user_entries_path)
    assert "degas" in p.EMList0
    assert p.nEM == 10
    assert p.ResetDegasDICContribution == 1
    assert p.DegasDICContributionMin == -2.5
    assert math.isinf(p.MaxFractionalContribution[1])  # carb max = inf


def test_endmember_parse(user_entries_path):
    tbl = load_endmember_group(user_entries_path, "PCK23_Alaska_SumCatSO4DIC_uniform")
    em = read_endmembers(tbl)
    assert em.disttype["Ca"] == "UNI"
    # carbonate min Ca/SumObs = 1/3
    assert em.get("carb", "Ca", "Min") == pytest.approx(1 / 3)
    assert em.get("slct_Ca", "Ca", "Min") == 1
    # degas d13C is a fractionation entry -> sentinel value + offset
    assert em.get("degas", "d13C", "Min") == FRAC_MIN
    assert em.get_addvalue("degas", "d13C", "Min") == pytest.approx(-10.1)
    assert em.get("degas", "d13C", "Max") == FRAC_MAX
    assert em.get_addvalue("degas", "d13C", "Max") == pytest.approx(-7.6)
    # corg d13C is a plain numeric range
    assert em.get("corg", "d13C", "Min") == pytest.approx(-30)
    assert em.get("corg", "d13C", "Max") == pytest.approx(-24)


def test_engineer_creek_not_yet_registered(user_entries_path):
    with pytest.raises(KeyError):
        find_scenario_parameters("does_not_exist", user_entries_path)
