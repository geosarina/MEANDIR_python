"""Scenario parameters — port of ``MEANDIR_FindScenarioParameters.m``.

Defines the :class:`ScenarioParameters` dataclass (the ~44 inversion variables,
plus the group-6 degassing controls), the registry of the five published Alaska
scenarios, and :func:`derive` which computes the derived quantities the engine
needs (``IonCharges``, the cation/anion/neutral normalization sub-lists,
``NormalizationType``, ``carbonisotopematch``).

User-defined scenarios (e.g. Engineer Creek) are added by registering a new
:class:`ScenarioParameters` instance via :func:`register`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

from .io.user_entries import load_charges


@dataclass
class ScenarioParameters:
    name: str

    # --- group 1: river observations & normalization ---
    Riverdatasource: str = ""
    AdjustRiverObs: int = 0
    ObsList: list = field(default_factory=list)
    CostFunType: list = field(default_factory=list)
    WeightingList: list = field(default_factory=list)
    ErrorCutMinMB: list = field(default_factory=list)
    ErrorCutMaxMB: list = field(default_factory=list)
    nCFList: list = field(default_factory=list)
    ConvertDelta2RList: list = field(default_factory=list)
    AllIonsExplicitlyResolved: int = 0
    ObsInNormalization: list = field(default_factory=list)
    ImposeNormalizationCheck: int = 1

    # --- group 2: general settings ---
    Solver: str = "mldivide_optimize"
    IterateOver: str = "Samples"
    maxiterations: float = 5e5
    maxsuccess: int = 200
    maxzerohits: float = 1e4
    numberiterations: float = float("nan")
    MisfitCuts: float = float("nan")
    CullOn: str = "EachSample"
    saveuncutdata: int = 0

    # --- group 3: end-members ---
    EMdatasource: str = ""
    EMList0: list = field(default_factory=list)
    MinFractionalContribution: list = field(default_factory=list)
    MaxFractionalContribution: list = field(default_factory=list)
    ListNormClosure: list = field(default_factory=list)
    ListChargeClosure: list = field(default_factory=list)
    EMUnits: str = "equi"
    EMsources: list = field(default_factory=list)
    EMsinks: list = field(default_factory=list)
    EndMembersWithNegativeRatios: list = field(default_factory=list)
    CoupleFeS2SO4intoEM: list = field(default_factory=list)
    CoupleFeS2d34SintoEM: list = field(default_factory=list)
    RecordFullFeS2Distribution: int = 0
    BalanceEvaporite: int = 0

    # --- group 4: Cl correction ---
    PrecProcessing: str = "EndMember"
    ClCriticalValuesGiven: int = 1

    # --- group 5: R/Z/C/W/Y ---
    CalculateRZCWY: int = 0
    R_Numerator_EMList: list = field(default_factory=list)
    R_Numerator_IonList: list = field(default_factory=list)
    Z_NumeratorType: list = field(default_factory=list)
    Z_Numerator_EMList: list = field(default_factory=list)
    C_Numerator_EMList: list = field(default_factory=list)
    RZC_Denominator_EMList: list = field(default_factory=list)
    RZC_Denominator_IonList: list = field(default_factory=list)

    # --- group 6: degassing DIC contribution ---
    ResetDegasDICContribution: int = 0
    DegasDICContributionMin: float = float("nan")
    DegasDICContributionMax: float = float("nan")

    # --- originals kept for the degassing reset (MATLAB *0 variables) ---
    MinFractionalContribution0: list = field(default_factory=list)
    MaxFractionalContribution0: list = field(default_factory=list)

    # --- derived (filled by derive()) ---
    IonCharges: list = field(default_factory=list)
    CationsInNormalization: list = field(default_factory=list)
    AnionsInNormalization: list = field(default_factory=list)
    NeutralInNormalization: list = field(default_factory=list)
    NormalizationType: str = ""
    carbonisotopematch: str = "NaN"

    @property
    def nEM(self) -> int:
        return len(self.EMList0)

    @property
    def nOL(self) -> int:
        return len(self.ObsList)


def derive(params: ScenarioParameters, charges: dict[str, str]) -> ScenarioParameters:
    """Compute derived quantities (section 2 of MEANDIR_FindScenarioParameters)."""
    obs = params.ObsList

    # IonCharges, in ObsList order.
    params.IonCharges = [charges.get(o, "NaN") for o in obs]

    # cation/anion/neutral sub-lists of the normalization, in conc2equi order is
    # not required downstream (membership tests only), so ObsInNormalization order is fine.
    norm = params.ObsInNormalization
    params.CationsInNormalization = [o for o in norm if charges.get(o) == "+"]
    params.AnionsInNormalization = [o for o in norm if charges.get(o) == "-"]
    params.NeutralInNormalization = [o for o in norm if charges.get(o) == "0"]

    # NormalizationType: a single observation name, or "SumObs".
    if len(norm) == 1:
        params.NormalizationType = norm[0]
    else:
        params.NormalizationType = "SumObs"

    # carbonisotopematch (test 19): DIC if present, else HCO3.
    if "DIC" in obs and "HCO3" in obs:
        params.carbonisotopematch = "NaN"
    elif "DIC" in obs:
        params.carbonisotopematch = "DIC"
    else:
        params.carbonisotopematch = "HCO3"

    # Keep originals for the degassing reset if not already set.
    if not params.MinFractionalContribution0:
        params.MinFractionalContribution0 = list(params.MinFractionalContribution)
    if not params.MaxFractionalContribution0:
        params.MaxFractionalContribution0 = list(params.MaxFractionalContribution)
    return params


# ---------------------------------------------------------------------------
# Registry of the five published Alaska scenarios (Kemeny et al., 2023).
# ---------------------------------------------------------------------------
_INF = math.inf

_ALASKA_COMMON = dict(
    Riverdatasource="PCKAlaskaData",
    AdjustRiverObs=1,
    ObsList=["Ca", "Mg", "Na", "K", "Cl", "SO4", "DIC", "d34S", "d13C"],
    CostFunType=["rel"] * 9,
    WeightingList=[1] * 9,
    ErrorCutMinMB=[95, 95, 95, 95, 95, 95, 95, -1, -1],
    ErrorCutMaxMB=[105, 105, 105, 105, 105, 105, 105, 1, 1],
    nCFList=[],
    ConvertDelta2RList=["d13C", "d34S"],
    AllIonsExplicitlyResolved=0,
    ObsInNormalization=["Ca", "Mg", "Na", "K", "SO4", "DIC"],
    ImposeNormalizationCheck=1,
    Solver="mldivide_optimize",
    IterateOver="Samples",
    maxiterations=5e5,
    maxsuccess=200,
    maxzerohits=1e4,
    numberiterations=float("nan"),
    MisfitCuts=float("nan"),
    CullOn="EachSample",
    saveuncutdata=0,
    EMUnits="equi",
    EMsources=["prec", "carb", "slct_Ca", "slct_Mg", "slct_Na", "slct_K",
               "pyri", "evap", "corg"],
    EMsinks=[],
    EndMembersWithNegativeRatios=[],
    CoupleFeS2SO4intoEM=[],
    CoupleFeS2d34SintoEM=[],
    RecordFullFeS2Distribution=0,
    BalanceEvaporite=1,
    PrecProcessing="EndMember",
    ClCriticalValuesGiven=1,
    CalculateRZCWY=1,
    R_Numerator_EMList=["carb"],
    R_Numerator_IonList=["Na", "Ca", "Mg", "K"],
    Z_NumeratorType=["ZfromEM"],
    Z_Numerator_EMList=["pyri"],
    C_Numerator_EMList=["corg"],
    RZC_Denominator_EMList=["carb", "slct_Ca", "slct_Mg", "slct_Na", "slct_K", "pyri"],
    RZC_Denominator_IonList=["Na", "Ca", "Mg", "K"],
)

_EM9 = ["prec", "carb", "slct_Ca", "slct_Mg", "slct_Na", "slct_K", "pyri", "evap", "corg"]
_EM10 = _EM9 + ["degas"]
_NORM9 = ["Ca", "Ca", "Ca", "Mg", "Na", "K", "SO4", "Ca", "DIC"]
_NORM10 = _NORM9 + ["DIC"]


def _alaska_registry() -> dict[str, ScenarioParameters]:
    reg: dict[str, ScenarioParameters] = {}

    # Scenario 1: no degassing end-member.
    reg["AK_scenario1_carbonate_slctindi"] = ScenarioParameters(
        name="AK_scenario1_carbonate_slctindi",
        EMdatasource="PCK23_Alaska_SumCatSO4DIC_uniform",
        EMList0=list(_EM9),
        MinFractionalContribution=[0] * 9,
        MaxFractionalContribution=[1] * 9,
        ListNormClosure=list(_NORM9),
        ListChargeClosure=[],
        ResetDegasDICContribution=0,
        **_ALASKA_COMMON,
    )

    # Scenarios 2-5 add a degassing end-member; differ in EM group / degas range.
    degas_min_max = {
        "AK_scenario2_carbonate_slctindi_degas_2p5": ("PCK23_Alaska_SumCatSO4DIC_uniform", -2.5),
        "AK_scenario3_carbonate_slctindi_degas_25": ("PCK23_Alaska_SumCatSO4DIC_uniform", -25.0),
        "AK_scenario4_carbonate_slctindi_degas_2p5_highfrac": ("PCK23_Alaska_SumCatSO4DIC_uniform_highfrac", -2.5),
        "AK_scenario5_carbonate_slctindi_degas_2p5_lowfrac": ("PCK23_Alaska_SumCatSO4DIC_uniform_lowfrac", -2.5),
    }
    for sname, (emsrc, degas_min) in degas_min_max.items():
        reg[sname] = ScenarioParameters(
            name=sname,
            EMdatasource=emsrc,
            EMList0=list(_EM10),
            MinFractionalContribution=[0, 0, 0, 0, 0, 0, 0, 0, 0, -9.99],
            MaxFractionalContribution=[1, _INF, 1, 1, 1, 1, 1, 1, _INF, -9.99],
            ListNormClosure=list(_NORM10),
            ListChargeClosure=[],
            ResetDegasDICContribution=1,
            DegasDICContributionMin=degas_min,
            DegasDICContributionMax=0.0,
            **_ALASKA_COMMON,
        )
    return reg


_REGISTRY: dict[str, ScenarioParameters] = _alaska_registry()


def register(params: ScenarioParameters) -> None:
    """Register a user-defined scenario."""
    _REGISTRY[params.name] = params


def available_scenarios() -> list[str]:
    return list(_REGISTRY)


def find_scenario_parameters(name: str, user_entries_path) -> ScenarioParameters:
    """Return a fully-derived copy of the named scenario's parameters."""
    if name not in _REGISTRY:
        raise KeyError(
            f"Unknown scenario {name!r}. Available: {available_scenarios()}"
        )
    params = replace(_REGISTRY[name])  # shallow copy so derive() doesn't mutate the registry
    charges = load_charges(user_entries_path)
    return derive(params, charges)
