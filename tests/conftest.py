"""Shared test fixtures: paths to the vendored Alaska reference workbooks."""

from pathlib import Path

import pytest

REF_DATA = Path(__file__).resolve().parent.parent / "reference" / "data"


@pytest.fixture(scope="session")
def user_entries_path() -> Path:
    return REF_DATA / "MEANDIR_UserEntries.xlsx"


@pytest.fixture(scope="session")
def river_data_path() -> Path:
    return REF_DATA / "RiverDataSpreadsheet_Kemeny_etal_2023.xlsx"
