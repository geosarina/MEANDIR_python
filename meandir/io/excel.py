"""Low-level Excel helpers.

MEANDIR's MATLAB code reads each worksheet with ``xlsread`` into a cell array
called ``all`` (1-based, header in row 1). We mirror that by loading the sheet
into a 0-based grid (list of rows of cell values, ``None`` for empty cells),
then provide small helpers for header-based column lookup.
"""

from __future__ import annotations

from pathlib import Path

import openpyxl


def read_grid(path: str | Path, sheet: str) -> list[list]:
    """Read a worksheet into a list-of-rows grid of cell values.

    Empty cells are ``None``. Trailing fully-empty rows are dropped so the grid
    matches the populated region (like ``xlsread``'s ``all``).
    """
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    try:
        ws = wb[sheet]
        rows = [list(r) for r in ws.iter_rows(values_only=True)]
    finally:
        wb.close()
    while rows and all(c is None for c in rows[-1]):
        rows.pop()
    return rows


def header_index(header: list, name: str) -> int:
    """Return the column index of ``name`` in a header row, or -1 if absent."""
    for i, h in enumerate(header):
        if h == name:
            return i
    return -1


def has_column(header: list, name: str) -> bool:
    return header_index(header, name) != -1
