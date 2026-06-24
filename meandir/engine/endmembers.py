"""Parse raw end-member tables into a structured form — port of the numeric
semantics of ``MEANDIR_ReadEndMembersIntoMatlab.m``.

Rather than reproduce the MATLAB struct field-name gymnastics (e.g.
``EM.group.carb.NOR_MenCaNa``), we keep a clean indexed structure
``EMData[(em, numerator, nature)]`` that the distribution builder consumes.
The numeric behaviour — sentinel substitution for ``sample#``/``frac#``
entries, the NaN→0 fill of active end-members, and the per-numerator
distribution type — is reproduced exactly.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..io.user_entries import EndMemberTable

# Sentinel codes (identical to ReadEndMembersIntoMatlab). The "mean" sentinels
# are imaginary, matching MATLAB's ``123456789i`` / ``314159265i``.
SAMPLE_MIN = -123456789.0
SAMPLE_MAX = 123456789.0
SAMPLE_MEN = 123456789j
FRAC_MIN = -314159265.0
FRAC_MAX = 314159265.0
FRAC_MEN = 314159265j

_SAMPLE_WORDS = {"Sample", "sample", "Samples", "samples",
                 "sampls", "Sampls", "sampl", "Sampl"}
_FRAC_WORDS = {"frac", "Frac", "fractionation", "Fractionation",
               "Fractination", "fracs"}


@dataclass
class EMData:
    """Structured end-member chemistry for one group.

    value[(em, numerator, nature)] holds the (possibly sentinel-substituted)
    numeric entry; addvalue[...] holds the offset after a ``#`` (NaN if none).
    nature is the 3-char code 'Min'/'Max'/'Men'/'Sig'.
    """

    group: str
    em_codes: list = field(default_factory=list)
    disttype: dict = field(default_factory=dict)   # numerator -> 'UNI'/'NOR'/'LGU'
    value: dict = field(default_factory=dict)        # (em, num, nature) -> complex/float
    addvalue: dict = field(default_factory=dict)     # (em, num, nature) -> float

    def get(self, em, numerator, nature):
        return self.value[(em, numerator, nature)]

    def get_addvalue(self, em, numerator, nature):
        return self.addvalue.get((em, numerator, nature), np.nan)

    def has(self, em, numerator, nature):
        return (em, numerator, nature) in self.value

    def has_numerator(self, numerator):
        return numerator in self.disttype

    def all_have(self, ems, numerator):
        """True if every end-member carries a distribution for ``numerator``
        (port of the EMHaveSO4 / EMHaved34S checks in FindScenarioParameters)."""
        natures = ("Min", "Max", "Men", "Std", "Mod")
        return all(any(self.has(e, numerator, n) for n in natures) for e in ems)


def _is_number(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def read_endmembers(table: EndMemberTable) -> EMData:
    """Parse an :class:`EndMemberTable` into structured :class:`EMData`."""
    em_codes = list(table.em_codes)
    data = EMData(group=table.group, em_codes=em_codes)

    # Build dense (row, em) value/addvalue arrays so we can do the NaN->0 fill.
    n_rows = len(table.rows)
    n_em = len(em_codes)
    raw = np.full((n_rows, n_em), np.nan, dtype=complex)
    addv = np.full((n_rows, n_em), np.nan, dtype=float)
    row_meta = []  # (numerator, dist3, nature3)

    for ri, row in enumerate(table.rows):
        mname = row.matlab_name
        dist3 = mname[0:3]
        nature3 = mname[4:7]
        row_meta.append((row.numerator, dist3, nature3))
        for ci, em in enumerate(em_codes):
            cell = row.values.get(em)
            val, av = _decode_cell(cell, nature3)
            raw[ri, ci] = val
            addv[ri, ci] = av

    # NaN->0 fill: for each end-member column that has at least one non-NaN
    # entry, set its remaining NaN entries to 0 (matches MATLAB).
    for ci in range(n_em):
        col = raw[:, ci]
        if np.sum(~np.isnan(col)) > 0:
            col[np.isnan(col)] = 0
            raw[:, ci] = col

    # Populate the structured dicts.
    for ri, (numerator, dist3, nature3) in enumerate(row_meta):
        data.disttype[numerator] = dist3
        for ci, em in enumerate(em_codes):
            v = raw[ri, ci]
            # collapse exact-real complex values back to float for cleanliness
            if v.imag == 0:
                v = v.real
            data.value[(em, numerator, nature3)] = v
            data.addvalue[(em, numerator, nature3)] = addv[ri, ci]
    return data


def _decode_cell(cell, nature3):
    """Return (value, addvalue) for one spreadsheet cell, applying sentinels."""
    if cell is None:
        return (complex(np.nan), np.nan)
    if _is_number(cell):
        return (complex(float(cell)), np.nan)

    s = str(cell)
    if s.count("#") == 1:
        pre, post = s.split("#")
        try:
            offset = float(post)
        except ValueError:
            offset = 0.0
        if pre in _SAMPLE_WORDS:
            sentinel = {"Min": SAMPLE_MIN, "Max": SAMPLE_MAX,
                        "Men": SAMPLE_MEN}.get(nature3)
        elif pre in _FRAC_WORDS:
            sentinel = {"Min": FRAC_MIN, "Max": FRAC_MAX,
                        "Men": FRAC_MEN}.get(nature3)
        else:
            # unknown keyword -> treated as 0 (MATLAB warns and uses 0)
            return (complex(0.0), 0.0)
        if sentinel is None:  # e.g. Sig defined relative to sample (invalid)
            return (complex(0.0), 0.0)
        return (complex(sentinel), offset)

    # string without a single '#': MATLAB warns and uses 0
    return (complex(0.0), 0.0)
