"""Fractionation pairing — port of ``MEANDIR_FindFractionationPairs.m``.

Identifies, for each end-member whose isotope entry is defined as a
fractionation (distribution code ``EPS-*``), the row/column positions of the
isotope, its carrier ion, and the end-member, in both the full inversion matrix
(``*0``) and the reduced one (``*_r``). All indices are 0-based.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..constants import isotope_ion

_EPS = {"EPS-NOR", "EPS-UNI", "EPS-LGU"}


@dataclass
class _Side:
    n: int = 0
    iso: list = field(default_factory=list)
    EM: list = field(default_factory=list)
    isopos0: list = field(default_factory=list)
    empos0: list = field(default_factory=list)
    isopos_r: list = field(default_factory=list)
    empos_r: list = field(default_factory=list)
    ionpos0: list = field(default_factory=list)
    ionpos_r: list = field(default_factory=list)


@dataclass
class Fractionation:
    red: _Side = field(default_factory=_Side)
    all: _Side = field(default_factory=_Side)


def _index_or_none(seq, value):
    return seq.index(value) if value in seq else None


def find_fractionation_pairs(distcode_full, distcode_red, obs, obs_r, ems,
                             ems_r, carbonmatch) -> Fractionation:
    f = Fractionation()

    # (1) reduced matrix (what is actually inverted)
    for ii in range(len(obs_r)):
        for jj in range(len(ems_r)):
            if distcode_red[ii][jj] not in _EPS:
                continue
            iso = obs_r[ii]
            ion = isotope_ion(iso, carbonmatch)
            f.red.n += 1
            f.red.iso.append(iso)
            f.red.EM.append(ems_r[jj])
            f.red.isopos_r.append(ii)
            f.red.empos_r.append(jj)
            f.red.isopos0.append(obs.index(iso))
            f.red.empos0.append(ems.index(ems_r[jj]))
            f.red.ionpos0.append(obs.index(ion))
            f.red.ionpos_r.append(obs_r.index(ion))

    # (2) full matrix (every fractionation instance)
    for ii in range(len(obs)):
        for jj in range(len(ems)):
            if distcode_full[ii][jj] not in _EPS:
                continue
            iso = obs[ii]
            ion = isotope_ion(iso, carbonmatch)
            f.all.n += 1
            f.all.iso.append(iso)
            f.all.EM.append(ems[jj])
            f.all.isopos0.append(ii)
            f.all.empos0.append(jj)
            f.all.isopos_r.append(_index_or_none(obs_r, iso))
            f.all.empos_r.append(_index_or_none(ems_r, ems[jj]))
            f.all.ionpos0.append(obs.index(ion))
            f.all.ionpos_r.append(_index_or_none(obs_r, ion))
    return f
