"""Port of ``MEANDIR_Go.m``.

Given the list of observations in the inversion, return a mapping that flags
which variables are active. In MATLAB this is the ``go`` struct with a 0/1
field per variable; here it is a plain dict ``{name: bool}`` covering every
variable in :data:`meandir.constants.ALL_VARIABLES`.
"""

from __future__ import annotations

from .constants import ALL_VARIABLES


def go(obs_list) -> dict[str, bool]:
    """Return ``{variable: bool}`` indicating membership in ``obs_list``."""
    active = set(obs_list)
    return {name: (name in active) for name in ALL_VARIABLES}
