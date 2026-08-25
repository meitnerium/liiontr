"""Named simulation variable representations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .field import Field


@dataclass(slots=True)
class Variable:
    """Represent a named physical simulation variable.

    Parameters
    ----------
    name : str
        Name of the physical variable.
    unit : str
        Physical unit associated with the variable.
    field : Field, optional
        Numerical field containing the variable values.
    description : str, optional
        Human-readable description of the variable.

    Examples
    --------
    Typical variables include temperature, pressure, species
    concentration, and reaction progress.
    """

    name: str

    unit: str

    field: Field | None = None

    description: str = ""

    @property
    def values(self) -> np.ndarray | None:
        """Return the numerical values associated with the variable.

        Returns
        -------
        numpy.ndarray or None
            Numerical field values, or ``None`` if no field has been
            assigned.
        """
        if self.field is None:
            return None

        return self.field.values
