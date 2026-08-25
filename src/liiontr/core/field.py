"""Base field representations defined over computational domains."""

from __future__ import annotations

from abc import ABC
from typing import Any

import numpy as np


class Field(ABC):
    """Base class for simulation fields.

    A field stores numerical values as a NumPy array of floating-point
    values. Specialized scalar and vector fields inherit from this
    class.
    """

    def __init__(self, values: Any) -> None:
        """Initialize the field.

        Parameters
        ----------
        values : Any
            Values that can be converted to a NumPy floating-point
            array.
        """
        self.values = np.asarray(
            values,
            dtype=float,
        )

    @property
    def size(self) -> int:
        """Return the total number of values in the field."""
        return int(self.values.size)

    def copy(self) -> Field:
        """Return an independent copy of the field.

        Returns
        -------
        Field
            New field instance containing a copy of the numerical
            values.
        """
        return self.__class__(self.values.copy())
