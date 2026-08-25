"""Named simulation variable representations."""

from __future__ import annotations

from dataclasses import dataclass

from .field import Field


@dataclass(slots=True)
class Variable:
    """
    Represents a physical simulation variable.

    Examples:
        Temperature
        Pressure
        Species concentration
    """

    name: str

    unit: str

    field: Field | None = None

    description: str = ""

    @property
    def values(self):
        if self.field is None:
            return None

        return self.field.values
