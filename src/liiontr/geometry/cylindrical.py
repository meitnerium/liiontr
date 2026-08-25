"""Cylindrical battery cell geometry models."""

from __future__ import annotations

from dataclasses import dataclass
from math import pi

from .geometry import Geometry


@dataclass(slots=True)
class CylindricalGeometry(Geometry):
    """Represent the geometry of a cylindrical battery cell.

    Parameters
    ----------
    radius : float
        Cell radius in m.
    height : float
        Cell height in m.
    """

    radius: float
    height: float

    @property
    def volume(self) -> float:
        """Return the cylindrical cell volume in m³."""
        return pi * self.radius**2 * self.height

    @property
    def surface_area(self) -> float:
        """Return the total cylindrical cell surface area in m²."""
        return 2 * pi * self.radius * self.height + 2 * pi * self.radius**2
