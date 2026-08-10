from __future__ import annotations

from dataclasses import dataclass
from math import pi

from .geometry import Geometry


@dataclass(slots=True)
class CylindricalGeometry(Geometry):
    radius: float
    height: float

    @property
    def volume(self) -> float:
        return pi * self.radius**2 * self.height

    @property
    def surface_area(self) -> float:
        return 2 * pi * self.radius * self.height + 2 * pi * self.radius**2
