from __future__ import annotations

from abc import ABC
from dataclasses import dataclass

from liiontr.geometry.geometry import Geometry
from liiontr.materials.material import Material
from liiontr.chemistry.chemistry import Chemistry


@dataclass(slots=True)
class Cell(ABC):
    """
    Base class for every battery cell.
    """

    name: str

    geometry: Geometry

    material: Material

    chemistry: Chemistry

    @property
    def volume(self) -> float:
        return self.geometry.volume

    @property
    def surface_area(self) -> float:
        return self.geometry.surface_area

    @property
    def mass(self) -> float:
        return self.material.density(298.15) * self.volume

    @property
    def thermal_capacity(self) -> float:
        return self.mass * self.material.heat_capacity(298.15)
