"""Abstract interfaces for battery geometry models."""

from __future__ import annotations

from abc import ABC, abstractmethod


class Geometry(ABC):
    """Base class for battery geometries."""

    @property
    @abstractmethod
    def volume(self) -> float:
        """Return the cell volume in m³."""

    @property
    @abstractmethod
    def surface_area(self) -> float:
        """Return the external surface area in m²."""
