"""Abstract interfaces for temperature-dependent kinetic models."""

from __future__ import annotations

from abc import ABC, abstractmethod


class KineticModel(ABC):
    """Abstract interface for temperature-dependent kinetic models."""

    @abstractmethod
    def rate(self, temperature: float) -> float:
        """Return the reaction rate at the given temperature."""
        ...
