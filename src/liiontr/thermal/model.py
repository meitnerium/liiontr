"""Abstract interfaces for LiionTR thermal models."""

from __future__ import annotations

from abc import ABC, abstractmethod


class ThermalModel(ABC):
    """Abstract interface for battery thermal models."""

    @abstractmethod
    def temperature_derivative(
        self,
        temperature: float,
        heat_generation: float,
    ) -> float:
        """Return the cell temperature rate in K/s."""
        raise NotImplementedError
