from __future__ import annotations

from abc import ABC, abstractmethod


class ThermalModel(ABC):
    """
    Base class for thermal models.
    """

    @abstractmethod
    def temperature_derivative(
        self,
        temperature: float,
        heat_generation: float,
    ) -> float:
        """
        Compute dT/dt.
        """
        raise NotImplementedError
