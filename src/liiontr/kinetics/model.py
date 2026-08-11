from __future__ import annotations

from abc import ABC, abstractmethod


class KineticModel(ABC):
    """
    Abstract kinetic model.
    """

    @abstractmethod
    def rate(self, temperature: float) -> float:
        """
        Return the reaction rate at the given temperature.
        """
        ...
