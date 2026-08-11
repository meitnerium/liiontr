from __future__ import annotations

import math
from dataclasses import dataclass

from .model import KineticModel


R = 8.314462618


@dataclass(slots=True)
class Arrhenius(KineticModel):
    """
    Arrhenius kinetic model.

    rate = A * exp(-Ea / (R * T))
    """

    activation_energy: float
    pre_exponential_factor: float

    def rate(
        self,
        temperature: float,
    ) -> float:
        """
        Return the reaction rate at the given temperature.
        """

        return self.pre_exponential_factor * math.exp(
            -self.activation_energy / (R * temperature)
        )
