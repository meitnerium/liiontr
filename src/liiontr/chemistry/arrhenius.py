from __future__ import annotations

import math
from dataclasses import dataclass

from .backend import ChemistryBackend


R = 8.314462618


@dataclass(slots=True)
class ArrheniusBackend(ChemistryBackend):
    activation_energy: float
    pre_exponential_factor: float
    enthalpy: float

    def heat_generation(
        self,
        temperature: float,
    ) -> float:
        rate = self.pre_exponential_factor * math.exp(
            -self.activation_energy / (R * temperature)
        )

        return self.enthalpy * rate
