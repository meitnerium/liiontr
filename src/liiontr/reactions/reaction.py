from __future__ import annotations

from dataclasses import dataclass

from liiontr.kinetics.model import KineticModel


@dataclass(slots=True)
class Reaction:
    """
    Chemical reaction contributing to thermal runaway.
    """

    name: str
    kinetics: KineticModel
    enthalpy: float
    mass_fraction: float = 1.0

    def rate(
        self,
        temperature: float,
    ) -> float:
        """
        Return the reaction rate.
        """

        return self.kinetics.rate(temperature)

    def heat_generation(
        self,
        temperature: float,
    ) -> float:
        """
        Return the heat generation rate.
        """

        return self.enthalpy * self.mass_fraction * self.rate(temperature)
