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
        Return the kinetic rate constant.
        """

        return self.kinetics.rate(temperature)

    def progress_rate(
        self,
        temperature: float,
        conversion: float,
    ) -> float:
        """
        Return the reaction progress rate.

        The reaction stops once conversion reaches 1.
        """

        if conversion >= 1.0:
            return 0.0

        if conversion < 0.0:
            conversion = 0.0

        return self.rate(temperature) * (1.0 - conversion)

    def heat_generation(
        self,
        temperature: float,
        conversion: float = 0.0,
    ) -> float:
        """
        Return mass-specific heat generation rate [W/kg].

        The heat release decreases as the reaction is consumed.
        """

        return (
            self.enthalpy
            * self.mass_fraction
            * self.progress_rate(
                temperature,
                conversion,
            )
        )
