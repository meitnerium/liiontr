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
    reaction_order: float = 1.0

    def __post_init__(self) -> None:
        if self.reaction_order <= 0.0:
            raise ValueError("Reaction order must be greater than zero.")

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

        The reaction follows:

            d(alpha)/dt = k(T) * (1 - alpha)^n

        where alpha is the conversion and n is the
        reaction order.
        """

        if conversion >= 1.0:
            return 0.0

        if conversion < 0.0:
            conversion = 0.0

        remaining_fraction = 1.0 - conversion

        return self.rate(temperature) * remaining_fraction**self.reaction_order

    def heat_generation(
        self,
        temperature: float,
        conversion: float = 0.0,
    ) -> float:
        """
        Return mass-specific heat generation rate [W/kg].
        """

        return (
            self.enthalpy
            * self.mass_fraction
            * self.progress_rate(
                temperature,
                conversion,
            )
        )
