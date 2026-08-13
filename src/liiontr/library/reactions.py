from __future__ import annotations

from dataclasses import dataclass

from liiontr.kinetics import Arrhenius
from liiontr.reactions import Reaction


@dataclass(slots=True, frozen=True)
class ReactionParameters:
    """
    Parameter set used to construct a thermal runaway reaction.

    Conventions
    -----------
    activation_energy:
        Activation energy [J/mol].

    pre_exponential_factor:
        Arrhenius pre-exponential factor [1/s].

    enthalpy:
        Specific heat released by the reaction [J/kg].

    mass_fraction:
        Fraction of the total cell mass participating in the reaction.

    reaction_order:
        Order of the power-law progress model.

    reference:
        Literature or calibration source for the parameter set.
    """

    name: str

    activation_energy: float
    pre_exponential_factor: float

    enthalpy: float

    mass_fraction: float = 1.0
    reaction_order: float = 1.0

    reference: str | None = None

    def __post_init__(self) -> None:
        if self.activation_energy <= 0.0:
            raise ValueError("Activation energy must be greater than zero.")

        if self.pre_exponential_factor <= 0.0:
            raise ValueError("Pre-exponential factor must be greater than zero.")

        if not 0.0 <= self.mass_fraction <= 1.0:
            raise ValueError("Mass fraction must be between 0 and 1.")

        if self.reaction_order <= 0.0:
            raise ValueError("Reaction order must be greater than zero.")

    def build(self) -> Reaction:
        """
        Construct a Reaction from this parameter set.
        """

        return Reaction(
            name=self.name,
            kinetics=Arrhenius(
                activation_energy=self.activation_energy,
                pre_exponential_factor=self.pre_exponential_factor,
            ),
            enthalpy=self.enthalpy,
            mass_fraction=self.mass_fraction,
            reaction_order=self.reaction_order,
        )
