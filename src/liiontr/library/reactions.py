from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from liiontr.kinetics import Arrhenius
from liiontr.reactions import Reaction

if TYPE_CHECKING:
    from liiontr.cells.cell import Cell


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


@dataclass(slots=True, frozen=True)
class VolumetricReactionParameters:
    """
    Reaction parameters reported using reactant content per cell volume.

    This representation is useful for literature models where the
    amount of reacting material is given in kg/m3 rather than directly
    as a fraction of the total cell mass.

    Conventions
    -----------
    activation_energy:
        Activation energy [J/mol].

    pre_exponential_factor:
        Arrhenius pre-exponential factor [1/s].

    enthalpy:
        Specific heat released by the reacting material [J/kg].

    specific_content:
        Reacting-material content per unit cell volume [kg/m3].

    initial_remaining_fraction:
        Initial fraction of reacting material remaining.

    reaction_order:
        Order of the power-law progress model.
    """

    name: str

    activation_energy: float
    pre_exponential_factor: float

    enthalpy: float

    specific_content: float

    initial_remaining_fraction: float = 1.0
    reaction_order: float = 1.0

    reference: str | None = None

    def __post_init__(self) -> None:
        if self.activation_energy <= 0.0:
            raise ValueError("Activation energy must be greater than zero.")

        if self.pre_exponential_factor <= 0.0:
            raise ValueError("Pre-exponential factor must be greater than zero.")

        if self.specific_content <= 0.0:
            raise ValueError("Specific content must be greater than zero.")

        if not 0.0 <= self.initial_remaining_fraction <= 1.0:
            raise ValueError("Initial remaining fraction must be between 0 and 1.")

        if self.reaction_order <= 0.0:
            raise ValueError("Reaction order must be greater than zero.")

    @property
    def initial_conversion(self) -> float:
        """
        Return the initial conversion alpha.

        alpha = 1 - remaining_fraction
        """

        return 1.0 - self.initial_remaining_fraction

    def mass_fraction(
        self,
        cell: Cell,
    ) -> float:
        """
        Convert volumetric reactant content to cell mass fraction.
        """

        cell_density = cell.material.density(298.15)

        mass_fraction = self.specific_content / cell_density

        if mass_fraction > 1.0:
            raise ValueError("Volumetric reactant content exceeds total cell mass.")

        return mass_fraction

    def build(
        self,
        cell: Cell,
    ) -> Reaction:
        """
        Construct a Reaction for the specified cell.
        """

        return Reaction(
            name=self.name,
            kinetics=Arrhenius(
                activation_energy=self.activation_energy,
                pre_exponential_factor=self.pre_exponential_factor,
            ),
            enthalpy=self.enthalpy,
            mass_fraction=self.mass_fraction(cell),
            reaction_order=self.reaction_order,
        )
