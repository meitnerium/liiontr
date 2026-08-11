from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from liiontr.reactions import ReactionNetwork

from .backend import ChemistryBackend

if TYPE_CHECKING:
    from liiontr.cells.cell import Cell


@dataclass(slots=True)
class ReactionNetworkBackend(ChemistryBackend):
    """
    Chemistry backend based on a thermal reaction network.

    The reaction network returns heat generation per unit mass [W/kg].
    This backend converts it to total cell heat generation [W].
    """

    reaction_network: ReactionNetwork
    cell: Cell

    def heat_generation(
        self,
        temperature: float,
    ) -> float:
        heat_generation_per_mass = self.reaction_network.heat_generation(
            temperature
        )

        return heat_generation_per_mass * self.cell.mass