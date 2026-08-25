"""Chemistry backend based on thermal reaction networks."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from liiontr.reactions import ReactionNetwork

from .backend import ChemistryBackend

if TYPE_CHECKING:
    from liiontr.cells.cell import Cell


@dataclass(slots=True)
class ReactionNetworkBackend(ChemistryBackend):
    r"""Couple a thermal reaction network to a battery cell.

    Parameters
    ----------
    reaction_network : ReactionNetwork
        Thermal reaction network used to compute reaction progress and
        mass-specific heat generation.
    cell : Cell
        Battery cell whose mass is used to convert mass-specific heat
        generation into total heat generation.

    Notes
    -----
    The reaction network returns heat generation in W/kg. This backend
    converts it to total cell heat generation in W according to

    .. math::

        \dot{Q} =
        \dot{q}_{m} m_{\mathrm{cell}},

    where :math:`\dot{q}_{m}` is the mass-specific heat-generation rate
    and :math:`m_{\mathrm{cell}}` is the cell mass.
    """

    reaction_network: ReactionNetwork
    cell: Cell

    def heat_generation(
        self,
        temperature: float,
        conversions: list[float] | None = None,
    ) -> float:
        """Return the total heat generation rate of the cell.

        Parameters
        ----------
        temperature : float
            Cell temperature in K.
        conversions : list[float] or None, optional
            Reaction conversion values for the reaction network.

        Returns
        -------
        float
            Total cell heat generation rate in W.
        """
        heat_generation_per_mass = self.reaction_network.heat_generation(
            temperature=temperature,
            conversions=conversions,
        )

        return heat_generation_per_mass * self.cell.mass

    def progress_rates(
        self,
        temperature: float,
        conversions: list[float],
    ) -> list[float]:
        """Return reaction conversion rates.

        Parameters
        ----------
        temperature : float
            Cell temperature in K.
        conversions : list[float]
            Current dimensionless reaction conversion values.

        Returns
        -------
        list[float]
            Reaction conversion rates in s⁻¹.
        """
        return self.reaction_network.progress_rates(
            temperature=temperature,
            conversions=conversions,
        )
