"""Empirical gas yields associated with thermal reactions."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True, frozen=True)
class ReactionGasYield:
    """
    Gas yields associated with a reaction.

    Parameters
    ----------
    reaction_name:
        Name of the associated reaction.

    species_yields:
        Gas yield for each species in mol/kg of reacted material.
    """

    reaction_name: str
    species_yields: dict[str, float]

    def __post_init__(self) -> None:
        """Validate the reaction gas-yield definition."""
        if not self.reaction_name:
            raise ValueError("Reaction name must not be empty.")

        for species_name, gas_yield in self.species_yields.items():
            if not species_name:
                raise ValueError("Gas species name must not be empty.")

            if gas_yield < 0.0:
                raise ValueError("Gas yield must not be negative.")

    def yield_of(
        self,
        species_name: str,
    ) -> float:
        """Return a species yield in mol/kg reacted material."""
        return self.species_yields.get(
            species_name,
            0.0,
        )

    def generation_rates(
        self,
        cell_mass: float,
        reaction_mass_fraction: float,
        progress_rate: float,
    ) -> dict[str, float]:
        r"""Return gas generation rates in mol/s.

        The reacted material rate is:

            dm_reacted/dt =
                m_cell * w_reaction * d(alpha)/dt
        """
        if cell_mass <= 0.0:
            raise ValueError("Cell mass must be greater than zero.")

        if not 0.0 <= reaction_mass_fraction <= 1.0:
            raise ValueError("Reaction mass fraction must be between 0 and 1.")

        if progress_rate < 0.0:
            raise ValueError("Progress rate must not be negative.")

        reacted_mass_rate = cell_mass * reaction_mass_fraction * progress_rate

        return {
            species_name: (gas_yield * reacted_mass_rate)
            for species_name, gas_yield in self.species_yields.items()
        }
