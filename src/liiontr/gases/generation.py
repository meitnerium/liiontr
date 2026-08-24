from __future__ import annotations

from dataclasses import dataclass

from liiontr.reactions import ReactionNetwork

from .yield_model import ReactionGasYield


@dataclass(slots=True)
class GasGenerationModel:
    """
    Compute gas generation rates from a reaction network.

    Gas generation rates are returned in mol/s.
    """

    reaction_network: ReactionNetwork
    gas_yields: list[ReactionGasYield]

    @property
    def species_names(self) -> list[str]:
        """
        Return gas species names in deterministic state-vector order.

        Species are ordered by first appearance in the configured
        reaction gas yields.
        """

        names: list[str] = []

        for gas_yield in self.gas_yields:
            for species_name in gas_yield.species_yields:
                if species_name not in names:
                    names.append(
                        species_name
                    )

        return names

    def __post_init__(self) -> None:
        yield_names = [gas_yield.reaction_name for gas_yield in self.gas_yields]

        if len(yield_names) != len(set(yield_names)):
            raise ValueError("Duplicate gas yield for reaction.")

        reaction_names = {reaction.name for reaction in self.reaction_network.reactions}

        for reaction_name in yield_names:
            if reaction_name not in reaction_names:
                raise ValueError(f"Unknown reaction: {reaction_name}")

    def generation_rates(
        self,
        temperature: float,
        conversions: list[float],
        cell_mass: float,
    ) -> dict[str, float]:
        """
        Return total gas generation rates by species.

        Parameters
        ----------
        temperature:
            Cell temperature in K.

        conversions:
            Reaction conversions in network order.

        cell_mass:
            Total cell mass in kg.
        """

        if cell_mass <= 0.0:
            raise ValueError("Cell mass must be greater than zero.")

        progress_rates = self.reaction_network.progress_rates(
            temperature=temperature,
            conversions=conversions,
        )

        yields_by_reaction = {
            gas_yield.reaction_name: gas_yield for gas_yield in self.gas_yields
        }

        total_rates: dict[str, float] = {}

        for reaction, progress_rate in zip(
            self.reaction_network.reactions,
            progress_rates,
            strict=True,
        ):
            gas_yield = yields_by_reaction.get(reaction.name)

            if gas_yield is None:
                continue

            reaction_rates = gas_yield.generation_rates(
                cell_mass=cell_mass,
                reaction_mass_fraction=(reaction.mass_fraction),
                progress_rate=progress_rate,
            )

            for species_name, rate in reaction_rates.items():
                total_rates[species_name] = (
                    total_rates.get(
                        species_name,
                        0.0,
                    )
                    + rate
                )

        return total_rates
