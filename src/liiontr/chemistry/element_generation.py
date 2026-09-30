"""Element generation from thermal runaway reactions."""

from __future__ import annotations

from dataclasses import dataclass

from liiontr.chemistry.elements import (
    ElementInventory,
    ReactionElementYield,
)
from liiontr.reactions import ReactionNetwork


@dataclass(slots=True)
class ElementGenerationModel:
    """Generate elemental inventories from reaction progress."""

    reaction_network: ReactionNetwork
    element_yields: list[ReactionElementYield]

    @property
    def element_names(self) -> list[str]:
        """Return element names in deterministic order."""
        names: list[str] = []

        for element_yield in self.element_yields:
            for element_name in element_yield.element_names:
                if element_name not in names:
                    names.append(element_name)

        return names

    def __post_init__(self) -> None:
        reaction_names = [
            reaction.name
            for reaction in self.reaction_network.reactions
        ]

        yield_reaction_names = [
            element_yield.reaction_name
            for element_yield in self.element_yields
        ]

        if len(
            yield_reaction_names
        ) != len(
            set(yield_reaction_names)
        ):
            raise ValueError(
                "Element yields must have unique reaction names."
            )

        unknown_reactions = [
            reaction_name
            for reaction_name in yield_reaction_names
            if reaction_name not in reaction_names
        ]

        if unknown_reactions:
            names = ", ".join(
                unknown_reactions
            )

            raise ValueError(
                "Element yield refers to unknown reaction: "
                f"{names}"
            )

    def generation_rates(
        self,
        temperature: float,
        conversions: list[float],
        cell_mass: float,
    ) -> dict[str, float]:
        """Return elemental generation rates in mol-atoms/s."""
        if cell_mass <= 0.0:
            raise ValueError(
                "Cell mass must be greater than zero."
            )

        progress_rates = (
            self.reaction_network.progress_rates(
                temperature=temperature,
                conversions=conversions,
            )
        )

        reaction_indices = {
            reaction.name: index
            for index, reaction in enumerate(
                self.reaction_network.reactions
            )
        }

        rates = {
            element_name: 0.0
            for element_name in self.element_names
        }

        for element_yield in self.element_yields:
            reaction_index = reaction_indices[
                element_yield.reaction_name
            ]

            reaction = (
                self.reaction_network.reactions[
                    reaction_index
                ]
            )

            progress_rate = (
                progress_rates[
                    reaction_index
                ]
            )

            reacted_mass_rate = (
                cell_mass
                * reaction.mass_fraction
                * progress_rate
            )

            for (
                element_name,
                yield_value,
            ) in element_yield.element_yields.items():
                rates[element_name] += (
                    yield_value
                    * reacted_mass_rate
                )

        return rates

    def generated_moles(
        self,
        initial_conversions: list[float],
        conversions: list[float],
        cell_mass: float,
    ) -> dict[str, float]:
        """Return cumulative generated element amounts in mol-atoms."""
        if cell_mass <= 0.0:
            raise ValueError(
                "Cell mass must be greater than zero."
            )

        reaction_count = len(
            self.reaction_network.reactions
        )

        if len(initial_conversions) != reaction_count:
            raise ValueError(
                "Number of initial conversions must match "
                "number of reactions."
            )

        if len(conversions) != reaction_count:
            raise ValueError(
                "Number of conversions must match "
                "number of reactions."
            )

        if any(
            conversion < 0.0
            or conversion > 1.0
            for conversion in initial_conversions
        ):
            raise ValueError(
                "Initial conversions must be between 0 and 1."
            )

        if any(
            conversion < 0.0
            or conversion > 1.0
            for conversion in conversions
        ):
            raise ValueError(
                "Conversions must be between 0 and 1."
            )

        reaction_indices = {
            reaction.name: index
            for index, reaction in enumerate(
                self.reaction_network.reactions
            )
        }

        amounts = {
            element_name: 0.0
            for element_name in self.element_names
        }

        for element_yield in self.element_yields:
            reaction_index = reaction_indices[
                element_yield.reaction_name
            ]

            reaction = (
                self.reaction_network.reactions[
                    reaction_index
                ]
            )

            conversion_change = (
                conversions[reaction_index]
                - initial_conversions[
                    reaction_index
                ]
            )

            if conversion_change < 0.0:
                raise ValueError(
                    "Current conversion must not be less "
                    "than initial conversion."
                )

            reacted_mass = (
                cell_mass
                * reaction.mass_fraction
                * conversion_change
            )

            for (
                element_name,
                yield_value,
            ) in element_yield.element_yields.items():
                amounts[element_name] += (
                    yield_value
                    * reacted_mass
                )

        return amounts

    def generated_inventory(
        self,
        initial_conversions: list[float],
        conversions: list[float],
        cell_mass: float,
    ) -> ElementInventory:
        """Return cumulative generated elements as an inventory."""
        return ElementInventory(
            moles=self.generated_moles(
                initial_conversions=initial_conversions,
                conversions=conversions,
                cell_mass=cell_mass,
            )
        )