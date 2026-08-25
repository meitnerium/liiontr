"""Reaction-network orchestration for coupled thermal reactions."""

from __future__ import annotations

from dataclasses import dataclass, field

from .context import ReactionContext
from .context_variable import ContextVariable
from .model import ReactionModel


@dataclass(slots=True)
class ReactionNetwork:
    """
    Collection of chemical reaction models.

    Reaction names must be unique because they are used as keys
    in the shared ReactionContext.

    Both simple reactions and multichannel reactions can be used
    as long as they implement the ReactionModel protocol.
    """

    reactions: list[ReactionModel] = field(default_factory=list)

    context_variables: dict[str, ContextVariable] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate the initial reaction-network definition."""
        self._validate_unique_names()

    def _validate_unique_names(self) -> None:
        """Ensure that all reaction names are unique."""
        names = [reaction.name for reaction in self.reactions]

        if len(names) != len(set(names)):
            raise ValueError("Reaction names must be unique.")

    def add(
        self,
        reaction: ReactionModel,
    ) -> None:
        """Add a reaction model to the network."""
        if any(existing.name == reaction.name for existing in self.reactions):
            raise ValueError("Reaction names must be unique.")

        self.reactions.append(reaction)

    def _validate_conversions(
        self,
        conversions: list[float],
    ) -> None:
        """Validate the number of reaction conversions."""
        if len(conversions) != len(self.reactions):
            raise ValueError("Number of conversions must match number of reactions.")

    def context(
        self,
        conversions: list[float],
    ) -> ReactionContext:
        """Build the shared reaction context.

        The context is first populated with reaction conversions.
        Derived context variables are then evaluated from the
        current reaction state.
        """
        self._validate_conversions(conversions)

        context = ReactionContext(
            conversions={
                reaction.name: float(conversion)
                for reaction, conversion in zip(
                    self.reactions,
                    conversions,
                )
            }
        )

        for name, variable in self.context_variables.items():
            context.variables[name] = variable.evaluate(context)

        return context

    def progress_rates(
        self,
        temperature: float,
        conversions: list[float],
    ) -> list[float]:
        """Return the progress rate of each reaction."""
        context = self.context(conversions)

        return [
            reaction.progress_rate(
                temperature=temperature,
                conversion=conversion,
                context=context,
            )
            for reaction, conversion in zip(
                self.reactions,
                conversions,
            )
        ]

    def heat_generation(
        self,
        temperature: float,
        conversions: list[float] | None = None,
    ) -> float:
        """Return total mass-specific heat generation in W/kg."""
        if conversions is None:
            conversions = [0.0 for _ in self.reactions]

        context = self.context(conversions)

        return sum(
            reaction.heat_generation(
                temperature=temperature,
                conversion=conversion,
                context=context,
            )
            for reaction, conversion in zip(
                self.reactions,
                conversions,
            )
        )

    def heat_generation_by_reaction(
        self,
        temperature: float,
        conversions: list[float] | None = None,
    ) -> dict[str, float]:
        """Return heat generation for each reaction.

        Values are expressed in W/kg of total cell mass.
        """
        if conversions is None:
            conversions = [0.0 for _ in self.reactions]

        context = self.context(
            conversions=conversions,
        )

        return {
            reaction.name: reaction.heat_generation(
                temperature=temperature,
                conversion=conversion,
                context=context,
            )
            for reaction, conversion in zip(
                self.reactions,
                conversions,
                strict=True,
            )
        }
