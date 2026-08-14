from __future__ import annotations

from dataclasses import dataclass, field

from .context import ReactionContext
from .reaction import Reaction


@dataclass(slots=True)
class ReactionNetwork:
    """
    Collection of chemical reactions.

    Reaction names must be unique because they are used as keys
    in the shared ReactionContext.
    """

    reactions: list[Reaction] = field(default_factory=list)

    def __post_init__(self) -> None:
        self._validate_unique_names()

    def _validate_unique_names(self) -> None:
        """
        Ensure that all reaction names are unique.
        """

        names = [
            reaction.name
            for reaction in self.reactions
        ]

        if len(names) != len(set(names)):
            raise ValueError(
                "Reaction names must be unique."
            )

    def add(
        self,
        reaction: Reaction,
    ) -> None:
        """
        Add a reaction to the network.
        """

        if any(
            existing.name == reaction.name
            for existing in self.reactions
        ):
            raise ValueError(
                "Reaction names must be unique."
            )

        self.reactions.append(
            reaction
        )

    def _validate_conversions(
        self,
        conversions: list[float],
    ) -> None:
        """
        Validate the number of reaction conversions.
        """

        if len(conversions) != len(self.reactions):
            raise ValueError(
                "Number of conversions must match number of reactions."
            )

    def context(
        self,
        conversions: list[float],
    ) -> ReactionContext:
        """
        Build the shared reaction context.

        Reaction conversions are stored by reaction name.
        """

        self._validate_conversions(
            conversions
        )

        return ReactionContext(
            conversions={
                reaction.name: float(conversion)
                for reaction, conversion in zip(
                    self.reactions,
                    conversions,
                )
            }
        )

    def progress_rates(
        self,
        temperature: float,
        conversions: list[float],
    ) -> list[float]:
        """
        Return the progress rate of each reaction.
        """

        context = self.context(
            conversions
        )

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
        """
        Return total mass-specific heat generation [W/kg].
        """

        if conversions is None:
            conversions = [
                0.0
                for _ in self.reactions
            ]

        context = self.context(
            conversions
        )

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