from __future__ import annotations

from dataclasses import dataclass, field

from .reaction import Reaction


@dataclass(slots=True)
class ReactionNetwork:
    """
    Collection of chemical reactions.
    """

    reactions: list[Reaction] = field(default_factory=list)

    def add(
        self,
        reaction: Reaction,
    ) -> None:
        self.reactions.append(reaction)

    def progress_rates(
        self,
        temperature: float,
        conversions: list[float],
    ) -> list[float]:
        """
        Return the progress rate of each reaction.
        """

        if len(conversions) != len(self.reactions):
            raise ValueError("Number of conversions must match number of reactions.")

        return [
            reaction.progress_rate(
                temperature,
                conversion,
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
            conversions = [0.0] * len(self.reactions)

        if len(conversions) != len(self.reactions):
            raise ValueError("Number of conversions must match number of reactions.")

        return sum(
            reaction.heat_generation(
                temperature,
                conversion,
            )
            for reaction, conversion in zip(
                self.reactions,
                conversions,
            )
        )
