from __future__ import annotations

from dataclasses import dataclass, field

from .reaction import Reaction


@dataclass(slots=True)
class ReactionNetwork:
    """
    Collection of thermal runaway reactions.
    """

    reactions: list[Reaction] = field(default_factory=list)

    def add(
        self,
        reaction: Reaction,
    ) -> None:
        self.reactions.append(reaction)

    def heat_generation(
        self,
        temperature: float,
    ) -> float:
        return sum(reaction.heat_generation(temperature) for reaction in self.reactions)
