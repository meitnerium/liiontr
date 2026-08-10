from __future__ import annotations

from dataclasses import dataclass, field

from .reaction import Reaction
from .arrhenius import heat_generation


@dataclass(slots=True)
class ReactionNetwork:
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
        return sum(
            heat_generation(
                reaction,
                temperature,
            )
            for reaction in self.reactions
        )
