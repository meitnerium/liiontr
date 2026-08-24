from __future__ import annotations

from typing import TYPE_CHECKING, Protocol

if TYPE_CHECKING:
    from .context import ReactionContext


class ReactionModel(Protocol):
    """
    Interface implemented by reaction models used in a ReactionNetwork.

    A reaction model must expose a name and provide both a progress
    rate and a mass-specific heat generation rate.
    """

    name: str
    mass_fraction: float

    def progress_rate(
        self,
        temperature: float,
        conversion: float,
        context: ReactionContext | None = None,
    ) -> float:
        """
        Return the reaction progress rate.
        """

        ...

    def heat_generation(
        self,
        temperature: float,
        conversion: float = 0.0,
        context: ReactionContext | None = None,
    ) -> float:
        """
        Return mass-specific heat generation rate [W/kg].
        """

        ...
