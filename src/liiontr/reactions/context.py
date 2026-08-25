"""Shared reaction state and context utilities."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(slots=True)
class ReactionContext:
    """
    Shared state available to coupled reactions.

    Conversions are stored by reaction name.

    Additional variables can be used for quantities such as
    SEI thickness, pressure, gas concentration, or other
    reaction-dependent state variables.
    """

    conversions: dict[str, float] = field(default_factory=dict)

    variables: dict[str, float] = field(default_factory=dict)

    def conversion(
        self,
        reaction_name: str,
    ) -> float:
        """Return the conversion of a named reaction."""
        if reaction_name not in self.conversions:
            raise KeyError(f"Unknown reaction conversion: {reaction_name}")

        return self.conversions[reaction_name]

    def remaining_fraction(
        self,
        reaction_name: str,
    ) -> float:
        """Return the remaining fraction of a named reaction."""
        return 1.0 - self.conversion(reaction_name)

    def variable(
        self,
        name: str,
    ) -> float:
        """Return an additional reaction state variable."""
        if name not in self.variables:
            raise KeyError(f"Unknown reaction variable: {name}")

        return self.variables[name]
