"""Derived variables evaluated from shared reaction contexts."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

from .context import ReactionContext


class ContextVariable(ABC):
    """
    Base class for variables derived from a ReactionContext.

    Context variables are algebraic quantities computed from the
    current reaction state. They are not independent ODE states.
    """

    @abstractmethod
    def evaluate(
            self,
            context: ReactionContext,
    ) -> float:
        """Evaluate the variable from the current reaction context."""
        raise NotImplementedError


@dataclass(slots=True, frozen=True)
class LinearConversionVariable(ContextVariable):
    """
    Variable derived linearly from a reaction conversion.

    The relation is:

        x = x_ref + slope * (alpha - alpha_ref)
    """

    reaction_name: str

    reference_conversion: float
    reference_value: float

    slope: float

    def __post_init__(self) -> None:
        """Validate the reference reaction conversion."""
        if not 0.0 <= self.reference_conversion <= 1.0:
            raise ValueError(
                "Reference conversion must be between 0 and 1."
            )

    def evaluate(
            self,
            context: ReactionContext,
    ) -> float:
        """Evaluate the derived conversion variable."""
        conversion = context.conversion(
            self.reaction_name
        )

        return self.reference_value + self.slope * (
            conversion - self.reference_conversion
        )


@dataclass(slots=True, frozen=True)
class RemainingFractionRatioVariable(ContextVariable):
    """
    Ratio of the current remaining fraction to a reference value.

    The relation is:

        x = (1 - alpha) / c_ref

    where alpha is the reaction conversion and c_ref is the
    reference remaining fraction.
    """

    reaction_name: str

    reference_remaining_fraction: float

    def __post_init__(self) -> None:
        """Validate the reference remaining fraction."""
        if not (
                0.0
                < self.reference_remaining_fraction
                <= 1.0
        ):
            raise ValueError(
                "Reference remaining fraction must be "
                "greater than zero and at most one."
            )

    def evaluate(
            self,
            context: ReactionContext,
    ) -> float:
        """Evaluate the remaining-fraction ratio."""
        remaining_fraction = (
            context.remaining_fraction(
                self.reaction_name
            )
        )

        return remaining_fraction / self.reference_remaining_fraction
