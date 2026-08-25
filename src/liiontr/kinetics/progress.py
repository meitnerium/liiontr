"""Reaction-progress models for thermal runaway kinetics."""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from liiontr.reactions.context import ReactionContext


class ProgressModel(ABC):
    """Base class for reaction progress models."""

    @abstractmethod
    def factor(
        self,
        conversion: float,
        context: ReactionContext | None = None,
    ) -> float:
        """Return the conversion-dependent reaction factor."""
        raise NotImplementedError


@dataclass(slots=True)
class PowerLawProgress(ProgressModel):
    """
    Power-law reaction progress model.

    f(alpha) = (1 - alpha)^n
    """

    order: float = 1.0

    def __post_init__(self) -> None:
        """Validate the reaction-progress order."""
        if self.order <= 0.0:
            raise ValueError(
                "Progress model order must be greater than zero."
            )

    def factor(
        self,
        conversion: float,
        context: ReactionContext | None = None,
    ) -> float:
        """
        Return the reaction progress factor.

        The reaction context is accepted for interface compatibility
        but is not used by this model.
        """
        del context

        if conversion >= 1.0:
            return 0.0

        if conversion < 0.0:
            conversion = 0.0

        return (1.0 - conversion) ** self.order


@dataclass(slots=True)
class AutocatalyticProgress(ProgressModel):
    """
    Autocatalytic reaction progress model.

    f(alpha) = alpha^m * (1 - alpha)^n
    """

    autocatalytic_order: float = 1.0
    remaining_order: float = 1.0

    def __post_init__(self) -> None:
        """Validate the autocatalytic reaction orders."""
        if (
                self.autocatalytic_order <= 0.0
                or self.remaining_order <= 0.0
        ):
            raise ValueError(
                "Autocatalytic progress orders must be greater than zero."
            )

    def factor(
        self,
        conversion: float,
        context: ReactionContext | None = None,
    ) -> float:
        """Return the reaction progress factor.

        The reaction context is accepted for interface compatibility
        but is not used by this model.
        """
        del context

        if conversion <= 0.0:
            return 0.0

        if conversion >= 1.0:
            return 0.0

        return (
            conversion**self.autocatalytic_order
            * (1.0 - conversion) ** self.remaining_order
        )


@dataclass(slots=True)
class ThresholdProgress(ProgressModel):
    """
    Progress model activated by the state of another reaction.

    The wrapped progress model is active only when the remaining
    fraction of the specified reaction is below a threshold.
    """

    progress_model: ProgressModel
    reaction_name: str
    remaining_below: float

    def __post_init__(self) -> None:
        """Validate the remaining-fraction activation threshold."""
        if not 0.0 <= self.remaining_below <= 1.0:
            raise ValueError(
                "Remaining fraction threshold must be between 0 and 1."
            )

    def factor(
        self,
        conversion: float,
        context: ReactionContext | None = None,
    ) -> float:
        """Return the wrapped progress factor when the threshold is met."""
        if context is None:
            raise ValueError("Reaction context is required.")

        remaining_fraction = context.remaining_fraction(self.reaction_name)

        if remaining_fraction >= self.remaining_below:
            return 0.0

        return self.progress_model.factor(
            conversion=conversion,
            context=context,
        )


@dataclass(slots=True)
class ExponentialInhibitionProgress(ProgressModel):
    """
    Progress model with exponential inhibition.

    f(alpha, context) =
        f_base(alpha, context) * exp(-x)

    where x is obtained from a ReactionContext variable.
    """

    progress_model: ProgressModel
    variable_name: str

    def factor(
        self,
        conversion: float,
        context: ReactionContext | None = None,
    ) -> float:
        """Return the exponentially inhibited progress factor."""
        if context is None:
            raise ValueError("Reaction context is required.")

        inhibition_variable = context.variable(self.variable_name)

        return self.progress_model.factor(
            conversion=conversion,
            context=context,
        ) * math.exp(-inhibition_variable)
