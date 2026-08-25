"""Individual thermal runaway reaction models."""

from __future__ import annotations

from dataclasses import dataclass

from liiontr.kinetics.model import KineticModel
from liiontr.kinetics.progress import (
    PowerLawProgress,
    ProgressModel,
)

from .context import ReactionContext


@dataclass(slots=True)
class Reaction:
    """Represent a reaction contributing to thermal runaway."""

    name: str
    kinetics: KineticModel
    enthalpy: float
    mass_fraction: float = 1.0

    reaction_order: float = 1.0

    progress_model: ProgressModel | None = None

    def __post_init__(self) -> None:
        """Validate the reaction definition."""
        if self.reaction_order <= 0.0:
            raise ValueError("Reaction order must be greater than zero.")

        if self.progress_model is None:
            self.progress_model = PowerLawProgress(
                order=self.reaction_order,
            )

        elif self.reaction_order != 1.0:
            raise ValueError(
                "Specify either reaction_order or progress_model, not both."
            )

    def rate(
        self,
        temperature: float,
    ) -> float:
        """Return the temperature-dependent kinetic rate constant."""
        return self.kinetics.rate(temperature)

    def progress_rate(
        self,
        temperature: float,
        conversion: float,
        context: ReactionContext | None = None,
    ) -> float:
        r"""Return the reaction conversion rate.

        The reaction progress follows

        .. math::

        \frac{d\alpha}{dt}
        =
        k(T) f(\alpha, \mathrm{context}),

        where :math:`k(T)` is the temperature-dependent kinetic
        rate and :math:`f` is the reaction-progress model.
        """
        progress_model = self.progress_model

        if progress_model is None:
            raise RuntimeError("Reaction progress model is not initialized.")

        return self.rate(temperature) * progress_model.factor(
            conversion=conversion,
            context=context,
        )

    def heat_generation(
        self,
        temperature: float,
        conversion: float = 0.0,
        context: ReactionContext | None = None,
    ) -> float:
        """Return mass-specific heat generation rate in W/kg."""
        return (
            self.enthalpy
            * self.mass_fraction
            * self.progress_rate(
                temperature=temperature,
                conversion=conversion,
                context=context,
            )
        )
