"""Multi-channel thermal runaway reaction models."""

from __future__ import annotations

from dataclasses import dataclass

from liiontr.kinetics import KineticModel, ProgressModel

from .context import ReactionContext


@dataclass(slots=True, frozen=True)
class ReactionChannel:
    """
    One kinetic and thermal channel of a multichannel reaction.

    Each channel has its own kinetic model and reaction enthalpy,
    but all channels share the same reaction progress variable.
    """

    kinetics: KineticModel
    enthalpy: float


@dataclass(slots=True)
class MultiChannelReaction:
    """
    Chemical reaction composed of multiple parallel channels.

    All channels share one conversion alpha and one progress model.

    The reaction progress rate is:

        d(alpha)/dt =
            f(alpha, context) * sum(k_i(T))

    The mass-specific heat generation rate is:

        q_dot =
            mass_fraction
            * f(alpha, context)
            * sum(H_i * k_i(T))
    """

    name: str

    channels: list[ReactionChannel]

    mass_fraction: float

    progress_model: ProgressModel

    def __post_init__(self) -> None:
        """Validate the multi-channel reaction definition."""
        if not self.channels:
            raise ValueError("At least one reaction channel is required.")

        if not 0.0 <= self.mass_fraction <= 1.0:
            raise ValueError("Mass fraction must be between 0 and 1.")

    def progress_rate(
        self,
        temperature: float,
        conversion: float,
        context: ReactionContext | None = None,
    ) -> float:
        """Return the shared reaction progress rate."""
        progress_factor = self.progress_model.factor(
            conversion=conversion,
            context=context,
        )

        total_rate = sum(
            channel.kinetics.rate(temperature) for channel in self.channels
        )

        return total_rate * progress_factor

    def heat_generation(
        self,
        temperature: float,
        conversion: float = 0.0,
        context: ReactionContext | None = None,
    ) -> float:
        """Return mass-specific heat generation rate in W/kg."""
        progress_factor = self.progress_model.factor(
            conversion=conversion,
            context=context,
        )

        channel_heat = sum(
            channel.enthalpy * channel.kinetics.rate(temperature)
            for channel in self.channels
        )

        return self.mass_fraction * progress_factor * channel_heat
