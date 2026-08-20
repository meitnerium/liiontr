from __future__ import annotations

from dataclasses import dataclass

from .model import KineticModel


@dataclass(slots=True)
class TemperatureThresholdKinetics(KineticModel):
    """
    Kinetic model activated above a minimum temperature.

    The wrapped kinetic model is inactive below the specified
    temperature threshold.
    """

    kinetics: KineticModel
    minimum_temperature: float

    def __post_init__(self) -> None:
        if self.minimum_temperature <= 0.0:
            raise ValueError("Minimum temperature must be greater than zero.")

    def rate(
        self,
        temperature: float,
    ) -> float:
        """
        Return the temperature-dependent kinetic rate.

        The rate is zero below the activation threshold.
        """

        if temperature < self.minimum_temperature:
            return 0.0

        return self.kinetics.rate(temperature)
