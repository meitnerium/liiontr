from __future__ import annotations

from dataclasses import dataclass


GAS_CONSTANT = 8.314462618


@dataclass(slots=True, frozen=True)
class IdealGasPressureModel:
    """
    Ideal-gas pressure model for the free internal volume of a cell.

    The model accounts for:
    - gas initially present in the free volume;
    - thermal expansion of that initial gas;
    - additional gas generated during thermal runaway.

    Parameters
    ----------
    free_volume:
        Internal free gas volume in m3.

    initial_pressure:
        Initial absolute pressure in Pa.

    initial_temperature:
        Initial gas temperature in K.
    """

    free_volume: float
    initial_pressure: float = 101325.0
    initial_temperature: float = 298.15

    def __post_init__(self) -> None:
        if self.free_volume <= 0.0:
            raise ValueError("Free volume must be greater than zero.")

        if self.initial_pressure < 0.0:
            raise ValueError("Initial pressure must not be negative.")

        if self.initial_temperature <= 0.0:
            raise ValueError("Initial temperature must be greater than zero.")

    @property
    def initial_moles(self) -> float:
        """
        Return the initial amount of gas in the free volume.

        Uses:

            n = P V / (R T)
        """

        return (
            self.initial_pressure
            * self.free_volume
            / (GAS_CONSTANT * self.initial_temperature)
        )

    def pressure(
        self,
        temperature: float,
        generated_moles: float = 0.0,
    ) -> float:
        """
        Return absolute internal pressure in Pa.

        The initial gas inventory is combined with the amount
        of gas generated since the initial state.

        Parameters
        ----------
        temperature:
            Gas temperature in K.

        generated_moles:
            Additional gas generated since the initial state, in mol.
        """

        if generated_moles < 0.0:
            raise ValueError("Generated moles must not be negative.")

        return self.pressure_from_total_moles(
            temperature=temperature,
            total_moles=(self.initial_moles + generated_moles),
        )

    def pressure_from_total_moles(
        self,
        temperature: float,
        total_moles: float,
    ) -> float:
        """
        Return absolute pressure from the total gas amount.

        Parameters
        ----------
        temperature:
            Gas temperature in K.

        total_moles:
            Total amount of gas currently contained in the
            free volume, in mol.
        """

        if temperature <= 0.0:
            raise ValueError("Temperature must be greater than zero.")

        if total_moles < 0.0:
            raise ValueError("Total moles must not be negative.")

        return total_moles * GAS_CONSTANT * temperature / self.free_volume
