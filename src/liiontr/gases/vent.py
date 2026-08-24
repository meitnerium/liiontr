from __future__ import annotations

from dataclasses import dataclass
from math import sqrt

from .ideal import GAS_CONSTANT


@dataclass(slots=True, frozen=True)
class CompressibleVentFlowModel:
    """
    Compressible ideal-gas flow through a vent or orifice.

    The model automatically distinguishes between choked and
    unchoked flow according to the downstream-to-upstream
    pressure ratio.

    Parameters
    ----------
    vent_area:
        Effective geometric vent area in m2.

    discharge_coefficient:
        Dimensionless discharge coefficient.

    heat_capacity_ratio:
        Gas heat capacity ratio gamma = Cp / Cv.
    """

    vent_area: float

    discharge_coefficient: float = 1.0

    heat_capacity_ratio: float = 1.40

    def __post_init__(self) -> None:
        if self.vent_area <= 0.0:
            raise ValueError("Vent area must be greater than zero.")

        if not 0.0 < self.discharge_coefficient <= 1.0:
            raise ValueError("Discharge coefficient must be between 0 and 1.")

        if self.heat_capacity_ratio <= 1.0:
            raise ValueError("Heat capacity ratio must be greater than one.")

    @property
    def critical_pressure_ratio(self) -> float:
        """
        Return the critical downstream-to-upstream pressure ratio.

        Flow becomes choked when:

            P_downstream / P_upstream <= critical_pressure_ratio
        """

        gamma = self.heat_capacity_ratio

        return (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))

    def mass_flow_rate(
        self,
        upstream_pressure: float,
        downstream_pressure: float,
        temperature: float,
        molar_mass: float,
    ) -> float:
        """
        Return gas mass flow rate through the vent in kg/s.

        Positive flow is defined from the cell toward the environment.
        Reverse flow is not modeled.
        """

        if upstream_pressure <= 0.0:
            raise ValueError("Upstream pressure must be greater than zero.")

        if downstream_pressure < 0.0:
            raise ValueError("Downstream pressure must not be negative.")

        if temperature <= 0.0:
            raise ValueError("Temperature must be greater than zero.")

        if molar_mass <= 0.0:
            raise ValueError("Molar mass must be greater than zero.")

        if downstream_pressure >= upstream_pressure:
            return 0.0

        # gamma = self.heat_capacity_ratio

        specific_gas_constant = GAS_CONSTANT / molar_mass

        pressure_ratio = downstream_pressure / upstream_pressure

        if pressure_ratio <= self.critical_pressure_ratio:
            return self._choked_mass_flow_rate(
                upstream_pressure=upstream_pressure,
                temperature=temperature,
                specific_gas_constant=(specific_gas_constant),
            )

        return self._unchoked_mass_flow_rate(
            upstream_pressure=upstream_pressure,
            pressure_ratio=pressure_ratio,
            temperature=temperature,
            specific_gas_constant=(specific_gas_constant),
        )

    def _choked_mass_flow_rate(
        self,
        upstream_pressure: float,
        temperature: float,
        specific_gas_constant: float,
    ) -> float:
        gamma = self.heat_capacity_ratio

        sonic_factor = (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))

        return (
            self.discharge_coefficient
            * self.vent_area
            * upstream_pressure
            * sqrt(gamma / (specific_gas_constant * temperature))
            * sonic_factor
        )

    def _unchoked_mass_flow_rate(
        self,
        upstream_pressure: float,
        pressure_ratio: float,
        temperature: float,
        specific_gas_constant: float,
    ) -> float:
        gamma = self.heat_capacity_ratio

        pressure_term = pressure_ratio ** (2.0 / gamma) - pressure_ratio ** (
            (gamma + 1.0) / gamma
        )

        coefficient = (
            2.0 * gamma / (specific_gas_constant * temperature * (gamma - 1.0))
        )

        return (
            self.discharge_coefficient
            * self.vent_area
            * upstream_pressure
            * sqrt(coefficient * pressure_term)
        )

    def molar_flow_rate(
        self,
        upstream_pressure: float,
        downstream_pressure: float,
        temperature: float,
        molar_mass: float,
    ) -> float:
        """
        Return gas molar flow rate through the vent in mol/s.
        """

        mass_flow = self.mass_flow_rate(
            upstream_pressure=upstream_pressure,
            downstream_pressure=downstream_pressure,
            temperature=temperature,
            molar_mass=molar_mass,
        )

        return mass_flow / molar_mass
