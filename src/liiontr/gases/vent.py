"""Compressible gas venting models for battery thermal runaway.

This module provides models for gas discharge through a battery
safety vent. Both choked and subsonic compressible flow regimes are
supported for ideal-gas mixtures.

Code written by François Dion (francois.dion@numericatech.ca, Numerica Techologies by LDV)
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt

from .ideal import GAS_CONSTANT

from .inventory import GasInventory


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
        """Validate the compressible vent-flow parameters."""
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

        Parameters
        ----------
            upstream_pressure : float
                Absolute cell internal pressure in Pa.
            downstream_pressure : float
                Absolute ambient pressure in Pa.
            temperature : float
                Upstream gas temperature in K.
            molar_mass : float
                Mean gas molar mass in kg/mol.

        Returns
        -------
            float
                Vent mass flow rate in kg/s.

        Notes
        -----
            The flow is treated as choked when the downstream-to-upstream
            pressure ratio is below the critical pressure ratio.

            Reverse flow into the cell is not modeled.
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
        """Return the choked-flow mass discharge rate."""
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
        """Return the unchoked mass discharge rate."""
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
        """Return gas molar flow rate through the vent in mol/s."""
        mass_flow = self.mass_flow_rate(
            upstream_pressure=upstream_pressure,
            downstream_pressure=downstream_pressure,
            temperature=temperature,
            molar_mass=molar_mass,
        )

        return mass_flow / molar_mass


@dataclass(slots=True, frozen=True)
class MixtureVentFlowModel:
    """
    Vent flow model for a gas mixture.

    Total compressible flow is calculated using the
    mole-weighted mean molar mass of the mixture.

    The total molar flow is then distributed among
    species according to their mole fractions.
    """

    flow_model: CompressibleVentFlowModel

    downstream_pressure: float = 101325.0

    def __post_init__(self) -> None:
        """Validate the gas-mixture vent model."""
        if self.downstream_pressure < 0.0:
            raise ValueError("Downstream pressure must not be negative.")

    def total_molar_flow_rate(
        self,
        inventory: GasInventory,
        upstream_pressure: float,
        temperature: float,
    ) -> float:
        """Return total vent molar flow rate in mol/s."""
        if inventory.total_moles == 0.0:
            return 0.0

        return self.flow_model.molar_flow_rate(
            upstream_pressure=upstream_pressure,
            downstream_pressure=self.downstream_pressure,
            temperature=temperature,
            molar_mass=inventory.mean_molar_mass,
        )

    def species_molar_flow_rates(
        self,
        inventory: GasInventory,
        upstream_pressure: float,
        temperature: float,
    ) -> dict[str, float]:
        """Return vent molar flow rate for each gas species."""
        if inventory.total_moles == 0.0:
            return {species.name: 0.0 for species in inventory.species}

        total_flow_rate = self.total_molar_flow_rate(
            inventory=inventory,
            upstream_pressure=upstream_pressure,
            temperature=temperature,
        )

        return {
            species.name: (inventory.mole_fraction(species.name) * total_flow_rate)
            for species in inventory.species
        }
