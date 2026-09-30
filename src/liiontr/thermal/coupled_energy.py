"""Coupled cell and equilibrium-gas thermal energy."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

from scipy.optimize import brentq

from liiontr.cells.cell import Cell
from liiontr.chemistry.elements import ElementInventory
from liiontr.chemistry.equilibrium import (
    GasEquilibriumState,
)
from liiontr.thermal.equilibrium_energy import (
    EquilibriumGasEnergyModel,
)


@dataclass(slots=True)
class CoupledThermalEnergyModel:
    """Represent cell and equilibrium-gas sensible thermal energy."""

    cell: Cell
    gas_energy_model: EquilibriumGasEnergyModel

    minimum_temperature: float = 200.0
    maximum_temperature: float = 4000.0

    def __post_init__(self) -> None:
        """Validate temperature bounds."""
        if (
            not isfinite(self.minimum_temperature)
            or self.minimum_temperature <= 0.0
        ):
            raise ValueError(
                "Minimum temperature must be finite "
                "and greater than zero."
            )

        if (
            not isfinite(self.maximum_temperature)
            or self.maximum_temperature
            <= self.minimum_temperature
        ):
            raise ValueError(
                "Maximum temperature must be finite "
                "and greater than minimum temperature."
            )

    @property
    def reference_temperature(self) -> float:
        """Return the common thermal reference temperature."""
        return (
            self.gas_energy_model.reference_temperature
        )

    def cell_sensible_energy(
        self,
        temperature: float,
    ) -> float:
        """Return cell sensible energy relative to the reference."""
        if not isfinite(temperature) or temperature <= 0.0:
            raise ValueError(
                "Temperature must be finite and greater than zero."
            )

        return (
            self.cell.thermal_capacity
            * (
                temperature
                - self.reference_temperature
            )
        )

    def total_energy(
        self,
        temperature: float,
        volume: float,
        element_inventory: ElementInventory,
    ) -> float:
        """Return coupled cell and gas thermal energy."""
        cell_energy = self.cell_sensible_energy(
            temperature
        )

        gas_energy = (
            self.gas_energy_model.relative_internal_energy(
                temperature=temperature,
                volume=volume,
                element_inventory=element_inventory,
            )
        )

        return (
            cell_energy
            + gas_energy
        )

    def temperature_and_state_from_energy(
        self,
        energy: float,
        volume: float,
        element_inventory: ElementInventory,
    ) -> tuple[float, GasEquilibriumState]:
        """Recover temperature and equilibrium state from thermal energy."""
        if not isfinite(energy):
            raise ValueError("Energy must be finite.")

        reference_state = self.gas_energy_model.reference_state(
            volume=volume,
            element_inventory=element_inventory,
        )

        reference_internal_energy = reference_state.internal_energy

        state_cache: dict[
            float,
            GasEquilibriumState,
        ] = {}

        def state_at_temperature(
            temperature: float,
        ) -> GasEquilibriumState:
            cached_state = state_cache.get(temperature)

            if cached_state is not None:
                return cached_state

            equilibrium_state = self.gas_energy_model.state(
                temperature=temperature,
                volume=volume,
                element_inventory=element_inventory,
            )

            state_cache[temperature] = equilibrium_state

            return equilibrium_state

        def residual(
            temperature: float,
        ) -> float:
            equilibrium_state = state_at_temperature(temperature)

            cell_energy = self.cell_sensible_energy(temperature)

            gas_energy = equilibrium_state.internal_energy - reference_internal_energy

            return cell_energy + gas_energy - energy

        lower_residual = residual(self.minimum_temperature)

        upper_residual = residual(self.maximum_temperature)

        if lower_residual == 0.0:
            temperature = self.minimum_temperature

        elif upper_residual == 0.0:
            temperature = self.maximum_temperature

        elif lower_residual * upper_residual > 0.0:
            raise ValueError(
                "Energy is outside the configured temperature inversion range."
            )

        else:
            temperature = float(
                brentq(
                    residual,
                    self.minimum_temperature,
                    self.maximum_temperature,
                )
            )

        equilibrium_state = state_at_temperature(temperature)

        return (
            temperature,
            equilibrium_state,
        )

    def temperature_from_energy(
        self,
        energy: float,
        volume: float,
        element_inventory: ElementInventory,
    ) -> float:
        """Recover temperature from coupled thermal energy."""
        temperature, _ = self.temperature_and_state_from_energy(
            energy=energy,
            volume=volume,
            element_inventory=element_inventory,
        )

        return temperature