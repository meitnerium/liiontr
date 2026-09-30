"""Thermochemical gas-energy models."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

from liiontr.chemistry.elements import ElementInventory
from liiontr.chemistry.equilibrium import (
    GasEquilibriumBackend,
    GasEquilibriumState,
)


@dataclass(slots=True)
class EquilibriumGasEnergyModel:
    """Evaluate equilibrium gas energy relative to a reference state."""

    equilibrium_backend: GasEquilibriumBackend
    reference_temperature: float = 298.15

    def __post_init__(self) -> None:
        """Validate the energy-reference configuration."""
        if (
            not isfinite(self.reference_temperature)
            or self.reference_temperature <= 0.0
        ):
            raise ValueError(
                "Reference temperature must be finite "
                "and greater than zero."
            )

    def reference_state(
        self,
        volume: float,
        element_inventory: ElementInventory,
    ) -> GasEquilibriumState:
        """Return the equilibrium reference state."""
        if not isfinite(volume) or volume <= 0.0:
            raise ValueError(
                "Volume must be finite and greater than zero."
            )

        return self.equilibrium_backend.equilibrate_tv(
            temperature=self.reference_temperature,
            volume=volume,
            element_inventory=element_inventory,
        )

    def state(
        self,
        temperature: float,
        volume: float,
        element_inventory: ElementInventory,
    ) -> GasEquilibriumState:
        """Return the equilibrium state at the requested temperature."""
        if not isfinite(temperature) or temperature <= 0.0:
            raise ValueError(
                "Temperature must be finite and greater than zero."
            )

        if not isfinite(volume) or volume <= 0.0:
            raise ValueError(
                "Volume must be finite and greater than zero."
            )

        return self.equilibrium_backend.equilibrate_tv(
            temperature=temperature,
            volume=volume,
            element_inventory=element_inventory,
        )

    def relative_internal_energy(
        self,
        temperature: float,
        volume: float,
        element_inventory: ElementInventory,
    ) -> float:
        """Return gas internal energy relative to the reference state."""
        current_state = self.state(
            temperature=temperature,
            volume=volume,
            element_inventory=element_inventory,
        )

        reference_state = self.reference_state(
            volume=volume,
            element_inventory=element_inventory,
        )

        return (
            current_state.internal_energy
            - reference_state.internal_energy
        )

    def relative_enthalpy(
        self,
        temperature: float,
        volume: float,
        element_inventory: ElementInventory,
    ) -> float:
        """Return gas enthalpy relative to the reference state."""
        current_state = self.state(
            temperature=temperature,
            volume=volume,
            element_inventory=element_inventory,
        )

        reference_state = self.reference_state(
            volume=volume,
            element_inventory=element_inventory,
        )

        current_enthalpy = (
            current_state.internal_energy
            + current_state.pressure
            * volume
        )

        reference_enthalpy = (
            reference_state.internal_energy
            + reference_state.pressure
            * volume
        )

        return (
            current_enthalpy
            - reference_enthalpy
        )

    def relative_specific_enthalpy(
        self,
        temperature: float,
        volume: float,
        element_inventory: ElementInventory,
    ) -> float:
        """Return relative gas enthalpy per unit mass in J/kg."""
        state = self.state(
            temperature=temperature,
            volume=volume,
            element_inventory=element_inventory,
        )

        return self.relative_specific_enthalpy_from_state(
            state=state,
            element_inventory=element_inventory,
        )

    def relative_specific_enthalpy_from_state(
        self,
        state: GasEquilibriumState,
        element_inventory: ElementInventory,
    ) -> float:
        """Return relative specific enthalpy from an existing state."""
        reference_state = self.reference_state(
            volume=state.volume,
            element_inventory=element_inventory,
        )

        current_enthalpy = state.internal_energy + state.pressure * state.volume

        reference_enthalpy = (
            reference_state.internal_energy + reference_state.pressure * state.volume
        )

        return (current_enthalpy - reference_enthalpy) / state.total_mass