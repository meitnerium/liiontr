"""Energy carried by vented equilibrium gas."""

from __future__ import annotations

from dataclasses import dataclass

from liiontr.chemistry.elements import ElementInventory
from liiontr.chemistry.equilibrium import GasEquilibriumState
from liiontr.gases.element_vent import ElementVentFlowModel
from liiontr.thermal.equilibrium_energy import (
    EquilibriumGasEnergyModel,
)


@dataclass(slots=True)
class ElementVentEnergyModel:
    """Evaluate thermal energy carried by vented gas."""

    element_vent_model: ElementVentFlowModel
    gas_energy_model: EquilibriumGasEnergyModel

    def __post_init__(self) -> None:
        """Validate backend consistency."""
        if (
            self.element_vent_model.equilibrium_backend
            is not self.gas_energy_model.equilibrium_backend
        ):
            raise ValueError(
                "Vent flow and gas energy models must use "
                "the same equilibrium backend."
            )

    def energy_flow_rate(
        self,
        state: GasEquilibriumState,
        element_inventory: ElementInventory,
    ) -> float:
        """Return outward thermal-energy flow rate in W."""
        mass_flow_rate = (
            self.element_vent_model.total_mass_flow_rate(
                state
            )
        )

        if mass_flow_rate <= 0.0:
            return 0.0

        specific_enthalpy = self.gas_energy_model.relative_specific_enthalpy_from_state(
            state=state,
            element_inventory=element_inventory,
        )

        return (
            mass_flow_rate
            * specific_enthalpy
        )