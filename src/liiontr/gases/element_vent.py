"""Element-conserving vent flow from equilibrium gas states."""

from __future__ import annotations

from dataclasses import dataclass

from liiontr.chemistry.equilibrium import (
    GasEquilibriumBackend,
    GasEquilibriumState,
)
from liiontr.gases.inventory import (
    GasInventory,
    GasSpecies,
)
from liiontr.gases.vent import (
    CompressibleVentFlowModel,
    MixtureVentFlowModel,
)


@dataclass(slots=True)
class ElementVentFlowModel:
    """Convert equilibrium-gas vent flow to elemental flow."""

    equilibrium_backend: GasEquilibriumBackend

    vent_area: float

    discharge_coefficient: float = 1.0

    downstream_pressure: float = 101325.0

    def __post_init__(self) -> None:
        """Validate vent-flow configuration."""
        if self.vent_area <= 0.0:
            raise ValueError(
                "Vent area must be greater than zero."
            )

        if (
            self.discharge_coefficient <= 0.0
            or self.discharge_coefficient > 1.0
        ):
            raise ValueError(
                "Discharge coefficient must be greater "
                "than zero and at most one."
            )

        if self.downstream_pressure <= 0.0:
            raise ValueError(
                "Downstream pressure must be greater "
                "than zero."
            )

    def gas_inventory(
        self,
        state: GasEquilibriumState,
    ) -> GasInventory:
        """Build a gas inventory from an equilibrium state."""
        species = [
            GasSpecies(
                name=species_name,
                molar_mass=(
                    self.equilibrium_backend.species_molar_mass(
                        species_name
                    )
                ),
            )
            for species_name
            in state.species_moles
        ]

        return GasInventory(
            species=species,
            moles=dict(
                state.species_moles
            ),
        )

    def species_molar_flow_rates(
        self,
        state: GasEquilibriumState,
    ) -> dict[str, float]:
        """Return outward species molar vent rates in mol/s."""
        if (
            state.pressure
            <= self.downstream_pressure
        ):
            return {
                species_name: 0.0
                for species_name
                in state.species_moles
            }

        flow_model = CompressibleVentFlowModel(
            vent_area=self.vent_area,
            discharge_coefficient=(
                self.discharge_coefficient
            ),
            heat_capacity_ratio=(
                state.heat_capacity_ratio
            ),
        )

        mixture_model = MixtureVentFlowModel(
            flow_model=flow_model,
            downstream_pressure=(
                self.downstream_pressure
            ),
        )

        inventory = self.gas_inventory(
            state
        )

        return (
            mixture_model.species_molar_flow_rates(
                inventory=inventory,
                upstream_pressure=state.pressure,
                temperature=state.temperature,
            )
        )

    def element_molar_flow_rates(
        self,
        state: GasEquilibriumState,
    ) -> dict[str, float]:
        """Return outward elemental vent rates in mol-atoms/s."""
        species_rates = (
            self.species_molar_flow_rates(
                state
            )
        )

        return (
            self.equilibrium_backend.element_flow_rates(
                species_rates
            )
        )

    def total_mass_flow_rate(
        self,
        state: GasEquilibriumState,
    ) -> float:
        """Return total vent mass flow rate in kg/s."""
        species_rates = (
            self.species_molar_flow_rates(
                state
            )
        )

        return sum(
            molar_rate
            * self.equilibrium_backend.species_molar_mass(
                species_name
            )
            for species_name, molar_rate
            in species_rates.items()
        )