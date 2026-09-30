"""Thermal runaway problem definitions."""

from dataclasses import dataclass

from liiontr.cells.cell import Cell
from liiontr.chemistry import ChemistryBackend, ReactionNetworkBackend
from liiontr.chemistry.element_generation import (
    ElementGenerationModel,
)
from liiontr.chemistry.element_release import (
    ElementReleaseModel,
)
from liiontr.chemistry.elements import (
    ElementInventory,
)
from liiontr.chemistry.energy_scaled_generation import (
    EnergyScaledElementGenerationModel,
)
from liiontr.chemistry.equilibrium import (
    GasEquilibriumBackend,
)
from liiontr.core.problem import Problem
from liiontr.gases import (
    GasGenerationModel,
    GasInventory,
    IdealGasPressureModel,
    MixtureVentFlowModel,
)
from liiontr.gases.element_vent import (
    ElementVentFlowModel,
)


@dataclass(slots=True)
class ThermalProblem(Problem):
    """Lumped thermal problem definition."""

    cell: Cell

    chemistry_backend: ChemistryBackend | None = None

    initial_temperature: float = 298.15

    initial_conversions: list[float] | None = None

    ambient_temperature: float = 298.15

    convection_coefficient: float = 10.0

    duration: float = 3600.0

    maximum_temperature: float | None = None

    gas_generation_model: GasGenerationModel | None = None

    pressure_model: IdealGasPressureModel | None = None

    maximum_pressure: float | None = None

    initial_gas_inventory: GasInventory | None = None

    vent_model: MixtureVentFlowModel | None = None
    vent_open_pressure: float | None = None

    element_generation_model: (
        ElementGenerationModel | EnergyScaledElementGenerationModel | None
    ) = None

    element_release_model: ElementReleaseModel | None = None

    gas_equilibrium_backend: GasEquilibriumBackend | None = None

    initial_element_inventory: ElementInventory | None = None

    gas_volume: float | None = None

    element_vent_model: ElementVentFlowModel | None = None

    def __post_init__(self) -> None:
        """Validate elemental thermochemistry configuration."""
        if not self.uses_elemental_thermochemistry:
            return

        legacy_models = (
            self.gas_generation_model,
            self.pressure_model,
            self.initial_gas_inventory,
            self.vent_model,
        )

        if any(
            model is not None
            for model in legacy_models
        ):
            raise ValueError(
                "Elemental thermochemistry cannot be combined "
                "with legacy gas models."
            )

        if self.element_generation_model is None:
            raise ValueError(
                "Elemental thermochemistry requires an "
                "element generation model."
            )
        if (
            isinstance(
                self.element_generation_model,
                ElementGenerationModel,
            )
            and isinstance(
                self.chemistry_backend,
                ReactionNetworkBackend,
            )
            and self.element_generation_model.reaction_network
            is not self.chemistry_backend.reaction_network
        ):
            raise ValueError(
                "Element generation model must use the "
                "same reaction network as the chemistry backend."
            )

        if self.element_vent_model is not None:
            if self.gas_equilibrium_backend is None:
                raise ValueError("Element venting requires a gas equilibrium backend.")

            if (
                self.element_vent_model.equilibrium_backend
                is not self.gas_equilibrium_backend
            ):
                raise ValueError(
                    "Element vent model must use the same gas equilibrium backend."
                )

        if self.initial_element_inventory is None:
            raise ValueError(
                "Elemental thermochemistry requires an initial "
                "element inventory."
            )

        if (
            self.initial_element_inventory.total_atom_moles
            <= 0.0
        ):
            raise ValueError(
                "Initial element inventory must contain a "
                "positive amount of atoms."
            )

        if self.gas_volume is None:
            raise ValueError(
                "Elemental thermochemistry requires a gas volume."
            )

        if self.gas_volume <= 0.0:
            raise ValueError(
                "Gas volume must be greater than zero."
            )

        if (
            self.element_vent_model is not None
            and self.vent_open_pressure is None
        ):
            raise ValueError(
                "Element vent model requires a vent opening "
                "pressure."
            )

        if (
            self.vent_open_pressure is not None
            and self.element_vent_model is None
        ):
            raise ValueError(
                "Vent opening pressure requires an element "
                "vent model in elemental thermochemistry mode."
            )

        if (
            isinstance(
                self.element_generation_model,
                ElementGenerationModel,
            )
            and isinstance(
                self.chemistry_backend,
                ReactionNetworkBackend,
            )
            and self.element_generation_model.reaction_network
            is not self.chemistry_backend.reaction_network
        ):
            raise ValueError(
                "Element generation model must use the "
                "same reaction network as the chemistry backend."
            )

        if (
            self.element_release_model is not None
            and self.element_generation_model is None
        ):
            raise ValueError(
                "Element release model requires an element generation model."
            )

    @property
    def uses_elemental_thermochemistry(self) -> bool:
        """Return whether elemental thermochemistry is configured."""
        return any(
            value is not None
            for value in (
                self.element_generation_model,
                self.gas_equilibrium_backend,
                self.initial_element_inventory,
                self.gas_volume,
                self.element_vent_model,
            )
        )