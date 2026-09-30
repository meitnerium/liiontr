from __future__ import annotations

import pytest

from liiontr.chemistry.element_generation import (
    ElementGenerationModel,
)
from liiontr.chemistry.elements import (
    ElementInventory,
    ReactionElementYield,
)
from liiontr.chemistry.equilibrium import (
    GasEquilibriumState,
)
from liiontr.gases.element_vent import (
    ElementVentFlowModel,
)
from liiontr.kinetics import Arrhenius
from liiontr.library import cell_21700_generic
from liiontr.problems import ThermalProblem
from liiontr.reactions import (
    Reaction,
    ReactionNetwork,
)


class DummyEquilibriumBackend:
    """Minimal thermochemical backend used for configuration tests."""

    @property
    def element_names(self) -> tuple[str, ...]:
        """Return supported elements."""
        return (
            "C",
            "H",
            "O",
            "N",
        )

    @property
    def species_names(self) -> tuple[str, ...]:
        """Return supported species."""
        return (
            "CO2",
            "H2O",
            "N2",
        )

    def species_molar_mass(
        self,
        species_name: str,
    ) -> float:
        """Return an arbitrary positive molar mass."""
        del species_name
        return 0.028

    def element_flow_rates(
        self,
        species_molar_flow_rates: dict[str, float],
    ) -> dict[str, float]:
        """Return dummy element rates."""
        del species_molar_flow_rates
        return {}

    def equilibrate_tv(
        self,
        temperature: float,
        volume: float,
        element_inventory: ElementInventory,
    ) -> GasEquilibriumState:
        """Not used by these tests."""
        raise NotImplementedError

    def equilibrate_uv(
        self,
        internal_energy: float,
        volume: float,
        element_inventory: ElementInventory,
    ) -> GasEquilibriumState:
        """Not used by these tests."""
        raise NotImplementedError


def build_element_generation_model() -> ElementGenerationModel:
    reaction = Reaction(
        name="Reaction A",
        kinetics=Arrhenius(
            activation_energy=1.0,
            pre_exponential_factor=2.0,
        ),
        enthalpy=10000.0,
        mass_fraction=0.10,
    )

    network = ReactionNetwork(
        reactions=[
            reaction,
        ]
    )

    return ElementGenerationModel(
        reaction_network=network,
        element_yields=[
            ReactionElementYield(
                reaction_name="Reaction A",
                element_yields={
                    "C": 1.0,
                    "H": 2.0,
                },
            )
        ],
    )


def build_inventory() -> ElementInventory:
    return ElementInventory(
        moles={
            "N": 1.0e-3,
        }
    )

def test_thermal_problem_has_no_elemental_mode_by_default():
    problem = ThermalProblem(
        cell=cell_21700_generic(),
    )

    assert not problem.uses_elemental_thermochemistry


def test_thermal_problem_accepts_elemental_configuration():
    backend = DummyEquilibriumBackend()

    problem = ThermalProblem(
        cell=cell_21700_generic(),
        element_generation_model=(
            build_element_generation_model()
        ),
        gas_equilibrium_backend=backend,
        initial_element_inventory=build_inventory(),
        gas_volume=1.0e-6,
    )

    assert problem.uses_elemental_thermochemistry

    assert problem.gas_equilibrium_backend is backend

    assert problem.gas_volume == pytest.approx(
        1.0e-6
    )

def test_elemental_mode_rejects_missing_volume():
    backend = DummyEquilibriumBackend()

    with pytest.raises(
        ValueError,
        match="gas volume",
    ):
        ThermalProblem(
            cell=cell_21700_generic(),
            element_generation_model=(
                build_element_generation_model()
            ),
            gas_equilibrium_backend=backend,
            initial_element_inventory=build_inventory(),
        )


def test_elemental_mode_rejects_nonpositive_volume():
    backend = DummyEquilibriumBackend()

    with pytest.raises(
        ValueError,
        match="Gas volume",
    ):
        ThermalProblem(
            cell=cell_21700_generic(),
            element_generation_model=(
                build_element_generation_model()
            ),
            gas_equilibrium_backend=backend,
            initial_element_inventory=build_inventory(),
            gas_volume=0.0,
        )

def test_elemental_mode_rejects_legacy_pressure_model():
    from liiontr.gases import (
        IdealGasPressureModel,
    )

    backend = DummyEquilibriumBackend()

    with pytest.raises(
        ValueError,
        match="legacy gas models",
    ):
        ThermalProblem(
            cell=cell_21700_generic(),
            element_generation_model=(
                build_element_generation_model()
            ),
            gas_equilibrium_backend=backend,
            initial_element_inventory=build_inventory(),
            gas_volume=1.0e-6,
            pressure_model=IdealGasPressureModel(
                free_volume=1.0e-6,
            ),
        )

def test_elemental_problem_accepts_element_vent():
    backend = DummyEquilibriumBackend()

    vent_model = ElementVentFlowModel(
        equilibrium_backend=backend,
        vent_area=1.0e-6,
    )

    problem = ThermalProblem(
        cell=cell_21700_generic(),
        element_generation_model=(
            build_element_generation_model()
        ),
        gas_equilibrium_backend=backend,
        initial_element_inventory=build_inventory(),
        gas_volume=1.0e-6,
        element_vent_model=vent_model,
        vent_open_pressure=2.0e5,
    )

    assert problem.element_vent_model is vent_model
    assert problem.vent_open_pressure == pytest.approx(
        2.0e5
    )

def test_element_vent_requires_same_backend():
    backend_1 = DummyEquilibriumBackend()
    backend_2 = DummyEquilibriumBackend()

    vent_model = ElementVentFlowModel(
        equilibrium_backend=backend_2,
        vent_area=1.0e-6,
    )

    with pytest.raises(
        ValueError,
        match="same gas equilibrium backend",
    ):
        ThermalProblem(
            cell=cell_21700_generic(),
            element_generation_model=(
                build_element_generation_model()
            ),
            gas_equilibrium_backend=backend_1,
            initial_element_inventory=build_inventory(),
            gas_volume=1.0e-6,
            element_vent_model=vent_model,
            vent_open_pressure=2.0e5,
        )