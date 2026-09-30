from __future__ import annotations

import pytest

pytest.importorskip(
    "cantera"
)

from liiontr.chemistry.cantera import (
    CanteraEquilibriumBackend,
)
from liiontr.chemistry.elements import (
    ElementInventory,
)
from liiontr.gases.element_vent import (
    ElementVentFlowModel,
)
from liiontr.thermal.equilibrium_energy import (
    EquilibriumGasEnergyModel,
)
from liiontr.thermal.vent_energy import (
    ElementVentEnergyModel,
)


def build_models():
    backend = CanteraEquilibriumBackend()

    vent_model = ElementVentFlowModel(
        equilibrium_backend=backend,
        vent_area=1.0e-6,
        discharge_coefficient=0.8,
        downstream_pressure=101325.0,
    )

    gas_energy_model = EquilibriumGasEnergyModel(
        equilibrium_backend=backend,
        reference_temperature=298.15,
    )

    energy_model = ElementVentEnergyModel(
        element_vent_model=vent_model,
        gas_energy_model=gas_energy_model,
    )

    return (
        backend,
        vent_model,
        gas_energy_model,
        energy_model,
    )

def test_hot_vented_gas_carries_energy():
    (
        backend,
        _,
        _,
        energy_model,
    ) = build_models()

    inventory = ElementInventory(
        moles={
            "N": 2.0e-2,
        }
    )

    state = backend.equilibrate_tv(
        temperature=1000.0,
        volume=1.0e-4,
        element_inventory=inventory,
    )

    energy_rate = (
        energy_model.energy_flow_rate(
            state=state,
            element_inventory=inventory,
        )
    )

    assert energy_rate > 0.0

def test_vent_energy_matches_mass_flow_times_enthalpy():
    (
        backend,
        vent_model,
        gas_energy_model,
        energy_model,
    ) = build_models()

    inventory = ElementInventory(
        moles={
            "N": 2.0e-2,
        }
    )

    state = backend.equilibrate_tv(
        temperature=1000.0,
        volume=1.0e-4,
        element_inventory=inventory,
    )

    mass_flow_rate = (
        vent_model.total_mass_flow_rate(
            state
        )
    )

    specific_enthalpy = (
        gas_energy_model.relative_specific_enthalpy(
            temperature=state.temperature,
            volume=state.volume,
            element_inventory=inventory,
        )
    )

    expected = (
        mass_flow_rate
        * specific_enthalpy
    )

    actual = (
        energy_model.energy_flow_rate(
            state=state,
            element_inventory=inventory,
        )
    )

    assert actual == pytest.approx(
        expected,
        rel=1.0e-12,
    )

def test_no_energy_flow_below_downstream_pressure():
    (
        backend,
        _,
        _,
        energy_model,
    ) = build_models()

    inventory = ElementInventory(
        moles={
            "N": 2.0e-3,
        }
    )

    state = backend.equilibrate_tv(
        temperature=1000.0,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    assert state.pressure < 101325.0

    energy_rate = (
        energy_model.energy_flow_rate(
            state=state,
            element_inventory=inventory,
        )
    )

    assert energy_rate == pytest.approx(
        0.0,
        abs=1.0e-15,
    )

def test_reference_temperature_has_zero_relative_energy_flow():
    (
        backend,
        vent_model,
        _,
        energy_model,
    ) = build_models()

    inventory = ElementInventory(
        moles={
            "N": 2.0e-2,
        }
    )

    state = backend.equilibrate_tv(
        temperature=298.15,
        volume=1.0e-4,
        element_inventory=inventory,
    )

    assert (
        vent_model.total_mass_flow_rate(
            state
        )
        > 0.0
    )

    energy_rate = (
        energy_model.energy_flow_rate(
            state=state,
            element_inventory=inventory,
        )
    )

    assert energy_rate == pytest.approx(
        0.0,
        abs=1.0e-10,
    )

def test_rejects_different_equilibrium_backends():
    backend_1 = CanteraEquilibriumBackend()
    backend_2 = CanteraEquilibriumBackend()

    vent_model = ElementVentFlowModel(
        equilibrium_backend=backend_1,
        vent_area=1.0e-6,
    )

    gas_energy_model = EquilibriumGasEnergyModel(
        equilibrium_backend=backend_2,
    )

    with pytest.raises(
        ValueError,
        match="same equilibrium backend",
    ):
        ElementVentEnergyModel(
            element_vent_model=vent_model,
            gas_energy_model=gas_energy_model,
        )