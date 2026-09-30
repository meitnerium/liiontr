from __future__ import annotations

import numpy as np
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
from liiontr.thermal.equilibrium_energy import (
    EquilibriumGasEnergyModel,
)


def build_inventory() -> ElementInventory:
    return ElementInventory(
        moles={
            "C": 1.0e-3,
            "H": 4.0e-3,
            "O": 4.0e-3,
            "N": 2.0e-3,
        }
    )

def test_relative_energy_is_zero_at_reference_temperature():
    backend = CanteraEquilibriumBackend()

    model = EquilibriumGasEnergyModel(
        equilibrium_backend=backend,
        reference_temperature=298.15,
    )

    energy = model.relative_internal_energy(
        temperature=298.15,
        volume=1.0e-3,
        element_inventory=build_inventory(),
    )

    assert energy == pytest.approx(
        0.0,
        abs=1.0e-10,
    )

def test_relative_energy_increases_with_temperature():
    backend = CanteraEquilibriumBackend()

    model = EquilibriumGasEnergyModel(
        equilibrium_backend=backend,
    )

    inventory = ElementInventory(
        moles={
            "N": 2.0e-3,
        }
    )

    energy = model.relative_internal_energy(
        temperature=1000.0,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    assert energy > 0.0

def test_relative_energy_matches_state_difference():
    backend = CanteraEquilibriumBackend()

    model = EquilibriumGasEnergyModel(
        equilibrium_backend=backend,
        reference_temperature=298.15,
    )

    inventory = build_inventory()

    current_state = backend.equilibrate_tv(
        temperature=1500.0,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    reference_state = backend.equilibrate_tv(
        temperature=298.15,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    expected = (
        current_state.internal_energy
        - reference_state.internal_energy
    )

    actual = model.relative_internal_energy(
        temperature=1500.0,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    assert actual == pytest.approx(
        expected,
        rel=1.0e-12,
        abs=1.0e-10,
    )

def test_relative_energy_is_finite():
    backend = CanteraEquilibriumBackend()

    model = EquilibriumGasEnergyModel(
        equilibrium_backend=backend,
    )

    energy = model.relative_internal_energy(
        temperature=1500.0,
        volume=1.0e-3,
        element_inventory=build_inventory(),
    )

    assert np.isfinite(
        energy
    )

@pytest.mark.parametrize(
    "reference_temperature",
    [
        0.0,
        -1.0,
        float("nan"),
        float("inf"),
    ],
)
def test_rejects_invalid_reference_temperature(
    reference_temperature: float,
):
    backend = CanteraEquilibriumBackend()

    with pytest.raises(
        ValueError,
        match="Reference temperature",
    ):
        EquilibriumGasEnergyModel(
            equilibrium_backend=backend,
            reference_temperature=reference_temperature,
        )

@pytest.mark.parametrize(
    "volume",
    [
        0.0,
        -1.0,
        float("nan"),
    ],
)
def test_relative_energy_rejects_invalid_volume(
    volume: float,
):
    backend = CanteraEquilibriumBackend()

    model = EquilibriumGasEnergyModel(
        equilibrium_backend=backend,
    )

    with pytest.raises(
        ValueError,
        match="Volume",
    ):
        model.relative_internal_energy(
            temperature=1000.0,
            volume=volume,
            element_inventory=build_inventory(),
        )

def test_relative_enthalpy_is_zero_at_reference_temperature():
    backend = CanteraEquilibriumBackend()

    model = EquilibriumGasEnergyModel(
        equilibrium_backend=backend,
        reference_temperature=298.15,
    )

    enthalpy = model.relative_enthalpy(
        temperature=298.15,
        volume=1.0e-3,
        element_inventory=build_inventory(),
    )

    assert enthalpy == pytest.approx(
        0.0,
        abs=1.0e-10,
    )

def test_relative_enthalpy_matches_energy_plus_pv():
    backend = CanteraEquilibriumBackend()

    model = EquilibriumGasEnergyModel(
        equilibrium_backend=backend,
        reference_temperature=298.15,
    )

    inventory = build_inventory()

    volume = 1.0e-3

    current_state = backend.equilibrate_tv(
        temperature=1500.0,
        volume=volume,
        element_inventory=inventory,
    )

    reference_state = backend.equilibrate_tv(
        temperature=298.15,
        volume=volume,
        element_inventory=inventory,
    )

    expected = (
        current_state.internal_energy
        + current_state.pressure * volume
        - reference_state.internal_energy
        - reference_state.pressure * volume
    )

    actual = model.relative_enthalpy(
        temperature=1500.0,
        volume=volume,
        element_inventory=inventory,
    )

    assert actual == pytest.approx(
        expected,
        rel=1.0e-12,
        abs=1.0e-10,
    )

def test_relative_specific_enthalpy_is_positive_when_hot():
    backend = CanteraEquilibriumBackend()

    model = EquilibriumGasEnergyModel(
        equilibrium_backend=backend,
    )

    inventory = ElementInventory(
        moles={
            "N": 2.0e-3,
        }
    )

    specific_enthalpy = (
        model.relative_specific_enthalpy(
            temperature=1000.0,
            volume=1.0e-3,
            element_inventory=inventory,
        )
    )

    assert specific_enthalpy > 0.0

def test_relative_specific_enthalpy_is_extensive_consistent():
    backend = CanteraEquilibriumBackend()

    model = EquilibriumGasEnergyModel(
        equilibrium_backend=backend,
    )

    inventory_1 = ElementInventory(
        moles={
            "N": 2.0e-3,
        }
    )

    inventory_2 = ElementInventory(
        moles={
            "N": 4.0e-3,
        }
    )

    enthalpy_1 = (
        model.relative_specific_enthalpy(
            temperature=1000.0,
            volume=1.0e-3,
            element_inventory=inventory_1,
        )
    )

    enthalpy_2 = (
        model.relative_specific_enthalpy(
            temperature=1000.0,
            volume=2.0e-3,
            element_inventory=inventory_2,
        )
    )

    assert enthalpy_2 == pytest.approx(
        enthalpy_1,
        rel=1.0e-10,
    )

def test_specific_enthalpy_from_state_matches_direct_calculation():
    backend = CanteraEquilibriumBackend()

    model = EquilibriumGasEnergyModel(
        equilibrium_backend=backend,
    )

    inventory = build_inventory()

    volume = 1.0e-3
    temperature = 1000.0

    state = backend.equilibrate_tv(
        temperature=temperature,
        volume=volume,
        element_inventory=inventory,
    )

    direct = model.relative_specific_enthalpy(
        temperature=temperature,
        volume=volume,
        element_inventory=inventory,
    )

    reused = (
        model.relative_specific_enthalpy_from_state(
            state=state,
            element_inventory=inventory,
        )
    )

    assert reused == pytest.approx(
        direct,
        rel=1.0e-12,
        abs=1.0e-10,
    )