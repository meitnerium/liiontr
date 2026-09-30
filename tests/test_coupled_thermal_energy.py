from __future__ import annotations

import pytest

pytest.importorskip(
    "cantera"
)

from itertools import pairwise

from liiontr.chemistry.cantera import (
    CanteraEquilibriumBackend,
)
from liiontr.chemistry.elements import (
    ElementInventory,
)
from liiontr.library import (
    cell_21700_generic,
)
from liiontr.thermal.coupled_energy import (
    CoupledThermalEnergyModel,
)
from liiontr.thermal.equilibrium_energy import (
    EquilibriumGasEnergyModel,
)


def build_model() -> CoupledThermalEnergyModel:
    backend = CanteraEquilibriumBackend()

    gas_model = EquilibriumGasEnergyModel(
        equilibrium_backend=backend,
        reference_temperature=298.15,
    )

    return CoupledThermalEnergyModel(
        cell=cell_21700_generic(),
        gas_energy_model=gas_model,
        minimum_temperature=250.0,
        maximum_temperature=2500.0,
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

def test_total_energy_is_zero_at_reference_temperature():
    model = build_model()

    energy = model.total_energy(
        temperature=298.15,
        volume=1.0e-3,
        element_inventory=build_inventory(),
    )

    assert energy == pytest.approx(
        0.0,
        abs=1.0e-10,
    )

@pytest.mark.parametrize(
    "temperature",
    [
        300.0,
        500.0,
        1000.0,
        1500.0,
        2000.0,
    ],
)
def test_temperature_energy_roundtrip(
    temperature: float,
):
    model = build_model()

    inventory = build_inventory()

    energy = model.total_energy(
        temperature=temperature,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    recovered_temperature = (
        model.temperature_from_energy(
            energy=energy,
            volume=1.0e-3,
            element_inventory=inventory,
        )
    )

    assert recovered_temperature == pytest.approx(
        temperature,
        rel=1.0e-8,
        abs=1.0e-6,
    )

def test_total_energy_increases_with_temperature():
    model = build_model()

    inventory = build_inventory()

    temperatures = [
        300.0,
        500.0,
        1000.0,
        1500.0,
        2000.0,
    ]

    energies = [
        model.total_energy(
            temperature=temperature,
            volume=1.0e-3,
            element_inventory=inventory,
        )
        for temperature in temperatures
    ]

    assert all(later > earlier for earlier, later in pairwise(energies))

def test_gas_changes_total_thermal_energy():
    model = build_model()

    inventory = build_inventory()

    temperature = 1500.0

    total_energy = model.total_energy(
        temperature=temperature,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    cell_energy = model.cell_sensible_energy(
        temperature
    )

    assert total_energy != pytest.approx(
        cell_energy
    )

def test_temperature_inversion_rejects_energy_above_range():
    model = build_model()

    inventory = build_inventory()

    maximum_energy = model.total_energy(
        temperature=model.maximum_temperature,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    with pytest.raises(
        ValueError,
        match="outside",
    ):
        model.temperature_from_energy(
            energy=maximum_energy + 1.0e9,
            volume=1.0e-3,
            element_inventory=inventory,
        )

def test_temperature_inversion_rejects_energy_below_range():
    model = build_model()

    inventory = build_inventory()

    minimum_energy = model.total_energy(
        temperature=model.minimum_temperature,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    with pytest.raises(
        ValueError,
        match="outside",
    ):
        model.temperature_from_energy(
            energy=minimum_energy - 1.0e9,
            volume=1.0e-3,
            element_inventory=inventory,
        )