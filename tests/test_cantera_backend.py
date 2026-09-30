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


def test_cantera_backend_loads_gri30():
    backend = CanteraEquilibriumBackend()

    assert "C" in backend.element_names
    assert "H" in backend.element_names
    assert "O" in backend.element_names
    assert "N" in backend.element_names

    assert "C" in backend.species_names
    assert "H2" in backend.species_names
    assert "O2" in backend.species_names
    assert "N2" in backend.species_names


def test_carrier_moles_preserve_atom_inventory():
    backend = CanteraEquilibriumBackend()

    inventory = ElementInventory(
        moles={
            "C": 1.0,
            "H": 4.0,
            "O": 2.0,
            "N": 2.0,
        }
    )

    carrier_moles = backend.carrier_moles(
        inventory
    )

    assert carrier_moles["C"] == pytest.approx(
        1.0
    )

    assert carrier_moles["H2"] == pytest.approx(
        2.0
    )

    assert carrier_moles["O2"] == pytest.approx(
        1.0
    )

    assert carrier_moles["N2"] == pytest.approx(
        1.0
    )


def test_inventory_mass():
    backend = CanteraEquilibriumBackend()

    inventory = ElementInventory(
        moles={
            "C": 1.0,
        }
    )

    mass = backend.inventory_mass(
        inventory
    )

    assert mass == pytest.approx(
        0.012011,
        rel=1.0e-4,
    )


def test_unsupported_element_is_rejected():
    backend = CanteraEquilibriumBackend()

    inventory = ElementInventory(
        moles={
            "F": 1.0,
        }
    )

    with pytest.raises(
        ValueError,
        match="not supported",
    ):
        backend.carrier_moles(
            inventory
        )


def test_missing_carrier_is_rejected():
    backend = CanteraEquilibriumBackend(
        carrier_species={
            "H": "H2",
            "O": "O2",
            "N": "N2",
        }
    )

    inventory = ElementInventory(
        moles={
            "C": 1.0,
        }
    )

    with pytest.raises(
        ValueError,
        match="No carrier species",
    ):
        backend.carrier_moles(
            inventory
        )


def test_invalid_mixed_element_carrier_is_rejected():
    with pytest.raises(
        ValueError,
        match="other elements",
    ):
        CanteraEquilibriumBackend(
            carrier_species={
                "C": "CO2",
            }
        )


def test_empty_inventory_is_rejected():
    backend = CanteraEquilibriumBackend()

    inventory = ElementInventory(
        moles={
            "C": 0.0,
        }
    )

    with pytest.raises(
        ValueError,
        match="positive amount",
    ):
        backend.carrier_moles(
            inventory
        )

def test_equilibrate_tv_returns_state():
    backend = CanteraEquilibriumBackend()

    inventory = ElementInventory(
        moles={
            "C": 1.0e-3,
            "H": 4.0e-3,
            "O": 4.0e-3,
            "N": 2.0e-3,
        }
    )

    state = backend.equilibrate_tv(
        temperature=1500.0,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    assert state.temperature == pytest.approx(
        1500.0
    )

    assert state.volume == pytest.approx(
        1.0e-3
    )

    assert state.pressure > 0.0

    assert state.total_moles > 0.0

    assert state.mean_molar_mass > 0.0

    assert state.cp_mass > 0.0
    assert state.cv_mass > 0.0

    assert state.heat_capacity_ratio > 1.0


def test_equilibrate_tv_conserves_mass():
    backend = CanteraEquilibriumBackend()

    inventory = ElementInventory(
        moles={
            "C": 1.0e-3,
            "H": 4.0e-3,
            "O": 4.0e-3,
            "N": 2.0e-3,
        }
    )

    expected_mass = backend.inventory_mass(
        inventory
    )

    state = backend.equilibrate_tv(
        temperature=1500.0,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    assert state.total_mass == pytest.approx(
        expected_mass,
        rel=1.0e-10,
    )


def test_equilibrate_tv_conserves_elements():
    import cantera as ct

    backend = CanteraEquilibriumBackend()

    inventory = ElementInventory(
        moles={
            "C": 1.0e-3,
            "H": 4.0e-3,
            "O": 4.0e-3,
            "N": 2.0e-3,
        }
    )

    state = backend.equilibrate_tv(
        temperature=1500.0,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    gas = ct.Solution(
        "gri30.yaml"
    )

    for element_name, expected_moles in (
        inventory.moles.items()
    ):
        calculated_moles = sum(
            species_moles
            * gas.n_atoms(
                species_name,
                element_name,
            )
            for (
                species_name,
                species_moles,
            ) in state.species_moles.items()
        )

        assert calculated_moles == pytest.approx(
            expected_moles,
            rel=1.0e-7,
            abs=1.0e-12,
        )


def test_equilibrate_tv_pressure_obeys_ideal_gas_law():
    backend = CanteraEquilibriumBackend()

    inventory = ElementInventory(
        moles={
            "C": 1.0e-3,
            "H": 4.0e-3,
            "O": 4.0e-3,
            "N": 2.0e-3,
        }
    )

    state = backend.equilibrate_tv(
        temperature=1500.0,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    universal_gas_constant = (
        8.31446261815324
    )

    expected_pressure = (
        state.total_moles
        * universal_gas_constant
        * state.temperature
        / state.volume
    )

    assert state.pressure == pytest.approx(
        expected_pressure,
        rel=1.0e-7,
    )


@pytest.mark.parametrize(
    ("temperature", "volume"),
    [
        (0.0, 1.0e-3),
        (-1.0, 1.0e-3),
        (1500.0, 0.0),
        (1500.0, -1.0),
    ],
)
def test_equilibrate_tv_rejects_invalid_state(
    temperature: float,
    volume: float,
):
    backend = CanteraEquilibriumBackend()

    inventory = ElementInventory(
        moles={
            "C": 1.0e-3,
            "H": 4.0e-3,
            "O": 4.0e-3,
        }
    )

    with pytest.raises(
        ValueError,
    ):
        backend.equilibrate_tv(
            temperature=temperature,
            volume=volume,
            element_inventory=inventory,
        )

def test_equilibrate_uv_recovers_tv_state():
    backend = CanteraEquilibriumBackend()

    inventory = ElementInventory(
        moles={
            "C": 1.0e-3,
            "H": 4.0e-3,
            "O": 4.0e-3,
            "N": 2.0e-3,
        }
    )

    tv_state = backend.equilibrate_tv(
        temperature=1500.0,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    uv_state = backend.equilibrate_uv(
        internal_energy=tv_state.internal_energy,
        volume=tv_state.volume,
        element_inventory=inventory,
    )

    assert uv_state.temperature == pytest.approx(
        tv_state.temperature,
        rel=1.0e-6,
    )

    assert uv_state.pressure == pytest.approx(
        tv_state.pressure,
        rel=1.0e-6,
    )

    assert uv_state.internal_energy == pytest.approx(
        tv_state.internal_energy,
        rel=1.0e-7,
        abs=1.0e-9,
    )

def test_equilibrate_uv_conserves_internal_energy():
    backend = CanteraEquilibriumBackend()

    inventory = ElementInventory(
        moles={
            "C": 1.0e-3,
            "H": 4.0e-3,
            "O": 4.0e-3,
            "N": 2.0e-3,
        }
    )

    reference = backend.equilibrate_tv(
        temperature=1200.0,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    result = backend.equilibrate_uv(
        internal_energy=reference.internal_energy,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    assert result.internal_energy == pytest.approx(
        reference.internal_energy,
        rel=1.0e-7,
        abs=1.0e-9,
    )

def test_equilibrate_uv_accepts_negative_internal_energy():
    backend = CanteraEquilibriumBackend()

    inventory = ElementInventory(
        moles={
            "H": 2.0e-3,
            "O": 1.0e-3,
        }
    )

    reference = backend.equilibrate_tv(
        temperature=500.0,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    result = backend.equilibrate_uv(
        internal_energy=reference.internal_energy,
        volume=reference.volume,
        element_inventory=inventory,
    )

    assert result.internal_energy == pytest.approx(
        reference.internal_energy,
        rel=1.0e-7,
    )

def test_equilibrate_uv_rejects_nonfinite_energy():
    backend = CanteraEquilibriumBackend()

    inventory = ElementInventory(
        moles={
            "H": 2.0e-3,
            "O": 1.0e-3,
        }
    )

    with pytest.raises(
        ValueError,
        match="must be finite",
    ):
        backend.equilibrate_uv(
            internal_energy=float("nan"),
            volume=1.0e-3,
            element_inventory=inventory,
        )


def test_equilibrate_uv_rejects_invalid_volume():
    backend = CanteraEquilibriumBackend()

    inventory = ElementInventory(
        moles={
            "H": 2.0e-3,
            "O": 1.0e-3,
        }
    )

    with pytest.raises(
        ValueError,
        match="Volume",
    ):
        backend.equilibrate_uv(
            internal_energy=0.0,
            volume=0.0,
            element_inventory=inventory,
        )


def test_rejects_invalid_uv_seed_temperature():
    with pytest.raises(
        ValueError,
        match="seed temperature",
    ):
        CanteraEquilibriumBackend(
            uv_seed_temperature=0.0,
        )

def test_equilibrium_state_recovers_element_inventory():
    backend = CanteraEquilibriumBackend()

    inventory = ElementInventory(
        moles={
            "C": 1.0e-3,
            "H": 4.0e-3,
            "O": 4.0e-3,
            "N": 2.0e-3,
        }
    )

    state = backend.equilibrate_tv(
        temperature=1500.0,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    recovered = (
        backend.element_inventory_from_state(
            state
        )
    )

    for element_name, expected in inventory.moles.items():
        assert recovered.moles_of(
            element_name
        ) == pytest.approx(
            expected,
            rel=1.0e-7,
            abs=1.0e-12,
        )

def test_species_flow_converts_to_element_flow():
    backend = CanteraEquilibriumBackend()

    species_rates = {
        "CO2": 1.0e-3,
        "H2O": 2.0e-3,
        "N2": 0.5e-3,
    }

    element_rates = backend.element_flow_rates(
        species_rates
    )

    assert element_rates["C"] == pytest.approx(
        1.0e-3
    )

    assert element_rates["H"] == pytest.approx(
        4.0e-3
    )

    assert element_rates["O"] == pytest.approx(
        4.0e-3
    )

    assert element_rates["N"] == pytest.approx(
        1.0e-3
    )

def test_element_flow_rejects_unknown_species():
    backend = CanteraEquilibriumBackend()

    with pytest.raises(
        ValueError,
        match="not present",
    ):
        backend.element_flow_rates(
            {
                "UNKNOWN": 1.0,
            }
        )


def test_element_flow_rejects_negative_rate():
    backend = CanteraEquilibriumBackend()

    with pytest.raises(
        ValueError,
        match="must not be negative",
    ):
        backend.element_flow_rates(
            {
                "CO2": -1.0,
            }
        )