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


def build_state():
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
        volume=1.0e-5,
        element_inventory=inventory,
    )

    return (
        backend,
        inventory,
        state,
    )


def test_element_vent_generates_positive_flow():
    backend, _, state = build_state()

    vent = ElementVentFlowModel(
        equilibrium_backend=backend,
        vent_area=1.0e-6,
        discharge_coefficient=0.8,
        downstream_pressure=101325.0,
    )

    assert state.pressure > 101325.0

    rates = vent.element_molar_flow_rates(
        state
    )

    assert rates["C"] > 0.0
    assert rates["H"] > 0.0
    assert rates["O"] > 0.0
    assert rates["N"] > 0.0


def test_element_vent_species_flow_is_positive():
    backend, _, state = build_state()

    vent = ElementVentFlowModel(
        equilibrium_backend=backend,
        vent_area=1.0e-6,
        discharge_coefficient=0.8,
    )

    rates = vent.species_molar_flow_rates(
        state
    )

    assert sum(
        rates.values()
    ) > 0.0

    assert all(
        rate >= 0.0
        for rate in rates.values()
    )


def test_element_vent_mass_flow_is_positive():
    backend, _, state = build_state()

    vent = ElementVentFlowModel(
        equilibrium_backend=backend,
        vent_area=1.0e-6,
        discharge_coefficient=0.8,
    )

    assert vent.total_mass_flow_rate(
        state
    ) > 0.0


def test_element_vent_preserves_element_ratios():
    backend, inventory, state = build_state()

    vent = ElementVentFlowModel(
        equilibrium_backend=backend,
        vent_area=1.0e-6,
        discharge_coefficient=0.8,
    )

    rates = vent.element_molar_flow_rates(
        state
    )

    carbon_rate = rates["C"]

    assert (
        rates["H"]
        / carbon_rate
    ) == pytest.approx(
        inventory.moles_of("H")
        / inventory.moles_of("C"),
        rel=1.0e-7,
    )

    assert (
        rates["O"]
        / carbon_rate
    ) == pytest.approx(
        inventory.moles_of("O")
        / inventory.moles_of("C"),
        rel=1.0e-7,
    )

    assert (
        rates["N"]
        / carbon_rate
    ) == pytest.approx(
        inventory.moles_of("N")
        / inventory.moles_of("C"),
        rel=1.0e-7,
    )


def test_element_vent_has_zero_flow_below_ambient_pressure():
    backend, _, state = build_state()

    vent = ElementVentFlowModel(
        equilibrium_backend=backend,
        vent_area=1.0e-6,
        downstream_pressure=(
            state.pressure
            + 1000.0
        ),
    )

    species_rates = (
        vent.species_molar_flow_rates(
            state
        )
    )

    assert all(
        rate == pytest.approx(0.0)
        for rate in species_rates.values()
    )


@pytest.mark.parametrize(
    "vent_area",
    [
        0.0,
        -1.0,
    ],
)
def test_element_vent_rejects_invalid_area(
    vent_area: float,
):
    backend = CanteraEquilibriumBackend()

    with pytest.raises(
        ValueError,
        match="Vent area",
    ):
        ElementVentFlowModel(
            equilibrium_backend=backend,
            vent_area=vent_area,
        )