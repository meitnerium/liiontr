import pytest

from liiontr.chemistry import (
    GasEquilibriumState,
)


def build_state() -> GasEquilibriumState:
    return GasEquilibriumState(
        temperature=800.0,
        pressure=200000.0,
        volume=1.0e-6,
        species_moles={
            "CO2": 2.0e-3,
            "CO": 1.0e-3,
            "H2": 1.0e-3,
        },
        mean_molar_mass=0.030,
        density=1.20,
        cp_mass=1200.0,
        cv_mass=900.0,
        internal_energy=100.0,
    )


def test_total_moles():
    state = build_state()

    assert state.total_moles == pytest.approx(
        4.0e-3
    )


def test_species_moles():
    state = build_state()

    assert state.moles_of(
        "CO2"
    ) == pytest.approx(
        2.0e-3
    )

    assert state.moles_of(
        "CH4"
    ) == pytest.approx(
        0.0
    )


def test_mole_fraction():
    state = build_state()

    assert state.mole_fraction(
        "CO2"
    ) == pytest.approx(
        0.5
    )

    assert state.mole_fraction(
        "CO"
    ) == pytest.approx(
        0.25
    )


def test_total_mass():
    state = build_state()

    assert state.total_mass == pytest.approx(
        1.20e-6
    )


def test_heat_capacity_ratio():
    state = build_state()

    assert state.heat_capacity_ratio == pytest.approx(
        1200.0 / 900.0
    )


def test_negative_internal_energy_is_allowed():
    state = GasEquilibriumState(
        temperature=500.0,
        pressure=101325.0,
        volume=1.0e-6,
        species_moles={
            "CO2": 1.0e-3,
        },
        mean_molar_mass=44.0e-3,
        density=1.0,
        cp_mass=1000.0,
        cv_mass=750.0,
        internal_energy=-100.0,
    )

    assert state.internal_energy == pytest.approx(
        -100.0
    )


@pytest.mark.parametrize(
    ("field_name", "value"),
    [
        ("temperature", 0.0),
        ("volume", 0.0),
        ("mean_molar_mass", 0.0),
        ("density", 0.0),
        ("cp_mass", 0.0),
        ("cv_mass", 0.0),
    ],
)
def test_rejects_nonpositive_required_properties(
    field_name: str,
    value: float,
):
    arguments = {
        "temperature": 500.0,
        "pressure": 101325.0,
        "volume": 1.0e-6,
        "species_moles": {
            "CO2": 1.0e-3,
        },
        "mean_molar_mass": 44.0e-3,
        "density": 1.0,
        "cp_mass": 1000.0,
        "cv_mass": 750.0,
        "internal_energy": 100.0,
    }

    arguments[field_name] = value

    with pytest.raises(
        ValueError,
    ):
        GasEquilibriumState(
            **arguments,
        )


def test_rejects_negative_species_moles():
    with pytest.raises(
        ValueError,
        match="must not be negative",
    ):
        GasEquilibriumState(
            temperature=500.0,
            pressure=101325.0,
            volume=1.0e-6,
            species_moles={
                "CO2": -1.0,
            },
            mean_molar_mass=44.0e-3,
            density=1.0,
            cp_mass=1000.0,
            cv_mass=750.0,
            internal_energy=100.0,
        )


def test_rejects_empty_gas_state():
    with pytest.raises(
        ValueError,
        match="positive amount",
    ):
        GasEquilibriumState(
            temperature=500.0,
            pressure=101325.0,
            volume=1.0e-6,
            species_moles={
                "CO2": 0.0,
            },
            mean_molar_mass=44.0e-3,
            density=1.0,
            cp_mass=1000.0,
            cv_mass=750.0,
            internal_energy=100.0,
        )