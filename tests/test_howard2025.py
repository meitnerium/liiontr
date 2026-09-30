"""Tests for the Howard et al. (2025) gas dataset."""

import pytest

from liiontr.library.howard2025 import (
    howard2025_nmc811_21700_100soc_inert,
)


def test_howard2025_cell_properties():
    """Preserve the reported cell properties."""
    dataset = (
        howard2025_nmc811_21700_100soc_inert()
    )

    assert dataset.cell_format == "21700"
    assert dataset.chemistry == "NMC811"

    assert dataset.capacity_ah == pytest.approx(
        5.0
    )

    assert dataset.cell_mass == pytest.approx(
        0.068
    )

    assert dataset.state_of_charge == pytest.approx(
        1.0
    )

    assert dataset.gas_volume == pytest.approx(
        8.3e-3
    )


def test_howard2025_species_composition():
    """Preserve the reported vent-gas composition."""
    dataset = (
        howard2025_nmc811_21700_100soc_inert()
    )

    expected = {
        "H2": 0.220,
        "CO2": 0.278,
        "CO": 0.378,
        "C2H6": 0.009,
        "C2H4": 0.039,
        "C3H8": 0.023,
        "C3H6": 0.005,
        "CH4": 0.049,
    }

    assert dataset.species_fractions == expected

    assert sum(
        dataset.species_fractions.values()
    ) == pytest.approx(
        1.001
    )

def test_howard2025_species_moles():
    """Convert the reported gas volume to species moles."""
    dataset = (
        howard2025_nmc811_21700_100soc_inert()
    )

    species_moles = dataset.species_moles(
        temperature=298.15,
        pressure=101325.0,
    )

    assert sum(
        species_moles.values()
    ) == pytest.approx(
        0.33925457,
        rel=1.0e-7,
    )

    assert species_moles["H2"] == pytest.approx(
        0.07456144,
        rel=1.0e-6,
    )

    assert species_moles["CO2"] == pytest.approx(
        0.09421855,
        rel=1.0e-6,
    )

    assert species_moles["CO"] == pytest.approx(
        0.12811012,
        rel=1.0e-6,
    )

def test_howard2025_species_moles_rejects_invalid_temperature():
    """Reject invalid gas conversion temperature."""
    dataset = (
        howard2025_nmc811_21700_100soc_inert()
    )

    with pytest.raises(
        ValueError,
        match="Temperature",
    ):
        dataset.species_moles(
            temperature=0.0,
            pressure=101325.0,
        )


def test_howard2025_species_moles_rejects_invalid_pressure():
    """Reject invalid gas conversion pressure."""
    dataset = (
        howard2025_nmc811_21700_100soc_inert()
    )

    with pytest.raises(
        ValueError,
        match="Pressure",
    ):
        dataset.species_moles(
            temperature=298.15,
            pressure=0.0,
        )