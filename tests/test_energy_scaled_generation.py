"""Tests for energy-scaled elemental gas generation."""

import pytest

from liiontr.chemistry.energy_scaled_generation import (
    EnergyScaledElementGenerationModel,
)


def build_model():
    """Build a simple energy-scaled generation model."""
    return EnergyScaledElementGenerationModel(
        element_yields_per_cell_mass={
            "C": 4.0,
            "H": 5.0,
            "O": 3.0,
        },
        reference_reaction_energy_per_cell_mass=400000.0,
    )


def test_energy_scaled_generation_rates():
    """Scale elemental rates with reaction heat generation."""
    model = build_model()

    rates = model.generation_rates(
        heat_generation=400000.0,
    )

    assert rates["C"] == pytest.approx(
        4.0
    )

    assert rates["H"] == pytest.approx(
        5.0
    )

    assert rates["O"] == pytest.approx(
        3.0
    )


def test_energy_scaled_generation_half_inventory():
    """Generate half the final inventory at half reference energy."""
    model = build_model()

    inventory = model.generated_inventory(
        released_reaction_energy=20000.0,
        cell_mass=0.1,
    )

    assert inventory.moles["C"] == pytest.approx(
        0.2
    )

    assert inventory.moles["H"] == pytest.approx(
        0.25
    )

    assert inventory.moles["O"] == pytest.approx(
        0.15
    )


def test_energy_scaled_generation_caps_inventory():
    """Cap cumulative gas generation at the literature yield."""
    model = build_model()

    inventory = model.generated_inventory(
        released_reaction_energy=80000.0,
        cell_mass=0.1,
    )

    assert inventory.moles["C"] == pytest.approx(
        0.4
    )

    assert inventory.moles["H"] == pytest.approx(
        0.5
    )

    assert inventory.moles["O"] == pytest.approx(
        0.3
    )


def test_energy_scaled_generation_rejects_negative_heat():
    """Reject negative reaction heat generation."""
    model = build_model()

    with pytest.raises(
        ValueError,
        match="Heat generation",
    ):
        model.generation_rates(
            heat_generation=-1.0,
        )

def test_generation_progress_is_bounded():
    model = EnergyScaledElementGenerationModel(
        element_yields_per_cell_mass={
            "C": 2.0,
        },
        reference_reaction_energy_per_cell_mass=(
            1000.0
        ),
    )

    assert model.generation_progress(
        released_reaction_energy=0.0,
        cell_mass=2.0,
    ) == pytest.approx(
        0.0
    )

    assert model.generation_progress(
        released_reaction_energy=1000.0,
        cell_mass=2.0,
    ) == pytest.approx(
        0.5
    )

    assert model.generation_progress(
        released_reaction_energy=2000.0,
        cell_mass=2.0,
    ) == pytest.approx(
        1.0
    )

    assert model.generation_progress(
        released_reaction_energy=4000.0,
        cell_mass=2.0,
    ) == pytest.approx(
        1.0
    )

def test_generated_inventory_does_not_exceed_yield():
    model = EnergyScaledElementGenerationModel(
        element_yields_per_cell_mass={
            "C": 2.0,
            "H": 3.0,
        },
        reference_reaction_energy_per_cell_mass=(
            1000.0
        ),
    )

    inventory = model.generated_inventory(
        released_reaction_energy=1.0e9,
        cell_mass=2.0,
    )

    assert inventory.moles_of(
        "C"
    ) == pytest.approx(
        4.0
    )

    assert inventory.moles_of(
        "H"
    ) == pytest.approx(
        6.0
    )