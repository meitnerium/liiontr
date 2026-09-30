"""Tests for the Hu 2020 and Howard 2025 coupling."""

import pytest

from liiontr.library.hu2020_howard2025 import (
    HU2020_REFERENCE_REACTION_ENERGY_PER_CELL_MASS,
    hu2020_howard2025_element_generation_model,
)


def test_hu2020_reference_reaction_energy():
    """Preserve the Hu baseline energy-balance reference."""
    assert (
        HU2020_REFERENCE_REACTION_ENERGY_PER_CELL_MASS
        == pytest.approx(
            386405.782042
        )
    )


def test_hu2020_howard2025_element_yields():
    """Preserve the Howard elemental yields in the coupling."""
    model = (
        hu2020_howard2025_element_generation_model()
    )

    yields = (
        model.element_yields_per_cell_mass
    )

    assert yields["C"] == pytest.approx(
        4.4108875,
        rel=1.0e-6,
    )

    assert yields["H"] == pytest.approx(
        5.2830969,
        rel=1.0e-6,
    )

    assert yields["O"] == pytest.approx(
        4.6551062,
        rel=1.0e-6,
    )


def test_hu2020_howard2025_full_energy_generates_full_inventory():
    """Generate the full Howard inventory at full Hu energy."""
    model = (
        hu2020_howard2025_element_generation_model()
    )

    cell_mass = 0.068

    released_energy = (
        HU2020_REFERENCE_REACTION_ENERGY_PER_CELL_MASS
        * cell_mass
    )

    inventory = model.generated_inventory(
        released_reaction_energy=released_energy,
        cell_mass=cell_mass,
    )

    assert inventory.moles["C"] == pytest.approx(
        0.29994035,
        rel=1.0e-6,
    )

    assert inventory.moles["H"] == pytest.approx(
        0.35925059,
        rel=1.0e-6,
    )

    assert inventory.moles["O"] == pytest.approx(
        0.31654722,
        rel=1.0e-6,
    )