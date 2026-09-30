"""Cantera integration tests for the Howard et al. (2025) dataset."""

from __future__ import annotations

import pytest

pytest.importorskip("cantera")

from liiontr.chemistry.cantera import (
    CanteraEquilibriumBackend,
)
from liiontr.chemistry.elements import (
    ElementInventory,
)
from liiontr.library.howard2025 import (
    howard2025_nmc811_21700_100soc_inert,
)


def test_howard2025_element_inventory():
    """Convert the Howard gas composition to elemental inventory."""
    dataset = (
        howard2025_nmc811_21700_100soc_inert()
    )

    species_moles = dataset.species_moles(
        temperature=298.15,
        pressure=101325.0,
    )

    inventory = ElementInventory.from_species_formulas(species_moles)

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

    assert inventory.moles.get(
        "N",
        0.0,
    ) == pytest.approx(
        0.0,
        abs=1.0e-12,
    )

    backend = CanteraEquilibriumBackend()

    equilibrium_state = backend.equilibrate_tv(
        temperature=1000.0,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    assert equilibrium_state.pressure > 0.0