import pytest

from liiontr.gases import (
    GasInventory,
    GasSpecies,
)


def test_gas_species_properties():
    species = GasSpecies(
        name="CO2",
        molar_mass=44.0095e-3,
    )

    assert species.name == "CO2"
    assert species.molar_mass == pytest.approx(44.0095e-3)


def test_gas_species_requires_positive_molar_mass():
    with pytest.raises(
        ValueError,
        match="Molar mass",
    ):
        GasSpecies(
            name="CO2",
            molar_mass=0.0,
        )


def test_empty_inventory_has_zero_total_moles():
    inventory = GasInventory(
        species=[
            GasSpecies(
                name="CO2",
                molar_mass=44.0095e-3,
            ),
            GasSpecies(
                name="H2",
                molar_mass=2.01588e-3,
            ),
        ]
    )

    assert inventory.total_moles == pytest.approx(0.0)


def test_inventory_tracks_moles_by_species():
    inventory = GasInventory(
        species=[
            GasSpecies(
                name="CO2",
                molar_mass=44.0095e-3,
            ),
            GasSpecies(
                name="H2",
                molar_mass=2.01588e-3,
            ),
        ],
        moles={
            "CO2": 2.0e-3,
            "H2": 3.0e-3,
        },
    )

    assert inventory.moles_of("CO2") == pytest.approx(2.0e-3)

    assert inventory.moles_of("H2") == pytest.approx(3.0e-3)

    assert inventory.total_moles == pytest.approx(5.0e-3)


def test_inventory_computes_total_gas_mass():
    inventory = GasInventory(
        species=[
            GasSpecies(
                name="CO2",
                molar_mass=44.0e-3,
            ),
            GasSpecies(
                name="H2",
                molar_mass=2.0e-3,
            ),
        ],
        moles={
            "CO2": 2.0,
            "H2": 3.0,
        },
    )

    expected_mass = 2.0 * 44.0e-3 + 3.0 * 2.0e-3

    assert inventory.total_mass == pytest.approx(expected_mass)


def test_inventory_rejects_unknown_species():
    inventory = GasInventory(
        species=[
            GasSpecies(
                name="CO2",
                molar_mass=44.0095e-3,
            ),
        ]
    )

    with pytest.raises(
        KeyError,
        match="Unknown gas species",
    ):
        inventory.moles_of("H2")


def test_inventory_rejects_negative_moles():
    with pytest.raises(
        ValueError,
        match="Gas moles",
    ):
        GasInventory(
            species=[
                GasSpecies(
                    name="CO2",
                    molar_mass=44.0095e-3,
                ),
            ],
            moles={
                "CO2": -1.0,
            },
        )
