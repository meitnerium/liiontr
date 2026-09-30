import pytest

from liiontr.chemistry.elements import (
    ElementInventory,
    ReactionElementYield,
)


def test_element_inventory_returns_moles():
    inventory = ElementInventory(
        moles={
            "C": 1.0,
            "H": 4.0,
            "O": 2.0,
        }
    )

    assert inventory.moles_of("C") == pytest.approx(
        1.0
    )

    assert inventory.moles_of("H") == pytest.approx(
        4.0
    )

    assert inventory.moles_of("O") == pytest.approx(
        2.0
    )


def test_element_inventory_missing_element_is_zero():
    inventory = ElementInventory(
        moles={
            "C": 1.0,
        }
    )

    assert inventory.moles_of("N") == pytest.approx(
        0.0
    )


def test_element_inventory_total_atom_moles():
    inventory = ElementInventory(
        moles={
            "C": 1.0,
            "H": 4.0,
            "O": 2.0,
        }
    )

    assert inventory.total_atom_moles == pytest.approx(
        7.0
    )


def test_element_inventory_rejects_negative_moles():
    with pytest.raises(
        ValueError,
        match="must not be negative",
    ):
        ElementInventory(
            moles={
                "C": -1.0,
            }
        )


def test_element_inventory_rejects_nonfinite_moles():
    with pytest.raises(
        ValueError,
        match="must be finite",
    ):
        ElementInventory(
            moles={
                "C": float("inf"),
            }
        )


def test_reaction_element_yield_returns_yield():
    element_yield = ReactionElementYield(
        reaction_name="Reaction A",
        element_yields={
            "C": 2.0,
            "H": 5.0,
        },
    )

    assert element_yield.yield_of(
        "C"
    ) == pytest.approx(2.0)

    assert element_yield.yield_of(
        "H"
    ) == pytest.approx(5.0)

    assert element_yield.yield_of(
        "O"
    ) == pytest.approx(0.0)


def test_reaction_element_yield_rejects_negative_yield():
    with pytest.raises(
        ValueError,
        match="must not be negative",
    ):
        ReactionElementYield(
            reaction_name="Reaction A",
            element_yields={
                "C": -1.0,
            },
        )