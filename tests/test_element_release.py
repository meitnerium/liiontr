"""Tests for delayed elemental release models."""

import pytest

from liiontr.chemistry.element_release import (
    FirstOrderElementReleaseModel,
)
from liiontr.chemistry.elements import (
    ElementInventory,
)


def test_first_order_element_release_rates():
    """Release stored elements using the configured time constant."""
    model = FirstOrderElementReleaseModel(
        time_constant=2.0,
    )

    inventory = ElementInventory(
        moles={
            "C": 0.4,
            "H": 0.6,
            "O": 0.2,
        }
    )

    rates = model.release_rates(
        inventory
    )

    assert rates["C"] == pytest.approx(
        0.2
    )

    assert rates["H"] == pytest.approx(
        0.3
    )

    assert rates["O"] == pytest.approx(
        0.1
    )


def test_first_order_element_release_zero_inventory():
    """Return zero release for an empty elemental inventory."""
    model = FirstOrderElementReleaseModel(
        time_constant=1.0,
    )

    inventory = ElementInventory(
        moles={}
    )

    assert model.release_rates(
        inventory
    ) == {}


@pytest.mark.parametrize(
    "time_constant",
    [
        0.0,
        -1.0,
        float("inf"),
        float("nan"),
    ],
)
def test_first_order_element_release_rejects_invalid_time_constant(
    time_constant,
):
    """Reject invalid release time constants."""
    with pytest.raises(
        ValueError,
        match="Release time constant",
    ):
        FirstOrderElementReleaseModel(
            time_constant=time_constant,
        )