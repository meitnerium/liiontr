import pytest

from liiontr.kinetics.progress import (
    AutocatalyticProgress,
    PowerLawProgress,
)


def test_power_law_progress_at_zero_conversion():
    model = PowerLawProgress(
        order=1.0,
    )

    factor = model.factor(
        conversion=0.0,
    )

    assert factor == pytest.approx(1.0)


def test_power_law_progress_at_half_conversion():
    model = PowerLawProgress(
        order=2.0,
    )

    factor = model.factor(
        conversion=0.5,
    )

    assert factor == pytest.approx(0.25)


def test_power_law_progress_at_complete_conversion():
    model = PowerLawProgress(
        order=1.0,
    )

    factor = model.factor(
        conversion=1.0,
    )

    assert factor == pytest.approx(0.0)


@pytest.mark.parametrize(
    "order",
    [
        0.0,
        -1.0,
    ],
)
def test_power_law_progress_order_must_be_positive(
    order: float,
):
    with pytest.raises(
        ValueError,
        match="Progress model order must be greater than zero",
    ):
        PowerLawProgress(
            order=order,
        )


def test_autocatalytic_progress_at_half_conversion():
    model = AutocatalyticProgress(
        autocatalytic_order=1.0,
        remaining_order=1.0,
    )

    factor = model.factor(
        conversion=0.5,
    )

    assert factor == pytest.approx(0.25)


def test_autocatalytic_progress_is_zero_at_zero_conversion():
    model = AutocatalyticProgress(
        autocatalytic_order=1.0,
        remaining_order=1.0,
    )

    factor = model.factor(
        conversion=0.0,
    )

    assert factor == pytest.approx(0.0)


def test_autocatalytic_progress_is_zero_at_complete_conversion():
    model = AutocatalyticProgress(
        autocatalytic_order=1.0,
        remaining_order=1.0,
    )

    factor = model.factor(
        conversion=1.0,
    )

    assert factor == pytest.approx(0.0)


@pytest.mark.parametrize(
    ("autocatalytic_order", "remaining_order"),
    [
        (0.0, 1.0),
        (-1.0, 1.0),
        (1.0, 0.0),
        (1.0, -1.0),
    ],
)
def test_autocatalytic_orders_must_be_positive(
    autocatalytic_order: float,
    remaining_order: float,
):
    with pytest.raises(
        ValueError,
        match="Autocatalytic progress orders must be greater than zero",
    ):
        AutocatalyticProgress(
            autocatalytic_order=autocatalytic_order,
            remaining_order=remaining_order,
        )
