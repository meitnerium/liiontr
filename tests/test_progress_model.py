import pytest

from liiontr.kinetics.progress import PowerLawProgress


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
