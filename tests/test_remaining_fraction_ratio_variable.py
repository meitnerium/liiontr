import pytest

from liiontr.reactions import (
    ReactionContext,
    RemainingFractionRatioVariable,
)


def test_remaining_fraction_ratio_at_reference_state():
    variable = RemainingFractionRatioVariable(
        reaction_name="SEI decomposition",
        reference_remaining_fraction=0.15,
    )

    context = ReactionContext(
        conversions={
            "SEI decomposition": 0.85,
        }
    )

    value = variable.evaluate(context)

    assert value == pytest.approx(1.0)


def test_remaining_fraction_ratio_decreases_with_conversion():
    variable = RemainingFractionRatioVariable(
        reaction_name="SEI decomposition",
        reference_remaining_fraction=0.15,
    )

    context = ReactionContext(
        conversions={
            "SEI decomposition": 0.90,
        }
    )

    value = variable.evaluate(context)

    expected = 0.10 / 0.15

    assert value == pytest.approx(expected)


def test_remaining_fraction_ratio_is_zero_at_complete_conversion():
    variable = RemainingFractionRatioVariable(
        reaction_name="SEI decomposition",
        reference_remaining_fraction=0.15,
    )

    context = ReactionContext(
        conversions={
            "SEI decomposition": 1.0,
        }
    )

    value = variable.evaluate(context)

    assert value == pytest.approx(0.0)


@pytest.mark.parametrize(
    "reference_remaining_fraction",
    [
        0.0,
        -0.1,
        1.1,
    ],
)
def test_reference_remaining_fraction_must_be_valid(
    reference_remaining_fraction: float,
):
    with pytest.raises(
        ValueError,
        match=("Reference remaining fraction must be greater than 0 and at most 1"),
    ):
        RemainingFractionRatioVariable(
            reaction_name="SEI decomposition",
            reference_remaining_fraction=reference_remaining_fraction,
        )
