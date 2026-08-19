import pytest

from liiontr.reactions import (
    LinearConversionVariable,
    ReactionContext,
)


def test_linear_conversion_variable_at_reference_state():
    variable = LinearConversionVariable(
        reaction_name="Anode-electrolyte",
        reference_conversion=0.25,
        reference_value=1.0,
        slope=2.0,
    )

    context = ReactionContext(
        conversions={
            "Anode-electrolyte": 0.25,
        }
    )

    value = variable.evaluate(context)

    assert value == pytest.approx(1.0)


def test_linear_conversion_variable_changes_with_conversion():
    variable = LinearConversionVariable(
        reaction_name="Anode-electrolyte",
        reference_conversion=0.25,
        reference_value=1.0,
        slope=2.0,
    )

    context = ReactionContext(
        conversions={
            "Anode-electrolyte": 0.35,
        }
    )

    value = variable.evaluate(context)

    expected = 1.0 + 2.0 * (0.35 - 0.25)

    assert value == pytest.approx(expected)


def test_linear_conversion_variable_requires_reaction():
    variable = LinearConversionVariable(
        reaction_name="Anode-electrolyte",
        reference_conversion=0.25,
        reference_value=1.0,
        slope=2.0,
    )

    context = ReactionContext()

    with pytest.raises(
        KeyError,
        match="Unknown reaction conversion",
    ):
        variable.evaluate(context)


@pytest.mark.parametrize(
    "reference_conversion",
    [
        -0.1,
        1.1,
    ],
)
def test_reference_conversion_must_be_bounded(
    reference_conversion: float,
):
    with pytest.raises(
        ValueError,
        match="Reference conversion must be between 0 and 1",
    ):
        LinearConversionVariable(
            reaction_name="Anode-electrolyte",
            reference_conversion=reference_conversion,
            reference_value=1.0,
            slope=2.0,
        )
