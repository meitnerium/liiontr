import pytest

from liiontr.reactions import ReactionContext


def test_reaction_context_stores_conversions():
    context = ReactionContext(
        conversions={
            "SEI decomposition": 0.90,
            "Anode-electrolyte": 0.25,
        }
    )

    assert context.conversion("SEI decomposition") == pytest.approx(0.90)

    assert context.conversion("Anode-electrolyte") == pytest.approx(0.25)


def test_reaction_context_returns_remaining_fraction():
    context = ReactionContext(
        conversions={
            "SEI decomposition": 0.90,
        }
    )

    remaining = context.remaining_fraction("SEI decomposition")

    assert remaining == pytest.approx(0.10)


def test_reaction_context_stores_additional_variables():
    context = ReactionContext(
        variables={
            "sei_thickness_ratio": 0.033,
        }
    )

    assert context.variable("sei_thickness_ratio") == pytest.approx(0.033)


def test_missing_conversion_raises_key_error():
    context = ReactionContext()

    with pytest.raises(
        KeyError,
        match="Unknown reaction conversion",
    ):
        context.conversion("Unknown reaction")


def test_missing_variable_raises_key_error():
    context = ReactionContext()

    with pytest.raises(
        KeyError,
        match="Unknown reaction variable",
    ):
        context.variable("unknown")
