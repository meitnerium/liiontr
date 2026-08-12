import pytest

from liiontr.kinetics import Arrhenius
from liiontr.reactions import Reaction


def make_reaction(reaction_order: float) -> Reaction:
    return Reaction(
        name="Synthetic reaction",
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=500000.0,
        reaction_order=reaction_order,
    )


def test_reaction_order_changes_progress_rate():
    first_order = make_reaction(
        reaction_order=1.0,
    )

    second_order = make_reaction(
        reaction_order=2.0,
    )

    rate_first_order = first_order.progress_rate(
        temperature=500.0,
        conversion=0.5,
    )

    rate_second_order = second_order.progress_rate(
        temperature=500.0,
        conversion=0.5,
    )

    assert rate_second_order < rate_first_order


def test_default_reaction_order_is_one():
    reaction = Reaction(
        name="Synthetic reaction",
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=500000.0,
    )

    rate = reaction.progress_rate(
        temperature=500.0,
        conversion=0.5,
    )

    expected = reaction.rate(500.0) * 0.5

    assert rate == pytest.approx(expected)

    def test_reaction_order_must_be_positive():
        with pytest.raises(
            ValueError,
            match="Reaction order must be greater than zero",
        ):
            make_reaction(
                reaction_order=0.0,
            )
