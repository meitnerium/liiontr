import pytest

from liiontr.kinetics import Arrhenius, PowerLawProgress
from liiontr.reactions import Reaction


def test_reaction_progress_rate_decreases_with_conversion():
    reaction = Reaction(
        name="SEI",
        kinetics=Arrhenius(
            activation_energy=120000.0,
            pre_exponential_factor=1e5,
        ),
        enthalpy=500000.0,
    )

    rate_initial = reaction.progress_rate(
        temperature=500.0,
        conversion=0.0,
    )

    rate_half = reaction.progress_rate(
        temperature=500.0,
        conversion=0.5,
    )

    rate_complete = reaction.progress_rate(
        temperature=500.0,
        conversion=1.0,
    )

    assert rate_initial > rate_half > rate_complete
    assert rate_complete == 0.0


def test_reaction_uses_explicit_progress_model():
    reaction = Reaction(
        name="Synthetic reaction",
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=500000.0,
        progress_model=PowerLawProgress(
            order=2.0,
        ),
    )

    rate = reaction.progress_rate(
        temperature=500.0,
        conversion=0.5,
    )

    expected = reaction.rate(500.0) * 0.25

    assert rate == pytest.approx(expected)
