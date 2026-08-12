from liiontr.kinetics import Arrhenius
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
