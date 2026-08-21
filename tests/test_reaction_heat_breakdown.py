import pytest

from liiontr.kinetics import Arrhenius
from liiontr.reactions import (
    Reaction,
    ReactionNetwork,
)


def test_heat_generation_by_reaction():
    reaction_1 = Reaction(
        name="Reaction 1",
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=200000.0,
        mass_fraction=0.10,
    )

    reaction_2 = Reaction(
        name="Reaction 2",
        kinetics=Arrhenius(
            activation_energy=90000.0,
            pre_exponential_factor=2.0e5,
        ),
        enthalpy=300000.0,
        mass_fraction=0.20,
    )

    network = ReactionNetwork(
        reactions=[
            reaction_1,
            reaction_2,
        ]
    )

    conversions = [
        0.25,
        0.50,
    ]

    heat = network.heat_generation_by_reaction(
        temperature=500.0,
        conversions=conversions,
    )

    assert set(heat) == {
        "Reaction 1",
        "Reaction 2",
    }

    assert heat["Reaction 1"] > 0.0
    assert heat["Reaction 2"] > 0.0

    total = network.heat_generation(
        temperature=500.0,
        conversions=conversions,
    )

    assert sum(heat.values()) == pytest.approx(total)
