from liiontr.kinetics import Arrhenius
from liiontr.reactions import (
    Reaction,
    ReactionNetwork,
)


def test_single_reaction():
    reaction = Reaction(
        name="SEI",
        kinetics=Arrhenius(
            activation_energy=120000,
            pre_exponential_factor=1e8,
        ),
        enthalpy=500000,
    )

    q = reaction.heat_generation(500.0)

    assert q > 0


def test_network():
    network = ReactionNetwork()

    network.add(
        Reaction(
            name="SEI",
            kinetics=Arrhenius(
                activation_energy=120000,
                pre_exponential_factor=1e8,
            ),
            enthalpy=500000,
        )
    )

    network.add(
        Reaction(
            name="Electrolyte",
            kinetics=Arrhenius(
                activation_energy=100000,
                pre_exponential_factor=1e7,
            ),
            enthalpy=700000,
        )
    )

    assert network.heat_generation(500) > 0
