from liiontr.kinetics import Arrhenius
from liiontr.reactions import Reaction, ReactionNetwork


def test_arrhenius_heat_generation():
    reaction = Reaction(
        name="SEI decomposition",
        kinetics=Arrhenius(
            activation_energy=120000.0,
            pre_exponential_factor=1e5,
        ),
        enthalpy=500000.0,
    )

    network = ReactionNetwork()

    network.add(reaction)

    q = network.heat_generation(500.0)

    assert q > 0
