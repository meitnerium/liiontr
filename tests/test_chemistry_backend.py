from liiontr.chemistry import ReactionNetworkBackend
from liiontr.kinetics import Arrhenius
from liiontr.library import cell_21700_generic
from liiontr.reactions import Reaction, ReactionNetwork


def test_reaction_network_backend():
    network = ReactionNetwork()

    network.add(
        Reaction(
            name="SEI",
            kinetics=Arrhenius(
                activation_energy=120000.0,
                pre_exponential_factor=1e5,
            ),
            enthalpy=500000.0,
        )
    )

    backend = ReactionNetworkBackend(
        reaction_network=network,
        cell=cell_21700_generic(),
    )

    assert backend.heat_generation(500.0) > 0.0


def test_reaction_network_backend_returns_cell_power():
    cell = cell_21700_generic()

    network = ReactionNetwork()

    network.add(
        Reaction(
            name="SEI",
            kinetics=Arrhenius(
                activation_energy=120000.0,
                pre_exponential_factor=1e5,
            ),
            enthalpy=500000.0,
        )
    )

    backend = ReactionNetworkBackend(
        reaction_network=network,
        cell=cell,
    )

    power = backend.heat_generation(500.0)

    assert power > 0.0
