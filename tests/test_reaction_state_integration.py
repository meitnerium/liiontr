from liiontr.chemistry import ReactionNetworkBackend
from liiontr.kinetics import Arrhenius
from liiontr.library import cell_21700_generic
from liiontr.problems import ThermalProblem
from liiontr.reactions import Reaction, ReactionNetwork
from liiontr.solver import ScipySolver


def test_reaction_conversion_is_integrated():
    cell = cell_21700_generic()

    network = ReactionNetwork()

    network.add(
        Reaction(
            name="Synthetic reaction",
            kinetics=Arrhenius(
                activation_energy=80000.0,
                pre_exponential_factor=1e5,
            ),
            enthalpy=500000.0,
        )
    )

    backend = ReactionNetworkBackend(
        reaction_network=network,
        cell=cell,
    )

    problem = ThermalProblem(
        cell=cell,
        chemistry_backend=backend,
        initial_temperature=500.0,
        duration=10.0,
    )

    results = ScipySolver().solve(problem)

    conversion = results.get("conversion_0")

    assert conversion[0] == 0.0
    assert conversion[-1] > conversion[0]
    assert conversion[-1] <= 1.0
