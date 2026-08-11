from liiontr.chemistry import ReactionNetworkBackend
from liiontr.kinetics import Arrhenius
from liiontr.library import cell_21700_generic
from liiontr.problems import ThermalProblem
from liiontr.reactions import Reaction, ReactionNetwork
from liiontr.solver import ScipySolver


def test_thermal_runaway_with_reaction():
    cell = cell_21700_generic()

    network = ReactionNetwork()

    network.add(
        Reaction(
            name="Synthetic exothermic reaction",
            kinetics=Arrhenius(
                activation_energy=80000.0,
                pre_exponential_factor=1e8,
            ),
            enthalpy=5.0e6,
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
        duration=100.0,
    )

    results = ScipySolver().solve(problem)

    temperatures = results.get("temperature")

    assert temperatures[0] == 500.0
    assert temperatures[-1] > temperatures[0]