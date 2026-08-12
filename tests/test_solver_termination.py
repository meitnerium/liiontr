from liiontr.chemistry import ReactionNetworkBackend
from liiontr.kinetics import Arrhenius
from liiontr.library import cell_21700_generic
from liiontr.problems import ThermalProblem
from liiontr.reactions import Reaction, ReactionNetwork
from liiontr.solver import ScipySolver


def test_solver_stops_at_maximum_temperature():
    cell = cell_21700_generic()

    network = ReactionNetwork()

    network.add(
        Reaction(
            name="Synthetic runaway reaction",
            kinetics=Arrhenius(
                activation_energy=80000.0,
                pre_exponential_factor=1.0e8,
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
        maximum_temperature=1000.0,
    )

    results = ScipySolver().solve(problem)

    temperature = results.get("temperature")

    assert temperature[-1] <= 1000.0 + 1.0e-6
    assert temperature[-1] >= 999.0
