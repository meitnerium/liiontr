from liiontr.chemistry import ReactionNetworkBackend
from liiontr.kinetics import Arrhenius, AutocatalyticProgress
from liiontr.library import cell_21700_generic
from liiontr.problems import ThermalProblem
from liiontr.reactions import Reaction, ReactionNetwork
from liiontr.solver import ScipySolver


def test_autocatalytic_reaction_is_integrated():
    cell = cell_21700_generic()

    network = ReactionNetwork()

    network.add(
        Reaction(
            name="Synthetic autocatalytic reaction",
            kinetics=Arrhenius(
                activation_energy=80000.0,
                pre_exponential_factor=1.0e5,
            ),
            enthalpy=500000.0,
            progress_model=AutocatalyticProgress(
                autocatalytic_order=1.0,
                remaining_order=1.0,
            ),
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
        initial_conversions=[0.10],
        duration=10.0,
    )

    results = ScipySolver().solve(problem)

    conversion = results.get("conversion_0")

    assert conversion[0] == 0.10
    assert conversion[-1] > conversion[0]
    assert conversion.min() >= 0.0
    assert conversion.max() <= 1.0
