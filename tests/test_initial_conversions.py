import pytest

from liiontr.chemistry import ReactionNetworkBackend
from liiontr.kinetics import Arrhenius
from liiontr.library import cell_21700_generic
from liiontr.problems import ThermalProblem
from liiontr.reactions import Reaction, ReactionNetwork
from liiontr.solver import ScipySolver


def make_backend(number_of_reactions: int) -> ReactionNetworkBackend:
    cell = cell_21700_generic()

    network = ReactionNetwork()

    for index in range(number_of_reactions):
        network.add(
            Reaction(
                name=f"Synthetic reaction {index}",
                kinetics=Arrhenius(
                    activation_energy=80000.0,
                    pre_exponential_factor=1.0e5,
                ),
                enthalpy=500000.0,
            )
        )

    return ReactionNetworkBackend(
        reaction_network=network,
        cell=cell,
    )


def test_reaction_initial_conversion():
    backend = make_backend(1)

    problem = ThermalProblem(
        cell=backend.cell,
        chemistry_backend=backend,
        initial_temperature=500.0,
        initial_conversions=[0.25],
        duration=1.0e-6,
    )

    results = ScipySolver().solve(problem)

    conversion = results.get("conversion_0")

    assert conversion[0] == 0.25


def test_initial_conversion_count_must_match_reactions():
    backend = make_backend(2)

    problem = ThermalProblem(
        cell=backend.cell,
        chemistry_backend=backend,
        initial_temperature=500.0,
        initial_conversions=[0.25],
        duration=1.0,
    )

    with pytest.raises(
        ValueError,
        match="Number of initial conversions must match number of reactions",
    ):
        ScipySolver().solve(problem)


@pytest.mark.parametrize(
    "conversion",
    [
        -0.1,
        1.1,
    ],
)
def test_initial_conversion_must_be_between_zero_and_one(
    conversion: float,
):
    backend = make_backend(1)

    problem = ThermalProblem(
        cell=backend.cell,
        chemistry_backend=backend,
        initial_temperature=500.0,
        initial_conversions=[conversion],
        duration=1.0,
    )

    with pytest.raises(
        ValueError,
        match="Initial conversions must be between 0 and 1",
    ):
        ScipySolver().solve(problem)
