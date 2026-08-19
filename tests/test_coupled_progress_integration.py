import pytest

from liiontr.chemistry import ReactionNetworkBackend
from liiontr.kinetics import (
    Arrhenius,
    ExponentialInhibitionProgress,
    PowerLawProgress,
    ThresholdProgress,
)
from liiontr.library import cell_21700_generic
from liiontr.problems import ThermalProblem
from liiontr.reactions import (
    LinearConversionVariable,
    Reaction,
    ReactionNetwork,
)
from liiontr.solver import ScipySolver


def test_coupled_progress_models_are_integrated():
    cell = cell_21700_generic()

    sei = Reaction(
        name="SEI decomposition",
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=200000.0,
    )

    anode = Reaction(
        name="Anode-electrolyte",
        kinetics=Arrhenius(
            activation_energy=90000.0,
            pre_exponential_factor=2.0e5,
        ),
        enthalpy=500000.0,
        progress_model=ThresholdProgress(
            progress_model=ExponentialInhibitionProgress(
                progress_model=PowerLawProgress(
                    order=1.0,
                ),
                variable_name="sei_thickness_ratio",
            ),
            reaction_name="SEI decomposition",
            remaining_below=0.10,
        ),
    )

    network = ReactionNetwork(
        reactions=[
            sei,
            anode,
        ],
        context_variables={
            "sei_thickness_ratio": LinearConversionVariable(
                reaction_name="Anode-electrolyte",
                reference_conversion=0.0,
                reference_value=0.0,
                slope=2.0,
            )
        },
    )

    backend = ReactionNetworkBackend(
        reaction_network=network,
        cell=cell,
    )

    problem = ThermalProblem(
        cell=cell,
        chemistry_backend=backend,
        initial_temperature=500.0,
        initial_conversions=[
            0.95,
            0.10,
        ],
        duration=1.0,
    )

    results = ScipySolver().solve(problem)

    sei_conversion = results.get("conversion_0")

    anode_conversion = results.get("conversion_1")

    assert sei_conversion[0] == pytest.approx(0.95)

    assert anode_conversion[0] == pytest.approx(0.10)

    assert sei_conversion[-1] > sei_conversion[0]

    assert anode_conversion[-1] > anode_conversion[0]

    assert sei_conversion.max() <= 1.0
    assert anode_conversion.max() <= 1.0
