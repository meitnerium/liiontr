import pytest

from liiontr.kinetics import (
    Arrhenius,
    ExponentialInhibitionProgress,
    PowerLawProgress,
)
from liiontr.reactions import (
    LinearConversionVariable,
    Reaction,
    ReactionNetwork,
)


def test_network_builds_derived_context_variables():
    reaction = Reaction(
        name="Anode-electrolyte",
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=500000.0,
    )

    network = ReactionNetwork(
        reactions=[
            reaction,
        ],
        context_variables={
            "sei_thickness_ratio": LinearConversionVariable(
                reaction_name="Anode-electrolyte",
                reference_conversion=0.0,
                reference_value=1.0,
                slope=2.0,
            )
        },
    )

    context = network.context(
        conversions=[
            0.25,
        ]
    )

    assert context.variable("sei_thickness_ratio") == pytest.approx(1.5)


def test_derived_variable_is_available_to_progress_model():
    reaction = Reaction(
        name="Anode-electrolyte",
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=500000.0,
        progress_model=ExponentialInhibitionProgress(
            progress_model=PowerLawProgress(
                order=1.0,
            ),
            variable_name="sei_thickness_ratio",
        ),
    )

    network = ReactionNetwork(
        reactions=[
            reaction,
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

    rates = network.progress_rates(
        temperature=500.0,
        conversions=[
            0.25,
        ],
    )

    expected = reaction.rate(500.0) * 0.75 * 0.6065306597

    assert rates[0] == pytest.approx(expected)
