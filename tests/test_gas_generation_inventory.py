import pytest

from liiontr.gases import (
    GasGenerationModel,
    ReactionGasYield,
)
from liiontr.kinetics import Arrhenius
from liiontr.reactions import (
    Reaction,
    ReactionNetwork,
)


def test_generated_moles_follow_conversion_change():
    reaction = Reaction(
        name="Gas reaction",
        kinetics=Arrhenius(
            activation_energy=1.0,
            pre_exponential_factor=1.0,
        ),
        enthalpy=10000.0,
        mass_fraction=0.10,
    )

    network = ReactionNetwork(
        reactions=[
            reaction,
        ]
    )

    model = GasGenerationModel(
        reaction_network=network,
        gas_yields=[
            ReactionGasYield(
                reaction_name="Gas reaction",
                species_yields={
                    "CO2": 2.0,
                    "H2": 0.5,
                },
            )
        ],
    )

    generated = model.generated_moles(
        initial_conversions=[
            0.20,
        ],
        conversions=[
            0.70,
        ],
        cell_mass=0.060,
    )

    reacted_mass = 0.060 * 0.10 * (0.70 - 0.20)

    assert generated["CO2"] == pytest.approx(2.0 * reacted_mass)

    assert generated["H2"] == pytest.approx(0.5 * reacted_mass)


def test_zero_conversion_change_generates_no_gas():
    reaction = Reaction(
        name="Gas reaction",
        kinetics=Arrhenius(
            activation_energy=1.0,
            pre_exponential_factor=1.0,
        ),
        enthalpy=10000.0,
        mass_fraction=0.10,
    )

    network = ReactionNetwork(
        reactions=[
            reaction,
        ]
    )

    model = GasGenerationModel(
        reaction_network=network,
        gas_yields=[
            ReactionGasYield(
                reaction_name="Gas reaction",
                species_yields={
                    "CO2": 2.0,
                },
            )
        ],
    )

    generated = model.generated_moles(
        initial_conversions=[
            0.30,
        ],
        conversions=[
            0.30,
        ],
        cell_mass=0.060,
    )

    assert generated["CO2"] == pytest.approx(0.0)


def test_decreasing_conversion_is_rejected():
    reaction = Reaction(
        name="Gas reaction",
        kinetics=Arrhenius(
            activation_energy=1.0,
            pre_exponential_factor=1.0,
        ),
        enthalpy=10000.0,
        mass_fraction=0.10,
    )

    network = ReactionNetwork(
        reactions=[
            reaction,
        ]
    )

    model = GasGenerationModel(
        reaction_network=network,
        gas_yields=[
            ReactionGasYield(
                reaction_name="Gas reaction",
                species_yields={
                    "CO2": 2.0,
                },
            )
        ],
    )

    with pytest.raises(
        ValueError,
        match="conversion must not decrease",
    ):
        model.generated_moles(
            initial_conversions=[
                0.50,
            ],
            conversions=[
                0.40,
            ],
            cell_mass=0.060,
        )
