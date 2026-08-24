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


def test_gas_generation_model_combines_reactions():
    reaction_1 = Reaction(
        name="Reaction 1",
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=200000.0,
        mass_fraction=0.10,
    )

    reaction_2 = Reaction(
        name="Reaction 2",
        kinetics=Arrhenius(
            activation_energy=90000.0,
            pre_exponential_factor=2.0e5,
        ),
        enthalpy=300000.0,
        mass_fraction=0.20,
    )

    network = ReactionNetwork(
        reactions=[
            reaction_1,
            reaction_2,
        ]
    )

    model = GasGenerationModel(
        reaction_network=network,
        gas_yields=[
            ReactionGasYield(
                reaction_name="Reaction 1",
                species_yields={
                    "CO2": 2.0,
                    "H2": 0.5,
                },
            ),
            ReactionGasYield(
                reaction_name="Reaction 2",
                species_yields={
                    "CO2": 1.0,
                    "CO": 0.25,
                },
            ),
        ],
    )

    temperature = 500.0
    conversions = [
        0.25,
        0.50,
    ]
    cell_mass = 0.060

    rates = model.generation_rates(
        temperature=temperature,
        conversions=conversions,
        cell_mass=cell_mass,
    )

    progress_rates = network.progress_rates(
        temperature=temperature,
        conversions=conversions,
    )

    expected_co2 = (
        2.0 * cell_mass * 0.10 * progress_rates[0]
        + 1.0 * cell_mass * 0.20 * progress_rates[1]
    )

    expected_h2 = 0.5 * cell_mass * 0.10 * progress_rates[0]

    expected_co = 0.25 * cell_mass * 0.20 * progress_rates[1]

    assert rates["CO2"] == pytest.approx(expected_co2)

    assert rates["H2"] == pytest.approx(expected_h2)

    assert rates["CO"] == pytest.approx(expected_co)


def test_reaction_without_gas_yield_produces_no_gas():
    reaction = Reaction(
        name="Reaction",
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=200000.0,
        mass_fraction=0.10,
    )

    network = ReactionNetwork(
        reactions=[
            reaction,
        ]
    )

    model = GasGenerationModel(
        reaction_network=network,
        gas_yields=[],
    )

    rates = model.generation_rates(
        temperature=500.0,
        conversions=[0.0],
        cell_mass=0.060,
    )

    assert rates == {}


def test_unknown_reaction_yield_is_rejected():
    network = ReactionNetwork(reactions=[])

    with pytest.raises(
        ValueError,
        match="Unknown reaction",
    ):
        GasGenerationModel(
            reaction_network=network,
            gas_yields=[
                ReactionGasYield(
                    reaction_name="Missing reaction",
                    species_yields={
                        "CO2": 1.0,
                    },
                ),
            ],
        )


def test_duplicate_reaction_yields_are_rejected():
    reaction = Reaction(
        name="Reaction",
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=200000.0,
    )

    network = ReactionNetwork(
        reactions=[
            reaction,
        ]
    )

    with pytest.raises(
        ValueError,
        match="Duplicate gas yield",
    ):
        GasGenerationModel(
            reaction_network=network,
            gas_yields=[
                ReactionGasYield(
                    reaction_name="Reaction",
                    species_yields={
                        "CO2": 1.0,
                    },
                ),
                ReactionGasYield(
                    reaction_name="Reaction",
                    species_yields={
                        "H2": 1.0,
                    },
                ),
            ],
        )

def test_gas_generation_model_has_deterministic_species_order():
    reaction_1 = Reaction(
        name="Reaction 1",
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=200000.0,
    )

    reaction_2 = Reaction(
        name="Reaction 2",
        kinetics=Arrhenius(
            activation_energy=90000.0,
            pre_exponential_factor=2.0e5,
        ),
        enthalpy=300000.0,
    )

    network = ReactionNetwork(
        reactions=[
            reaction_1,
            reaction_2,
        ]
    )

    model = GasGenerationModel(
        reaction_network=network,
        gas_yields=[
            ReactionGasYield(
                reaction_name="Reaction 1",
                species_yields={
                    "CO2": 1.0,
                    "H2": 2.0,
                },
            ),
            ReactionGasYield(
                reaction_name="Reaction 2",
                species_yields={
                    "H2": 1.0,
                    "CO": 3.0,
                },
            ),
        ],
    )

    assert model.species_names == [
        "CO2",
        "H2",
        "CO",
    ]