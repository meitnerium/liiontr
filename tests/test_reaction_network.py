from liiontr.kinetics import Arrhenius
from liiontr.reactions import (
    Reaction,
    ReactionNetwork,
)
import pytest


def test_single_reaction():
    reaction = Reaction(
        name="SEI",
        kinetics=Arrhenius(
            activation_energy=120000,
            pre_exponential_factor=1e8,
        ),
        enthalpy=500000,
    )

    q = reaction.heat_generation(500.0)

    assert q > 0


def test_network_progress_rates():
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

    network.add(
        Reaction(
            name="Electrolyte",
            kinetics=Arrhenius(
                activation_energy=100000.0,
                pre_exponential_factor=1e4,
            ),
            enthalpy=300000.0,
        )
    )

    rates = network.progress_rates(
        temperature=500.0,
        conversions=[0.0, 0.5],
    )

    assert len(rates) == 2
    assert rates[0] > 0.0
    assert rates[1] > 0.0


def test_network():
    network = ReactionNetwork()

    network.add(
        Reaction(
            name="SEI",
            kinetics=Arrhenius(
                activation_energy=120000,
                pre_exponential_factor=1e8,
            ),
            enthalpy=500000,
        )
    )

    network.add(
        Reaction(
            name="Electrolyte",
            kinetics=Arrhenius(
                activation_energy=100000,
                pre_exponential_factor=1e7,
            ),
            enthalpy=700000,
        )
    )

    assert network.heat_generation(500) > 0


def test_heat_generation_decreases_with_conversion():
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

    q_initial = network.heat_generation(
        temperature=500.0,
        conversions=[0.0],
    )

    q_half = network.heat_generation(
        temperature=500.0,
        conversions=[0.5],
    )

    q_complete = network.heat_generation(
        temperature=500.0,
        conversions=[1.0],
    )

    assert q_initial > q_half > q_complete
    assert q_complete == 0.0


def test_reaction_names_must_be_unique():
    network = ReactionNetwork()

    reaction_1 = Reaction(
        name="SEI decomposition",
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=200000.0,
    )

    reaction_2 = Reaction(
        name="SEI decomposition",
        kinetics=Arrhenius(
            activation_energy=90000.0,
            pre_exponential_factor=2.0e5,
        ),
        enthalpy=300000.0,
    )

    network.add(reaction_1)

    with pytest.raises(
        ValueError,
        match="Reaction names must be unique",
    ):
        network.add(reaction_2)


def test_duplicate_reaction_names_are_rejected_at_construction():
    reaction_1 = Reaction(
        name="SEI decomposition",
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=200000.0,
    )

    reaction_2 = Reaction(
        name="SEI decomposition",
        kinetics=Arrhenius(
            activation_energy=90000.0,
            pre_exponential_factor=2.0e5,
        ),
        enthalpy=300000.0,
    )

    with pytest.raises(
        ValueError,
        match="Reaction names must be unique",
    ):
        ReactionNetwork(
            reactions=[
                reaction_1,
                reaction_2,
            ]
        )
