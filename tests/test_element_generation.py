import pytest

from liiontr.chemistry import (
    ElementGenerationModel,
    ReactionElementYield,
)
from liiontr.kinetics import Arrhenius
from liiontr.reactions import (
    Reaction,
    ReactionNetwork,
)


def build_network() -> ReactionNetwork:
    reaction = Reaction(
        name="Reaction A",
        kinetics=Arrhenius(
            activation_energy=1.0,
            pre_exponential_factor=2.0,
        ),
        enthalpy=10000.0,
        mass_fraction=0.10,
    )

    return ReactionNetwork(
        reactions=[
            reaction,
        ]
    )


def test_element_generation_names():
    network = build_network()

    model = ElementGenerationModel(
        reaction_network=network,
        element_yields=[
            ReactionElementYield(
                reaction_name="Reaction A",
                element_yields={
                    "C": 2.0,
                    "H": 4.0,
                    "O": 1.0,
                },
            )
        ],
    )

    assert model.element_names == [
        "C",
        "H",
        "O",
    ]


def test_element_generation_rates():
    network = build_network()

    model = ElementGenerationModel(
        reaction_network=network,
        element_yields=[
            ReactionElementYield(
                reaction_name="Reaction A",
                element_yields={
                    "C": 2.0,
                    "H": 4.0,
                },
            )
        ],
    )

    temperature = 400.0
    conversions = [0.25]
    cell_mass = 0.060

    progress_rate = (
        network.progress_rates(
            temperature=temperature,
            conversions=conversions,
        )[0]
    )

    reacted_mass_rate = (
        cell_mass
        * 0.10
        * progress_rate
    )

    rates = model.generation_rates(
        temperature=temperature,
        conversions=conversions,
        cell_mass=cell_mass,
    )

    assert rates["C"] == pytest.approx(
        2.0 * reacted_mass_rate
    )

    assert rates["H"] == pytest.approx(
        4.0 * reacted_mass_rate
    )


def test_element_generation_cumulative_moles():
    network = build_network()

    model = ElementGenerationModel(
        reaction_network=network,
        element_yields=[
            ReactionElementYield(
                reaction_name="Reaction A",
                element_yields={
                    "C": 2.0,
                    "H": 4.0,
                },
            )
        ],
    )

    amounts = model.generated_moles(
        initial_conversions=[
            0.20,
        ],
        conversions=[
            0.70,
        ],
        cell_mass=0.060,
    )

    reacted_mass = (
        0.060
        * 0.10
        * (
            0.70
            - 0.20
        )
    )

    assert amounts["C"] == pytest.approx(
        2.0 * reacted_mass
    )

    assert amounts["H"] == pytest.approx(
        4.0 * reacted_mass
    )


def test_element_generation_inventory():
    network = build_network()

    model = ElementGenerationModel(
        reaction_network=network,
        element_yields=[
            ReactionElementYield(
                reaction_name="Reaction A",
                element_yields={
                    "C": 2.0,
                    "O": 1.0,
                },
            )
        ],
    )

    inventory = model.generated_inventory(
        initial_conversions=[
            0.0,
        ],
        conversions=[
            1.0,
        ],
        cell_mass=0.060,
    )

    assert inventory.moles_of(
        "C"
    ) == pytest.approx(
        0.060
        * 0.10
        * 2.0
    )

    assert inventory.moles_of(
        "O"
    ) == pytest.approx(
        0.060
        * 0.10
        * 1.0
    )


def test_element_generation_rejects_unknown_reaction():
    network = build_network()

    with pytest.raises(
        ValueError,
        match="unknown reaction",
    ):
        ElementGenerationModel(
            reaction_network=network,
            element_yields=[
                ReactionElementYield(
                    reaction_name="Unknown reaction",
                    element_yields={
                        "C": 1.0,
                    },
                )
            ],
        )


def test_element_generation_rejects_duplicate_reaction_yields():
    network = build_network()

    with pytest.raises(
        ValueError,
        match="unique reaction names",
    ):
        ElementGenerationModel(
            reaction_network=network,
            element_yields=[
                ReactionElementYield(
                    reaction_name="Reaction A",
                    element_yields={
                        "C": 1.0,
                    },
                ),
                ReactionElementYield(
                    reaction_name="Reaction A",
                    element_yields={
                        "H": 2.0,
                    },
                ),
            ],
        )


def test_element_generation_rejects_nonpositive_mass():
    network = build_network()

    model = ElementGenerationModel(
        reaction_network=network,
        element_yields=[
            ReactionElementYield(
                reaction_name="Reaction A",
                element_yields={
                    "C": 1.0,
                },
            )
        ],
    )

    with pytest.raises(
        ValueError,
        match="Cell mass",
    ):
        model.generation_rates(
            temperature=400.0,
            conversions=[0.0],
            cell_mass=0.0,
        )


def test_generated_moles_rejects_decreasing_conversion():
    network = build_network()

    model = ElementGenerationModel(
        reaction_network=network,
        element_yields=[
            ReactionElementYield(
                reaction_name="Reaction A",
                element_yields={
                    "C": 1.0,
                },
            )
        ],
    )

    with pytest.raises(
        ValueError,
        match="must not be less",
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