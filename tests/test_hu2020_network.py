import pytest

from liiontr.library import cell_21700_generic
from liiontr.library.hu2020 import (
    hu2020_initial_conversions,
    hu2020_reaction_network,
)
from liiontr.reactions import (
    MultiChannelReaction,
    RemainingFractionRatioVariable,
)


def test_hu2020_reaction_network_contains_four_reactions():
    cell = cell_21700_generic()

    network = hu2020_reaction_network(
        cell=cell,
    )

    assert len(network.reactions) == 4

    assert network.reactions[0].name == "SEI decomposition"
    assert network.reactions[1].name == "Anode-electrolyte"
    assert network.reactions[2].name == "Cathode decomposition"
    assert network.reactions[3].name == "Electrolyte decomposition"


def test_hu2020_cathode_is_multichannel():
    cell = cell_21700_generic()

    network = hu2020_reaction_network(
        cell=cell,
    )

    cathode = network.reactions[2]

    assert isinstance(
        cathode,
        MultiChannelReaction,
    )


def test_hu2020_network_defines_sei_thickness_ratio():
    cell = cell_21700_generic()

    network = hu2020_reaction_network(
        cell=cell,
    )

    variable = network.context_variables["sei_thickness_ratio"]

    assert isinstance(
        variable,
        RemainingFractionRatioVariable,
    )

    assert variable.reaction_name == "SEI decomposition"

    assert variable.reference_remaining_fraction == pytest.approx(0.15)


def test_hu2020_initial_conversions():
    conversions = hu2020_initial_conversions()

    assert conversions == pytest.approx(
        [
            0.85,
            0.25,
            0.04,
            0.00,
        ]
    )


def test_hu2020_initial_context_has_unit_sei_thickness_ratio():
    cell = cell_21700_generic()

    network = hu2020_reaction_network(
        cell=cell,
    )

    context = network.context(conversions=hu2020_initial_conversions())

    assert context.variable("sei_thickness_ratio") == pytest.approx(1.0)
