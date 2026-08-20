import pytest

from liiontr.library import cell_21700_generic
from liiontr.library.hu2020 import (
    hu2020_initial_conversions,
    hu2020_reaction_network,
)
from liiontr.reactions import RemainingFractionRatioVariable


def test_hu2020_reaction_network_contains_sei_and_anode():
    cell = cell_21700_generic()

    network = hu2020_reaction_network(
        cell=cell,
    )

    assert len(network.reactions) == 2

    assert network.reactions[0].name == "SEI decomposition"
    assert network.reactions[1].name == "Anode-electrolyte"


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
        ]
    )


def test_hu2020_initial_context_has_unit_sei_thickness_ratio():
    cell = cell_21700_generic()

    network = hu2020_reaction_network(
        cell=cell,
    )

    context = network.context(conversions=hu2020_initial_conversions())

    assert context.variable("sei_thickness_ratio") == pytest.approx(1.0)
