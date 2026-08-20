import math

import pytest

from liiontr.kinetics import (
    ExponentialInhibitionProgress,
    PowerLawProgress,
    ThresholdProgress,
)
from liiontr.library import cell_21700_generic
from liiontr.library.hu2020 import (
    hu2020_anode_electrolyte,
    hu2020_anode_progress_model,
)
from liiontr.reactions import ReactionContext


def test_hu2020_anode_parameters():
    parameters = hu2020_anode_electrolyte()

    assert parameters.name == "Anode-electrolyte"

    assert parameters.activation_energy == pytest.approx(1.3508e5)

    assert parameters.pre_exponential_factor == pytest.approx(2.5e13)

    assert parameters.enthalpy == pytest.approx(1.714e6)

    assert parameters.specific_content == pytest.approx(610.4)

    assert parameters.initial_remaining_fraction == pytest.approx(0.75)

    assert parameters.initial_conversion == pytest.approx(0.25)

    assert parameters.reaction_order == pytest.approx(1.0)


def test_hu2020_anode_progress_model_structure():
    model = hu2020_anode_progress_model()

    assert isinstance(
        model,
        ThresholdProgress,
    )

    assert model.reaction_name == "SEI decomposition"

    assert model.remaining_below == pytest.approx(0.10)

    assert isinstance(
        model.progress_model,
        ExponentialInhibitionProgress,
    )

    assert model.progress_model.variable_name == "sei_thickness_ratio"

    assert isinstance(
        model.progress_model.progress_model,
        PowerLawProgress,
    )


def test_hu2020_anode_is_inactive_before_sei_threshold():
    model = hu2020_anode_progress_model()

    context = ReactionContext(
        conversions={
            "SEI decomposition": 0.80,
        },
        variables={
            "sei_thickness_ratio": 1.0,
        },
    )

    factor = model.factor(
        conversion=0.25,
        context=context,
    )

    assert factor == pytest.approx(0.0)


def test_hu2020_anode_is_active_after_sei_threshold():
    model = hu2020_anode_progress_model()

    context = ReactionContext(
        conversions={
            "SEI decomposition": 0.95,
        },
        variables={
            "sei_thickness_ratio": 1.0,
        },
    )

    factor = model.factor(
        conversion=0.25,
        context=context,
    )

    expected = 0.75 * math.exp(-1.0)

    assert factor == pytest.approx(expected)


def test_hu2020_anode_builds_reaction():
    cell = cell_21700_generic()

    parameters = hu2020_anode_electrolyte()

    progress_model = hu2020_anode_progress_model()

    reaction = parameters.build(
        cell=cell,
        progress_model=progress_model,
    )

    expected_mass_fraction = 610.4 / cell.material.density(298.15)

    assert reaction.name == "Anode-electrolyte"

    assert reaction.mass_fraction == pytest.approx(expected_mass_fraction)

    assert reaction.progress_model is progress_model
