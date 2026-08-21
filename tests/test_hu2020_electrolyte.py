import pytest

from liiontr.kinetics import TemperatureThresholdKinetics
from liiontr.library import cell_21700_generic
from liiontr.library.hu2020 import (
    hu2020_electrolyte_decomposition,
    hu2020_electrolyte_reaction,
)


def test_hu2020_electrolyte_parameters():
    parameters = hu2020_electrolyte_decomposition()

    assert parameters.name == "Electrolyte decomposition"

    assert parameters.activation_energy == pytest.approx(2.74e5)

    assert parameters.pre_exponential_factor == pytest.approx(5.14e25)

    assert parameters.enthalpy == pytest.approx(155000.0)

    assert parameters.specific_content == pytest.approx(406.9)

    assert parameters.initial_remaining_fraction == pytest.approx(1.0)

    assert parameters.initial_conversion == pytest.approx(0.0)

    assert parameters.reaction_order == pytest.approx(1.0)


def test_hu2020_electrolyte_builds_reaction():
    cell = cell_21700_generic()

    reaction = hu2020_electrolyte_reaction(
        cell=cell,
    )

    expected_mass_fraction = 406.9 / cell.material.density(298.15)

    assert reaction.name == "Electrolyte decomposition"

    assert reaction.mass_fraction == pytest.approx(expected_mass_fraction)

    assert isinstance(
        reaction.kinetics,
        TemperatureThresholdKinetics,
    )

    assert reaction.kinetics.minimum_temperature == pytest.approx(473.15)

    base = reaction.kinetics.kinetics

    assert base.activation_energy == pytest.approx(2.74e5)

    assert base.pre_exponential_factor == pytest.approx(5.14e25)


def test_hu2020_electrolyte_is_inactive_below_threshold():
    cell = cell_21700_generic()

    reaction = hu2020_electrolyte_reaction(
        cell=cell,
    )

    rate = reaction.progress_rate(
        temperature=470.0,
        conversion=0.0,
    )

    assert rate == pytest.approx(0.0)


def test_hu2020_electrolyte_is_active_above_threshold():
    cell = cell_21700_generic()

    reaction = hu2020_electrolyte_reaction(
        cell=cell,
    )

    rate = reaction.progress_rate(
        temperature=480.0,
        conversion=0.0,
    )

    assert rate > 0.0
