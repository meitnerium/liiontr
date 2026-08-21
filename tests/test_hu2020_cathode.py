import pytest

from liiontr.kinetics import (
    AutocatalyticProgress,
    TemperatureThresholdKinetics,
)
from liiontr.library import cell_21700_generic
from liiontr.library.hu2020 import (
    hu2020_cathode_decomposition,
    hu2020_cathode_initial_conversion,
)
from liiontr.reactions import MultiChannelReaction


def test_hu2020_cathode_initial_conversion():
    assert hu2020_cathode_initial_conversion() == pytest.approx(0.04)


def test_hu2020_cathode_is_multichannel():
    cell = cell_21700_generic()

    reaction = hu2020_cathode_decomposition(
        cell=cell,
    )

    assert isinstance(
        reaction,
        MultiChannelReaction,
    )

    assert reaction.name == "Cathode decomposition"

    assert len(reaction.channels) == 2


def test_hu2020_cathode_mass_fraction():
    cell = cell_21700_generic()

    reaction = hu2020_cathode_decomposition(
        cell=cell,
    )

    expected_mass_fraction = 1221.0 / cell.material.density(298.15)

    assert reaction.mass_fraction == pytest.approx(expected_mass_fraction)


def test_hu2020_cathode_progress_model():
    cell = cell_21700_generic()

    reaction = hu2020_cathode_decomposition(
        cell=cell,
    )

    assert isinstance(
        reaction.progress_model,
        AutocatalyticProgress,
    )

    assert reaction.progress_model.autocatalytic_order == pytest.approx(1.0)

    assert reaction.progress_model.remaining_order == pytest.approx(1.0)


def test_hu2020_cathode_channel_one():
    cell = cell_21700_generic()

    reaction = hu2020_cathode_decomposition(
        cell=cell,
    )

    channel = reaction.channels[0]

    assert channel.enthalpy == pytest.approx(77000.0)

    assert isinstance(
        channel.kinetics,
        TemperatureThresholdKinetics,
    )

    assert channel.kinetics.minimum_temperature == pytest.approx(393.15)

    base = channel.kinetics.kinetics

    assert base.activation_energy == pytest.approx(1.1495e5)

    assert base.pre_exponential_factor == pytest.approx(1.75e9)


def test_hu2020_cathode_channel_two():
    cell = cell_21700_generic()

    reaction = hu2020_cathode_decomposition(
        cell=cell,
    )

    channel = reaction.channels[1]

    assert channel.enthalpy == pytest.approx(84000.0)

    assert isinstance(
        channel.kinetics,
        TemperatureThresholdKinetics,
    )

    assert channel.kinetics.minimum_temperature == pytest.approx(393.15)

    base = channel.kinetics.kinetics

    assert base.activation_energy == pytest.approx(1.5888e5)

    assert base.pre_exponential_factor == pytest.approx(1.077e12)


def test_hu2020_cathode_is_inactive_below_temperature_threshold():
    cell = cell_21700_generic()

    reaction = hu2020_cathode_decomposition(
        cell=cell,
    )

    rate = reaction.progress_rate(
        temperature=390.0,
        conversion=0.04,
    )

    assert rate == pytest.approx(0.0)


def test_hu2020_cathode_is_active_above_temperature_threshold():
    cell = cell_21700_generic()

    reaction = hu2020_cathode_decomposition(
        cell=cell,
    )

    rate = reaction.progress_rate(
        temperature=400.0,
        conversion=0.04,
    )

    assert rate > 0.0
