import pytest

from liiontr.library import cell_21700_generic
from liiontr.library.hu2020 import hu2020_sei_decomposition


def test_hu2020_sei_parameters():
    parameters = hu2020_sei_decomposition()

    assert parameters.name == "SEI decomposition"

    assert parameters.activation_energy == pytest.approx(1.3508e5)

    assert parameters.pre_exponential_factor == pytest.approx(1.667e15)

    assert parameters.enthalpy == pytest.approx(257000.0)

    assert parameters.specific_content == pytest.approx(610.4)

    assert parameters.initial_remaining_fraction == pytest.approx(0.15)

    assert parameters.initial_conversion == pytest.approx(0.85)

    assert parameters.reaction_order == pytest.approx(1.0)


def test_hu2020_sei_builds_reaction():
    cell = cell_21700_generic()

    parameters = hu2020_sei_decomposition()

    reaction = parameters.build(
        cell=cell,
    )

    expected_mass_fraction = 610.4 / cell.material.density(298.15)

    assert reaction.name == "SEI decomposition"

    assert reaction.mass_fraction == pytest.approx(expected_mass_fraction)

    assert reaction.kinetics.activation_energy == pytest.approx(1.3508e5)

    assert reaction.kinetics.pre_exponential_factor == pytest.approx(1.667e15)


def test_hu2020_sei_has_reference():
    parameters = hu2020_sei_decomposition()

    assert parameters.reference is not None
    assert "10.1021/acsomega.0c01862" in parameters.reference
