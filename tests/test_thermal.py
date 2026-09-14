import pytest

from liiontr.library.cells import cell_21700_generic
from liiontr.thermal.lumped import LumpedThermalModel
from tests.helpers import create_test_cell


def test_lumped_temperature():
    cell = create_test_cell()

    model = LumpedThermalModel(cell=cell)

    rate = model.temperature_derivative(
        temperature=300,
        heat_generation=0,
    )

    assert rate < 0

def test_lumped_thermal_heat_loss():
    cell = cell_21700_generic()

    model = LumpedThermalModel(
        cell=cell,
        convection_coefficient=10.0,
        ambient_temperature=300.0,
    )

    expected = (
        10.0
        * cell.surface_area
        * (400.0 - 300.0)
    )

    assert model.heat_loss(400.0) == pytest.approx(
        expected
    )

    assert model.heat_loss(300.0) == pytest.approx(
        0.0
    )

    assert model.heat_loss(250.0) < 0.0