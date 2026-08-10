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
