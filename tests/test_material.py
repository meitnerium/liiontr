from liiontr.materials.material import Material
from liiontr.materials.properties import ConstantProperty


def test_material():
    aluminum = Material(
        name="Aluminum",
        density=ConstantProperty(2700.0, "kg/m^3"),
        heat_capacity=ConstantProperty(900.0, "J/kg/K"),
        thermal_conductivity=ConstantProperty(237.0, "W/m/K"),
    )

    assert aluminum.density(300.0) == 2700.0
    assert aluminum.heat_capacity(500.0) == 900.0
