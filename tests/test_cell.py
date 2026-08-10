from liiontr.cells.cylindrical import CylindricalCell
from liiontr.geometry.cylindrical import CylindricalGeometry
from liiontr.materials.material import Material
from liiontr.materials.properties import ConstantProperty
from liiontr.chemistry.nmc import NMC811


def test_cell():
    geometry = CylindricalGeometry(
        radius=0.0105,
        height=0.070,
    )

    material = Material(
        name="Generic",
        density=ConstantProperty(2500.0, "kg/m^3"),
        heat_capacity=ConstantProperty(1000.0, "J/kg/K"),
        thermal_conductivity=ConstantProperty(1.0, "W/m/K"),
    )

    cell = CylindricalCell(
        name="21700",
        geometry=geometry,
        material=material,
        chemistry=NMC811(),
    )

    assert cell.volume > 0.0
    assert cell.mass > 0.0
    assert cell.thermal_capacity > 0.0
