from liiontr import (
    CylindricalCell,
    CylindricalGeometry,
    ConstantProperty,
    Material,
    NMC811,
)


def create_test_cell():
    geometry = CylindricalGeometry(
        radius=0.0105,
        height=0.070,
    )

    material = Material(
        name="Generic",
        density=ConstantProperty(
            2500,
            "kg/m3",
        ),
        heat_capacity=ConstantProperty(
            1000,
            "J/kg/K",
        ),
        thermal_conductivity=ConstantProperty(
            1,
            "W/m/K",
        ),
    )

    return CylindricalCell(
        name="21700",
        geometry=geometry,
        material=material,
        chemistry=NMC811(),
    )
