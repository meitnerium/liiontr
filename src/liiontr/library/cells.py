from liiontr import (
    CylindricalCell,
    CylindricalGeometry,
    ConstantProperty,
    Material,
    NMC811,
)


def cell_21700_generic() -> CylindricalCell:
    """
    Generic 21700 NMC811 cell.

    Parameters are approximate and intended
    for software testing.
    """

    geometry = CylindricalGeometry(
        radius=0.0105,
        height=0.070,
    )

    material = Material(
        name="Generic Li-ion cell",
        density=ConstantProperty(
            2500.0,
            "kg/m3",
        ),
        heat_capacity=ConstantProperty(
            1000.0,
            "J/kg/K",
        ),
        thermal_conductivity=ConstantProperty(
            1.0,
            "W/m/K",
        ),
    )

    return CylindricalCell(
        name="21700-NMC811",
        geometry=geometry,
        material=material,
        chemistry=NMC811(),
    )
