"""Reference battery-cell definitions provided by LiionTR."""

from liiontr import (
    LCO,
    NMC811,
    ConstantProperty,
    CylindricalCell,
    CylindricalGeometry,
    Material,
)


def cell_21700_generic() -> CylindricalCell:
    """Return a generic 21700 NMC811 cell.

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


def cell_18650_lco_hatchard2001() -> CylindricalCell:
    """Return the 18650 LiCoO2/graphite cell used by Hatchard et al. (2001).

    The effective homogeneous properties follow the standard cylindrical
    model parameters reported by Hatchard et al. for E-One/Moli Energy
    18650 cells: 9 mm radius, 65 mm length, 2.58 g/cm3 average jelly-roll
    density, 0.83 J/(g K) average specific heat, and 0.034 W/(cm K)
    average radial thermal conductivity.

    Notes
    -----
    LiionTR represents the cell as a homogeneous cylinder. The original
    Hatchard model also treated the steel can separately in its outer
    radial ring. The present cell definition is therefore intended for
    the lumped/uniform-temperature version of the benchmark.
    """
    geometry = CylindricalGeometry(
        radius=0.009,
        height=0.065,
    )

    material = Material(
        name="Hatchard 2001 effective 18650 cell",
        density=ConstantProperty(
            2580.0,
            "kg/m3",
        ),
        heat_capacity=ConstantProperty(
            830.0,
            "J/kg/K",
        ),
        thermal_conductivity=ConstantProperty(
            3.4,
            "W/m/K",
        ),
    )

    return CylindricalCell(
        name="18650-LCO-Hatchard2001",
        geometry=geometry,
        material=material,
        chemistry=LCO(),
    )
