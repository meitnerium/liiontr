from math import pi

from liiontr.geometry.cylindrical import CylindricalGeometry


def test_cylindrical_geometry():
    g = CylindricalGeometry(
        radius=0.01,
        height=0.07,
    )

    assert abs(g.volume - pi * 0.01**2 * 0.07) < 1e-12

    assert g.surface_area > 0
