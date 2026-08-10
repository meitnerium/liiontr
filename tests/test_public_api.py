from liiontr import (
    CylindricalCell,
    CylindricalGeometry,
    ConstantProperty,
    Material,
    NMC811,
)


def test_public_api():
    assert isinstance(CylindricalCell.__name__, str)
    assert isinstance(CylindricalGeometry.__name__, str)
    assert isinstance(ConstantProperty.__name__, str)
    assert isinstance(Material.__name__, str)
    assert isinstance(NMC811.__name__, str)
