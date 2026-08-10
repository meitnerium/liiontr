from liiontr.core.scalar_field import ScalarField


def test_scalar_field():
    field = ScalarField([1, 2, 3])

    assert field.size == 3

    assert field.values[1] == 2
