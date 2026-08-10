from liiontr.core.variable import Variable
from liiontr.core.scalar_field import ScalarField


def test_variable():
    field = ScalarField([300, 310, 350])

    temperature = Variable(
        name="temperature",
        unit="K",
        field=field,
    )

    assert temperature.name == "temperature"

    assert temperature.unit == "K"

    assert temperature.values[1] == 310
