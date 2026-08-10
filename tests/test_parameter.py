from liiontr.core.parameter import Parameter
from liiontr.core.parameter_set import ParameterSet


def test_add_parameter():
    parameters = ParameterSet()

    parameters.add(
        Parameter(
            name="sei.A",
            value=2.0e13,
            unit="1/s",
            description="SEI pre-exponential factor",
        )
    )

    assert len(parameters) == 1

    assert parameters.value("sei.A") == 2.0e13
