from liiontr.chemistry import (
    ArrheniusBackend,
    CanteraBackend,
)


def test_arrhenius_backend():
    model = ArrheniusBackend(
        activation_energy=120000,
        pre_exponential_factor=1e5,
        enthalpy=500000,
    )

    assert model.heat_generation(500) > 0


def test_cantera_backend():
    model = CanteraBackend("nmc811.yaml")

    assert model.mechanism == "nmc811.yaml"
