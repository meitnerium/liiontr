import pytest

from liiontr.kinetics import (
    Arrhenius,
    TemperatureThresholdKinetics,
)


def test_temperature_threshold_is_inactive_below_threshold():
    kinetics = TemperatureThresholdKinetics(
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        minimum_temperature=400.0,
    )

    rate = kinetics.rate(
        temperature=399.0,
    )

    assert rate == pytest.approx(0.0)


def test_temperature_threshold_is_active_above_threshold():
    base_kinetics = Arrhenius(
        activation_energy=80000.0,
        pre_exponential_factor=1.0e5,
    )

    kinetics = TemperatureThresholdKinetics(
        kinetics=base_kinetics,
        minimum_temperature=400.0,
    )

    rate = kinetics.rate(
        temperature=401.0,
    )

    assert rate == pytest.approx(base_kinetics.rate(401.0))


def test_temperature_threshold_is_active_at_threshold():
    base_kinetics = Arrhenius(
        activation_energy=80000.0,
        pre_exponential_factor=1.0e5,
    )

    kinetics = TemperatureThresholdKinetics(
        kinetics=base_kinetics,
        minimum_temperature=400.0,
    )

    rate = kinetics.rate(
        temperature=400.0,
    )

    assert rate == pytest.approx(base_kinetics.rate(400.0))


def test_temperature_threshold_must_be_positive():
    with pytest.raises(
        ValueError,
        match="Minimum temperature must be greater than zero",
    ):
        TemperatureThresholdKinetics(
            kinetics=Arrhenius(
                activation_energy=80000.0,
                pre_exponential_factor=1.0e5,
            ),
            minimum_temperature=0.0,
        )
