import pytest

from liiontr.gases import IdealGasPressureModel


def test_initial_pressure_is_recovered():
    model = IdealGasPressureModel(
        free_volume=1.0e-6,
        initial_pressure=101325.0,
        initial_temperature=298.15,
    )

    pressure = model.pressure(
        temperature=298.15,
        generated_moles=0.0,
    )

    assert pressure == pytest.approx(101325.0)


def test_pressure_increases_with_temperature():
    model = IdealGasPressureModel(
        free_volume=1.0e-6,
        initial_pressure=101325.0,
        initial_temperature=300.0,
    )

    pressure = model.pressure(
        temperature=600.0,
        generated_moles=0.0,
    )

    assert pressure == pytest.approx(2.0 * 101325.0)


def test_generated_gas_increases_pressure():
    model = IdealGasPressureModel(
        free_volume=1.0e-6,
        initial_pressure=101325.0,
        initial_temperature=298.15,
    )

    initial_pressure = model.pressure(
        temperature=298.15,
        generated_moles=0.0,
    )

    final_pressure = model.pressure(
        temperature=298.15,
        generated_moles=1.0e-4,
    )

    assert final_pressure > initial_pressure


def test_initial_moles_follow_ideal_gas_law():
    model = IdealGasPressureModel(
        free_volume=2.0e-6,
        initial_pressure=101325.0,
        initial_temperature=298.15,
    )

    expected = 101325.0 * 2.0e-6 / (8.314462618 * 298.15)

    assert model.initial_moles == pytest.approx(expected)


@pytest.mark.parametrize(
    "free_volume",
    [
        0.0,
        -1.0e-6,
    ],
)
def test_free_volume_must_be_positive(
    free_volume: float,
):
    with pytest.raises(
        ValueError,
        match="Free volume",
    ):
        IdealGasPressureModel(
            free_volume=free_volume,
        )


def test_pressure_from_total_moles():
    model = IdealGasPressureModel(
        free_volume=1.0e-6,
        initial_pressure=101325.0,
        initial_temperature=298.15,
    )

    pressure = model.pressure_from_total_moles(
        temperature=298.15,
        total_moles=model.initial_moles,
    )

    assert pressure == pytest.approx(101325.0)


def test_pressure_and_total_moles_formulations_are_equivalent():
    model = IdealGasPressureModel(
        free_volume=1.0e-6,
        initial_pressure=101325.0,
        initial_temperature=298.15,
    )

    generated_moles = 2.0e-4

    pressure_from_generated = model.pressure(
        temperature=500.0,
        generated_moles=generated_moles,
    )

    pressure_from_total = model.pressure_from_total_moles(
        temperature=500.0,
        total_moles=(model.initial_moles + generated_moles),
    )

    assert pressure_from_total == pytest.approx(pressure_from_generated)


def test_zero_total_moles_gives_zero_pressure():
    model = IdealGasPressureModel(
        free_volume=1.0e-6,
    )

    pressure = model.pressure_from_total_moles(
        temperature=500.0,
        total_moles=0.0,
    )

    assert pressure == pytest.approx(0.0)


def test_negative_total_moles_are_rejected():
    model = IdealGasPressureModel(
        free_volume=1.0e-6,
    )

    with pytest.raises(
        ValueError,
        match="Total moles",
    ):
        model.pressure_from_total_moles(
            temperature=500.0,
            total_moles=-1.0,
        )
