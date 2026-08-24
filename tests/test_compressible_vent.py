import pytest

from liiontr.gases import (
    GAS_CONSTANT,
    CompressibleVentFlowModel,
)


def test_vent_has_no_flow_without_pressure_difference():
    model = CompressibleVentFlowModel(
        vent_area=1.0e-6,
        discharge_coefficient=0.8,
        heat_capacity_ratio=1.30,
    )

    mass_flow = model.mass_flow_rate(
        upstream_pressure=101325.0,
        downstream_pressure=101325.0,
        temperature=500.0,
        molar_mass=44.0e-3,
    )

    assert mass_flow == pytest.approx(0.0)


def test_vent_has_no_reverse_flow():
    model = CompressibleVentFlowModel(
        vent_area=1.0e-6,
        discharge_coefficient=0.8,
        heat_capacity_ratio=1.30,
    )

    mass_flow = model.mass_flow_rate(
        upstream_pressure=90000.0,
        downstream_pressure=101325.0,
        temperature=500.0,
        molar_mass=44.0e-3,
    )

    assert mass_flow == pytest.approx(0.0)


def test_critical_pressure_ratio():
    model = CompressibleVentFlowModel(
        vent_area=1.0e-6,
        discharge_coefficient=0.8,
        heat_capacity_ratio=1.40,
    )

    expected = (2.0 / (1.40 + 1.0)) ** (1.40 / (1.40 - 1.0))

    assert model.critical_pressure_ratio == pytest.approx(expected)


def test_choked_mass_flow_rate():
    model = CompressibleVentFlowModel(
        vent_area=1.0e-6,
        discharge_coefficient=0.8,
        heat_capacity_ratio=1.30,
    )

    upstream_pressure = 2.0e6
    downstream_pressure = 101325.0
    temperature = 600.0
    molar_mass = 44.0e-3

    specific_gas_constant = GAS_CONSTANT / molar_mass

    gamma = 1.30

    expected = (
        0.8
        * 1.0e-6
        * upstream_pressure
        * (gamma / (specific_gas_constant * temperature)) ** 0.5
        * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    )

    mass_flow = model.mass_flow_rate(
        upstream_pressure=upstream_pressure,
        downstream_pressure=downstream_pressure,
        temperature=temperature,
        molar_mass=molar_mass,
    )

    assert mass_flow == pytest.approx(expected)


def test_unchoked_flow_is_positive():
    model = CompressibleVentFlowModel(
        vent_area=1.0e-6,
        discharge_coefficient=0.8,
        heat_capacity_ratio=1.30,
    )

    mass_flow = model.mass_flow_rate(
        upstream_pressure=120000.0,
        downstream_pressure=101325.0,
        temperature=500.0,
        molar_mass=44.0e-3,
    )

    assert mass_flow > 0.0


def test_molar_flow_rate_matches_mass_flow_rate():
    model = CompressibleVentFlowModel(
        vent_area=1.0e-6,
        discharge_coefficient=0.8,
        heat_capacity_ratio=1.30,
    )

    molar_mass = 44.0e-3

    mass_flow = model.mass_flow_rate(
        upstream_pressure=2.0e6,
        downstream_pressure=101325.0,
        temperature=600.0,
        molar_mass=molar_mass,
    )

    molar_flow = model.molar_flow_rate(
        upstream_pressure=2.0e6,
        downstream_pressure=101325.0,
        temperature=600.0,
        molar_mass=molar_mass,
    )

    assert molar_flow == pytest.approx(mass_flow / molar_mass)


@pytest.mark.parametrize(
    "vent_area",
    [
        0.0,
        -1.0e-6,
    ],
)
def test_vent_area_must_be_positive(
    vent_area: float,
):
    with pytest.raises(
        ValueError,
        match="Vent area",
    ):
        CompressibleVentFlowModel(
            vent_area=vent_area,
        )


def test_heat_capacity_ratio_must_be_greater_than_one():
    with pytest.raises(
        ValueError,
        match="Heat capacity ratio",
    ):
        CompressibleVentFlowModel(
            vent_area=1.0e-6,
            heat_capacity_ratio=1.0,
        )
