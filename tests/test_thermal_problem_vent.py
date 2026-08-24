import pytest

from liiontr.solver import ScipySolver

from liiontr.gases import (
    CompressibleVentFlowModel,
    GasInventory,
    GasSpecies,
    IdealGasPressureModel,
    MixtureVentFlowModel,
)
from liiontr.library import cell_21700_generic
from liiontr.problems import ThermalProblem


def test_thermal_problem_accepts_vent_configuration():
    cell = cell_21700_generic()

    pressure_model = IdealGasPressureModel(
        free_volume=1.0e-6,
        initial_pressure=101325.0,
        initial_temperature=298.15,
    )

    initial_gas_inventory = GasInventory(
        species=[
            GasSpecies(
                name="N2",
                molar_mass=28.0134e-3,
            ),
        ],
        moles={
            "N2": pressure_model.initial_moles,
        },
    )

    vent_model = MixtureVentFlowModel(
        flow_model=CompressibleVentFlowModel(
            vent_area=1.0e-6,
            discharge_coefficient=0.8,
            heat_capacity_ratio=1.40,
        ),
        downstream_pressure=101325.0,
    )

    problem = ThermalProblem(
        cell=cell,
        pressure_model=pressure_model,
        initial_gas_inventory=initial_gas_inventory,
        vent_model=vent_model,
        vent_open_pressure=1.5e6,
    )

    assert problem.vent_model is vent_model
    assert problem.vent_open_pressure == 1.5e6


def test_thermal_problem_has_no_vent_by_default():
    cell = cell_21700_generic()

    problem = ThermalProblem(
        cell=cell,
    )

    assert problem.vent_model is None
    assert problem.vent_open_pressure is None


def test_vent_model_requires_opening_pressure():
    cell = cell_21700_generic()

    pressure_model = IdealGasPressureModel(
        free_volume=1.0e-6,
    )

    inventory = GasInventory(
        species=[
            GasSpecies(
                name="N2",
                molar_mass=28.0134e-3,
            ),
        ],
        moles={
            "N2": pressure_model.initial_moles,
        },
    )

    vent_model = MixtureVentFlowModel(
        flow_model=CompressibleVentFlowModel(
            vent_area=1.0e-6,
        )
    )

    problem = ThermalProblem(
        cell=cell,
        pressure_model=pressure_model,
        initial_gas_inventory=inventory,
        vent_model=vent_model,
    )

    with pytest.raises(
        ValueError,
        match="opening pressure",
    ):
        ScipySolver().solve(problem)


def test_vent_requires_initial_gas_inventory():
    cell = cell_21700_generic()

    pressure_model = IdealGasPressureModel(
        free_volume=1.0e-6,
    )

    vent_model = MixtureVentFlowModel(
        flow_model=CompressibleVentFlowModel(
            vent_area=1.0e-6,
        )
    )

    problem = ThermalProblem(
        cell=cell,
        pressure_model=pressure_model,
        vent_model=vent_model,
        vent_open_pressure=1.5e6,
    )

    with pytest.raises(
        ValueError,
        match="initial gas inventory",
    ):
        ScipySolver().solve(problem)
