import pytest

from liiontr.chemistry import ReactionNetworkBackend
from liiontr.gases import (
    GasGenerationModel,
    IdealGasPressureModel,
    ReactionGasYield,
)
from liiontr.kinetics import Arrhenius
from liiontr.library import cell_21700_generic
from liiontr.problems import ThermalProblem
from liiontr.reactions import (
    Reaction,
    ReactionNetwork,
)
from liiontr.solver import ScipySolver


def test_solver_stops_at_maximum_pressure():
    cell = cell_21700_generic()

    reaction = Reaction(
        name="Gas reaction",
        kinetics=Arrhenius(
            activation_energy=1.0,
            pre_exponential_factor=2.0,
        ),
        enthalpy=10000.0,
        mass_fraction=0.10,
    )

    network = ReactionNetwork(
        reactions=[
            reaction,
        ]
    )

    backend = ReactionNetworkBackend(
        reaction_network=network,
        cell=cell,
    )

    gas_generation_model = GasGenerationModel(
        reaction_network=network,
        gas_yields=[
            ReactionGasYield(
                reaction_name="Gas reaction",
                species_yields={
                    "CO2": 2.0,
                },
            )
        ],
    )

    pressure_model = IdealGasPressureModel(
        free_volume=1.0e-6,
        initial_pressure=101325.0,
        initial_temperature=400.0,
    )

    problem = ThermalProblem(
        cell=cell,
        chemistry_backend=backend,
        gas_generation_model=gas_generation_model,
        pressure_model=pressure_model,
        initial_temperature=400.0,
        initial_conversions=[
            0.0,
        ],
        ambient_temperature=400.0,
        convection_coefficient=0.0,
        duration=10.0,
        maximum_pressure=200000.0,
    )

    results = ScipySolver().solve(problem)

    pressure = results.get("pressure")

    assert pressure is not None

    assert results.time is not None

    assert results.time[-1] < problem.duration

    assert pressure[-1] == pytest.approx(
        200000.0,
        rel=1.0e-6,
    )


def test_maximum_pressure_requires_pressure_model():
    cell = cell_21700_generic()

    problem = ThermalProblem(
        cell=cell,
        maximum_pressure=200000.0,
    )

    with pytest.raises(
        ValueError,
        match="Maximum pressure requires",
    ):
        ScipySolver().solve(problem)
