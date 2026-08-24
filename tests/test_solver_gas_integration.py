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


def test_solver_integrates_generated_gas_and_pressure():
    cell = cell_21700_generic()

    reaction = Reaction(
        name="Gas reaction",
        kinetics=Arrhenius(
            activation_energy=1.0,
            pre_exponential_factor=1.0,
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
        duration=1.0,
    )

    results = ScipySolver().solve(problem)

    temperature = results.get("temperature")
    conversion = results.get("conversion_0")
    co2 = results.get("gas_CO2")
    pressure = results.get("pressure")

    assert temperature is not None
    assert conversion is not None
    assert co2 is not None
    assert pressure is not None

    assert co2[0] == pytest.approx(0.0)

    assert co2[-1] > 0.0

    expected_moles = 2.0 * cell.mass * reaction.mass_fraction * conversion[-1]

    assert co2[-1] == pytest.approx(
        expected_moles,
        rel=1.0e-5,
    )

    assert pressure[0] == pytest.approx(101325.0)

    expected_pressure = pressure_model.pressure(
        temperature=float(temperature[-1]),
        generated_moles=float(co2[-1]),
    )

    assert pressure[-1] == pytest.approx(expected_pressure)

    assert pressure[-1] > pressure[0]
