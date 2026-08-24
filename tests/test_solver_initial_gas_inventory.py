import pytest

from liiontr.chemistry import ReactionNetworkBackend
from liiontr.gases import (
    GasGenerationModel,
    GasInventory,
    GasSpecies,
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


def test_solver_tracks_initial_and_generated_gases():
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

    initial_gas_inventory = GasInventory(
        species=[
            GasSpecies(
                name="N2",
                molar_mass=28.0134e-3,
            ),
            GasSpecies(
                name="CO2",
                molar_mass=44.0095e-3,
            ),
        ],
        moles={
            "N2": 2.0e-5,
        },
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
        initial_gas_inventory=initial_gas_inventory,
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

    nitrogen = results.get("gas_N2")

    co2 = results.get("gas_CO2")

    pressure = results.get("pressure")

    assert nitrogen is not None
    assert co2 is not None
    assert pressure is not None

    assert nitrogen[0] == pytest.approx(2.0e-5)

    assert nitrogen[-1] == pytest.approx(nitrogen[0])

    assert co2[0] == pytest.approx(0.0)

    assert co2[-1] > 0.0

    expected_initial_pressure = pressure_model.pressure_from_total_moles(
        temperature=400.0,
        total_moles=2.0e-5,
    )

    assert pressure[0] == pytest.approx(expected_initial_pressure)
