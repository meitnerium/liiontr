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


def test_thermal_problem_accepts_gas_models():
    cell = cell_21700_generic()

    reaction = Reaction(
        name="Gas reaction",
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=200000.0,
        mass_fraction=0.10,
    )

    network = ReactionNetwork(
        reactions=[
            reaction,
        ]
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
        initial_temperature=298.15,
    )

    problem = ThermalProblem(
        cell=cell,
        gas_generation_model=gas_generation_model,
        pressure_model=pressure_model,
    )

    assert problem.gas_generation_model is gas_generation_model
    assert problem.pressure_model is pressure_model


def test_thermal_problem_has_no_gas_model_by_default():
    cell = cell_21700_generic()

    problem = ThermalProblem(
        cell=cell,
    )

    assert problem.gas_generation_model is None
    assert problem.pressure_model is None
