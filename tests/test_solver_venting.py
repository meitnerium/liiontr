import numpy as np

from liiontr.chemistry import ReactionNetworkBackend
from liiontr.gases import (
    CompressibleVentFlowModel,
    GasGenerationModel,
    GasInventory,
    GasSpecies,
    IdealGasPressureModel,
    MixtureVentFlowModel,
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


def build_problem(
    with_vent: bool,
) -> ThermalProblem:
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
            "N2": pressure_model.initial_moles,
        },
    )

    vent_model = None
    vent_open_pressure = None

    if with_vent:
        vent_model = MixtureVentFlowModel(
            flow_model=CompressibleVentFlowModel(
                vent_area=1.0e-5,
                discharge_coefficient=0.8,
                heat_capacity_ratio=1.30,
            ),
            downstream_pressure=101325.0,
        )

        vent_open_pressure = 200000.0

    return ThermalProblem(
        cell=cell,
        chemistry_backend=backend,
        gas_generation_model=gas_generation_model,
        initial_gas_inventory=initial_gas_inventory,
        pressure_model=pressure_model,
        vent_model=vent_model,
        vent_open_pressure=vent_open_pressure,
        initial_temperature=400.0,
        initial_conversions=[
            0.0,
        ],
        ambient_temperature=400.0,
        convection_coefficient=0.0,
        duration=1.0,
    )


def test_open_vent_reduces_internal_pressure():
    closed_results = ScipySolver().solve(
        build_problem(
            with_vent=False,
        )
    )

    vented_results = ScipySolver().solve(
        build_problem(
            with_vent=True,
        )
    )

    closed_pressure = closed_results.get("pressure")

    vented_pressure = vented_results.get("pressure")

    assert closed_pressure is not None
    assert vented_pressure is not None

    assert vented_pressure[-1] < closed_pressure[-1]


def test_vent_opens_when_threshold_is_reached():
    results = ScipySolver().solve(
        build_problem(
            with_vent=True,
        )
    )

    pressure = results.get("pressure")

    vent_open = results.get("vent_open")

    assert pressure is not None
    assert vent_open is not None

    assert max(pressure) >= 200000.0 * 0.999

    assert vent_open[0] == 0.0
    assert vent_open[-1] == 1.0


def test_vent_never_closes_after_opening():
    results = ScipySolver().solve(
        build_problem(
            with_vent=True,
        )
    )

    vent_open = results.get("vent_open")

    assert vent_open is not None

    differences = np.diff(vent_open)

    assert np.all(differences >= 0.0)
