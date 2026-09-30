from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from liiontr.chemistry.energy_scaled_generation import (
    EnergyScaledElementGenerationModel,
)

pytest.importorskip(
    "cantera"
)

from liiontr.chemistry.cantera import (
    CanteraEquilibriumBackend,
)
from liiontr.chemistry.element_generation import (
    ElementGenerationModel,
)
from liiontr.chemistry.element_release import (
    FirstOrderElementReleaseModel,
)
from liiontr.chemistry.elements import (
    ElementInventory,
    ReactionElementYield,
)
from liiontr.chemistry.reaction_backend import (
    ReactionNetworkBackend,
)
from liiontr.gases.element_vent import (
    ElementVentFlowModel,
)
from liiontr.kinetics import Arrhenius
from liiontr.library import cell_21700_generic
from liiontr.problems import ThermalProblem
from liiontr.reactions import (
    Reaction,
    ReactionNetwork,
)
from liiontr.solver import ScipySolver


@pytest.fixture(scope="module")
def closed_solution():
    problem = build_problem(
        with_vent=False,
    )

    results = ScipySolver().solve(
        problem
    )

    return problem, results

@pytest.fixture(scope="module")
def vented_solution():
    problem = build_problem(
        with_vent=True,
    )

    results = ScipySolver().solve(
        problem
    )

    return problem, results

def build_problem(
    with_vent: bool,
    convection_coefficient: float = 0.0,
    element_release_model=None,
) -> ThermalProblem:
    cell = cell_21700_generic()

    reaction = Reaction(
        name="Element reaction",
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

    reaction_backend = ReactionNetworkBackend(
        reaction_network=network,
        cell=cell,
    )

    generation_model = ElementGenerationModel(
        reaction_network=network,
        element_yields=[
            ReactionElementYield(
                reaction_name="Element reaction",
                element_yields={
                    "C": 1.0,
                    "H": 4.0,
                    "O": 4.0,
                },
            )
        ],
    )

    equilibrium_backend = (
        CanteraEquilibriumBackend()
    )

    gas_volume = 1.0e-3

    temperature = 1000.0

    pressure = 101325.0

    gas_constant = 8.31446261815324

    initial_n2_moles = (
        pressure
        * gas_volume
        / (
            gas_constant
            * temperature
        )
    )

    initial_inventory = ElementInventory(
        moles={
            "N": (
                2.0
                * initial_n2_moles
            ),
        }
    )

    element_vent_model = None
    vent_open_pressure = None

    if with_vent:
        element_vent_model = (
            ElementVentFlowModel(
                equilibrium_backend=(
                    equilibrium_backend
                ),
                vent_area=1.0e-5,
                discharge_coefficient=0.8,
                downstream_pressure=101325.0,
            )
        )

        vent_open_pressure = 180000.0

    return ThermalProblem(
        cell=cell,
        chemistry_backend=(reaction_backend),
        initial_temperature=temperature,
        initial_conversions=[
            0.0,
        ],
        ambient_temperature=temperature,
        convection_coefficient=convection_coefficient,
        duration=0.5,
        element_generation_model=(generation_model),
        gas_equilibrium_backend=(equilibrium_backend),
        initial_element_inventory=(initial_inventory),
        gas_volume=gas_volume,
        element_vent_model=(element_vent_model),
        vent_open_pressure=(vent_open_pressure),
        element_release_model=(element_release_model),
    )

def test_elemental_solver_generates_elements():
    results = ScipySolver().solve(
        build_problem(
            with_vent=False,
        )
    )

    assert (
        results.get_variable(
            "element_C"
        )[-1]
        > 0.0
    )

    assert (
        results.get_variable(
            "element_H"
        )[-1]
        > 0.0
    )

    assert (
        results.get_variable(
            "element_O"
        )[-1]
        > 0.0
    )

    assert np.all(
        results.get_variable(
            "pressure"
        )
        > 0.0
    )

def test_elemental_vent_opens(
    vented_solution,
):
    _, results = vented_solution

    vent_open = results.get_variable(
        "vent_open"
    )

    assert vent_open[0] == pytest.approx(
        0.0
    )

    assert vent_open[-1] == pytest.approx(
        1.0
    )

    assert np.all(
        np.diff(
            vent_open
        )
        >= 0.0
    )

def test_elemental_vent_reduces_pressure():
    closed = ScipySolver().solve(
        build_problem(
            with_vent=False,
        )
    )

    vented = ScipySolver().solve(
        build_problem(
            with_vent=True,
        )
    )

    assert (
        vented.get_variable(
            "pressure"
        )[-1]
        <
        closed.get_variable(
            "pressure"
        )[-1]
    )

def test_elemental_state_remains_nonnegative(
    vented_solution,
):
    _, results = vented_solution

    for element_name in (
        "C",
        "H",
        "O",
        "N",
    ):
        assert np.all(
            results.get_variable(
                f"element_{element_name}"
            )
            >= 0.0
        )

def test_elemental_vent_conserves_elements(
    vented_solution,
):
    problem, results = vented_solution

    # Keep the rest of the existing test unchanged.

    conversion = results.get_variable(
        "conversion_0"
    )[-1]

    reaction = (
        problem
        .element_generation_model
        .reaction_network
        .reactions[0]
    )

    reacted_mass = (
        problem.cell.mass
        * reaction.mass_fraction
        * conversion
    )

    generated_yields = {
        "C": 1.0,
        "H": 4.0,
        "O": 4.0,
        "N": 0.0,
    }

    assert (
        problem.initial_element_inventory
        is not None
    )

    for element_name in (
        "C",
        "H",
        "O",
        "N",
    ):
        initial_amount = (
            problem.initial_element_inventory.moles_of(
                element_name
            )
        )

        generated_amount = (
            generated_yields[
                element_name
            ]
            * reacted_mass
        )

        remaining_amount = (
            results.get_variable(
                f"element_{element_name}"
            )[-1]
        )

        vented_amount = (
            results.get_variable(
                f"vented_element_{element_name}"
            )[-1]
        )

        expected_total = (
            initial_amount
            + generated_amount
        )

        actual_total = (
            remaining_amount
            + vented_amount
        )

        assert actual_total == pytest.approx(
            expected_total,
            rel=1.0e-5,
            abs=1.0e-10,
        )

def test_elemental_vent_records_vented_inventory():
    results = ScipySolver().solve(
        build_problem(
            with_vent=True,
        )
    )

    assert (
        results.get_variable(
            "vented_element_N"
        )[-1]
        > 0.0
    )

    assert (
        results.get_variable(
            "vented_element_C"
        )[-1]
        > 0.0
    )

def test_elemental_solver_reports_thermochemical_properties():
    results = ScipySolver().solve(
        build_problem(
            with_vent=False,
        )
    )

    pressure = results.get_variable(
        "pressure"
    )

    density = results.get_variable(
        "gas_density"
    )

    mean_molar_mass = results.get_variable(
        "gas_mean_molar_mass"
    )

    cp_mass = results.get_variable(
        "gas_cp_mass"
    )

    cv_mass = results.get_variable(
        "gas_cv_mass"
    )

    gamma = results.get_variable(
        "gas_heat_capacity_ratio"
    )

    assert np.all(
        pressure > 0.0
    )

    assert np.all(
        density > 0.0
    )

    assert np.all(
        mean_molar_mass > 0.0
    )

    assert np.all(
        cp_mass > 0.0
    )

    assert np.all(
        cv_mass > 0.0
    )

    assert np.all(
        cp_mass > cv_mass
    )

    assert np.all(
        gamma > 1.0
    )

    assert np.allclose(
        gamma,
        cp_mass / cv_mass,
        rtol=1.0e-10,
        atol=1.0e-12,
    )

def test_elemental_solver_reports_equilibrium_species():
    results = ScipySolver().solve(
        build_problem(
            with_vent=False,
        )
    )

    co2 = results.get_variable(
        "equilibrium_species_CO2"
    )

    h2o = results.get_variable(
        "equilibrium_species_H2O"
    )

    n2 = results.get_variable(
        "equilibrium_species_N2"
    )

    assert n2[0] > 0.0

    assert co2[-1] > 0.0

    assert h2o[-1] > 0.0

    assert n2[-1] > 0.0

def test_elemental_internal_energy_is_finite():
    results = ScipySolver().solve(
        build_problem(
            with_vent=False,
        )
    )

    internal_energy = results.get_variable(
        "gas_internal_energy"
    )

    assert np.all(
        np.isfinite(
            internal_energy
        )
    )

def test_elemental_solver_integrates_thermal_energy():
    problem = build_problem(
        with_vent=False,
    )

    results = ScipySolver().solve(
        problem
    )

    energy = results.get_variable(
        "thermal_energy"
    )

    assert np.all(
        np.isfinite(
            energy
        )
    )

    assert energy[-1] > energy[0]

def test_closed_elemental_energy_matches_reaction_heat(
    closed_solution,
):
    problem, results = (
        closed_solution
    )

    results = ScipySolver().solve(
        problem
    )

    conversion = results.get_variable(
        "conversion_0"
    )

    energy = results.get_variable(
        "thermal_energy"
    )

    reaction = (
        problem
        .element_generation_model
        .reaction_network
        .reactions[0]
    )

    expected_energy_change = (
        problem.cell.mass
        * reaction.mass_fraction
        * reaction.enthalpy
        * (
            conversion[-1]
            - conversion[0]
        )
    )

    actual_energy_change = (
        energy[-1]
        - energy[0]
    )

    assert actual_energy_change == pytest.approx(
        expected_energy_change,
        rel=1.0e-5,
        abs=1.0e-8,
    )

def test_elemental_vent_records_thermal_energy(
    vented_solution,
):
    _, results = vented_solution

    vented_energy = results.get_variable(
        "vented_thermal_energy"
    )

    assert vented_energy[0] == pytest.approx(
        0.0
    )

    assert vented_energy[-1] > 0.0

    assert np.all(
        np.diff(
            vented_energy
        )
        >= -1.0e-10
    )

def test_elemental_vent_energy_balance(
    vented_solution,
):
    problem, results = (
        vented_solution
    )

    results = ScipySolver().solve(
        problem
    )

    conversion = results.get_variable(
        "conversion_0"
    )

    thermal_energy = results.get_variable(
        "thermal_energy"
    )

    vented_energy = results.get_variable(
        "vented_thermal_energy"
    )

    reaction = (
        problem
        .element_generation_model
        .reaction_network
        .reactions[0]
    )

    reaction_heat = (
        problem.cell.mass
        * reaction.mass_fraction
        * reaction.enthalpy
        * (
            conversion[-1]
            - conversion[0]
        )
    )

    stored_energy_change = (
        thermal_energy[-1]
        - thermal_energy[0]
    )

    accounted_energy = (
        stored_energy_change
        + vented_energy[-1]
    )

    assert accounted_energy == pytest.approx(
        reaction_heat,
        rel=2.0e-5,
        abs=1.0e-8,
    )

def test_closed_elemental_case_has_no_vented_energy():
    results = ScipySolver().solve(
        build_problem(
            with_vent=False,
        )
    )

    vented_energy = results.get_variable(
        "vented_thermal_energy"
    )

    assert np.allclose(
        vented_energy,
        0.0,
        atol=1.0e-12,
    )

def test_zero_convection_has_zero_thermal_loss_energy():
    results = ScipySolver().solve(
        build_problem(
            with_vent=False,
        )
    )

    loss_energy = results.get_variable(
        "thermal_loss_energy"
    )

    assert np.allclose(
        loss_energy,
        0.0,
        atol=1.0e-12,
    )

def test_elemental_energy_balance_with_convection():
    problem = build_problem(
        with_vent=False,
        convection_coefficient=25.0,
    )

    results = ScipySolver().solve(
        problem
    )

    conversion = results.get_variable(
        "conversion_0"
    )

    thermal_energy = results.get_variable(
        "thermal_energy"
    )

    loss_energy = results.get_variable(
        "thermal_loss_energy"
    )

    reaction = (
        problem
        .element_generation_model
        .reaction_network
        .reactions[0]
    )

    reaction_heat = (
        problem.cell.mass
        * reaction.mass_fraction
        * reaction.enthalpy
        * (
            conversion[-1]
            - conversion[0]
        )
    )

    stored_energy_change = (
        thermal_energy[-1]
        - thermal_energy[0]
    )

    accounted_energy = (
        stored_energy_change
        + loss_energy[-1]
    )

    assert accounted_energy == pytest.approx(
        reaction_heat,
        rel=2.0e-5,
        abs=1.0e-8,
    )

    assert not np.isclose(
        loss_energy[-1],
        0.0,
        atol=1.0e-10,
    )

def test_elemental_full_energy_balance():
    problem = build_problem(
        with_vent=True,
        convection_coefficient=25.0,
    )

    results = ScipySolver().solve(
        problem
    )

    conversion = results.get_variable(
        "conversion_0"
    )

    thermal_energy = results.get_variable(
        "thermal_energy"
    )

    vented_energy = results.get_variable(
        "vented_thermal_energy"
    )

    loss_energy = results.get_variable(
        "thermal_loss_energy"
    )

    reaction = (
        problem
        .element_generation_model
        .reaction_network
        .reactions[0]
    )

    reaction_heat = (
        problem.cell.mass
        * reaction.mass_fraction
        * reaction.enthalpy
        * (
            conversion[-1]
            - conversion[0]
        )
    )

    accounted_energy = (
        thermal_energy[-1]
        - thermal_energy[0]
        + vented_energy[-1]
        + loss_energy[-1]
    )

    assert accounted_energy == pytest.approx(
        reaction_heat,
        rel=3.0e-5,
        abs=1.0e-8,
    )

def test_elemental_solver_supports_delayed_release():
    """Store generated elements before releasing them to the gas phase."""
    problem = build_problem(
        with_vent=False,
        element_release_model=(
            FirstOrderElementReleaseModel(
                time_constant=10.0,
            )
        ),
    )

    results = ScipySolver().solve(
        problem
    )

    pending_carbon = results.get_variable(
        "pending_element_C"
    )

    gas_carbon = results.get_variable(
        "element_C"
    )

    assert np.max(
        pending_carbon
    ) > 0.0

    assert np.all(
        pending_carbon >= -1.0e-10
    )

    assert np.all(
        gas_carbon >= -1.0e-10
    )

    generation_model = (
        problem.element_generation_model
    )

    assert generation_model is not None

    temperature = results.get_variable(
        "temperature"
    )

    time = results.time

    reaction_count = len(
        problem.chemistry_backend
        .reaction_network
        .reactions
    )

    conversion_histories = [
        results.get_variable(
            f"conversion_{index}"
        )
        for index in range(
            reaction_count
        )
    ]

    carbon_generation_rate = np.empty(
        time.size,
        dtype=float,
    )

    for time_index in range(
        time.size
    ):
        conversions = [
            float(
                conversion_history[
                    time_index
                ]
            )
            for conversion_history
            in conversion_histories
        ]

        generation_rates = (
            generation_model.generation_rates(
                temperature=float(
                    temperature[
                        time_index
                    ]
                ),
                conversions=conversions,
                cell_mass=problem.cell.mass,
            )
        )

        carbon_generation_rate[
            time_index
        ] = generation_rates.get(
            "C",
            0.0,
        )

    generated_carbon = float(
        np.trapezoid(
            carbon_generation_rate,
            time,
        )
    )

    initial_inventory = (
        problem.initial_element_inventory
    )

    assert initial_inventory is not None

    initial_carbon = (
        initial_inventory.moles_of(
            "C"
        )
    )

    final_accounted_carbon = (
        pending_carbon[-1]
        + gas_carbon[-1]
    )

    expected_final_carbon = (
        initial_carbon
        + generated_carbon
    )

    assert final_accounted_carbon == pytest.approx(
        expected_final_carbon,
        rel=1.0e-3,
        abs=1.0e-8,
    )

def test_energy_scaled_generation_is_capped_in_solver():
    """Limit energy-scaled elemental generation to its full yield."""
    element_generation_model = (
        EnergyScaledElementGenerationModel(
            element_yields_per_cell_mass={
                "C": 2.0,
            },
            reference_reaction_energy_per_cell_mass=(
                100.0
            ),
        )
    )

    problem = replace(
        build_problem(
            with_vent=False,
        ),
        element_generation_model=(
            element_generation_model
        ),
    )

    results = ScipySolver().solve(
        problem
    )

    reaction_energy = results.get_variable(
        "reaction_energy"
    )

    carbon = results.get_variable(
        "element_C"
    )

    initial_inventory = (
        problem.initial_element_inventory
    )

    assert initial_inventory is not None

    initial_carbon = (
        initial_inventory.moles_of(
            "C"
        )
    )

    generated_carbon = (
        carbon[-1]
        - initial_carbon
    )

    maximum_carbon = (
        element_generation_model
        .element_yields_per_cell_mass["C"]
        * problem.cell.mass
    )

    maximum_generation_energy = (
        element_generation_model
        .reference_reaction_energy_per_cell_mass
        * problem.cell.mass
    )

    assert reaction_energy[-1] > (
        maximum_generation_energy
    )

    assert generated_carbon > 0.0

    assert generated_carbon <= (maximum_carbon * (1.0 + 1.0e-3) + 1.0e-8)

    assert generated_carbon >= (maximum_carbon * (1.0 - 1.0e-3))