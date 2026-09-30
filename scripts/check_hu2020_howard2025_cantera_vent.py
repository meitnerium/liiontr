"""Run the coupled Hu-Howard-Cantera model with venting."""

from __future__ import annotations

import os

import numpy as np

from liiontr.chemistry.cantera import (
    CanteraEquilibriumBackend,
)
from liiontr.chemistry.element_release import (
    FirstOrderElementReleaseModel,
)
from liiontr.chemistry.elements import (
    ElementInventory,
)
from liiontr.chemistry.reaction_backend import (
    ReactionNetworkBackend,
)
from liiontr.gases.element_vent import (
    ElementVentFlowModel,
)
from liiontr.library.cells import (
    cell_21700_generic,
)
from liiontr.library.hu2020 import (
    hu2020_initial_conversions,
    hu2020_reaction_network,
)
from liiontr.library.hu2020_howard2025 import (
    hu2020_howard2025_element_generation_model,
)
from liiontr.problems.thermal import (
    ThermalProblem,
)
from liiontr.solver.scipy_solver import (
    ScipySolver,
)

IDEAL_GAS_CONSTANT = 8.31446261815324

RELEASE_TIME_CONSTANT = float(
    os.environ.get(
        "LIIONTR_RELEASE_TIME_CONSTANT",
        "1.0",
    )
)


def main() -> None:
    """Run the vented Hu-Howard-Cantera simulation."""
    cell = cell_21700_generic()

    reaction_network = (
        hu2020_reaction_network(
            cell
        )
    )

    chemistry_backend = (
        ReactionNetworkBackend(
            reaction_network=reaction_network,
            cell=cell,
        )
    )

    element_generation_model = (
        hu2020_howard2025_element_generation_model()
    )

    FirstOrderElementReleaseModel(
        time_constant=(RELEASE_TIME_CONSTANT),
    )
    
    equilibrium_backend = (
        CanteraEquilibriumBackend()
    )

    initial_temperature = 480.0
    initial_pressure = 101325.0

    gas_volume = 5.0e-6

    vent_open_pressure = 5.0e5

    initial_n2_moles = (
        initial_pressure
        * gas_volume
        / (
            IDEAL_GAS_CONSTANT
            * initial_temperature
        )
    )

    initial_element_inventory = (
        ElementInventory(
            moles={
                "N": 2.0 * initial_n2_moles,
            }
        )
    )

    element_vent_model = (
        ElementVentFlowModel(
            equilibrium_backend=equilibrium_backend,
            vent_area=1.0e-6,
            discharge_coefficient=0.8,
            downstream_pressure=101325.0,
        )
    )

    element_release_model = FirstOrderElementReleaseModel(
        time_constant=(RELEASE_TIME_CONSTANT),
    )

    problem = ThermalProblem(
        cell=cell,
        chemistry_backend=chemistry_backend,
        element_generation_model=(element_generation_model),
        element_release_model=(element_release_model),
        gas_equilibrium_backend=(equilibrium_backend),
        initial_element_inventory=(initial_element_inventory),
        gas_volume=gas_volume,
        element_vent_model=(element_vent_model),
        vent_open_pressure=(vent_open_pressure),
        initial_temperature=(initial_temperature),
        initial_conversions=(hu2020_initial_conversions()),
        ambient_temperature=298.15,
        convection_coefficient=10.0,
        duration=20.0,
        maximum_temperature=2000.0,
    )

    results = ScipySolver().solve(
        problem
    )

    temperature = results.get_variable(
        "temperature"
    )

    pressure = results.get_variable(
        "pressure"
    )

    reaction_energy = results.get_variable("reaction_energy")

    time = results.time

    reaction_network = chemistry_backend.reaction_network

    reaction_names = [reaction.name for reaction in reaction_network.reactions]

    conversion_histories = [
        results.get_variable(f"conversion_{index}")
        for index in range(len(reaction_names))
    ]

    reaction_heat_rates = {
        reaction_name: np.empty(
            time.size,
            dtype=float,
        )
        for reaction_name in reaction_names
    }

    for time_index in range(time.size):
        conversions = [
            float(conversion_history[time_index])
            for conversion_history in conversion_histories
        ]

        heat_by_reaction = reaction_network.heat_generation_by_reaction(
            temperature=float(temperature[time_index]),
            conversions=conversions,
        )

        for reaction_name in reaction_names:
            reaction_heat_rates[reaction_name][time_index] = (
                heat_by_reaction[reaction_name] * cell.mass
            )

    contained_carbon = results.get_variable("element_C")

    vented_carbon = results.get_variable("vented_element_C")

    contained_carbon = results.get_variable("element_C")

    pending_carbon = results.get_variable("pending_element_C")

    vented_carbon = results.get_variable("vented_element_C")

    total_generated_carbon = pending_carbon + contained_carbon + vented_carbon

    maximum_generation_energy = (
        element_generation_model.reference_reaction_energy_per_cell_mass * cell.mass
    )

    final_reaction_energy = float(reaction_energy[-1])

    reaction_energy_per_cell_mass = final_reaction_energy / cell.mass

    reaction_energies = {
        reaction_name: float(
            np.trapezoid(
                reaction_heat_rates[reaction_name],
                time,
            )
        )
        for reaction_name in reaction_names
    }

    postprocessed_reaction_energy = sum(reaction_energies.values())

    postprocessed_reaction_energy_per_mass = postprocessed_reaction_energy / cell.mass

    quadrature_difference = postprocessed_reaction_energy - final_reaction_energy

    quadrature_difference_per_mass = quadrature_difference / cell.mass

    gas_generation_progress = element_generation_model.generation_progress(
        released_reaction_energy=(final_reaction_energy),
        cell_mass=cell.mass,
    )

    maximum_carbon_yield = (
        element_generation_model.element_yields_per_cell_mass["C"] * cell.mass
    )

    final_generated_carbon = float(total_generated_carbon[-1])

    carbon_yield_fraction = final_generated_carbon / maximum_carbon_yield

    total_released_carbon = contained_carbon + vented_carbon

    reaction_energy = np.asarray(
        results.get_variable("reaction_energy"),
        dtype=float,
    )

    maximum_generation_energy = (
        element_generation_model.reference_reaction_energy_per_cell_mass * cell.mass
    )

    carbon_generation_rate = np.zeros_like(
        time,
        dtype=float,
    )

    for time_index in range(time.size):
        conversions = [
            float(results.get_variable(f"conversion_{reaction_index}")[time_index])
            for reaction_index in range(4)
        ]

        heat_generation = chemistry_backend.heat_generation(
            temperature=float(temperature[time_index]),
            conversions=conversions,
        )

        released_reaction_energy = max(
            float(reaction_energy[time_index]),
            0.0,
        )

        if released_reaction_energy >= maximum_generation_energy:
            carbon_generation_rate[time_index] = 0.0

        else:
            generation_rates = element_generation_model.generation_rates(
                heat_generation=(heat_generation),
            )

            carbon_generation_rate[time_index] = generation_rates.get(
                "C",
                0.0,
            )

    carbon_release_rate = np.gradient(
        total_released_carbon,
        time,
    )

    maximum_pending_index = int(np.argmax(pending_carbon))

    maximum_release_index = int(np.argmax(carbon_release_rate))

    carbon_vent_rate = np.gradient(
        vented_carbon,
        time,
    )

    carbon_accumulation_rate = np.gradient(
        contained_carbon,
        time,
    )

    # carbon_pending_rate = np.gradient(
    #    pending_carbon,
    #    time,
    # )
    
    maximum_pressure_index = int(np.argmax(pressure))

    maximum_pressure_time = float(time[maximum_pressure_index])

    final_conversions = [
        float(results.get_variable(f"conversion_{index}")[-1]) for index in range(4)
    ]

    print(
        "Hu 2020 + Howard 2025 + Cantera + venting"
    )
    print(
        "------------------------------------------"
    )

    print(
        f"Cell mass:            "
        f"{cell.mass:.6f} kg"
    )

    print(
        f"Gas volume:           "
        f"{gas_volume * 1.0e6:.3f} mL"
    )

    print(
        f"Vent opening pressure:"
        f" {vent_open_pressure / 1.0e5:.3f} bar"
    )

    print(
        f"Vent area:            "
        f"{1.0e-6 * 1.0e6:.3f} mm2"
    )

    print(f"Release time constant: {RELEASE_TIME_CONSTANT:.3f} s")

    print()

    print(
        f"Initial temperature:  "
        f"{temperature[0]:.6f} K"
    )

    print(
        f"Final temperature:    "
        f"{temperature[-1]:.6f} K"
    )

    print(
        f"Maximum temperature:  "
        f"{temperature.max():.6f} K"
    )

    print()

    print(
        f"Initial pressure:     "
        f"{pressure[0] / 1.0e5:.6f} bar"
    )

    print(
        f"Final pressure:       "
        f"{pressure[-1] / 1.0e5:.6f} bar"
    )

    print(
        f"Maximum pressure:     "
        f"{pressure.max() / 1.0e5:.6f} bar"
    )

    print(f"Maximum pressure time: {maximum_pressure_time:.9f} s")

    print(
        f"Maximum release:    "
        f"{carbon_release_rate[maximum_release_index]:.8e} mol/s "
        f"at {time[maximum_release_index]:.9f} s"
    )

    print(
        f"Maximum pending C:  "
        f"{pending_carbon[maximum_pending_index]:.8e} mol "
        f"at {time[maximum_pending_index]:.9f} s"
    )

    if "vent_open" in results.variables:
        vent_open = results.get_variable("vent_open")

        open_indices = np.flatnonzero(vent_open > 0.5)

        if open_indices.size > 0:
            first_open_index = int(open_indices[0])

            print(f"Vent opening time:    {time[first_open_index]:.9f} s")

            print(f"Pressure at opening:  {pressure[first_open_index] / 1.0e5:.6f} bar")

    print()
    print("Pressure peak neighborhood")
    print("--------------------------")
    print("time [s]       T [K]       P [bar]       vent")

    start_index = max(
        maximum_pressure_index - 4,
        0,
    )

    end_index = min(
        maximum_pressure_index + 5,
        len(time),
    )

    for index in range(
        start_index,
        end_index,
    ):
        if "vent_open" in results.variables:
            vent_state = int(results.get_variable("vent_open")[index] > 0.5)
        else:
            vent_state = -1

        print(
            f"{time[index]:12.8f} "
            f"{temperature[index]:10.3f} "
            f"{pressure[index] / 1.0e5:12.3f} "
            f"{vent_state:4d}"
        )

    print()
    print("Carbon rate neighborhood")
    print("------------------------")
    print(
        "time [s]       "
        "generation       "
        "release          "
        "vent             "
        "gas accumulation"
    )

    for index in range(
        start_index,
        end_index,
    ):
        print(
            f"{time[index]:12.8f} "
            f"{carbon_generation_rate[index]:14.6e} "
            f"{carbon_release_rate[index]:14.6e} "
            f"{carbon_vent_rate[index]:14.6e} "
            f"{carbon_accumulation_rate[index]:14.6e}"
        )

    maximum_generation_index = int(np.argmax(carbon_generation_rate))

    maximum_vent_index = int(np.argmax(carbon_vent_rate))

    print()
    print("Carbon rate maxima")
    print("------------------")

    print(
        f"Maximum generation: "
        f"{carbon_generation_rate[maximum_generation_index]:.8e} mol/s "
        f"at {time[maximum_generation_index]:.9f} s"
    )

    print(
        f"Maximum vent rate:  "
        f"{carbon_vent_rate[maximum_vent_index]:.8e} mol/s "
        f"at {time[maximum_vent_index]:.9f} s"
    )

    print()
    print("Energy-scaled gas generation")
    print("----------------------------")

    print(f"Reaction energy:       {final_reaction_energy:.6f} J")

    print(f"Reaction energy/mass:  {reaction_energy_per_cell_mass:.6f} J/kg")

    print(f"Reference energy:      {maximum_generation_energy:.6f} J")

    print(
        f"Reference energy/mass: "
        f"{element_generation_model.reference_reaction_energy_per_cell_mass:.6f} J/kg"
    )

    print(f"Generation progress:   {gas_generation_progress:.9f}")

    print(f"Generated carbon:      {final_generated_carbon:.9f} mol-atoms")

    print(f"Maximum carbon yield:  {maximum_carbon_yield:.9f} mol-atoms")

    print(f"Carbon yield fraction: {carbon_yield_fraction:.9f}")

    print()
    print("Final Hu conversions")
    print("--------------------")

    for index, conversion in enumerate(final_conversions):
        print(f"Reaction {index}: {conversion:.12f}")

    print()

    print("Reaction-energy contributions")
    print("-----------------------------")

    for reaction_name in reaction_names:
        reaction_energy_value = reaction_energies[reaction_name]

        print(
            f"{reaction_name:<30s} "
            f"{reaction_energy_value / cell.mass / 1000.0:.6f} "
            f"kJ/kg"
        )

    print()

    print(
        f"Trapezoidal total:     "
        f"{postprocessed_reaction_energy_per_mass / 1000.0:.6f} "
        f"kJ/kg"
    )

    print(f"ODE-state total:       {reaction_energy_per_cell_mass / 1000.0:.6f} kJ/kg")

    print(f"Quadrature difference: {quadrature_difference_per_mass:.6f} J/kg")

    print()
    print("Final pending elemental inventory")
    print("---------------------------------")

    for element_name in (
        "C",
        "H",
        "O",
        "N",
    ):
        variable_name = f"pending_element_{element_name}"

        if variable_name in results.variables:
            values = results.get_variable(variable_name)

            print(f"{element_name:<3} {values[-1]:.8f} mol-atoms")

    print()
    print("Final contained elemental inventory")
    print("-----------------------------------")

    for element_name in (
        "C",
        "H",
        "O",
        "N",
    ):
        values = results.get_variable(
            f"element_{element_name}"
        )

        print(
            f"{element_name:<3} "
            f"{values[-1]:.8f} mol-atoms"
        )

    print()
    print("Cumulative vented elemental inventory")
    print("-------------------------------------")

    for element_name in (
        "C",
        "H",
        "O",
        "N",
    ):
        variable_name = (
            f"vented_element_{element_name}"
        )

        if variable_name in results.variables:
            values = results.get_variable(
                variable_name
            )

            print(
                f"{element_name:<3} "
                f"{values[-1]:.8f} mol-atoms"
            )

    if (
        "vented_thermal_energy"
        in results.variables
    ):
        vented_energy = results.get_variable(
            "vented_thermal_energy"
        )

        print()
        print(
            f"Vented thermal energy: "
            f"{vented_energy[-1]:.6f} J"
        )

    print()
    print("Final equilibrium species")
    print("-------------------------")

    for species_name in (
        "H2",
        "H2O",
        "CO",
        "CO2",
        "CH4",
        "N2",
    ):
        variable_name = (
            f"equilibrium_species_"
            f"{species_name}"
        )

        if variable_name in results.variables:
            values = results.get_variable(
                variable_name
            )

            print(
                f"{species_name:<5} "
                f"{values[-1]:.8f} mol"
            )


if __name__ == "__main__":
    main()