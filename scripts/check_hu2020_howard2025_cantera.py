"""Run the coupled Hu 2020, Howard 2025, and Cantera model."""

from __future__ import annotations

from liiontr.chemistry.cantera import (
    CanteraEquilibriumBackend,
)
from liiontr.chemistry.elements import (
    ElementInventory,
)
from liiontr.chemistry.reaction_backend import (
    ReactionNetworkBackend,
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


def main() -> None:
    """Run the closed Hu-Howard-Cantera simulation."""
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

    equilibrium_backend = (
        CanteraEquilibriumBackend()
    )

    initial_temperature = 480.0
    initial_pressure = 101325.0

    # Diagnostic control volume.
    #
    # This is currently a modeling parameter rather than
    # a literature-derived internal free volume.
    gas_volume = 1.0e-4

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

    problem = ThermalProblem(
        cell=cell,
        chemistry_backend=chemistry_backend,
        element_generation_model=(
            element_generation_model
        ),
        gas_equilibrium_backend=(
            equilibrium_backend
        ),
        initial_element_inventory=(
            initial_element_inventory
        ),
        gas_volume=gas_volume,
        initial_temperature=initial_temperature,
        initial_conversions=(
            hu2020_initial_conversions()
        ),
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

    print(
        "Hu 2020 + Howard 2025 + Cantera"
    )
    print(
        "--------------------------------"
    )

    print(
        f"Cell mass:            "
        f"{cell.mass:.6f} kg"
    )

    print(
        f"Gas volume:           "
        f"{gas_volume:.6e} m3"
    )

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

    print()
    print("Final elemental inventory")
    print("-------------------------")

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
        values = results.get_variable(
            f"equilibrium_species_"
            f"{species_name}"
        )

        print(
            f"{species_name:<5} "
            f"{values[-1]:.8f} mol"
        )


if __name__ == "__main__":
    main()