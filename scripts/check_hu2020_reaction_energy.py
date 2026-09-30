"""Check integrated reaction energy for the Hu 2020 reference case."""

from __future__ import annotations

import numpy as np

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
from liiontr.problems.thermal import (
    ThermalProblem,
)
from liiontr.solver.scipy_solver import (
    ScipySolver,
)
from liiontr.thermal.lumped import (
    LumpedThermalModel,
)


def main() -> None:
    """Run the Hu baseline and integrate its reaction heat."""
    cell = cell_21700_generic()

    reaction_network = (
        hu2020_reaction_network(
            cell
        )
    )

    backend = ReactionNetworkBackend(
        reaction_network=reaction_network,
        cell=cell,
    )

    problem = ThermalProblem(
        cell=cell,
        chemistry_backend=backend,
        initial_temperature=480.0,
        initial_conversions=(
            hu2020_initial_conversions()
        ),
        ambient_temperature=298.15,
        convection_coefficient=10.0,
        duration=20.0,
    )

    results = ScipySolver().solve(
        problem
    )

    time = results.time

    temperature = results.get_variable(
        "temperature"
    )

    reaction_count = len(
        reaction_network.reactions
    )

    conversion_histories = [
        results.get_variable(
            f"conversion_{index}"
        )
        for index in range(
            reaction_count
        )
    ]

    heat_generation = np.empty(
        time.size,
        dtype=float,
    )

    for index in range(
        time.size
    ):
        conversions = [
            float(
                conversion_history[index]
            )
            for conversion_history
            in conversion_histories
        ]

        heat_generation[index] = (
            backend.heat_generation(
                temperature=float(
                    temperature[index]
                ),
                conversions=conversions,
            )
        )

    reaction_energy = float(
        np.trapezoid(
            heat_generation,
            time,
        )
    )

    specific_reaction_energy = (
        reaction_energy
        / cell.mass
    )

    thermal_model = (
        LumpedThermalModel(
            cell=cell,
            convection_coefficient=(
                problem.convection_coefficient
            ),
            ambient_temperature=(
                problem.ambient_temperature
            ),
        )
    )

    heat_loss = np.array(
        [
            thermal_model.heat_loss(
                float(value)
            )
            for value in temperature
        ],
        dtype=float,
    )

    thermal_loss_energy = float(
        np.trapezoid(
            heat_loss,
            time,
        )
    )

    sensible_energy = (
        cell.thermal_capacity
        * (
            float(temperature[-1])
            - problem.initial_temperature
        )
    )

    balance_residual = (
        reaction_energy
        - sensible_energy
        - thermal_loss_energy
    )

    relative_residual = (
        balance_residual
        / reaction_energy
    )

    print(
        "Hu 2020 reaction-energy check"
    )

    print(
        "-----------------------------"
    )

    print(
        f"Cell mass:             "
        f"{cell.mass:.9f} kg"
    )

    print(
        f"Initial temperature:   "
        f"{temperature[0]:.6f} K"
    )

    print(
        f"Final temperature:     "
        f"{temperature[-1]:.6f} K"
    )

    print(
        f"Maximum temperature:   "
        f"{np.max(temperature):.6f} K"
    )

    print()

    print(
        f"Reaction energy:        "
        f"{reaction_energy:.6f} J"
    )

    print(
        f"Reaction energy / mass: "
        f"{specific_reaction_energy:.6f} J/kg"
    )

    print(
        f"                         "
        f"{specific_reaction_energy / 1000.0:.6f} kJ/kg"
    )

    print()

    print(
        f"Sensible energy:        "
        f"{sensible_energy:.6f} J"
    )

    print(
        f"Thermal loss energy:    "
        f"{thermal_loss_energy:.6f} J"
    )

    print(
        f"Balance residual:       "
        f"{balance_residual:.6f} J"
    )

    print(
        f"Relative residual:      "
        f"{100.0 * relative_residual:.6f} %"
    )

    print()

    print(
        "Final conversions"
    )

    print(
        "-----------------"
    )

    for index, reaction in enumerate(
        reaction_network.reactions
    ):
        print(
            f"{reaction.name:<30} "
            f"{conversion_histories[index][-1]:.8f}"
        )


if __name__ == "__main__":
    main()