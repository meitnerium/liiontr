"""Run the baseline Hu et al. (2020) thermal runaway model."""

from __future__ import annotations

import numpy as np
from scipy.interpolate import PchipInterpolator

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
    """Run the Hu 2020 baseline simulation."""
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

    problem = ThermalProblem(
        cell=cell,
        chemistry_backend=chemistry_backend,
        initial_temperature=480.0,
        initial_conversions=(
            hu2020_initial_conversions()
        ),
        ambient_temperature=298.15,
        convection_coefficient=10.0,
        duration=20.0,
        maximum_temperature=1200.0,
    )

    results = ScipySolver().solve(
        problem
    )

    temperature = results.get_variable(
        "temperature"
    )

    time = np.asarray(
        results.time,
        dtype=float,
    )

    temperature = np.asarray(
        results.get_variable("temperature"),
        dtype=float,
    )

    thermal_model = LumpedThermalModel(
        cell=cell,
        convection_coefficient=10.0,
        ambient_temperature=298.15,
    )

    heat_loss = np.asarray(
        [thermal_model.heat_loss(float(value)) for value in temperature],
        dtype=float,
    )

    heat_loss_interpolator = PchipInterpolator(
        time,
        heat_loss,
    )

    thermal_loss_energy = float(
        heat_loss_interpolator.integrate(
            float(time[0]),
            float(time[-1]),
        )
    )

    sensible_energy = cell.thermal_capacity * (
        float(temperature[-1]) - float(temperature[0])
    )

    reaction_energy = sensible_energy + thermal_loss_energy

    print("Hu 2020 baseline")
    print("----------------")
    print(
        f"Initial temperature: "
        f"{temperature[0]:.6f} K"
    )
    print(
        f"Final temperature:   "
        f"{temperature[-1]:.6f} K"
    )
    print(
        f"Maximum temperature: "
        f"{temperature.max():.6f} K"
    )
    print(
        f"Final time:          "
        f"{results.time[-1]:.6f} s"
    )

    print()
    print("Final conversions")
    print("-----------------")

    for index, reaction in enumerate(
        reaction_network.reactions
    ):
        conversion = results.get_variable(
            f"conversion_{index}"
        )

        print(
            f"{reaction.name:<30} "
            f"{conversion[-1]:.8f}"
        )
        


    print()
    print("Energy balance")
    print("--------------")

    print(f"Sensible energy:      {sensible_energy:.6f} J")

    print(f"Thermal loss energy:  {thermal_loss_energy:.6f} J")

    print(f"Reaction energy:      {reaction_energy:.6f} J")

    print(f"Reaction energy/mass: {reaction_energy / cell.mass:.6f} J/kg")

if __name__ == "__main__":
    main()