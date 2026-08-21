from __future__ import annotations

from liiontr.chemistry import ReactionNetworkBackend
from liiontr.library import cell_21700_generic
from liiontr.library.hu2020 import (
    hu2020_initial_conversions,
    hu2020_reaction_network,
)
from liiontr.problems import ThermalProblem
from liiontr.solver import ScipySolver


def integrate_trapezoidal(
    time: list[float],
    values: list[float],
) -> float:
    """
    Integrate values over time using the trapezoidal rule.
    """

    return sum(
        0.5
        * (value_0 + value_1)
        * (time_1 - time_0)
        for time_0, time_1, value_0, value_1 in zip(
            time[:-1],
            time[1:],
            values[:-1],
            values[1:],
            strict=True,
        )
    )


def main() -> None:
    cell = cell_21700_generic()

    network = hu2020_reaction_network(
        cell=cell,
    )

    backend = ReactionNetworkBackend(
        reaction_network=network,
        cell=cell,
    )

    problem = ThermalProblem(
        cell=cell,
        chemistry_backend=backend,
        initial_temperature=480.0,
        initial_conversions=hu2020_initial_conversions(),
        ambient_temperature=298.15,
        convection_coefficient=10.0,
        duration=20.0,
        maximum_temperature=1200.0,
    )

    results = ScipySolver().solve(
        problem
    )

    temperature = results.get(
        "temperature"
    )

    sei_conversion = results.get(
        "conversion_0"
    )

    anode_conversion = results.get(
        "conversion_1"
    )

    cathode_conversion = results.get(
        "conversion_2"
    )

    electrolyte_conversion = results.get(
        "conversion_3"
    )

    reaction_names = [
        reaction.name
        for reaction in network.reactions
    ]

    heat_history = {
        name: []
        for name in reaction_names
    }

    for index, current_temperature in enumerate(
        temperature
    ):
        conversions = [
            float(sei_conversion[index]),
            float(anode_conversion[index]),
            float(cathode_conversion[index]),
            float(electrolyte_conversion[index]),
        ]

        heat = network.heat_generation_by_reaction(
            temperature=float(current_temperature),
            conversions=conversions,
        )

        for name in reaction_names:
            heat_history[name].append(
                heat[name]
            )

    time = [
        float(value)
        for value in results.time
    ]

    print()
    print("Hu 2020 thermal runaway simulation")
    print("----------------------------------")

    print(
        f"Simulation time: "
        f"{results.time[-1]:.6g} s"
    )

    print(
        f"Temperature: "
        f"{temperature[0]:.3f} -> "
        f"{temperature[-1]:.3f} K"
    )

    print()

    print(
        f"SEI conversion: "
        f"{sei_conversion[0]:.6f} -> "
        f"{sei_conversion[-1]:.6f}"
    )

    print(
        f"Anode conversion: "
        f"{anode_conversion[0]:.6f} -> "
        f"{anode_conversion[-1]:.6f}"
    )

    print(
        f"Cathode conversion: "
        f"{cathode_conversion[0]:.6f} -> "
        f"{cathode_conversion[-1]:.6f}"
    )

    print(
        f"Electrolyte conversion: "
        f"{electrolyte_conversion[0]:.6f} -> "
        f"{electrolyte_conversion[-1]:.6f}"
    )

    print()
    print("Heat generation by reaction")
    print("---------------------------")

    total_energy = 0.0

    for name in reaction_names:
        heat = heat_history[name]

        peak_heat = max(
            heat
        )

        energy = integrate_trapezoidal(
            time=time,
            values=heat,
        )

        total_energy += energy

        print()
        print(name)

        print(
            f"  Peak heat generation: "
            f"{peak_heat:.3f} W/kg"
        )

        print(
            f"  Released energy: "
            f"{energy / 1000.0:.3f} kJ/kg"
        )

    print()
    print(
        f"Total reaction energy: "
        f"{total_energy / 1000.0:.3f} kJ/kg"
    )

    print()

    print(
        f"Number of solver points: "
        f"{len(results.time)}"
    )


if __name__ == "__main__":
    main()