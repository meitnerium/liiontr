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
        0.5 * (value_0 + value_1) * (time_1 - time_0)
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

    results = ScipySolver().solve(problem)

    temperature = results.get("temperature")
    sei_conversion = results.get("conversion_0")
    anode_conversion = results.get("conversion_1")
    cathode_conversion = results.get("conversion_2")
    electrolyte_conversion = results.get("conversion_3")

    if temperature is None:
        raise RuntimeError("Temperature results are missing.")

    if sei_conversion is None:
        raise RuntimeError("SEI conversion results are missing.")

    if anode_conversion is None:
        raise RuntimeError("Anode conversion results are missing.")

    if cathode_conversion is None:
        raise RuntimeError("Cathode conversion results are missing.")

    if electrolyte_conversion is None:
        raise RuntimeError("Electrolyte conversion results are missing.")

    if results.time is None:
        raise RuntimeError("Simulation time results are missing.")

    time = [float(value) for value in results.time]

    reaction_names = [reaction.name for reaction in network.reactions]

    heat_history: dict[str, list[float]] = {name: [] for name in reaction_names}

    for index, current_temperature in enumerate(temperature):
        conversions = [
            float(sei_conversion[index]),
            float(anode_conversion[index]),
            float(cathode_conversion[index]),
            float(electrolyte_conversion[index]),
        ]

        heat_breakdown = network.heat_generation_by_reaction(
            temperature=float(current_temperature),
            conversions=conversions,
        )

        for name in reaction_names:
            heat_history[name].append(heat_breakdown[name])

    print()
    print("Hu 2020 thermal runaway simulation")
    print("----------------------------------")

    print(f"Simulation time: {time[-1]:.6g} s")

    print(f"Temperature: {temperature[0]:.3f} -> {temperature[-1]:.3f} K")

    print()

    print(f"SEI conversion: {sei_conversion[0]:.6f} -> {sei_conversion[-1]:.6f}")

    print(f"Anode conversion: {anode_conversion[0]:.6f} -> {anode_conversion[-1]:.6f}")

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
        heat_values = heat_history[name]

        peak_heat = max(heat_values)

        energy = integrate_trapezoidal(
            time=time,
            values=heat_values,
        )

        total_energy += energy

        print()
        print(name)

        print(f"  Peak heat generation: {peak_heat:.3f} W/kg")

        print(f"  Released energy: {energy / 1000.0:.3f} kJ/kg")

    print()

    print(f"Total reaction energy: {total_energy / 1000.0:.3f} kJ/kg")

    print()

    print(f"Number of solver points: {len(time)}")


if __name__ == "__main__":
    main()
