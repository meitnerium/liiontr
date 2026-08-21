import numpy as np
import pytest

from liiontr.chemistry import ReactionNetworkBackend
from liiontr.library import cell_21700_generic
from liiontr.library.hu2020 import (
    hu2020_initial_conversions,
    hu2020_reaction_network,
)
from liiontr.problems import ThermalProblem
from liiontr.solver import ScipySolver


def test_hu2020_energy_balance():
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

    conversions = [
        results.get("conversion_0"),
        results.get("conversion_1"),
        results.get("conversion_2"),
        results.get("conversion_3"),
    ]

    total_heat_generation = []

    for index, current_temperature in enumerate(temperature):
        current_conversions = [float(conversion[index]) for conversion in conversions]

        total_heat_generation.append(
            network.heat_generation(
                temperature=float(current_temperature),
                conversions=current_conversions,
            )
        )

    generated_energy = np.trapezoid(
        total_heat_generation,
        results.time,
    )

    thermal_capacity = cell.thermal_capacity

    stored_energy = (
        thermal_capacity * (float(temperature[-1]) - float(temperature[0])) / cell.mass
    )

    assert generated_energy > stored_energy

    assert generated_energy == pytest.approx(
        386000.0,
        rel=0.02,
    )

    assert stored_energy == pytest.approx(
        377000.0,
        rel=0.02,
    )
