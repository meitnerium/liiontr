import pytest

from liiontr.chemistry import ReactionNetworkBackend
from liiontr.library import cell_21700_generic
from liiontr.library.hu2020 import (
    hu2020_initial_conversions,
    hu2020_reaction_network,
)
from liiontr.problems import ThermalProblem
from liiontr.solver import ScipySolver


def test_hu2020_network_runs_in_thermal_solver():
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
        initial_temperature=450.0,
        initial_conversions=hu2020_initial_conversions(),
        ambient_temperature=298.15,
        convection_coefficient=10.0,
        duration=5.0,
        maximum_temperature=1200.0,
    )

    results = ScipySolver().solve(problem)

    temperature = results.get("temperature")

    sei_conversion = results.get("conversion_0")

    anode_conversion = results.get("conversion_1")

    cathode_conversion = results.get("conversion_2")

    electrolyte_conversion = results.get("conversion_3")

    assert temperature[0] == pytest.approx(450.0)

    assert sei_conversion[0] == pytest.approx(0.85)

    assert anode_conversion[0] == pytest.approx(0.25)

    assert cathode_conversion[0] == pytest.approx(0.04)

    assert electrolyte_conversion[0] == pytest.approx(0.0)

    assert temperature[-1] > temperature[0]

    assert sei_conversion[-1] > sei_conversion[0]

    assert anode_conversion[-1] > anode_conversion[0]

    assert cathode_conversion[-1] > cathode_conversion[0]

    assert sei_conversion.min() >= 0.0
    assert sei_conversion.max() <= 1.0

    assert anode_conversion.min() >= 0.0
    assert anode_conversion.max() <= 1.0

    assert cathode_conversion.min() >= 0.0
    assert cathode_conversion.max() <= 1.0

    assert electrolyte_conversion.min() >= 0.0
    assert electrolyte_conversion.max() <= 1.0


def test_hu2020_electrolyte_remains_inactive_below_threshold():
    cell = cell_21700_generic()

    network = hu2020_reaction_network(
        cell=cell,
    )

    conversions = hu2020_initial_conversions()

    rates = network.progress_rates(
        temperature=470.0,
        conversions=conversions,
    )

    electrolyte_rate = rates[3]

    assert electrolyte_rate == pytest.approx(0.0)


def test_hu2020_electrolyte_is_active_above_threshold():
    cell = cell_21700_generic()

    network = hu2020_reaction_network(
        cell=cell,
    )

    conversions = hu2020_initial_conversions()

    rates = network.progress_rates(
        temperature=480.0,
        conversions=conversions,
    )

    electrolyte_rate = rates[3]

    assert electrolyte_rate > 0.0
