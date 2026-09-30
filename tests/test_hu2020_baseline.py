"""Regression test for the Hu et al. (2020) thermal model."""

from __future__ import annotations

import pytest

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


def test_hu2020_baseline_regression():
    """Preserve the Hu 2020 baseline thermal response."""
    cell = cell_21700_generic()

    network = hu2020_reaction_network(
        cell
    )

    backend = ReactionNetworkBackend(
        reaction_network=network,
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
        maximum_temperature=1200.0,
    )

    results = ScipySolver().solve(
        problem
    )

    temperature = results.get_variable(
        "temperature"
    )

    assert temperature[0] == pytest.approx(
        480.0
    )

    assert temperature[-1] == pytest.approx(
        856.786897,
        rel=1.0e-6,
    )

    assert temperature.max() == pytest.approx(
        866.243743,
        rel=1.0e-6,
    )

    for index in range(4):
        conversion = results.get_variable(
            f"conversion_{index}"
        )

        assert conversion[-1] == pytest.approx(
            1.0,
            abs=1.0e-6,
        )