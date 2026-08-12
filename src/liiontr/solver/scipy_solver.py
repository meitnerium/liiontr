from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from liiontr.chemistry import ReactionNetworkBackend
from liiontr.core.results import Results
from liiontr.problems.thermal import ThermalProblem
from liiontr.thermal.lumped import LumpedThermalModel


@dataclass(slots=True)
class ScipySolver:
    """
    SciPy-based ODE solver.

    BDF is used by default because thermal runaway kinetics
    can become stiff at elevated temperatures.
    """

    method: str = "BDF"

    relative_tolerance: float = 1.0e-6

    absolute_tolerance: float = 1.0e-9

    def solve(
        self,
        problem: ThermalProblem,
    ) -> Results:
        model = LumpedThermalModel(
            cell=problem.cell,
            convection_coefficient=problem.convection_coefficient,
            ambient_temperature=problem.ambient_temperature,
        )

        backend = problem.chemistry_backend

        n_reactions = 0

        if isinstance(
            backend,
            ReactionNetworkBackend,
        ):
            n_reactions = len(backend.reaction_network.reactions)

        if problem.initial_conversions is None:
            initial_conversions = [0.0] * n_reactions

        else:
            initial_conversions = list(problem.initial_conversions)

            if len(initial_conversions) != n_reactions:
                raise ValueError(
                    "Number of initial conversions must match number of reactions."
                )

            if any(
                conversion < 0.0 or conversion > 1.0
                for conversion in initial_conversions
            ):
                raise ValueError("Initial conversions must be between 0 and 1.")

        initial_state = [
            problem.initial_temperature,
            *initial_conversions,
        ]

        def rhs(
            time: float,
            state: np.ndarray,
        ) -> list[float]:
            del time

            temperature = float(state[0])

            conversions = [
                min(
                    1.0,
                    max(
                        0.0,
                        float(value),
                    ),
                )
                for value in state[1 : 1 + n_reactions]
            ]

            if backend is None:
                heat_generation = 0.0
                progress_rates: list[float] = []

            elif isinstance(
                backend,
                ReactionNetworkBackend,
            ):
                heat_generation = backend.heat_generation(
                    temperature=temperature,
                    conversions=conversions,
                )

                progress_rates = backend.progress_rates(
                    temperature=temperature,
                    conversions=conversions,
                )

            else:
                heat_generation = backend.heat_generation(temperature)

                progress_rates = []

            temperature_rate = model.temperature_derivative(
                temperature=temperature,
                heat_generation=heat_generation,
            )

            return [
                temperature_rate,
                *progress_rates,
            ]

        events: Any = None

        maximum_temperature = problem.maximum_temperature

        if maximum_temperature is not None:

            def maximum_temperature_event(
                time: float,
                state: np.ndarray,
            ) -> float:
                del time

                return maximum_temperature - float(state[0])

            maximum_temperature_event.terminal = True  # type: ignore[attr-defined]
            maximum_temperature_event.direction = -1.0  # type: ignore[attr-defined]

            events = maximum_temperature_event

        solution = solve_ivp(
            rhs,
            (0.0, problem.duration),
            initial_state,
            method=self.method,
            rtol=self.relative_tolerance,
            atol=self.absolute_tolerance,
            events=events,
        )

        if not solution.success:
            raise RuntimeError(f"ODE integration failed: {solution.message}")

        results = Results(
            time=solution.t,
        )

        results.add_variable(
            "temperature",
            solution.y[0],
        )

        for index in range(n_reactions):
            conversion = np.clip(
                solution.y[index + 1],
                0.0,
                1.0,
            )

            results.add_variable(
                f"conversion_{index}",
                conversion,
            )

        return results
