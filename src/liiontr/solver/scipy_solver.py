from dataclasses import dataclass

from scipy.integrate import solve_ivp

from liiontr.core.results import Results
from liiontr.problems.thermal import ThermalProblem
from liiontr.thermal.lumped import LumpedThermalModel


@dataclass(slots=True)
class ScipySolver:
    def solve(
        self,
        problem: ThermalProblem,
    ) -> Results:
        model = LumpedThermalModel(
            cell=problem.cell,
            convection_coefficient=problem.convection_coefficient,
            ambient_temperature=problem.ambient_temperature,
        )

        def rhs(t, y):
            if problem.chemistry_backend is None:
                heat_generation = 0.0
            else:
                heat_generation = problem.chemistry_backend.heat_generation(y[0])

            return [
                model.temperature_derivative(
                    y[0],
                    heat_generation,
                )
            ]

        solution = solve_ivp(
            rhs,
            (0.0, problem.duration),
            [problem.initial_temperature],
        )

        results = Results(
            time=solution.t,
        )

        results.add_variable(
            "temperature",
            solution.y[0],
        )

        return results
