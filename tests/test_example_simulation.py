from liiontr.library import cell_21700_generic
from liiontr.problems import ThermalProblem
from liiontr.solver.scipy_solver import ScipySolver


def test_first_simulation():
    cell = cell_21700_generic()

    problem = ThermalProblem(
        cell=cell,
        initial_temperature=350.0,
        duration=100.0,
    )

    solution = ScipySolver().solve(problem)

    assert solution.y[0][-1] < 350.0
