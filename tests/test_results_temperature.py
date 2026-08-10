from liiontr.library import cell_21700_generic
from liiontr.problems import ThermalProblem
from liiontr.solver.scipy_solver import ScipySolver


def test_temperature_results():
    problem = ThermalProblem(
        cell=cell_21700_generic(),
        initial_temperature=350.0,
        duration=100.0,
    )

    results = ScipySolver().solve(problem)

    assert results.temperature[0] == 350.0
    assert results.temperature[-1] < 350.0
