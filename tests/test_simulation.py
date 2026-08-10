from liiontr.core.problem import Problem
from liiontr.core.results import Results
from liiontr.core.simulation import Simulation
from liiontr.core.solver import AbstractSolver


class DummyProblem(Problem):
    pass


class DummySolver(AbstractSolver):
    def solve(self, problem: Problem) -> Results:
        return Results()


def test_simulation():
    simulation = Simulation(
        problem=DummyProblem(),
        solver=DummySolver(),
    )

    results = simulation.run()

    assert isinstance(results, Results)
