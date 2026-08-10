from dataclasses import dataclass

from .problem import Problem
from .solver import AbstractSolver


@dataclass(slots=True)
class Simulation:
    problem: Problem

    solver: AbstractSolver

    def run(self):
        self.problem.validate()

        return self.solver.solve(self.problem)
