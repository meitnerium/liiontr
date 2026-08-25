"""High-level simulation orchestration."""

from dataclasses import dataclass

from .problem import Problem
from .solver import AbstractSolver


@dataclass(slots=True)
class Simulation:
    """Couple a physical problem with a numerical solver.

    Parameters
    ----------
    problem : Problem
        Physical problem to solve.
    solver : AbstractSolver
        Numerical solver used to integrate or otherwise solve the
        problem.
    """

    problem: Problem

    solver: AbstractSolver

    def run(self):
        """Validate and solve the configured problem.

        Returns
        -------
        Any
            Result returned by the configured solver.

        Raises
        ------
        Exception
            Propagates any exception raised during problem validation
            or numerical solution.
        """
        self.problem.validate()

        return self.solver.solve(self.problem)
