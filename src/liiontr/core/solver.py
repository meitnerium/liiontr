"""Abstract numerical solver interfaces."""

from __future__ import annotations

from abc import ABC, abstractmethod

from .problem import Problem
from .results import Results


class AbstractSolver(ABC):
    """
    Base class for numerical solvers.
    """

    @abstractmethod
    def solve(
        self,
        problem: Problem,
    ) -> Results:
        """
        Solve a physical problem.
        """
        raise NotImplementedError
