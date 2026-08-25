"""Base problem definitions for LiionTR simulations."""

from __future__ import annotations

from abc import ABC
from typing import Any

from .domain import Domain


class Problem(ABC):
    """Base class for physical problems.

    A problem defines what must be solved, including the physical
    domain, model parameters, physics components, and initial state.

    The problem does not define how the governing equations are
    solved. Numerical solution is delegated to a solver.

    Parameters
    ----------
    parameters : Any, optional
        Parameter collection associated with the problem.
    domain : Domain, optional
        Computational domain.
    physics : Any, optional
        Iterable of physics models initially associated with the
        problem.
    initial_state : Any, optional
        Initial state of the physical system.
    """

    def __init__(
        self,
        *,
        parameters: Any = None,
        domain: Domain | None = None,
        physics: Any = None,
        initial_state: Any = None,
    ) -> None:
        """Initialize a physical problem.

        Parameters
        ----------
        parameters : Any, optional
            Parameter collection associated with the problem.
        domain : Domain, optional
            Computational domain.
        physics : Any, optional
            Iterable of physics models.
        initial_state : Any, optional
            Initial state of the physical system.
        """
        self.parameters = parameters

        self.domain = domain

        self.physics = physics or []

        self.initial_state = initial_state

    def add_physics(
        self,
        physics_model: Any,
    ) -> None:
        """Add a physics model to the problem.

        Parameters
        ----------
        physics_model : Any
            Physics model to append to the problem.
        """
        self.physics.append(physics_model)

    def validate(self) -> bool:
        """Validate the problem definition.

        Returns
        -------
        bool
            ``True`` when the problem definition is valid.

        Notes
        -----
        The base implementation performs no validation and always
        returns ``True``. Specialized problem classes may override
        this method.
        """
        return True
