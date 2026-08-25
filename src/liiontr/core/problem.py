"""Base problem definitions for LiionTR simulations."""

from __future__ import annotations
from .domain import Domain
from abc import ABC  # , abstractmethod


class Problem(ABC):
    """
    Base class for all physical problems.

    A Problem defines WHAT must be solved.
    It does not define HOW it is solved.
    """

    def __init__(
        self,
        *,
        parameters=None,
        domain: Domain | None = None,
        physics=None,
        initial_state=None,
    ):
        self.parameters = parameters

        self.domain = domain

        self.physics = physics or []

        self.initial_state = initial_state

    def add_physics(self, physics_model):
        self.physics.append(physics_model)

    def validate(self):
        """
        Validate problem definition.
        """

        return True
