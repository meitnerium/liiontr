"""Simulation result containers and result access utilities."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(slots=True)
class Results:
    """
    Generic simulation results container.

    Stores time evolution of physical variables.
    """

    time: np.ndarray | None = None

    variables: dict[str, np.ndarray] = field(default_factory=dict)

    def add_variable(
        self,
        name: str,
        values: np.ndarray,
    ) -> None:
        """
        Add a simulation variable.
        """

        self.variables[name] = values

    def get_variable(
        self,
        name: str,
    ) -> np.ndarray:
        """
        Retrieve a stored variable.
        """

        return self.variables[name]

    def get(
        self,
        name: str,
    ) -> np.ndarray:
        """
        Short alias for get_variable().
        """

        return self.get_variable(name)

    @property
    def temperature(self) -> np.ndarray:
        """
        Return temperature history.
        """

        return self.variables["temperature"]

    @property
    def y(self) -> np.ndarray:
        """
        Compatibility with scipy.integrate.solve_ivp output.
        """

        if not self.variables:
            raise ValueError("No variables stored in Results")

        return np.vstack(list(self.variables.values()))
