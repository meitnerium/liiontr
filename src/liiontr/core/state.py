"""Simulation state containers and state access utilities."""

from __future__ import annotations

from dataclasses import dataclass, field

from .variable import Variable


@dataclass(slots=True)
class State:
    """Represent the complete state of a physical system.

    Parameters
    ----------
    variables : dict[str, Variable]
        Mapping from variable names to simulation variables.
    """

    variables: dict[str, Variable] = field(default_factory=dict)

    def add(self, variable: Variable) -> None:
        """Add or replace a variable in the state.

        Parameters
        ----------
        variable : Variable
            Variable to store. Its ``name`` attribute is used as the
            dictionary key.
        """
        self.variables[variable.name] = variable

    def get(self, name: str) -> Variable:
        """Return a variable by name.

        Parameters
        ----------
        name : str
            Name of the requested variable.

        Returns
        -------
        Variable
            Stored simulation variable.

        Raises
        ------
        KeyError
            If no variable with the requested name exists.
        """
        return self.variables[name]

    def __contains__(self, name: str) -> bool:
        """Return whether a variable is present in the state.

        Parameters
        ----------
        name : str
            Variable name.

        Returns
        -------
        bool
            ``True`` if the variable exists, otherwise ``False``.
        """
        return name in self.variables

    def __len__(self) -> int:
        """Return the number of variables stored in the state."""
        return len(self.variables)
