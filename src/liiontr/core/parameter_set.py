"""Collections of named model parameters."""

from __future__ import annotations

from dataclasses import dataclass, field

from .parameter import Parameter


@dataclass(slots=True)
class ParameterSet:
    """Store and retrieve named model parameters.

    Parameters
    ----------
    _parameters : dict[str, Parameter]
        Mapping from parameter names to parameter objects.
    """

    _parameters: dict[str, Parameter] = field(default_factory=dict)

    def add(self, parameter: Parameter) -> None:
        """Add or replace a parameter.

        Parameters
        ----------
        parameter : Parameter
            Parameter to store. Its ``name`` attribute is used as
            the dictionary key.
        """
        self._parameters[parameter.name] = parameter

    def get(self, name: str) -> Parameter:
        """Return a parameter by name.

        Parameters
        ----------
        name : str
            Name of the requested parameter.

        Returns
        -------
        Parameter
            Stored parameter.

        Raises
        ------
        KeyError
            If no parameter with the requested name exists.
        """
        return self._parameters[name]

    def value(self, name: str) -> float:
        """Return the numerical value of a named parameter.

        Parameters
        ----------
        name : str
            Name of the requested parameter.

        Returns
        -------
        float
            Numerical parameter value.

        Raises
        ------
        KeyError
            If no parameter with the requested name exists.
        """
        return self.get(name).value

    def __contains__(self, name: str) -> bool:
        """Return whether a parameter is present.

        Parameters
        ----------
        name : str
            Parameter name.

        Returns
        -------
        bool
            ``True`` if the parameter exists, otherwise ``False``.
        """
        return name in self._parameters

    def __len__(self) -> int:
        """Return the number of stored parameters."""
        return len(self._parameters)
