from __future__ import annotations

from dataclasses import dataclass, field

from .variable import Variable


@dataclass(slots=True)
class State:
    """
    Complete state of the physical system.
    """

    variables: dict[str, Variable] = field(default_factory=dict)

    def add(self, variable: Variable) -> None:
        """
        Add a variable to the state.
        """
        self.variables[variable.name] = variable

    def get(self, name: str) -> Variable:
        """
        Retrieve a variable.
        """
        return self.variables[name]

    def __contains__(self, name: str) -> bool:
        return name in self.variables

    def __len__(self) -> int:
        return len(self.variables)
