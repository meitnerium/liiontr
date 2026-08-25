"""Collections of named model parameters."""

from __future__ import annotations

from dataclasses import dataclass, field

from .parameter import Parameter


@dataclass(slots=True)
class ParameterSet:
    _parameters: dict[str, Parameter] = field(default_factory=dict)

    def add(self, parameter: Parameter) -> None:
        self._parameters[parameter.name] = parameter

    def get(self, name: str) -> Parameter:
        return self._parameters[name]

    def value(self, name: str) -> float:
        return self.get(name).value

    def __contains__(self, name: str) -> bool:
        return name in self._parameters

    def __len__(self):
        return len(self._parameters)
