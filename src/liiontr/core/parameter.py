"""Physical and numerical parameter representations."""

from dataclasses import dataclass


@dataclass(slots=True, frozen=True)
class Parameter:
    """
    Represents a single model parameter.
    """

    name: str
    value: float
    unit: str
    description: str = ""
