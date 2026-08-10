from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class Reaction:
    """
    Chemical reaction contributing to thermal runaway.
    """

    name: str

    activation_energy: float
    pre_exponential_factor: float

    enthalpy: float

    mass_fraction: float = 1.0
