from __future__ import annotations

from abc import ABC
from dataclasses import dataclass


@dataclass(slots=True)
class Chemistry(ABC):
    """
    Base class describing a Li-ion chemistry.
    """

    name: str

    nominal_voltage: float

    specific_capacity: float

    density: float
