"""Battery chemistry definitions used by LiionTR."""

from __future__ import annotations

from abc import ABC
from dataclasses import dataclass


@dataclass(slots=True)
class Chemistry(ABC):
    """Describe a lithium-ion battery chemistry."""

    name: str

    nominal_voltage: float

    specific_capacity: float

    density: float
