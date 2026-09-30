"""Models for delayed release of generated elements."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Protocol

from liiontr.chemistry.elements import (
    ElementInventory,
)


class ElementReleaseModel(Protocol):
    """Define the interface for elemental release models."""

    def release_rates(
        self,
        stored_inventory: ElementInventory,
    ) -> dict[str, float]:
        """Return elemental release rates in mol/s."""
        ...

@dataclass(frozen=True, slots=True)
class FirstOrderElementReleaseModel:
    """Release stored elemental inventory with first-order kinetics."""

    time_constant: float

    def __post_init__(self) -> None:
        """Validate the release model."""
        if (
            not isfinite(self.time_constant)
            or self.time_constant <= 0.0
        ):
            raise ValueError(
                "Release time constant must be finite "
                "and greater than zero."
            )

    def release_rates(
        self,
        stored_inventory: ElementInventory,
    ) -> dict[str, float]:
        """Return elemental release rates in mol/s."""
        return {
            element_name: (
                amount / self.time_constant
            )
            for element_name, amount
            in stored_inventory.moles.items()
        }