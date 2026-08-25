"""Abstract interfaces for LiionTR chemistry backends."""

from __future__ import annotations

from abc import ABC, abstractmethod


class ChemistryBackend(ABC):
    """Backend used to compute thermal chemistry."""

    @abstractmethod
    def heat_generation(
        self,
        temperature: float,
        conversions: list[float] | None = None,
    ) -> float:
        """Return the total cell heat generation rate in W."""
        raise NotImplementedError
