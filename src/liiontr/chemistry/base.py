from __future__ import annotations

from abc import ABC, abstractmethod


class ChemistryBackend(ABC):
    """
    Abstract interface for battery chemistry models.

    A chemistry backend provides:
    - heat generation
    - reaction information
    - optional gas generation
    """

    @abstractmethod
    def heat_generation(
        self,
        temperature: float,
    ) -> float:
        """
        Return volumetric or mass specific
        heat generation rate.
        """

        raise NotImplementedError
