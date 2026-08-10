from __future__ import annotations

from abc import ABC, abstractmethod


class ChemistryBackend(ABC):
    """
    Backend used to compute thermal chemistry.
    """

    @abstractmethod
    def heat_generation(
        self,
        temperature: float,
    ) -> float:
        """
        Heat generation rate.
        """

        raise NotImplementedError
