from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass


class ProgressModel(ABC):
    """
    Base class for reaction progress models.
    """

    @abstractmethod
    def factor(
        self,
        conversion: float,
    ) -> float:
        """
        Return the conversion-dependent reaction factor.
        """

        raise NotImplementedError


@dataclass(slots=True)
class PowerLawProgress(ProgressModel):
    """
    Power-law reaction progress model.

    f(alpha) = (1 - alpha)^n
    """

    order: float = 1.0

    def __post_init__(self) -> None:
        if self.order <= 0.0:
            raise ValueError("Progress model order must be greater than zero.")

    def factor(
        self,
        conversion: float,
    ) -> float:
        """
        Return the reaction progress factor.
        """

        if conversion >= 1.0:
            return 0.0

        if conversion < 0.0:
            conversion = 0.0

        return (1.0 - conversion) ** self.order


@dataclass(slots=True)
class AutocatalyticProgress(ProgressModel):
    """
    Autocatalytic reaction progress model.

    f(alpha) = alpha^m * (1 - alpha)^n
    """

    autocatalytic_order: float = 1.0
    remaining_order: float = 1.0

    def __post_init__(self) -> None:
        if self.autocatalytic_order <= 0.0 or self.remaining_order <= 0.0:
            raise ValueError("Autocatalytic progress orders must be greater than zero.")

    def factor(
        self,
        conversion: float,
    ) -> float:
        """
        Return the autocatalytic reaction progress factor.
        """

        if conversion <= 0.0:
            return 0.0

        if conversion >= 1.0:
            return 0.0

        return (
            conversion**self.autocatalytic_order
            * (1.0 - conversion) ** self.remaining_order
        )
