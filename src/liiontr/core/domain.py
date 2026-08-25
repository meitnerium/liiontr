"""Computational domain definitions used by LiionTR models."""

from __future__ import annotations

from abc import ABC


class Domain(ABC):
    """
    Base class for physical domains.

    A domain describes where physics takes place.
    """

    def __init__(
        self,
        name: str,
        dimension: int,
    ):
        self.name = name
        self.dimension = dimension

    def describe(self) -> str:
        """
        Return a human-readable description.
        """

        return f"{self.name} ({self.dimension}D)"
