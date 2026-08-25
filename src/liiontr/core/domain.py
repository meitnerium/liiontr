"""Computational domain definitions used by LiionTR models."""

from __future__ import annotations

from abc import ABC


class Domain(ABC):
    """Base class for physical domains.

    A domain describes the dimensional space over which a physical
    model is defined.

    Parameters
    ----------
    name : str
        Human-readable name of the domain.
    dimension : int
        Spatial dimension of the domain.
    """

    def __init__(
        self,
        name: str,
        dimension: int,
    ) -> None:
        """Initialize a physical domain.

        Parameters
        ----------
        name : str
            Human-readable name of the domain.
        dimension : int
            Spatial dimension of the domain.
        """
        self.name = name
        self.dimension = dimension

    def describe(self) -> str:
        """Return a human-readable description of the domain.

        Returns
        -------
        str
            Domain name followed by its spatial dimension.
        """
        return f"{self.name} ({self.dimension}D)"
