"""Cylindrical lithium-ion battery cell models."""

from dataclasses import dataclass

from .cell import Cell


@dataclass(slots=True)
class CylindricalCell(Cell):
    """Represent a cylindrical lithium-ion battery cell."""
