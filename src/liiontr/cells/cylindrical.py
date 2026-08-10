from dataclasses import dataclass

from .cell import Cell


@dataclass(slots=True)
class CylindricalCell(Cell):
    """
    Cylindrical Li-ion cell.
    """
