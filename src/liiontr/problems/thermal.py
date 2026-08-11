from dataclasses import dataclass

from liiontr.cells.cell import Cell
from liiontr.chemistry import ChemistryBackend
from liiontr.core.problem import Problem


@dataclass(slots=True)
class ThermalProblem(Problem):
    """
    Lumped thermal problem definition.
    """

    cell: Cell

    chemistry_backend: ChemistryBackend | None = None

    initial_temperature: float = 298.15

    ambient_temperature: float = 298.15

    convection_coefficient: float = 10.0

    duration: float = 3600.0
