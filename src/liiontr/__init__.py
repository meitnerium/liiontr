"""LiionTR.

Research framework for lithium-ion battery thermal runaway.
"""

from .cells import CylindricalCell
from .chemistry import NMC811
from .core import (
    Problem,
    Results,
    Simulation,
    State,
    Variable,
)
from .geometry import CylindricalGeometry
from .materials import ConstantProperty, Material

__version__ = "0.0.1"


__all__ = [
    "NMC811",
    "ConstantProperty",
    "CylindricalCell",
    "CylindricalGeometry",
    "Material",
    "Problem",
    "Results",
    "Simulation",
    "State",
    "Variable",
]
