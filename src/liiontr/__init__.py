"""
LiionTR

Research framework for lithium-ion battery thermal runaway.
"""

from .cells import CylindricalCell
from .chemistry import NMC811
from .geometry import CylindricalGeometry
from .materials import Material, ConstantProperty

from .core import (
    Problem,
    Results,
    Simulation,
    State,
    Variable,
)

__version__ = "0.0.1"


__all__ = [
    "CylindricalCell",
    "CylindricalGeometry",
    "Material",
    "ConstantProperty",
    "NMC811",
    "Problem",
    "Results",
    "Simulation",
    "State",
    "Variable",
]
