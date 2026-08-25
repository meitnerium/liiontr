"""Material models and thermophysical properties used by LiionTR."""

from .material import Material
from .properties import ConstantProperty

__all__ = [
    "Material",
    "ConstantProperty",
]
