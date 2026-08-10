from dataclasses import dataclass

from .properties import ConstantProperty


@dataclass(slots=True)
class Material:
    """
    Homogeneous material definition.
    """

    name: str

    density: ConstantProperty

    heat_capacity: ConstantProperty

    thermal_conductivity: ConstantProperty
