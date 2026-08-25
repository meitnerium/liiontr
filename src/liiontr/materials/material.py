"""Material definitions for battery thermal models."""

from dataclasses import dataclass

from .properties import ConstantProperty


@dataclass(slots=True)
class Material:
    """Represent a homogeneous material.

    Parameters
    ----------
    name : str
        Human-readable material name.
    density : ConstantProperty
        Mass density in kg/m³.
    heat_capacity : ConstantProperty
        Specific heat capacity in J/(kg·K).
    thermal_conductivity : ConstantProperty
        Thermal conductivity in W/(m·K).

    Notes
    -----
    The current implementation uses temperature-independent material
    properties through :class:`ConstantProperty`.
    """

    name: str

    density: ConstantProperty

    heat_capacity: ConstantProperty

    thermal_conductivity: ConstantProperty
