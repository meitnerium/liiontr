"""Base battery cell definitions for LiionTR."""

from __future__ import annotations

from abc import ABC
from dataclasses import dataclass

from liiontr.chemistry.chemistry import Chemistry
from liiontr.geometry.geometry import Geometry
from liiontr.materials.material import Material


@dataclass(slots=True)
class Cell(ABC):
    """Base class for battery cells.

    Parameters
    ----------
    name : str
        Human-readable cell name.
    geometry : Geometry
        Geometric representation of the cell.
    material : Material
        Effective homogeneous material properties used for thermal
        calculations.
    chemistry : Chemistry
        Electrochemical or thermal-runaway chemistry associated with
        the cell.
    """

    name: str
    geometry: Geometry
    material: Material
    chemistry: Chemistry

    @property
    def volume(self) -> float:
        """Return the cell volume in m³."""
        return self.geometry.volume

    @property
    def surface_area(self) -> float:
        """Return the external cell surface area in m²."""
        return self.geometry.surface_area

    @property
    def mass(self) -> float:
        r"""Return the estimated cell mass in kg.

        Notes
        -----
        The mass is computed from the effective material density
        evaluated at 298.15 K:

        .. math::

            m = \rho(298.15\\,\mathrm{K}) V

        where ``V`` is the cell volume.
        """
        return self.material.density(298.15) * self.volume

    @property
    def thermal_capacity(self) -> float:
        r"""Return the lumped cell thermal capacity in J/K.

        Notes
        -----
        The thermal capacity is computed from the cell mass and the
        effective specific heat capacity evaluated at 298.15 K:

        .. math::

            C_{th} =
            m c_p(298.15\,\mathrm{K})
        """
        return self.mass * self.material.heat_capacity(298.15)
