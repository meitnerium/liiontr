"""Lumped-capacitance thermal model for battery cells."""

from __future__ import annotations

from dataclasses import dataclass

from liiontr.cells.cell import Cell

from .model import ThermalModel


@dataclass(slots=True)
class LumpedThermalModel(ThermalModel):
    r"""Represent a zero-dimensional lumped thermal model.

    The cell is assumed to have a spatially uniform temperature.
    Heat generated internally by electrochemical or thermal-runaway
    reactions is balanced by convective heat loss to the surroundings.

    Parameters
    ----------
    cell : Cell
        Battery cell whose mass, surface area, and thermal capacity are
        used in the energy balance.
    convection_coefficient : float, optional
        Convective heat-transfer coefficient in W/(m²·K).
    ambient_temperature : float, optional
        Ambient temperature in K.

    Notes
    -----
    The governing energy balance is

    .. math::

        C_{\mathrm{th}}
        \frac{dT}{dt}
        =
        \dot{Q}_{\mathrm{gen}}
        -
        h A (T - T_{\infty}),

    where :math:`C_{\mathrm{th}}` is the lumped thermal capacity in J/K,
    :math:`\dot{Q}_{\mathrm{gen}}` is the total internal heat-generation
    rate in W, :math:`h` is the convective heat-transfer coefficient,
    and :math:`A` is the external cell surface area.

    This model neglects internal temperature gradients, radiative heat
    transfer, and conductive heat transfer to neighbouring objects.
    """

    cell: Cell

    convection_coefficient: float = 10.0

    ambient_temperature: float = 298.15

    def temperature_derivative(
        self,
        temperature: float,
        heat_generation: float,
    ) -> float:
        """Return the cell temperature rate from the energy balance.

        Parameters
        ----------
        temperature : float
            Uniform cell temperature in K.
        heat_generation : float
            Total internal cell heat-generation rate in W.

        Returns
        -------
        float
            Cell temperature rate in K/s.
        """
        heat_loss = (
            self.convection_coefficient
            * self.cell.surface_area
            * (temperature - self.ambient_temperature)
        )

        return (heat_generation - heat_loss) / self.cell.thermal_capacity
