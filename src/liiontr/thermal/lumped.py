from __future__ import annotations

from dataclasses import dataclass

from liiontr.cells.cell import Cell

from .model import ThermalModel


@dataclass(slots=True)
class LumpedThermalModel(ThermalModel):
    cell: Cell

    convection_coefficient: float = 10.0
    ambient_temperature: float = 298.15

    def temperature_derivative(
        self,
        temperature: float,
        heat_generation: float,
    ) -> float:
        heat_loss = (
            self.convection_coefficient
            * self.cell.surface_area
            * (temperature - self.ambient_temperature)
        )

        return (heat_generation - heat_loss) / self.cell.thermal_capacity
