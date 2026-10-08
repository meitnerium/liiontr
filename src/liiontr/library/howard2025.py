"""Literature gas-generation data from Howard et al. (2025)."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

from liiontr.chemistry.elements import (
    ElementInventory,
)


@dataclass(frozen=True, slots=True)
class VentGasDataset:
    """Store a literature vent-gas dataset."""

    name: str
    cell_format: str
    chemistry: str
    capacity_ah: float
    cell_mass: float
    state_of_charge: float
    gas_volume: float
    species_fractions: dict[str, float]
    reference: str

    def species_moles(
        self,
        temperature: float,
        pressure: float,
    ) -> dict[str, float]:
        """Convert reported gas volume and composition to species moles."""
        if not isfinite(temperature) or temperature <= 0.0:
            raise ValueError("Temperature must be finite and greater than zero.")

        if not isfinite(pressure) or pressure <= 0.0:
            raise ValueError("Pressure must be finite and greater than zero.")

        fraction_sum = sum(self.species_fractions.values())

        if fraction_sum <= 0.0:
            raise ValueError("Species fractions must have a positive sum.")

        total_moles = pressure * self.gas_volume / (IDEAL_GAS_CONSTANT * temperature)

        return {
            species_name: (total_moles * fraction / fraction_sum)
            for species_name, fraction in self.species_fractions.items()
        }

    def element_yields_per_cell_mass(
        self,
        temperature: float,
        pressure: float,
    ) -> dict[str, float]:
        """Return elemental gas yields in mol per kg of cell."""
        species_moles = self.species_moles(
            temperature=temperature,
            pressure=pressure,
        )

        inventory = ElementInventory.from_species_formulas(species_moles)

        return {
            element_name: (amount / self.cell_mass)
            for element_name, amount in inventory.moles.items()
        }


IDEAL_GAS_CONSTANT = 8.31446261815324


def howard2025_nmc811_21700_100soc_inert() -> VentGasDataset:
    """Return the Howard 2025 NMC811 21700 inert-atmosphere dataset."""
    return VentGasDataset(
        name="Howard 2025 NMC811 21700 100% SoC inert",
        cell_format="21700",
        chemistry="NMC811",
        capacity_ah=5.0,
        cell_mass=0.068,
        gas_volume=8.3e-3,
        state_of_charge=1.0,
        species_fractions={
            "H2": 0.220,
            "CO2": 0.278,
            "CO": 0.378,
            "C2H6": 0.009,
            "C2H4": 0.039,
            "C3H8": 0.023,
            "C3H6": 0.005,
            "CH4": 0.049,
        },
        reference=(
            "Howard et al. (2025), Batteries 11, 320, "
            "doi:10.3390/batteries11090320"
        ),
    )
