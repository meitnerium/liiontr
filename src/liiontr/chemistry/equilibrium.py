"""Thermochemical gas-equilibrium interfaces and states."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Protocol, runtime_checkable

from liiontr.chemistry.elements import ElementInventory


@dataclass(slots=True)
class GasEquilibriumState:
    """Thermodynamic state of an equilibrated gas mixture."""

    temperature: float
    pressure: float
    volume: float

    species_moles: dict[str, float]

    mean_molar_mass: float
    density: float

    cp_mass: float
    cv_mass: float

    internal_energy: float

    def __post_init__(self) -> None:
        scalar_values = {
            "temperature": self.temperature,
            "pressure": self.pressure,
            "volume": self.volume,
            "mean_molar_mass": self.mean_molar_mass,
            "density": self.density,
            "cp_mass": self.cp_mass,
            "cv_mass": self.cv_mass,
            "internal_energy": self.internal_energy,
        }

        for name, value in scalar_values.items():
            if not isfinite(value):
                raise ValueError(
                    f"{name} must be finite."
                )

        if self.temperature <= 0.0:
            raise ValueError(
                "Temperature must be greater than zero."
            )

        if self.pressure < 0.0:
            raise ValueError(
                "Pressure must not be negative."
            )

        if self.volume <= 0.0:
            raise ValueError(
                "Volume must be greater than zero."
            )

        if self.mean_molar_mass <= 0.0:
            raise ValueError(
                "Mean molar mass must be greater than zero."
            )

        if self.density <= 0.0:
            raise ValueError(
                "Density must be greater than zero."
            )

        if self.cp_mass <= 0.0:
            raise ValueError(
                "cp_mass must be greater than zero."
            )

        if self.cv_mass <= 0.0:
            raise ValueError(
                "cv_mass must be greater than zero."
            )

        normalized_species: dict[str, float] = {}

        for species_name, amount in self.species_moles.items():
            name = species_name.strip()

            if not name:
                raise ValueError(
                    "Species name must not be empty."
                )

            value = float(amount)

            if not isfinite(value):
                raise ValueError(
                    "Species mole amount must be finite."
                )

            if value < 0.0:
                raise ValueError(
                    "Species mole amount must not be negative."
                )

            normalized_species[name] = value

        if sum(normalized_species.values()) <= 0.0:
            raise ValueError(
                "Gas state must contain a positive amount "
                "of gas."
            )

        self.species_moles = normalized_species

    @property
    def total_moles(self) -> float:
        """Return total gas amount in mol."""
        return sum(
            self.species_moles.values()
        )

    @property
    def total_mass(self) -> float:
        """Return total gas mass in kg."""
        return (
            self.density
            * self.volume
        )

    @property
    def heat_capacity_ratio(self) -> float:
        """Return the heat-capacity ratio cp/cv."""
        return (
            self.cp_mass
            / self.cv_mass
        )

    def moles_of(
        self,
        species_name: str,
    ) -> float:
        """Return the mole amount of one species."""
        return self.species_moles.get(
            species_name,
            0.0,
        )

    def mole_fraction(
        self,
        species_name: str,
    ) -> float:
        """Return the mole fraction of one species."""
        return (
            self.moles_of(species_name)
            / self.total_moles
        )


@runtime_checkable
class GasEquilibriumBackend(Protocol):
    """Protocol for thermochemical gas-equilibrium backends."""

    @property
    def element_names(self) -> tuple[str, ...]:
        """Return elements supported by the backend."""
        ...

    @property
    def species_names(self) -> tuple[str, ...]:
        """Return species supported by the backend."""
        ...

    def species_molar_mass(
        self,
        species_name: str,
    ) -> float:
        """Return a species molar mass in kg/mol."""
        ...

    def element_flow_rates(
        self,
        species_molar_flow_rates: dict[str, float],
    ) -> dict[str, float]:
        """Convert species molar flow rates to elemental rates."""
        ...

    def equilibrate_tv(
        self,
        temperature: float,
        volume: float,
        element_inventory: ElementInventory,
    ) -> GasEquilibriumState:
        """Equilibrate at constant temperature and volume."""
        ...

    def equilibrate_uv(
        self,
        internal_energy: float,
        volume: float,
        element_inventory: ElementInventory,
    ) -> GasEquilibriumState:
        """Equilibrate at constant internal energy and volume."""
        ...