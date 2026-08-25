"""Gas species definitions and gas-mixture inventories."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(slots=True, frozen=True)
class GasSpecies:
    """
    Definition of a gas species.

    Parameters
    ----------
    name:
        Species name.

    molar_mass:
        Molar mass in kg/mol.
    """

    name: str
    molar_mass: float

    def __post_init__(self) -> None:
        """Validate the gas species definition."""
        if self.molar_mass <= 0.0:
            raise ValueError("Molar mass must be greater than zero.")


@dataclass(slots=True)
class GasInventory:
    """Gas composition represented by species mole amounts."""

    species: list[GasSpecies]
    moles: dict[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate gas species and mole amounts."""
        names = [species.name for species in self.species]

        if len(names) != len(set(names)):
            raise ValueError("Gas species names must be unique.")

        known_species = set(names)

        for name, amount in self.moles.items():
            if name not in known_species:
                raise KeyError(f"Unknown gas species: {name}")

            if amount < 0.0:
                raise ValueError("Gas moles must not be negative.")

    def moles_of(
        self,
        name: str,
    ) -> float:
        """Return the mole amount of one species."""
        if name not in {species.name for species in self.species}:
            raise KeyError(f"Unknown gas species: {name}")

        return self.moles.get(
            name,
            0.0,
        )

    @property
    def total_moles(self) -> float:
        """Return the total gas amount in mol."""
        return sum(self.moles_of(species.name) for species in self.species)

    @property
    def total_mass(self) -> float:
        """Return total gas mass in kg."""
        return sum(
            self.moles_of(species.name) * species.molar_mass for species in self.species
        )

    def mole_fraction(
        self,
        name: str,
    ) -> float:
        """Return the mole fraction of one gas species."""
        amount = self.moles_of(name)

        total = self.total_moles

        if total == 0.0:
            return 0.0

        return amount / total

    @property
    def mean_molar_mass(self) -> float:
        """Return the mole-weighted mean molar mass in kg/mol."""
        total = self.total_moles

        if total == 0.0:
            raise ValueError(
                "Cannot compute mean molar mass of an empty gas inventory."
            )

        return self.total_mass / total
