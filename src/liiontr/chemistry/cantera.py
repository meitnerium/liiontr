"""Cantera thermochemical equilibrium backend."""

from __future__ import annotations

from dataclasses import dataclass, field
from math import isfinite

import cantera as ct

from liiontr.chemistry.elements import ElementInventory
from liiontr.chemistry.equilibrium import (
    GasEquilibriumState,
)

DEFAULT_GRI30_CARRIER_SPECIES = {
    "C": "C",
    "H": "H2",
    "O": "O2",
    "N": "N2",
}


@dataclass(slots=True)
class CanteraEquilibriumBackend:
    """Prepare elemental inventories for Cantera equilibrium."""

    mechanism: str = "gri30.yaml"
    phase_name: str | None = None
    uv_seed_temperature: float = 1000.0
    carrier_species: dict[str, str] = field(
        default_factory=lambda: dict(
            DEFAULT_GRI30_CARRIER_SPECIES
        )
    )

    _gas: ct.Solution = field(
        init=False,
        repr=False,
    )

    def __post_init__(self) -> None:
        """Load the Cantera phase and validate carrier species."""
        if self.phase_name is None:
            self._gas = ct.Solution(
                self.mechanism
            )
        else:
            self._gas = ct.Solution(
                self.mechanism,
                self.phase_name,
            )

        self.carrier_species = {
            element.strip(): species.strip()
            for element, species
            in self.carrier_species.items()
        }

        if self.uv_seed_temperature <= 0.0:
            raise ValueError(
                "UV seed temperature must be greater than zero."
            )

        self._validate_carrier_species()

    @property
    def element_names(self) -> tuple[str, ...]:
        """Return elements supported by the Cantera phase."""
        return tuple(
            self._gas.element_names
        )

    @property
    def species_names(self) -> tuple[str, ...]:
        """Return species supported by the Cantera phase."""
        return tuple(
            self._gas.species_names
        )

    def species_molar_mass(
        self,
        species_name: str,
    ) -> float:
        """Return a species molar mass in kg/mol."""
        if species_name not in self._gas.species_names:
            raise ValueError(
                f"Species {species_name!r} is not present "
                f"in mechanism {self.mechanism!r}."
            )

        species_index = self._gas.species_index(
            species_name
        )

        return (
            float(
                self._gas.molecular_weights[
                    species_index
                ]
            )
            / 1000.0
        )

    def _validate_carrier_species(
        self,
    ) -> None:
        """Validate the elemental carrier-species mapping."""
        mechanism_elements = set(
            self._gas.element_names
        )

        mechanism_species = set(
            self._gas.species_names
        )

        for (
            element,
            species_name,
        ) in self.carrier_species.items():
            if not element:
                raise ValueError(
                    "Carrier element name must not be empty."
                )

            if not species_name:
                raise ValueError(
                    "Carrier species name must not be empty."
                )

            if element not in mechanism_elements:
                raise ValueError(
                    f"Element {element!r} is not present "
                    "in the Cantera mechanism."
                )

            if species_name not in mechanism_species:
                raise ValueError(
                    f"Carrier species {species_name!r} "
                    "is not present in the Cantera mechanism."
                )

            target_atoms = self._gas.n_atoms(
                species_name,
                element,
            )

            if target_atoms <= 0.0:
                raise ValueError(
                    f"Carrier species {species_name!r} "
                    f"does not contain element {element!r}."
                )

            other_elements = [
                other_element
                for other_element
                in self._gas.element_names
                if (
                    other_element != element
                    and self._gas.n_atoms(
                        species_name,
                        other_element,
                    )
                    > 0.0
                )
            ]

            if other_elements:
                names = ", ".join(
                    other_elements
                )

                raise ValueError(
                    f"Carrier species {species_name!r} "
                    f"for element {element!r} also contains "
                    f"other elements: {names}."
                )

    def validate_inventory(
        self,
        inventory: ElementInventory,
    ) -> None:
        """Validate an element inventory against the mechanism."""
        if inventory.total_atom_moles <= 0.0:
            raise ValueError(
                "Element inventory must contain a positive "
                "amount of atoms."
            )

        mechanism_elements = set(
            self._gas.element_names
        )

        for (
            element,
            amount,
        ) in inventory.moles.items():
            if amount == 0.0:
                continue

            if element not in mechanism_elements:
                raise ValueError(
                    f"Element {element!r} is not supported "
                    f"by mechanism {self.mechanism!r}."
                )

            if element not in self.carrier_species:
                raise ValueError(
                    f"No carrier species is configured "
                    f"for element {element!r}."
                )

    def carrier_moles(
        self,
        inventory: ElementInventory,
    ) -> dict[str, float]:
        """Convert an elemental inventory to carrier-species moles."""
        self.validate_inventory(
            inventory
        )

        species_moles: dict[str, float] = {}

        for (
            element,
            atom_moles,
        ) in inventory.moles.items():
            if atom_moles == 0.0:
                continue

            species_name = (
                self.carrier_species[
                    element
                ]
            )

            atoms_per_species = (
                self._gas.n_atoms(
                    species_name,
                    element,
                )
            )

            amount = (
                atom_moles
                / atoms_per_species
            )

            species_moles[
                species_name
            ] = (
                species_moles.get(
                    species_name,
                    0.0,
                )
                + amount
            )

        return species_moles

    def inventory_mass(
        self,
        inventory: ElementInventory,
    ) -> float:
        """Return the mass represented by an element inventory."""
        self.validate_inventory(
            inventory
        )

        mass = 0.0

        for (
            element,
            atom_moles,
        ) in inventory.moles.items():
            if atom_moles == 0.0:
                continue

            atomic_weight = (
                self._gas.atomic_weight(
                    element
                )
            )

            mass += (
                atom_moles
                * atomic_weight
                / 1000.0
            )

        return mass

    def _state_from_current_gas(
        self,
        volume: float,
        total_mass: float,
    ) -> GasEquilibriumState:
        """Build a LiionTR state from the current Cantera state."""
        mean_molar_mass = (
            self._gas.mean_molecular_weight
            / 1000.0
        )

        total_moles = (
            total_mass
            / mean_molar_mass
        )

        species_moles = {
            species_name: (
                float(mole_fraction)
                * total_moles
            )
            for species_name, mole_fraction
            in zip(
                self._gas.species_names,
                self._gas.X,
            )
            if mole_fraction > 0.0
        }

        return GasEquilibriumState(
            temperature=float(
                self._gas.T
            ),
            pressure=float(
                self._gas.P
            ),
            volume=volume,
            species_moles=species_moles,
            mean_molar_mass=mean_molar_mass,
            density=float(
                self._gas.density_mass
            ),
            cp_mass=float(
                self._gas.cp_mass
            ),
            cv_mass=float(
                self._gas.cv_mass
            ),
            internal_energy=(
                float(
                    self._gas.int_energy_mass
                )
                * total_mass
            ),
        )

    def equilibrate_tv(
        self,
        temperature: float,
        volume: float,
        element_inventory: ElementInventory,
    ) -> GasEquilibriumState:
        """Equilibrate an elemental inventory at constant T and V."""
        if temperature <= 0.0:
            raise ValueError(
                "Temperature must be greater than zero."
            )

        if volume <= 0.0:
            raise ValueError(
                "Volume must be greater than zero."
            )

        carrier_moles = self.carrier_moles(
            element_inventory
        )

        total_mass = self.inventory_mass(
            element_inventory
        )

        density = (
            total_mass
            / volume
        )

        self._gas.TDX = (
            temperature,
            density,
            carrier_moles,
        )

        self._gas.equilibrate(
            "TV"
        )

        return self._state_from_current_gas(
            volume=volume,
            total_mass=total_mass,
        )

    def equilibrate_uv(
        self,
        internal_energy: float,
        volume: float,
        element_inventory: ElementInventory,
    ) -> GasEquilibriumState:
        """Equilibrate an elemental inventory at constant U and V."""
        if not isfinite(internal_energy):
            raise ValueError(
                "Internal energy must be finite."
            )

        if volume <= 0.0:
            raise ValueError(
                "Volume must be greater than zero."
            )

        carrier_moles = self.carrier_moles(
            element_inventory
        )

        total_mass = self.inventory_mass(
            element_inventory
        )

        density = (
            total_mass
            / volume
        )

        self._gas.TDX = (
            self.uv_seed_temperature,
            density,
            carrier_moles,
        )

        # First obtain a chemically reasonable composition
        # with the correct elemental inventory.
        self._gas.equilibrate(
            "TV"
        )

        specific_internal_energy = (
            internal_energy
            / total_mass
        )

        specific_volume = (
            volume
            / total_mass
        )

        self._gas.UV = (
            specific_internal_energy,
            specific_volume,
        )

        self._gas.equilibrate(
            "UV"
        )

        return self._state_from_current_gas(
            volume=volume,
            total_mass=total_mass,
        )

    def element_inventory_from_species_moles(
        self,
        species_moles: dict[str, float],
    ) -> ElementInventory:
        """Convert species mole amounts to an elemental inventory."""
        mechanism_species = set(
            self._gas.species_names
        )

        element_moles = {
            element_name: 0.0
            for element_name in self._gas.element_names
        }

        for species_name, amount in species_moles.items():
            if species_name not in mechanism_species:
                raise ValueError(
                    f"Species {species_name!r} is not present "
                    f"in mechanism {self.mechanism!r}."
                )

            if amount < 0.0:
                raise ValueError(
                    "Species mole amount must not be negative."
                )

            for element_name in self._gas.element_names:
                atom_count = self._gas.n_atoms(
                    species_name,
                    element_name,
                )

                element_moles[element_name] += (
                    amount
                    * atom_count
                )

        return ElementInventory(
            moles={
                element_name: amount
                for element_name, amount
                in element_moles.items()
                if amount > 0.0
            }
        )

    def element_inventory_from_state(
        self,
        state: GasEquilibriumState,
    ) -> ElementInventory:
        """Return the elemental inventory represented by a gas state."""
        return self.element_inventory_from_species_moles(
            state.species_moles
        )

    def element_flow_rates(
        self,
        species_molar_flow_rates: dict[str, float],
    ) -> dict[str, float]:
        """Convert species molar flow rates to element flow rates."""
        mechanism_species = set(
            self._gas.species_names
        )

        element_rates = {
            element_name: 0.0
            for element_name in self._gas.element_names
        }

        for (
            species_name,
            species_rate,
        ) in species_molar_flow_rates.items():
            if species_name not in mechanism_species:
                raise ValueError(
                    f"Species {species_name!r} is not present "
                    f"in mechanism {self.mechanism!r}."
                )

            if species_rate < 0.0:
                raise ValueError(
                    "Species molar flow rate must not be negative."
                )

            for element_name in self._gas.element_names:
                atom_count = self._gas.n_atoms(
                    species_name,
                    element_name,
                )

                element_rates[element_name] += (
                    species_rate
                    * atom_count
                )

        return {
            element_name: rate
            for element_name, rate
            in element_rates.items()
            if rate > 0.0
        }