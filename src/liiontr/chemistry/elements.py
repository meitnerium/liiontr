"""Element inventories and reaction element yields."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from math import isfinite


def element_counts_from_formula(
    formula: str,
) -> dict[str, int]:
    """Return elemental atom counts from a chemical formula."""
    if not formula:
        raise ValueError(
            "Chemical formula must not be empty."
        )

    matches = list(
        re.finditer(
            r"([A-Z][a-z]?)(\d*)",
            formula,
        )
    )

    reconstructed = "".join(
        match.group(0)
        for match in matches
    )

    if reconstructed != formula:
        raise ValueError(
            f"Unsupported chemical formula: {formula!r}."
        )

    counts: dict[str, int] = {}

    for match in matches:
        element_name = match.group(1)
        count_text = match.group(2)

        count = (
            int(count_text)
            if count_text
            else 1
        )

        counts[element_name] = (
            counts.get(
                element_name,
                0,
            )
            + count
        )

    return counts

@dataclass(slots=True)
class ElementInventory:
    """Elemental inventory expressed as moles of atoms."""

    moles: dict[str, float] = field(
        default_factory=dict
    )

    def __post_init__(self) -> None:
        normalized: dict[str, float] = {}

        for element, amount in self.moles.items():
            name = element.strip()

            if not name:
                raise ValueError(
                    "Element name must not be empty."
                )

            value = float(amount)

            if not isfinite(value):
                raise ValueError(
                    "Element mole amount must be finite."
                )

            if value < 0.0:
                raise ValueError(
                    "Element mole amount must not be negative."
                )

            if name in normalized:
                raise ValueError(
                    f"Duplicate element name: {name}"
                )

            normalized[name] = value

        self.moles = normalized

    @classmethod
    def from_species_formulas(
        cls,
        species_moles: dict[str, float],
    ) -> ElementInventory:
        """Build an elemental inventory from molecular formulas."""
        element_moles: dict[str, float] = {}

        for formula, amount in species_moles.items():
            if not isfinite(amount) or amount < 0.0:
                raise ValueError("Species mole amounts must be finite and nonnegative.")

            counts = element_counts_from_formula(formula)

            for element_name, atom_count in counts.items():
                element_moles[element_name] = (
                    element_moles.get(
                        element_name,
                        0.0,
                    )
                    + amount * atom_count
                )

        return cls(moles=element_moles)

    @property
    def element_names(self) -> list[str]:
        """Return element names in deterministic order."""
        return list(self.moles)

    def moles_of(
        self,
        element: str,
    ) -> float:
        """Return the mole amount of one element."""
        return self.moles.get(
            element,
            0.0,
        )

    @property
    def total_atom_moles(self) -> float:
        """Return the total amount of elemental atoms in mol."""
        return sum(
            self.moles.values()
        )


@dataclass(slots=True)
class ReactionElementYield:
    """Element yields associated with one physical reaction."""

    reaction_name: str
    element_yields: dict[str, float]

    def __post_init__(self) -> None:
        self.reaction_name = (
            self.reaction_name.strip()
        )

        if not self.reaction_name:
            raise ValueError(
                "Reaction name must not be empty."
            )

        normalized: dict[str, float] = {}

        for element, element_yield in (
            self.element_yields.items()
        ):
            name = element.strip()

            if not name:
                raise ValueError(
                    "Element name must not be empty."
                )

            value = float(
                element_yield
            )

            if not isfinite(value):
                raise ValueError(
                    "Element yield must be finite."
                )

            if value < 0.0:
                raise ValueError(
                    "Element yield must not be negative."
                )

            if name in normalized:
                raise ValueError(
                    f"Duplicate element name: {name}"
                )

            normalized[name] = value

        self.element_yields = normalized

    @property
    def element_names(self) -> list[str]:
        """Return elements produced by this reaction."""
        return list(
            self.element_yields
        )

    def yield_of(
        self,
        element: str,
    ) -> float:
        """Return the yield of one element."""
        return self.element_yields.get(
            element,
            0.0,
        )