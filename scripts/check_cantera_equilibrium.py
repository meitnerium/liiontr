"""Inspect a synthetic Cantera equilibrium state."""

from liiontr.chemistry.cantera import (
    CanteraEquilibriumBackend,
)
from liiontr.chemistry.elements import (
    ElementInventory,
)


def main() -> None:
    """Run and display a synthetic equilibrium calculation."""
    backend = CanteraEquilibriumBackend()

    inventory = ElementInventory(
        moles={
            "C": 1.0e-3,
            "H": 4.0e-3,
            "O": 4.0e-3,
            "N": 2.0e-3,
        }
    )

    state = backend.equilibrate_tv(
        temperature=1500.0,
        volume=1.0e-3,
        element_inventory=inventory,
    )

    uv_state = backend.equilibrate_uv(
        internal_energy=state.internal_energy,
        volume=state.volume,
        element_inventory=inventory,
    )
    
    print()
    print("Cantera TV equilibrium")
    print("----------------------")

    print(
        f"Temperature:       "
        f"{state.temperature:.3f} K"
    )

    print(
        f"Pressure:          "
        f"{state.pressure / 1000.0:.3f} kPa"
    )

    print(
        f"Density:           "
        f"{state.density:.6f} kg/m3"
    )

    print(
        f"Mean molar mass:   "
        f"{state.mean_molar_mass * 1000.0:.6f} g/mol"
    )

    print(
        f"cp:                "
        f"{state.cp_mass:.3f} J/(kg K)"
    )

    print(
        f"cv:                "
        f"{state.cv_mass:.3f} J/(kg K)"
    )

    print(
        f"gamma:             "
        f"{state.heat_capacity_ratio:.6f}"
    )

    print(
        f"Internal energy:   "
        f"{state.internal_energy:.6f} J"
    )

    print()
    print("Dominant species")
    print("----------------")

    species = sorted(
        state.species_moles,
        key=state.mole_fraction,
        reverse=True,
    )

    for species_name in species:
        mole_fraction = state.mole_fraction(
            species_name
        )

        if mole_fraction < 1.0e-6:
            continue

        print(
            f"{species_name:8s} "
            f"{mole_fraction:12.6e} "
            f"{state.moles_of(species_name):12.6e} mol"
        )

    print()
    print("TV -> UV consistency")
    print("--------------------")

    print(f"TV temperature:    {state.temperature:.6f} K")

    print(f"UV temperature:    {uv_state.temperature:.6f} K")

    print(f"Delta T:           {uv_state.temperature - state.temperature:.6e} K")

    print(f"TV pressure:       {state.pressure / 1000.0:.6f} kPa")

    print(f"UV pressure:       {uv_state.pressure / 1000.0:.6f} kPa")

    print(f"Delta P:           {uv_state.pressure - state.pressure:.6e} Pa")

    print(f"TV energy:         {state.internal_energy:.9f} J")

    print(f"UV energy:         {uv_state.internal_energy:.9f} J")

    print(
        f"Delta U:           {uv_state.internal_energy - state.internal_energy:.6e} J"
    )

if __name__ == "__main__":
    main()