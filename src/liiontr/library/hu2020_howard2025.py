"""Coupling between Hu 2020 kinetics and Howard 2025 gas data."""

from __future__ import annotations

from liiontr.chemistry.energy_scaled_generation import (
    EnergyScaledElementGenerationModel,
)
from liiontr.library.howard2025 import (
    howard2025_nmc811_21700_100soc_inert,
)

# Specific reaction energy of the uncoupled Hu 2020
# baseline, obtained from the thermal energy balance:
#
#     Q_reaction = Delta U_sensible + Q_loss
#
# The earlier value of 386562 J/kg was obtained by
# trapezoidal post-processing of the sharply peaked
# reaction heat rate and slightly overestimated the
# baseline energy.

HU2020_REFERENCE_REACTION_ENERGY_PER_CELL_MASS = (
    386405.782042
)


def hu2020_howard2025_element_generation_model(
    gas_reference_temperature: float = 298.15,
    gas_reference_pressure: float = 101325.0,
) -> EnergyScaledElementGenerationModel:
    """Build the Hu-Howard elemental gas-generation coupling.

    The total elemental gas yield is taken from Howard et al.
    (2025) for an NMC811 21700 cell at 100% state of charge
    under inert conditions.

    Gas-generation timing is assumed to scale with the total
    reaction heat-release rate of the Hu et al. (2020) thermal
    runaway model.

    The reference reaction energy is obtained from the LiionTR
    Hu 2020 baseline simulation and is therefore a coupling
    calibration parameter, not a value reported by Howard et al.
    """
    dataset = (
        howard2025_nmc811_21700_100soc_inert()
    )

    element_yields = (
        dataset.element_yields_per_cell_mass(
            temperature=gas_reference_temperature,
            pressure=gas_reference_pressure,
        )
    )

    return EnergyScaledElementGenerationModel(
        element_yields_per_cell_mass=element_yields,
        reference_reaction_energy_per_cell_mass=(
            HU2020_REFERENCE_REACTION_ENERGY_PER_CELL_MASS
        ),
    )