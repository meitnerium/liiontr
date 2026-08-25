from __future__ import annotations

import numpy as np

from liiontr.chemistry import ReactionNetworkBackend
from liiontr.gases import (
    CompressibleVentFlowModel,
    GasGenerationModel,
    GasInventory,
    GasSpecies,
    IdealGasPressureModel,
    MixtureVentFlowModel,
    ReactionGasYield,
)
from liiontr.kinetics import Arrhenius
from liiontr.library import cell_21700_generic
from liiontr.problems import ThermalProblem
from liiontr.reactions import (
    Reaction,
    ReactionNetwork,
)
from liiontr.solver import ScipySolver


def main() -> None:
    cell = cell_21700_generic()

    reaction = Reaction(
        name="Gas-producing reaction",
        kinetics=Arrhenius(
            activation_energy=1.0,
            pre_exponential_factor=2.0,
        ),
        enthalpy=10000.0,
        mass_fraction=0.10,
    )

    network = ReactionNetwork(
        reactions=[
            reaction,
        ]
    )

    backend = ReactionNetworkBackend(
        reaction_network=network,
        cell=cell,
    )

    gas_generation_model = GasGenerationModel(
        reaction_network=network,
        gas_yields=[
            ReactionGasYield(
                reaction_name="Gas-producing reaction",
                species_yields={
                    "CO2": 2.0,
                },
            )
        ],
    )

    pressure_model = IdealGasPressureModel(
        free_volume=1.0e-6,
        initial_pressure=101325.0,
        initial_temperature=400.0,
    )

    initial_gas_inventory = GasInventory(
        species=[
            GasSpecies(
                name="N2",
                molar_mass=28.0134e-3,
            ),
            GasSpecies(
                name="CO2",
                molar_mass=44.0095e-3,
            ),
        ],
        moles={
            "N2": pressure_model.initial_moles,
        },
    )

    vent_model = MixtureVentFlowModel(
        flow_model=CompressibleVentFlowModel(
            vent_area=1.0e-5,
            discharge_coefficient=0.8,
            heat_capacity_ratio=1.30,
        ),
        downstream_pressure=101325.0,
    )

    problem = ThermalProblem(
        cell=cell,
        chemistry_backend=backend,
        gas_generation_model=gas_generation_model,
        initial_gas_inventory=initial_gas_inventory,
        pressure_model=pressure_model,
        vent_model=vent_model,
        vent_open_pressure=200000.0,
        initial_temperature=400.0,
        initial_conversions=[
            0.0,
        ],
        ambient_temperature=400.0,
        convection_coefficient=0.0,
        duration=1.0,
    )

    results = ScipySolver().solve(problem)

    time = results.time

    temperature = results.get("temperature")

    conversion = results.get("conversion_0")

    pressure = results.get("pressure")

    nitrogen = results.get("gas_N2")

    carbon_dioxide = results.get("gas_CO2")

    vent_open = results.get("vent_open")

    if time is None:
        raise RuntimeError("Time results are missing.")

    if temperature is None:
        raise RuntimeError("Temperature results are missing.")

    if conversion is None:
        raise RuntimeError("Conversion results are missing.")

    if pressure is None:
        raise RuntimeError("Pressure results are missing.")

    if nitrogen is None:
        raise RuntimeError("N2 results are missing.")

    if carbon_dioxide is None:
        raise RuntimeError("CO2 results are missing.")

    if vent_open is None:
        raise RuntimeError("Vent state results are missing.")

    vent_indices = np.flatnonzero(np.asarray(vent_open) > 0.5)

    print()
    print("Vented thermal runaway simulation")
    print("---------------------------------")

    print(f"Cell mass: {cell.mass * 1000.0:.3f} g")

    print(f"Initial temperature: {temperature[0]:.3f} K")

    print(f"Final temperature: {temperature[-1]:.3f} K")

    print()

    print(f"Initial pressure: {pressure[0] / 1.0e5:.3f} bar")

    print(f"Peak pressure: {max(pressure) / 1.0e5:.3f} bar")

    print(f"Final pressure: {pressure[-1] / 1.0e5:.3f} bar")

    if len(vent_indices) > 0:
        vent_index = int(vent_indices[0])

        print()
        print("Vent opening")

        print(f"  Time: {time[vent_index]:.6f} s")

        print(f"  Temperature: {temperature[vent_index]:.3f} K")

        print(f"  Pressure: {pressure[vent_index] / 1.0e5:.3f} bar")

    print()
    print("Gas inventory")

    print(f"  N2: {nitrogen[0] * 1000.0:.6f} -> {nitrogen[-1] * 1000.0:.6f} mmol")

    print(
        f"  CO2: "
        f"{carbon_dioxide[0] * 1000.0:.6f} -> "
        f"{carbon_dioxide[-1] * 1000.0:.6f} mmol"
    )

    print()

    print(f"Reaction conversion: {conversion[0]:.6f} -> {conversion[-1]:.6f}")

    print(f"Solver points: {len(time)}")


if __name__ == "__main__":
    main()
