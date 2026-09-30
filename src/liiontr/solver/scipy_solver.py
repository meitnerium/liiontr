"""SciPy-based numerical integration for LiionTR thermal problems."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from liiontr.chemistry import ReactionNetworkBackend
from liiontr.chemistry.element_generation import (
    ElementGenerationModel,
)
from liiontr.chemistry.elements import ElementInventory
from liiontr.chemistry.energy_scaled_generation import (
    EnergyScaledElementGenerationModel,
)
from liiontr.chemistry.equilibrium import (
    GasEquilibriumState,
)
from liiontr.core.results import Results
from liiontr.gases import GasInventory
from liiontr.problems.thermal import ThermalProblem
from liiontr.thermal.coupled_energy import (
    CoupledThermalEnergyModel,
)
from liiontr.thermal.equilibrium_energy import (
    EquilibriumGasEnergyModel,
)
from liiontr.thermal.lumped import LumpedThermalModel
from liiontr.thermal.vent_energy import (
    ElementVentEnergyModel,
)


@dataclass(slots=True)
class ScipySolver:
    """
    SciPy-based ODE solver.

    BDF is used by default because thermal runaway kinetics
    can become stiff at elevated temperatures.
    """

    method: str = "BDF"

    relative_tolerance: float = 1.0e-6

    absolute_tolerance: float = 1.0e-9

    def solve(
        self,
        problem: ThermalProblem,
    ) -> Results:
        """Solve a lumped thermal-runaway problem using SciPy.

        The solver integrates the cell temperature, reaction conversions,
        and gas-species mole amounts. When a vent model is configured,
        integration is split into a closed phase and an irreversibly open
        venting phase.

        Parameters
        ----------
        problem : ThermalProblem
            Thermal-runaway problem defining the cell, chemistry, thermal
            boundary conditions, gas models, pressure limits, and venting
            configuration.

        Returns
        -------
        Results
            Time-dependent simulation results containing temperature,
            reaction conversions, gas amounts, pressure, and vent state
            when applicable.

        Raises
        ------
        ValueError
            If the problem configuration is inconsistent.
        RuntimeError
            If numerical integration fails or a required model is missing.

        Notes
        -----
        The state vector is organized as

        ``[T, alpha_0, ..., alpha_n, n_gas_0, ..., n_gas_m]``,

        where ``T`` is the cell temperature in K, ``alpha`` values are
        dimensionless reaction conversions, and gas amounts are in mol.

        The default BDF integration method is suitable for the stiff
        kinetics commonly encountered during thermal runaway.
        """
        if problem.uses_elemental_thermochemistry:
            return self._solve_elemental(
                problem
            )
        model = LumpedThermalModel(
            cell=problem.cell,
            convection_coefficient=problem.convection_coefficient,
            ambient_temperature=problem.ambient_temperature,
        )

        backend = problem.chemistry_backend

        n_reactions = 0

        if isinstance(
            backend,
            ReactionNetworkBackend,
        ):
            n_reactions = len(backend.reaction_network.reactions)

        if problem.initial_conversions is None:
            initial_conversions = [0.0] * n_reactions

        else:
            initial_conversions = list(problem.initial_conversions)

            if len(initial_conversions) != n_reactions:
                raise ValueError(
                    "Number of initial conversions must match number of reactions."
                )

            if any(
                conversion < 0.0 or conversion > 1.0
                for conversion in initial_conversions
            ):
                raise ValueError("Initial conversions must be between 0 and 1.")

        gas_generation_model = problem.gas_generation_model

        pressure_model = problem.pressure_model

        maximum_pressure = problem.maximum_pressure

        if maximum_pressure is not None and pressure_model is None:
            raise ValueError("Maximum pressure requires a pressure model.")

        initial_gas_inventory = problem.initial_gas_inventory

        vent_model = problem.vent_model
        vent_open_pressure = problem.vent_open_pressure

        if vent_model is not None:
            if vent_open_pressure is None:
                raise ValueError("Vent model requires a vent opening pressure.")

            if pressure_model is None:
                raise ValueError("Vent model requires a pressure model.")

            if initial_gas_inventory is None:
                raise ValueError(
                    "Vent model requires an explicit initial gas inventory."
                )

        if vent_open_pressure is not None and vent_model is None:
            raise ValueError("Vent opening pressure requires a vent model.")

        gas_species_names: list[str] = []

        if initial_gas_inventory is not None:
            gas_species_names.extend(
                species.name for species in initial_gas_inventory.species
            )

        if gas_generation_model is not None:
            if not isinstance(
                backend,
                ReactionNetworkBackend,
            ):
                raise ValueError("Gas generation requires a ReactionNetworkBackend.")

            if gas_generation_model.reaction_network is not backend.reaction_network:
                raise ValueError(
                    "Gas generation model must use "
                    "the same reaction network as "
                    "the chemistry backend."
                )

            for species_name in gas_generation_model.species_names:
                if species_name not in gas_species_names:
                    gas_species_names.append(species_name)

        if vent_model is not None:
            if initial_gas_inventory is None:
                raise RuntimeError("Initial gas inventory is missing.")

            defined_species = {
                species.name for species in initial_gas_inventory.species
            }

            missing_species = [
                species_name
                for species_name in gas_species_names
                if species_name not in defined_species
            ]

            if missing_species:
                names = ", ".join(missing_species)

                raise ValueError(
                    f"Vent model requires gas species definitions for: {names}"
                )

        if initial_gas_inventory is None:
            initial_gas_moles = [0.0 for _ in gas_species_names]

        else:
            initial_gas_moles = [
                initial_gas_inventory.moles_of(species_name)
                for species_name in gas_species_names
            ]

        initial_state = [
            problem.initial_temperature,
            *initial_conversions,
            *initial_gas_moles,
        ]

        reaction_count = len(initial_conversions)

        conversion_start = 1

        conversion_end = conversion_start + reaction_count

        gas_start = conversion_end

        def pressure_from_state(
            state: np.ndarray,
        ) -> float:
            """Return the internal gas pressure for an ODE state in Pa."""
            if pressure_model is None:
                raise RuntimeError("Pressure model is not configured.")

            temperature = float(state[0])

            gas_moles = sum(
                max(
                    float(value),
                    0.0,
                )
                for value in state[gas_start:]
            )

            if initial_gas_inventory is not None:
                return pressure_model.pressure_from_total_moles(
                    temperature=temperature,
                    total_moles=gas_moles,
                )

            return pressure_model.pressure(
                temperature=temperature,
                generated_moles=gas_moles,
            )

        def base_rates(
            state: np.ndarray,
        ) -> tuple[
            float,
            list[float],
            list[float],
        ]:
            """Return thermal, reaction-progress, and gas-generation rates."""
            temperature = float(state[0])

            conversions = [
                min(
                    max(
                        float(value),
                        0.0,
                    ),
                    1.0,
                )
                for value in state[conversion_start:conversion_end]
            ]

            if backend is None:
                heat_generation = 0.0
                progress_rates: list[float] = []

            elif isinstance(
                backend,
                ReactionNetworkBackend,
            ):
                heat_generation = backend.heat_generation(
                    temperature=temperature,
                    conversions=conversions,
                )

                progress_rates = backend.progress_rates(
                    temperature=temperature,
                    conversions=conversions,
                )

            else:
                heat_generation = backend.heat_generation(temperature)

                progress_rates = []

            temperature_rate = model.temperature_derivative(
                temperature=temperature,
                heat_generation=(heat_generation),
            )

            gas_rates: list[float] = [0.0 for _ in gas_species_names]

            if gas_generation_model is not None:
                generation_rates = gas_generation_model.generation_rates(
                    temperature=temperature,
                    conversions=conversions,
                    cell_mass=problem.cell.mass,
                )

                gas_rates = [
                    generation_rates.get(
                        species_name,
                        0.0,
                    )
                    for species_name in gas_species_names
                ]

            return (
                temperature_rate,
                progress_rates,
                gas_rates,
            )

        def closed_rhs(
            time: float,
            state: np.ndarray,
        ) -> list[float]:
            """Return ODE rates while the cell vent remains closed."""
            del time

            (
                temperature_rate,
                progress_rates,
                gas_rates,
            ) = base_rates(state)

            return [
                temperature_rate,
                *progress_rates,
                *gas_rates,
            ]

        def open_rhs(
            time: float,
            state: np.ndarray,
        ) -> list[float]:
            """Return ODE rates after irreversible vent opening."""
            del time

            if vent_model is None:
                raise RuntimeError("Vent model is not configured.")

            if initial_gas_inventory is None:
                raise RuntimeError("Initial gas inventory is missing.")

            (
                temperature_rate,
                progress_rates,
                generation_rates,
            ) = base_rates(state)

            temperature = float(state[0])

            current_moles = {
                species_name: max(
                    float(state[gas_start + index]),
                    0.0,
                )
                for index, species_name in enumerate(gas_species_names)
            }

            inventory = GasInventory(
                species=list(initial_gas_inventory.species),
                moles=current_moles,
            )

            upstream_pressure = pressure_from_state(state)

            vent_rates = vent_model.species_molar_flow_rates(
                inventory=inventory,
                upstream_pressure=(upstream_pressure),
                temperature=temperature,
            )

            net_gas_rates: list[float] = []

            for index, species_name in enumerate(gas_species_names):
                net_rate = generation_rates[index] - vent_rates.get(
                    species_name,
                    0.0,
                )

                state_amount = float(state[gas_start + index])

                if state_amount <= 0.0 and net_rate < 0.0:
                    net_rate = 0.0

                net_gas_rates.append(net_rate)

            return [
                temperature_rate,
                *progress_rates,
                *net_gas_rates,
            ]

        def add_temperature_event(
            events: list[Any],
        ) -> None:
            """Add the configured maximum-temperature termination event."""
            maximum_temperature = problem.maximum_temperature

            if maximum_temperature is None:
                return

            def temperature_event(
                time: float,
                state: np.ndarray,
            ) -> float:
                """Return the maximum-temperature event residual."""
                del time

                return maximum_temperature - float(state[0])

            temperature_event.terminal = True  # type: ignore[attr-defined]
            temperature_event.direction = -1.0  # type: ignore[attr-defined]

            events.append(temperature_event)

        def add_maximum_pressure_event(
            events: list[Any],
        ) -> None:
            """Add the configured maximum-pressure termination event."""
            if maximum_pressure is None:
                return

            def maximum_pressure_event(
                time: float,
                state: np.ndarray,
            ) -> float:
                """Return the maximum-pressure event residual."""
                del time

                return maximum_pressure - pressure_from_state(state)

            maximum_pressure_event.terminal = True  # type: ignore[attr-defined]
            maximum_pressure_event.direction = -1.0  # type: ignore[attr-defined]

            events.append(maximum_pressure_event)

        phase_1_events: list[Any] = []

        add_temperature_event(phase_1_events)

        add_maximum_pressure_event(phase_1_events)

        vent_event_index: int | None = None

        vent_initially_open = False

        initial_state_array = np.asarray(
            initial_state,
            dtype=float,
        )

        if vent_model is not None and vent_open_pressure is not None:
            initial_pressure = pressure_from_state(initial_state_array)

            vent_initially_open = initial_pressure >= vent_open_pressure

            if not vent_initially_open:

                def vent_open_event(
                    time: float,
                    state: np.ndarray,
                ) -> float:
                    """Return the vent-opening pressure event residual."""
                    del time

                    return vent_open_pressure - pressure_from_state(state)

                vent_open_event.terminal = True  # type: ignore[attr-defined]
                vent_open_event.direction = -1.0  # type: ignore[attr-defined]

                vent_event_index = len(phase_1_events)

                phase_1_events.append(vent_open_event)

        phase_1_solution = None

        vent_triggered = vent_initially_open

        if not vent_initially_open:
            phase_1_solution = solve_ivp(
                closed_rhs,
                (
                    0.0,
                    problem.duration,
                ),
                initial_state,
                method=self.method,
                rtol=self.relative_tolerance,
                atol=self.absolute_tolerance,
                events=phase_1_events or None,
            )

            if not phase_1_solution.success:
                raise RuntimeError(
                    f"ODE integration failed: {phase_1_solution.message}"
                )

            if vent_event_index is not None:
                vent_triggered = len(phase_1_solution.t_events[vent_event_index]) > 0

        phase_2_solution = None

        if vent_triggered:
            if vent_initially_open:
                phase_2_start_time = 0.0

                phase_2_initial_state = initial_state_array

            else:
                if phase_1_solution is None:
                    raise RuntimeError("Closed-phase solution is missing.")

                phase_2_start_time = float(phase_1_solution.t[-1])

                phase_2_initial_state = phase_1_solution.y[
                    :,
                    -1,
                ]

            if phase_2_start_time < problem.duration:
                phase_2_events: list[Any] = []

                add_temperature_event(phase_2_events)

                add_maximum_pressure_event(phase_2_events)

                phase_2_solution = solve_ivp(
                    open_rhs,
                    (
                        phase_2_start_time,
                        problem.duration,
                    ),
                    phase_2_initial_state,
                    method=self.method,
                    rtol=self.relative_tolerance,
                    atol=self.absolute_tolerance,
                    events=(phase_2_events or None),
                )

                if not phase_2_solution.success:
                    raise RuntimeError(
                        "ODE integration failed "
                        "after vent opening: "
                        f"{phase_2_solution.message}"
                    )

        if vent_initially_open:
            if phase_2_solution is None:
                time_values = np.asarray(
                    [0.0],
                    dtype=float,
                )

                state_values = initial_state_array[
                    :,
                    np.newaxis,
                ]

            else:
                time_values = phase_2_solution.t

                state_values = phase_2_solution.y

            vent_open_values = np.ones(
                len(time_values),
                dtype=float,
            )

        elif vent_triggered and phase_1_solution is not None:
            if phase_2_solution is None:
                time_values = phase_1_solution.t

                state_values = phase_1_solution.y

                vent_open_values = np.zeros(
                    len(time_values),
                    dtype=float,
                )

                vent_open_values[-1] = 1.0

            else:
                time_values = np.concatenate(
                    (
                        phase_1_solution.t,
                        phase_2_solution.t[1:],
                    )
                )

                state_values = np.concatenate(
                    (
                        phase_1_solution.y,
                        phase_2_solution.y[:, 1:],
                    ),
                    axis=1,
                )

                closed_flags = np.zeros(
                    len(phase_1_solution.t),
                    dtype=float,
                )

                closed_flags[-1] = 1.0

                open_flags = np.ones(
                    max(
                        len(phase_2_solution.t) - 1,
                        0,
                    ),
                    dtype=float,
                )

                vent_open_values = np.concatenate(
                    (
                        closed_flags,
                        open_flags,
                    )
                )

        else:
            if phase_1_solution is None:
                raise RuntimeError("ODE solution is missing.")

            time_values = phase_1_solution.t

            state_values = phase_1_solution.y

            vent_open_values = np.zeros(
                len(time_values),
                dtype=float,
            )

        results = Results(
            time=time_values,
        )

        results.add_variable(
            "temperature",
            state_values[0],
        )

        for index in range(n_reactions):
            conversion = np.clip(
                state_values[index + 1],
                0.0,
                1.0,
            )

            results.add_variable(
                f"conversion_{index}",
                conversion,
            )

        gas_arrays: list[np.ndarray] = []

        for index, species_name in enumerate(gas_species_names):
            gas_values = np.clip(
                state_values[gas_start + index],
                0.0,
                None,
            )

            gas_arrays.append(gas_values)

            results.add_variable(
                f"gas_{species_name}",
                gas_values,
            )

        if pressure_model is not None:
            pressure_values: list[float] = []

            for time_index in range(len(time_values)):
                pressure_values.append(pressure_from_state(state_values[:, time_index]))

            results.add_variable(
                "pressure",
                pressure_values,
            )

        if vent_model is not None:
            results.add_variable(
                "vent_open",
                vent_open_values,
            )

        return results

    def _solve_elemental(
        self,
        problem: ThermalProblem,
    ) -> Results:
        """Solve using conserved elemental gas inventories."""
        thermal_model = LumpedThermalModel(
            cell=problem.cell,
            convection_coefficient=(
                problem.convection_coefficient
            ),
            ambient_temperature=(
                problem.ambient_temperature
            ),
        )

        reaction_backend = problem.chemistry_backend

        if not isinstance(
            reaction_backend,
            ReactionNetworkBackend,
        ):
            raise TypeError(
                "Elemental thermochemistry requires a "
                "ReactionNetworkBackend."
            )

        element_generation_model = (
            problem.element_generation_model
        )

        uses_energy_scaled_generation = isinstance(
            element_generation_model,
            EnergyScaledElementGenerationModel,
        )
        
        element_release_model = problem.element_release_model

        uses_delayed_release = element_release_model is not None

        equilibrium_backend = (
            problem.gas_equilibrium_backend
        )

        initial_inventory = (
            problem.initial_element_inventory
        )

        gas_volume = problem.gas_volume

        if (
            element_generation_model is None
            or equilibrium_backend is None
            or initial_inventory is None
            or gas_volume is None
        ):
            raise RuntimeError(
                "Incomplete elemental thermochemistry "
                "configuration."
            )

        if (
            isinstance(
                element_generation_model,
                ElementGenerationModel,
            )
            and element_generation_model.reaction_network
            is not reaction_backend.reaction_network
        ):
            raise ValueError(
                "Element generation model must use the "
                "same reaction network as the chemistry "
                "backend."
            )

        reaction_count = len(
            reaction_backend.reaction_network.reactions
        )

        if problem.initial_conversions is None:
            initial_conversions = [
                0.0
                for _ in range(reaction_count)
            ]
        else:
            initial_conversions = list(
                problem.initial_conversions
            )

        if len(initial_conversions) != reaction_count:
            raise ValueError(
                "Number of initial conversions must match "
                "number of reactions."
            )

        if any(
            conversion < 0.0
            or conversion > 1.0
            for conversion in initial_conversions
        ):
            raise ValueError(
                "Initial conversions must be between 0 and 1."
            )



        gas_energy_model = (
            EquilibriumGasEnergyModel(
                equilibrium_backend=equilibrium_backend,
            )
        )

        coupled_energy_model = (
            CoupledThermalEnergyModel(
                cell=problem.cell,
                gas_energy_model=gas_energy_model,
            )
        )

        vent_energy_model = None

        if problem.element_vent_model is not None:
            vent_energy_model = (
                ElementVentEnergyModel(
                    element_vent_model=(
                        problem.element_vent_model
                    ),
                    gas_energy_model=gas_energy_model,
                )
            )


        element_names = list(
            initial_inventory.element_names
        )

        for element_name in (
            element_generation_model.element_names
        ):
            if element_name not in element_names:
                element_names.append(
                    element_name
                )

        element_count = len(element_names)

        supported_elements = set(
            equilibrium_backend.element_names
        )

        unsupported_elements = [
            element_name
            for element_name in element_names
            if element_name not in supported_elements
        ]

        if unsupported_elements:
            names = ", ".join(
                unsupported_elements
            )

            raise ValueError(
                "Gas equilibrium backend does not support "
                f"elements: {names}"
            )

        initial_element_moles = [
            initial_inventory.moles_of(element_name) for element_name in element_names
        ]

        initial_vented_moles = [0.0 for _ in element_names]

        initial_vented_energy = 0.0

        initial_energy = coupled_energy_model.total_energy(
            temperature=(problem.initial_temperature),
            volume=gas_volume,
            element_inventory=initial_inventory,
        )

        initial_thermal_loss_energy = 0.0

        if uses_energy_scaled_generation:
            initial_reaction_energy = [0.0]
        else:
            initial_reaction_energy = []

        if uses_delayed_release:
            initial_pending_moles = [0.0 for _ in element_names]
        else:
            initial_pending_moles = []

        initial_state = [
            initial_energy,
            *initial_conversions,
            *initial_reaction_energy,
            *initial_pending_moles,
            *initial_element_moles,
            *initial_vented_moles,
            initial_vented_energy,
            initial_thermal_loss_energy,
        ]

        conversion_start = 1

        conversion_end = conversion_start + reaction_count

        if uses_energy_scaled_generation:
            reaction_energy_index = conversion_end

            pending_start = reaction_energy_index + 1
        else:
            reaction_energy_index = None

            pending_start = conversion_end

        if uses_delayed_release:
            pending_end = pending_start + element_count
        else:
            pending_end = pending_start

        element_start = pending_end

        element_end = element_start + element_count

        vented_start = element_end

        vented_end = vented_start + element_count

        vented_energy_index = vented_end

        thermal_loss_energy_index = vented_energy_index + 1

        if uses_energy_scaled_generation:
            if not isinstance(
                element_generation_model,
                EnergyScaledElementGenerationModel,
            ):
                raise RuntimeError("Missing energy-scaled element generation model.")

            maximum_generation_energy = (
                element_generation_model.reference_reaction_energy_per_cell_mass
                * problem.cell.mass
            )
        else:
            maximum_generation_energy = None
        
        def conversions_from_state(
            state: np.ndarray,
        ) -> list[float]:
            return [
                min(
                    max(float(value), 0.0),
                    1.0,
                )
                for value in state[
                    conversion_start:
                    conversion_end
                ]
            ]

        def pending_inventory_from_state(
            state: np.ndarray,
        ) -> ElementInventory:
            """Return the stored unreleased elemental inventory."""
            if not uses_delayed_release:
                return ElementInventory(moles={})

            return ElementInventory(
                moles={
                    element_name: max(
                        float(state[pending_start + index]),
                        0.0,
                    )
                    for index, element_name in enumerate(element_names)
                }
            )

        def inventory_from_state(
            state: np.ndarray,
        ) -> ElementInventory:
            return ElementInventory(
                moles={
                    element_name: max(
                        float(
                            state[
                                element_start
                                + index
                            ]
                        ),
                        0.0,
                    )
                    for index, element_name
                    in enumerate(element_names)
                }
            )

        def temperature_and_equilibrium_from_state(
            state: np.ndarray,
        ) -> tuple[
            float,
            ElementInventory,
            GasEquilibriumState,
        ]:
            inventory = inventory_from_state(state)

            (
                temperature,
                equilibrium_state,
            ) = coupled_energy_model.temperature_and_state_from_energy(
                energy=float(state[0]),
                volume=gas_volume,
                element_inventory=inventory,
            )

            return (
                temperature,
                inventory,
                equilibrium_state,
            )

        def temperature_from_state(
            state: np.ndarray,
        ) -> float:
            inventory = inventory_from_state(
                state
            )

            return (
                coupled_energy_model.temperature_from_energy(
                    energy=float(state[0]),
                    volume=gas_volume,
                    element_inventory=inventory,
                )
            )

        def equilibrium_from_state(
            state: np.ndarray,
            temperature: float | None = None,
            inventory: ElementInventory | None = None,
        ) -> GasEquilibriumState:
            if temperature is None and inventory is None:
                _, _, equilibrium_state = temperature_and_equilibrium_from_state(state)

                return equilibrium_state

            if inventory is None:
                inventory = inventory_from_state(state)

            if temperature is None:
                temperature = coupled_energy_model.temperature_from_energy(
                    energy=float(state[0]),
                    volume=gas_volume,
                    element_inventory=inventory,
                )

            return equilibrium_backend.equilibrate_tv(
                temperature=temperature,
                volume=gas_volume,
                element_inventory=inventory,
            )

        def pressure_from_state(
            state: np.ndarray,
        ) -> float:
            return equilibrium_from_state(
                state
            ).pressure

        def base_rates(
            state: np.ndarray,
            temperature: float | None = None,
        ) -> list[float]:
            if temperature is None:
                temperature = temperature_from_state(state)

            conversions = (
                conversions_from_state(
                    state
                )
            )

            heat_generation = (
                reaction_backend.heat_generation(
                    temperature=temperature,
                    conversions=conversions,
                )
            )

            progress_rates = (
                reaction_backend.progress_rates(
                    temperature=temperature,
                    conversions=conversions,
                )
            )

            heat_loss = thermal_model.heat_loss(temperature)

            energy_rate = heat_generation - heat_loss

            if uses_energy_scaled_generation:
                reaction_energy_rates = [
                    max(
                        heat_generation,
                        0.0,
                    )
                ]
            else:
                reaction_energy_rates = []

            if isinstance(
                element_generation_model,
                EnergyScaledElementGenerationModel,
            ):
                if reaction_energy_index is None or maximum_generation_energy is None:
                    raise RuntimeError(
                        "Missing reaction-energy state for "
                        "energy-scaled element generation."
                    )

                released_reaction_energy = max(
                    float(state[reaction_energy_index]),
                    0.0,
                )

                if released_reaction_energy >= maximum_generation_energy:
                    generated_element_rates = {
                        element_name: 0.0 for element_name in element_names
                    }

                else:
                    generated_element_rates = element_generation_model.generation_rates(
                        heat_generation=(heat_generation),
                    )

            else:
                generated_element_rates = element_generation_model.generation_rates(
                    temperature=temperature,
                    conversions=conversions,
                    cell_mass=problem.cell.mass,
                )

            pending_rates: list[float] = []

            if uses_delayed_release:
                if element_release_model is None:
                    raise RuntimeError("Missing element release model.")

                pending_inventory = pending_inventory_from_state(state)

                released_element_rates = element_release_model.release_rates(
                    pending_inventory
                )

                pending_rates = [
                    (
                        generated_element_rates.get(
                            element_name,
                            0.0,
                        )
                        - released_element_rates.get(
                            element_name,
                            0.0,
                        )
                    )
                    for element_name in element_names
                ]

                element_rates = [
                    released_element_rates.get(
                        element_name,
                        0.0,
                    )
                    for element_name in element_names
                ]

            else:
                element_rates = [
                    generated_element_rates.get(
                        element_name,
                        0.0,
                    )
                    for element_name in element_names
                ]

            vented_rates = [0.0 for _ in element_names]

            return [
                energy_rate,
                *progress_rates,
                *reaction_energy_rates,
                *pending_rates,
                *element_rates,
                *vented_rates,
                0.0,
                heat_loss,
            ]

        def closed_rhs(
            time: float,
            state: np.ndarray,
        ) -> list[float]:
            del time

            return base_rates(
                state
            )

        def open_rhs(
            time: float,
            state: np.ndarray,
        ) -> list[float]:
            del time

            (
                temperature,
                inventory,
                equilibrium_state,
            ) = temperature_and_equilibrium_from_state(state)

            rates = base_rates(
                state,
                temperature=temperature,
            )

            element_vent_model = problem.element_vent_model

            if element_vent_model is None:
                return rates

            vent_rates = (
                element_vent_model.element_molar_flow_rates(
                    equilibrium_state
                )
            )

            for index, element_name in enumerate(
                element_names
            ):
                element_state_index = element_start + index

                element_rate_index = element_start + index

                vented_rate_index = vented_start + index

                vent_rate = vent_rates.get(
                    element_name,
                    0.0,
                )

                net_rate = rates[element_rate_index] - vent_rate

                current_amount = max(
                    float(state[element_state_index]),
                    0.0,
                )

                if current_amount <= 0.0 and net_rate < 0.0:
                    vent_rate = max(
                        rates[element_rate_index],
                        0.0,
                    )

                    net_rate = rates[element_rate_index] - vent_rate

                rates[element_rate_index] = net_rate

                rates[vented_rate_index] = vent_rate

            if vent_energy_model is not None:
                vent_energy_rate = vent_energy_model.energy_flow_rate(
                    state=equilibrium_state,
                    element_inventory=inventory,
                )

                rates[0] -= vent_energy_rate

                rates[vented_energy_index] = vent_energy_rate

            return rates

        def temperature_event(
            time: float,
            state: np.ndarray,
        ) -> float:
            del time

            maximum_temperature = (
                problem.maximum_temperature
            )

            if maximum_temperature is None:
                return 1.0

            return maximum_temperature - temperature_from_state(state)

        temperature_event.terminal = True  # type: ignore[attr-defined]
        temperature_event.direction = -1.0  # type: ignore[attr-defined]

        def maximum_pressure_event(
            time: float,
            state: np.ndarray,
        ) -> float:
            del time

            maximum_pressure = (
                problem.maximum_pressure
            )

            if maximum_pressure is None:
                return 1.0

            return (
                maximum_pressure
                - pressure_from_state(
                    state
                )
            )

        maximum_pressure_event.terminal = True  # type: ignore[attr-defined]
        maximum_pressure_event.direction = -1.0  # type: ignore[attr-defined]

        closed_events = []

        if problem.maximum_temperature is not None:
            closed_events.append(
                temperature_event
            )

        if problem.maximum_pressure is not None:
            closed_events.append(
                maximum_pressure_event
            )

        vent_event_index: int | None = None

        if (
                problem.element_vent_model is not None
                and problem.vent_open_pressure is not None
        ):

            def vent_open_event(
                    time: float,
                    state: np.ndarray,
            ) -> float:
                del time

                return (
                        pressure_from_state(
                            state
                        )
                        - problem.vent_open_pressure
                )

            vent_open_event.terminal = True  # type: ignore[attr-defined]
            vent_open_event.direction = 1.0  # type: ignore[attr-defined]

            vent_event_index = len(
                closed_events
            )

            closed_events.append(
                vent_open_event
            )

        closed_solution = solve_ivp(
            closed_rhs,
            (
                0.0,
                problem.duration,
            ),
            initial_state,
            method=self.method,
            rtol=self.relative_tolerance,
            atol=self.absolute_tolerance,
            events=(
                    closed_events
                    or None
            ),
        )

        if not closed_solution.success:
            raise RuntimeError(
                "ODE integration failed: "
                f"{closed_solution.message}"
            )

        vent_open_time: float | None = None

        open_solution = None

        if vent_event_index is not None:
            vent_times = (
                closed_solution.t_events[
                    vent_event_index
                ]
            )

            if len(vent_times) > 0:
                vent_open_time = float(
                    vent_times[0]
                )

                if (
                        closed_solution.y_events
                        is None
                ):
                    raise RuntimeError(
                        "Missing vent-opening event state."
                    )

                open_initial_state = (
                    closed_solution.y_events[
                        vent_event_index
                    ][0]
                )

                if (
                        vent_open_time
                        < problem.duration
                ):
                    open_events = []

                    if (
                            problem.maximum_temperature
                            is not None
                    ):
                        open_events.append(
                            temperature_event
                        )

                    if (
                            problem.maximum_pressure
                            is not None
                    ):
                        open_events.append(
                            maximum_pressure_event
                        )

                    open_solution = solve_ivp(
                        open_rhs,
                        (
                            vent_open_time,
                            problem.duration,
                        ),
                        open_initial_state,
                        method=self.method,
                        rtol=self.relative_tolerance,
                        atol=self.absolute_tolerance,
                        events=(
                                open_events
                                or None
                        ),
                    )

                    if not open_solution.success:
                        raise RuntimeError(
                            "ODE integration failed after "
                            "vent opening: "
                            f"{open_solution.message}"
                        )

        time = closed_solution.t
        state_history = closed_solution.y

        if (
            open_solution is not None
            and open_solution.t.size > 1
        ):
            time = np.concatenate(
                [
                    closed_solution.t,
                    open_solution.t[1:],
                ]
            )

            state_history = np.concatenate(
                [
                    closed_solution.y,
                    open_solution.y[:, 1:],
                ],
                axis=1,
            )

        results = Results(
            time=time,
        )

        temperature_history = np.empty(
            time.size,
            dtype=float,
        )

        for index in range(
            time.size
        ):
            temperature_history[index] = (
                temperature_from_state(
                    state_history[:, index]
                )
            )

        results.add_variable(
            "temperature",
            temperature_history,
        )

        results.add_variable(
            "thermal_energy",
            state_history[0],
        )

        if reaction_energy_index is not None:
            results.add_variable(
                "reaction_energy",
                np.maximum(
                    state_history[reaction_energy_index],
                    0.0,
                ),
            )

        results.add_variable(
            "vented_thermal_energy",
            np.maximum(
                state_history[
                    vented_energy_index
                ],
                0.0,
            ),
        )

        for index in range(
            reaction_count
        ):
            results.add_variable(
                f"conversion_{index}",
                np.clip(
                    state_history[
                        1 + index
                    ],
                    0.0,
                    1.0,
                ),
            )

        if uses_delayed_release:
            for index, element_name in enumerate(element_names):
                results.add_variable(
                    f"pending_element_{element_name}",
                    state_history[pending_start + index],
                )

        for index, element_name in enumerate(
            element_names
        ):
            results.add_variable(
                f"element_{element_name}",
                np.maximum(
                    state_history[
                        element_start
                        + index
                    ],
                    0.0,
                ),
            )

        for index, element_name in enumerate(element_names):
            results.add_variable(
                f"vented_element_{element_name}",
                np.maximum(
                    state_history[vented_start + index],
                    0.0,
                ),
            )

        results.add_variable(
            "thermal_loss_energy",
            state_history[thermal_loss_energy_index],
        )

        pressure_history = np.empty(
            time.size,
            dtype=float,
        )

        density_history = np.empty(
            time.size,
            dtype=float,
        )

        mean_molar_mass_history = np.empty(
            time.size,
            dtype=float,
        )

        cp_mass_history = np.empty(
            time.size,
            dtype=float,
        )

        cv_mass_history = np.empty(
            time.size,
            dtype=float,
        )

        gamma_history = np.empty(
            time.size,
            dtype=float,
        )

        internal_energy_history = np.empty(
            time.size,
            dtype=float,
        )

        species_histories = {
            species_name: np.zeros(
                time.size,
                dtype=float,
            )
            for species_name in equilibrium_backend.species_names
        }

        for index in range(time.size):
            equilibrium_state = equilibrium_from_state(state_history[:, index])

            pressure_history[index] = equilibrium_state.pressure

            density_history[index] = equilibrium_state.density

            mean_molar_mass_history[index] = equilibrium_state.mean_molar_mass

            cp_mass_history[index] = equilibrium_state.cp_mass

            cv_mass_history[index] = equilibrium_state.cv_mass

            gamma_history[index] = equilibrium_state.heat_capacity_ratio

            internal_energy_history[index] = equilibrium_state.internal_energy

            for (
                species_name,
                species_moles,
            ) in equilibrium_state.species_moles.items():
                species_histories[species_name][index] = species_moles

        results.add_variable(
            "pressure",
            pressure_history,
        )

        results.add_variable(
            "gas_density",
            density_history,
        )

        results.add_variable(
            "gas_mean_molar_mass",
            mean_molar_mass_history,
        )

        results.add_variable(
            "gas_cp_mass",
            cp_mass_history,
        )

        results.add_variable(
            "gas_cv_mass",
            cv_mass_history,
        )

        results.add_variable(
            "gas_heat_capacity_ratio",
            gamma_history,
        )

        results.add_variable(
            "gas_internal_energy",
            internal_energy_history,
        )

        for (
            species_name,
            species_history,
        ) in species_histories.items():
            results.add_variable(
                f"equilibrium_species_{species_name}",
                species_history,
            )

        vent_open = np.zeros(
            time.size,
            dtype=float,
        )

        if vent_open_time is not None:
            vent_open[
                time >= vent_open_time
            ] = 1.0

        results.add_variable(
            "vent_open",
            vent_open,
        )

        return results