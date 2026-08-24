from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from liiontr.chemistry import ReactionNetworkBackend
from liiontr.core.results import Results
from liiontr.gases import GasInventory
from liiontr.problems.thermal import ThermalProblem
from liiontr.thermal.lumped import LumpedThermalModel


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
            maximum_temperature = problem.maximum_temperature

            if maximum_temperature is None:
                return

            def temperature_event(
                time: float,
                state: np.ndarray,
            ) -> float:
                del time

                return maximum_temperature - float(state[0])

            temperature_event.terminal = True  # type: ignore[attr-defined]
            temperature_event.direction = -1.0  # type: ignore[attr-defined]

            events.append(temperature_event)

        def add_maximum_pressure_event(
            events: list[Any],
        ) -> None:
            if maximum_pressure is None:
                return

            def maximum_pressure_event(
                time: float,
                state: np.ndarray,
            ) -> float:
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
