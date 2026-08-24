from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from liiontr.chemistry import ReactionNetworkBackend
from liiontr.core.results import Results
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

        gas_species_names: list[str] = []

        if gas_generation_model is not None:
            if not isinstance(
                problem.chemistry_backend,
                ReactionNetworkBackend,
            ):
                raise ValueError("Gas generation requires a ReactionNetworkBackend.")

            if (
                gas_generation_model.reaction_network
                is not problem.chemistry_backend.reaction_network
            ):
                raise ValueError(
                    "Gas generation model must use the "
                    "same reaction network as the "
                    "chemistry backend."
                )

            gas_species_names = gas_generation_model.species_names

        initial_gas_moles = [0.0 for _ in gas_species_names]

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

            generated_moles = sum(max(float(value), 0.0) for value in state[gas_start:])

            return pressure_model.pressure(
                temperature=temperature,
                generated_moles=generated_moles,
            )

        def rhs(
            time: float,
            state: np.ndarray,
        ) -> list[float]:
            del time

            temperature = float(state[0])

            conversions = [
                min(
                    max(float(value), 0.0),
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
                heat_generation=heat_generation,
            )

            gas_rates: list[float] = []

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

            return [
                temperature_rate,
                *progress_rates,
                *gas_rates,
            ]

        events: list[Any] = []

        maximum_temperature = problem.maximum_temperature

        if maximum_temperature is not None:

            def temperature_event(
                time: float,
                state: np.ndarray,
            ) -> float:
                del time

                return maximum_temperature - float(state[0])

            temperature_event.terminal = True  # type: ignore[attr-defined]
            temperature_event.direction = -1.0  # type: ignore[attr-defined]

            events.append(temperature_event)

        maximum_pressure = problem.maximum_pressure

        if maximum_pressure is not None:

            def pressure_event(
                time: float,
                state: np.ndarray,
            ) -> float:
                del time

                return maximum_pressure - pressure_from_state(state)

            pressure_event.terminal = True  # type: ignore[attr-defined]
            pressure_event.direction = -1.0  # type: ignore[attr-defined]

            events.append(pressure_event)

        solution = solve_ivp(
            rhs,
            (
                0.0,
                problem.duration,
            ),
            initial_state,
            method=self.method,
            rtol=self.relative_tolerance,
            atol=self.absolute_tolerance,
            events=events or None,
        )

        if not solution.success:
            raise RuntimeError(f"ODE integration failed: {solution.message}")

        results = Results(
            time=solution.t,
        )

        results.add_variable(
            "temperature",
            solution.y[0],
        )

        for index in range(n_reactions):
            conversion = np.clip(
                solution.y[index + 1],
                0.0,
                1.0,
            )

            results.add_variable(
                f"conversion_{index}",
                conversion,
            )

            gas_arrays = []

            for index, species_name in enumerate(gas_species_names):
                gas_values = solution.y[gas_start + index].clip(min=0.0)

                gas_arrays.append(gas_values)

                results.add_variable(
                    f"gas_{species_name}",
                    gas_values,
                )

                if pressure_model is not None:
                    pressure_values = []

                    for time_index, temperature in enumerate(solution.y[0]):
                        generated_moles = sum(
                            float(gas_values[time_index]) for gas_values in gas_arrays
                        )

                        pressure_values.append(
                            pressure_model.pressure(
                                temperature=float(temperature),
                                generated_moles=generated_moles,
                            )
                        )

                    results.add_variable(
                        "pressure",
                        pressure_values,
                    )
        return results
