"""Generate reproducible numerical figures for the LiionTR report."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid

from liiontr.chemistry.reaction_backend import ReactionNetworkBackend
from liiontr.core.results import Results
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
from liiontr.library.cells import cell_21700_generic
from liiontr.library.hu2020 import (
    hu2020_initial_conversions,
    hu2020_reaction_network,
)
from liiontr.problems.thermal import ThermalProblem
from liiontr.reactions import Reaction, ReactionNetwork
from liiontr.solver.scipy_solver import ScipySolver
from liiontr.thermal.lumped import LumpedThermalModel

PROJECT_ROOT = Path(__file__).resolve().parents[1]

REPORT_DIR = PROJECT_ROOT / "report"
FIGURE_DIR = REPORT_DIR / "figures"
DATA_DIR = REPORT_DIR / "data"
GENERATED_DIR = REPORT_DIR / "generated"

REACTION_LABELS = [
    "SEI",
    "Anode-electrolyte",
    "Cathode",
    "Electrolyte",
]

def build_synthetic_vent_problem(
    *,
    with_vent: bool,
) -> ThermalProblem:
    """Build the synthetic coupled gas-pressure-venting problem."""
    cell = cell_21700_generic()

    reaction = Reaction(
        name="Gas reaction",
        kinetics=Arrhenius(
            activation_energy=1.0,
            pre_exponential_factor=2.0,
        ),
        enthalpy=10000.0,
        mass_fraction=0.10,
    )

    reaction_network = ReactionNetwork(
        reactions=[reaction]
    )

    chemistry_backend = ReactionNetworkBackend(
        reaction_network=reaction_network,
        cell=cell,
    )

    gas_generation_model = GasGenerationModel(
        reaction_network=reaction_network,
        gas_yields=[
            ReactionGasYield(
                reaction_name="Gas reaction",
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

    vent_model = None
    vent_open_pressure = None

    if with_vent:
        vent_model = MixtureVentFlowModel(
            flow_model=CompressibleVentFlowModel(
                vent_area=1.0e-5,
                discharge_coefficient=0.8,
                heat_capacity_ratio=1.30,
            ),
            downstream_pressure=101325.0,
        )

        vent_open_pressure = 200000.0

    return ThermalProblem(
        cell=cell,
        chemistry_backend=chemistry_backend,
        gas_generation_model=gas_generation_model,
        initial_gas_inventory=initial_gas_inventory,
        pressure_model=pressure_model,
        vent_model=vent_model,
        vent_open_pressure=vent_open_pressure,
        initial_temperature=400.0,
        initial_conversions=[0.0],
        ambient_temperature=400.0,
        convection_coefficient=0.0,
        duration=1.0,
    )

def calculate_synthetic_vent_histories(
    problem: ThermalProblem,
    results: Results,
) -> dict[str, np.ndarray]:
    """Reconstruct gas-generation and vent-flow histories."""
    if problem.initial_gas_inventory is None:
        raise RuntimeError(
            "Synthetic vent problem requires an initial gas inventory."
        )

    if problem.vent_model is None:
        raise RuntimeError(
            "Synthetic vent problem requires a vent model."
        )

    if problem.gas_generation_model is None:
        raise RuntimeError(
            "Synthetic vent problem requires a gas-generation model."
        )

    time = np.asarray(
        results.time,
        dtype=float,
    )

    temperature = np.asarray(
        results.temperature,
        dtype=float,
    )

    pressure = np.asarray(
        results.get_variable("pressure"),
        dtype=float,
    )

    vent_open = np.asarray(
        results.get_variable("vent_open"),
        dtype=float,
    )

    conversion = np.asarray(
        results.get_variable("conversion_0"),
        dtype=float,
    )

    n2_moles = np.asarray(
        results.get_variable("gas_N2"),
        dtype=float,
    )

    co2_moles = np.asarray(
        results.get_variable("gas_CO2"),
        dtype=float,
    )

    total_moles = (
        n2_moles
        + co2_moles
    )

    n2_mole_fraction = np.divide(
        n2_moles,
        total_moles,
        out=np.zeros_like(total_moles),
        where=total_moles > 0.0,
    )

    co2_mole_fraction = np.divide(
        co2_moles,
        total_moles,
        out=np.zeros_like(total_moles),
        where=total_moles > 0.0,
    )

    downstream_pressure = problem.vent_model.downstream_pressure

    critical_pressure_ratio = problem.vent_model.flow_model.critical_pressure_ratio

    pressure_ratio = downstream_pressure / pressure

    choked_flow = (
        (vent_open >= 0.5)
        & (pressure > downstream_pressure)
        & (pressure_ratio <= critical_pressure_ratio)
    )

    n2_generation_rate = np.zeros_like(time)
    co2_generation_rate = np.zeros_like(time)
    total_generation_rate = np.zeros_like(time)

    n2_molar_flow = np.zeros_like(time)
    co2_molar_flow = np.zeros_like(time)
    total_molar_flow = np.zeros_like(time)

    n2_mass_flow = np.zeros_like(time)
    co2_mass_flow = np.zeros_like(time)
    total_mass_flow = np.zeros_like(time)

    species = (
        problem.initial_gas_inventory.species
    )

    species_by_name = {
        gas_species.name: gas_species
        for gas_species in species
    }

    for index in range(len(time)):
        generation_rates = (
            problem.gas_generation_model.generation_rates(
                temperature=float(
                    temperature[index]
                ),
                conversions=[
                    float(conversion[index])
                ],
                cell_mass=problem.cell.mass,
            )
        )

        n2_generation_rate[index] = (
            generation_rates.get(
                "N2",
                0.0,
            )
        )

        co2_generation_rate[index] = (
            generation_rates.get(
                "CO2",
                0.0,
            )
        )

        total_generation_rate[index] = (
            n2_generation_rate[index]
            + co2_generation_rate[index]
        )

        if vent_open[index] < 0.5:
            continue

        inventory = GasInventory(
            species=list(species),
            moles={
                "N2": float(
                    n2_moles[index]
                ),
                "CO2": float(
                    co2_moles[index]
                ),
            },
        )

        molar_rates = (
            problem.vent_model.species_molar_flow_rates(
                inventory=inventory,
                upstream_pressure=float(
                    pressure[index]
                ),
                temperature=float(
                    temperature[index]
                ),
            )
        )

        n2_molar_flow[index] = (
            molar_rates.get(
                "N2",
                0.0,
            )
        )

        co2_molar_flow[index] = (
            molar_rates.get(
                "CO2",
                0.0,
            )
        )

        total_molar_flow[index] = (
            n2_molar_flow[index]
            + co2_molar_flow[index]
        )

        n2_mass_flow[index] = (
            n2_molar_flow[index]
            * species_by_name[
                "N2"
            ].molar_mass
        )

        co2_mass_flow[index] = (
            co2_molar_flow[index]
            * species_by_name[
                "CO2"
            ].molar_mass
        )

        total_mass_flow[index] = (
            n2_mass_flow[index]
            + co2_mass_flow[index]
        )

    return {
        "time_s": time,
        "temperature_K": temperature,
        "pressure_Pa": pressure,
        "vent_open": vent_open,
        "conversion": conversion,
        "n2_moles": n2_moles,
        "co2_moles": co2_moles,
        "total_moles": total_moles,
        "n2_generation_rate_mol_per_s": n2_generation_rate,
        "co2_generation_rate_mol_per_s": co2_generation_rate,
        "total_generation_rate_mol_per_s": total_generation_rate,
        "n2_molar_flow_mol_per_s": n2_molar_flow,
        "co2_molar_flow_mol_per_s": co2_molar_flow,
        "total_molar_flow_mol_per_s": total_molar_flow,
        "n2_mass_flow_kg_per_s": n2_mass_flow,
        "co2_mass_flow_kg_per_s": co2_mass_flow,
        "total_mass_flow_kg_per_s": total_mass_flow,
        "n2_mole_fraction": n2_mole_fraction,
        "co2_mole_fraction": co2_mole_fraction,
        "pressure_ratio": pressure_ratio,
        "choked_flow": choked_flow,
    }

def calculate_synthetic_gas_balance(
    problem: ThermalProblem,
    histories: dict[str, np.ndarray],
) -> dict[str, float]:
    """Calculate the global molar balance of the synthetic vent case."""
    if problem.initial_gas_inventory is None:
        raise RuntimeError(
            "Gas balance requires an initial gas inventory."
        )

    time = histories["time_s"]

    initial_moles = float(
        problem.initial_gas_inventory.total_moles
    )

    generated_moles = float(
        np.trapezoid(
            histories[
                "total_generation_rate_mol_per_s"
            ],
            time,
        )
    )

    vent_open = histories["vent_open"]

    open_indices = np.flatnonzero(vent_open >= 0.5)

    if len(open_indices) == 0:
        vented_moles = 0.0
    else:
        opening_index = int(open_indices[0])

        vented_moles = float(
            np.trapezoid(
                histories["total_molar_flow_mol_per_s"][opening_index:],
                time[opening_index:],
            )
        )

    final_moles = float(
        histories["total_moles"][-1]
    )

    residual = (
        initial_moles
        + generated_moles
        - vented_moles
        - final_moles
    )

    total_available_moles = (
        initial_moles
        + generated_moles
    )

    if total_available_moles > 0.0:
        relative_residual = (
            residual
            / total_available_moles
        )
    else:
        relative_residual = 0.0

    return {
        "initial_moles": initial_moles,
        "generated_moles": generated_moles,
        "vented_moles": vented_moles,
        "final_moles": final_moles,
        "residual_moles": residual,
        "relative_residual": relative_residual,
    }

def calculate_synthetic_flow_regime_summary(
    problem: ThermalProblem,
    histories: dict[str, np.ndarray],
) -> dict[str, float]:
    """Calculate characteristic vent-flow-regime values."""
    if problem.vent_model is None:
        raise RuntimeError(
            "Flow-regime analysis requires a vent model."
        )

    time = histories["time_s"]
    pressure = histories["pressure_Pa"]
    vent_open = histories["vent_open"]
    choked_flow = histories["choked_flow"]

    flow_model = (
        problem.vent_model.flow_model
    )

    downstream_pressure = (
        problem.vent_model.downstream_pressure
    )

    critical_ratio = (
        flow_model.critical_pressure_ratio
    )

    critical_upstream_pressure = (
        downstream_pressure
        / critical_ratio
    )

    open_indices = np.flatnonzero(
        vent_open >= 0.5
    )

    choked_indices = np.flatnonzero(
        choked_flow
    )

    if len(open_indices) == 0:
        raise RuntimeError(
            "Vent never opened."
        )

    vent_open_index = int(
        open_indices[0]
    )

    if len(choked_indices) > 0:
        choked_start_time = float(
            time[choked_indices[0]]
        )

        choked_end_time = float(
            time[choked_indices[-1]]
        )
    else:
        choked_start_time = float("nan")
        choked_end_time = float("nan")

    

    crossing_indices = np.flatnonzero(
        (np.arange(len(time)) > vent_open_index)
        & (pressure < critical_upstream_pressure)
    )

    if len(crossing_indices) > 0:
        after_index = int(crossing_indices[0])

        before_index = after_index - 1

        time_before = float(time[before_index])

        time_after = float(time[after_index])

        pressure_before = float(pressure[before_index])

        pressure_after = float(pressure[after_index])

        if pressure_after != pressure_before:
            fraction = (critical_upstream_pressure - pressure_before) / (
                pressure_after - pressure_before
            )

            transition_time = time_before + fraction * (time_after - time_before)
        else:
            transition_time = time_after

        transition_pressure = critical_upstream_pressure
    else:
        transition_time = float("nan")
        transition_pressure = float("nan")

    if np.isfinite(transition_time):
        choked_end_time = transition_time

    if np.isfinite(transition_time):
        choked_duration = transition_time - float(time[vent_open_index])
    else:
        choked_duration = float("nan")

    return {
        "critical_pressure_ratio": float(critical_ratio),
        "critical_upstream_pressure_Pa": float(critical_upstream_pressure),
        "choked_start_time_s": choked_start_time,
        "choked_end_time_s": choked_end_time,
        "choked_duration_s": choked_duration,
        "transition_time_s": transition_time,
        "transition_pressure_Pa": transition_pressure,
        "initial_co2_mole_fraction": float(histories["co2_mole_fraction"][0]),
        "peak_co2_mole_fraction": float(np.max(histories["co2_mole_fraction"])),
        "final_co2_mole_fraction": float(histories["co2_mole_fraction"][-1]),
    }

def calculate_synthetic_vent_summary(
    histories: dict[str, np.ndarray],
    closed_results: Results,
) -> dict[str, float]:
    """Calculate characteristic values of the venting simulation."""
    time = histories["time_s"]
    temperature = histories["temperature_K"]
    pressure = histories["pressure_Pa"]
    vent_open = histories["vent_open"]
    total_moles = histories["total_moles"]
    mass_flow = histories[
        "total_mass_flow_kg_per_s"
    ]

    open_indices = np.flatnonzero(
        vent_open >= 0.5
    )

    if len(open_indices) == 0:
        raise RuntimeError(
            "Vent did not open in the synthetic reference case."
        )

    opening_index = int(open_indices[0])

    closed_pressure = np.asarray(
        closed_results.get_variable("pressure"),
        dtype=float,
    )

    return {
        "vent_open_time_s": float(
            time[opening_index]
        ),
        "vent_open_temperature_K": float(
            temperature[opening_index]
        ),
        "peak_pressure_Pa": float(
            np.max(pressure)
        ),
        "final_pressure_Pa": float(
            pressure[-1]
        ),
        "closed_final_pressure_Pa": float(
            closed_pressure[-1]
        ),
        "maximum_total_moles": float(
            np.max(total_moles)
        ),
        "final_total_moles": float(
            total_moles[-1]
        ),
        "peak_mass_flow_kg_per_s": float(
            np.max(mass_flow)
        ),
    }

def save_synthetic_vent_temperature_figure(
    histories: dict[str, np.ndarray],
    vent_open_time: float,
) -> None:
    """Save the synthetic venting temperature history."""
    figure, axis = plt.subplots(
        figsize=(6.4, 4.0)
    )

    axis.plot(
        histories["time_s"],
        histories["temperature_K"],
        linewidth=1.8,
    )

    axis.axvline(
        vent_open_time,
        linestyle="--",
        label="Vent opening",
    )

    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Temperature (K)")
    axis.grid(True, alpha=0.25)
    axis.legend()

    figure.tight_layout()

    figure.savefig(
        FIGURE_DIR
        / "synthetic_vent_temperature.pdf",
        bbox_inches="tight",
    )

    plt.close(figure)

def save_synthetic_vent_pressure_figure(
    histories: dict[str, np.ndarray],
    vent_open_time: float,
    vent_open_pressure: float,
) -> None:
    """Save the synthetic vented-system pressure history."""
    figure, axis = plt.subplots(
        figsize=(6.4, 4.2)
    )

    axis.plot(
        histories["time_s"],
        histories["pressure_Pa"] / 1000.0,
        label="Internal pressure",
        linewidth=2.0,
    )

    axis.axhline(
        vent_open_pressure / 1000.0,
        label="Vent-opening pressure",
        linestyle="--",
    )

    axis.axhline(
        101325.0 / 1000.0,
        label="Ambient pressure",
        linestyle=":",
    )

    axis.axvline(
        vent_open_time,
        label="Vent opening",
        linestyle="-.",
    )

    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Absolute pressure (kPa)")
    axis.grid(True, alpha=0.25)
    axis.legend()

    figure.tight_layout()

    figure.savefig(
        FIGURE_DIR
        / "synthetic_vent_pressure.pdf",
        bbox_inches="tight",
    )

    plt.close(figure)

def save_synthetic_vent_gas_inventory_figure(
    histories: dict[str, np.ndarray],
    vent_open_time: float,
) -> None:
    """Save the internal gas-inventory histories."""
    figure, axis = plt.subplots(
        figsize=(6.4, 4.2)
    )

    axis.plot(
        histories["time_s"],
        histories["n2_moles"] * 1000.0,
        label="N2",
        linewidth=1.6,
    )

    axis.plot(
        histories["time_s"],
        histories["co2_moles"] * 1000.0,
        label="CO2",
        linewidth=1.6,
    )

    axis.plot(
        histories["time_s"],
        histories["total_moles"] * 1000.0,
        label="Total",
        linewidth=2.0,
        linestyle="--",
    )

    axis.axvline(
        vent_open_time,
        label="Vent opening",
        linestyle=":",
    )

    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Gas amount (mmol)")
    axis.grid(True, alpha=0.25)
    axis.legend()

    figure.tight_layout()

    figure.savefig(
        FIGURE_DIR
        / "synthetic_vent_gas_inventory.pdf",
        bbox_inches="tight",
    )

    plt.close(figure)

def save_synthetic_vent_flow_figure(
    histories: dict[str, np.ndarray],
    vent_open_time: float,
    choked_end_time: float,
) -> None:
    """Save the reconstructed vent mass-flow history."""
    figure, axis = plt.subplots(
        figsize=(6.4, 4.2)
    )

    conversion = 1.0e6
    if np.isfinite(choked_end_time):
        axis.axvspan(
            vent_open_time,
            choked_end_time,
            alpha=0.12,
            label="Choked flow",
        )
    axis.plot(
        histories["time_s"],
        (
            histories[
                "total_mass_flow_kg_per_s"
            ]
            * conversion
        ),
        label="Total",
        linewidth=2.0,
    )

    axis.plot(
        histories["time_s"],
        (
            histories[
                "co2_mass_flow_kg_per_s"
            ]
            * conversion
        ),
        label="CO2",
        linewidth=1.5,
        linestyle="--",
    )

    axis.plot(
        histories["time_s"],
        (
            histories[
                "n2_mass_flow_kg_per_s"
            ]
            * conversion
        ),
        label="N2",
        linewidth=1.5,
        linestyle=":",
    )

    axis.axvline(
        vent_open_time,
        label="Vent opening",
        linestyle="-.",
    )

    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Vent mass flow (mg/s)")
    axis.grid(True, alpha=0.25)
    axis.legend()

    figure.tight_layout()

    figure.savefig(
        FIGURE_DIR
        / "synthetic_vent_flow.pdf",
        bbox_inches="tight",
    )

    plt.close(figure)

def save_synthetic_vent_csv(
    histories: dict[str, np.ndarray],
    closed_results: Results,
) -> None:
    """Save the synthetic venting reference data."""
    output_path = (
        DATA_DIR
        / "synthetic_vent_reference.csv"
    )

    closed_time = np.asarray(
        closed_results.time,
        dtype=float,
    )
    closed_pressure = np.asarray(
        closed_results.get_variable("pressure"),
        dtype=float,
    )

    closed_pressure_interpolated = np.interp(
        histories["time_s"],
        closed_time,
        closed_pressure,
    )

    with output_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as file:
        writer = csv.writer(file)

        writer.writerow(
            [
                "time_s",
                "temperature_K",
                "vented_pressure_Pa",
                "closed_pressure_Pa",
                "vent_open",
                "N2_mol",
                "CO2_mol",
                "total_mol",
                "generation_rate_mol_per_s",
                "vent_rate_mol_per_s",
                "N2_vent_kg_per_s",
                "CO2_vent_kg_per_s",
                "total_vent_kg_per_s",
            ]
        )

        for index in range(
            len(histories["time_s"])
        ):
            writer.writerow(
                [
                    histories["time_s"][index],
                    histories["temperature_K"][index],
                    histories["pressure_Pa"][index],
                    closed_pressure_interpolated[index],
                    histories["vent_open"][index],
                    histories["n2_moles"][index],
                    histories["co2_moles"][index],
                    histories["total_moles"][index],
                    histories["total_generation_rate_mol_per_s"][index],
                    histories["total_molar_flow_mol_per_s"][index],
                    histories["n2_mass_flow_kg_per_s"][index],
                    histories["co2_mass_flow_kg_per_s"][index],
                    histories["total_mass_flow_kg_per_s"][index],
                ]
            )

def save_synthetic_vent_results(
    summary: dict[str, float],
    synthetic_gas_balance: dict[str, float],
    synthetic_flow_summary: dict[str, float],
) -> None:
    """Write synthetic venting results as LaTeX macros."""
    output_path = (
        GENERATED_DIR
        / "synthetic_vent_results.tex"
    )

    lines = [
        "% Automatically generated by scripts/generate_report_figures.py.",
        "% Do not edit manually.",
        "",
        (
            "\\newcommand{\\SyntheticVentOpenTime}"
            f"{{{summary['vent_open_time_s']:.4f}}}"
        ),
        (
            "\\newcommand{\\SyntheticVentOpenTemperature}"
            f"{{{summary['vent_open_temperature_K']:.2f}}}"
        ),
        (
            "\\newcommand{\\SyntheticVentPeakPressure}"
            f"{{{summary['peak_pressure_Pa'] / 1000.0:.2f}}}"
        ),
        (
            "\\newcommand{\\SyntheticVentFinalPressure}"
            f"{{{summary['final_pressure_Pa'] / 1000.0:.2f}}}"
        ),
        (
            "\\newcommand{\\SyntheticClosedFinalPressure}"
            f"{{{summary['closed_final_pressure_Pa'] / 1000.0:.2f}}}"
        ),
        (
            "\\newcommand{\\SyntheticVentPeakMassFlow}"
            f"{{{summary['peak_mass_flow_kg_per_s'] * 1.0e6:.3f}}}"
        ),
        (
            "\\newcommand{\\SyntheticVentMaximumMoles}"
            f"{{{summary['maximum_total_moles'] * 1000.0:.4f}}}"
        ),
        (
            "\\newcommand{\\SyntheticVentFinalMoles}"
            f"{{{summary['final_total_moles'] * 1000.0:.4f}}}"
        ),
        (
            "\\newcommand{\\SyntheticGasInitialMoles}"
            f"{{{synthetic_gas_balance['initial_moles'] * 1000.0:.6f}}}"
        ),
        (
            "\\newcommand{\\SyntheticGasGeneratedMoles}"
            f"{{{synthetic_gas_balance['generated_moles'] * 1000.0:.6f}}}"
        ),
        (
            "\\newcommand{\\SyntheticGasVentedMoles}"
            f"{{{synthetic_gas_balance['vented_moles'] * 1000.0:.6f}}}"
        ),
        (
            "\\newcommand{\\SyntheticGasFinalMoles}"
            f"{{{synthetic_gas_balance['final_moles'] * 1000.0:.6f}}}"
        ),
        (
            "\\newcommand{\\SyntheticGasBalanceResidual}"
            f"{{{synthetic_gas_balance['residual_moles'] * 1.0e6:.6f}}}"
        ),
        (
            "\\newcommand{\\SyntheticGasRelativeResidual}"
            f"{{{100.0 * synthetic_gas_balance['relative_residual']:.6f}}}"
        ),
        (
            "\\newcommand{\\SyntheticCriticalPressureRatio}"
            f"{{{synthetic_flow_summary['critical_pressure_ratio']:.6f}}}"
        ),
        (
            "\\newcommand{\\SyntheticCriticalUpstreamPressure}"
            f"{{{synthetic_flow_summary['critical_upstream_pressure_Pa'] / 1000.0:.3f}}}"
        ),
        (
            "\\newcommand{\\SyntheticChokedFlowStart}"
            f"{{{synthetic_flow_summary['choked_start_time_s'] * 1000.0:.6f}}}"
        ),
        (
            "\\newcommand{\\SyntheticChokedFlowEnd}"
            f"{{{synthetic_flow_summary['choked_end_time_s'] * 1000.0:.6f}}}"
        ),
        (
            "\\newcommand{\\SyntheticChokedFlowDuration}"
            f"{{{synthetic_flow_summary['choked_duration_s'] * 1.0e6:.3f}}}"
        ),
        (
            "\\newcommand{\\SyntheticFlowTransitionPressure}"
            f"{{{synthetic_flow_summary['transition_pressure_Pa'] / 1000.0:.3f}}}"
        ),
        (
            "\\newcommand{\\SyntheticPeakCOtwoMoleFraction}"
            f"{{{synthetic_flow_summary['peak_co2_mole_fraction']:.6f}}}"
        ),
        (
            "\\newcommand{\\SyntheticFinalCOtwoMoleFraction}"
            f"{{{synthetic_flow_summary['final_co2_mole_fraction']:.6f}}}"
        ),
        "",
    ]

    output_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

def save_synthetic_gas_balance_latex(
    gas_balance: dict[str, float],
) -> None:
    """Write the synthetic gas molar-balance table."""
    output_path = (
        GENERATED_DIR
        / "synthetic_gas_balance.tex"
    )

    initial_mmol = (
        gas_balance["initial_moles"]
        * 1000.0
    )

    generated_mmol = (
        gas_balance["generated_moles"]
        * 1000.0
    )

    vented_mmol = (
        gas_balance["vented_moles"]
        * 1000.0
    )

    final_mmol = (
        gas_balance["final_moles"]
        * 1000.0
    )

    residual_umol = (
        gas_balance["residual_moles"]
        * 1.0e6
    )

    relative_residual_percent = (
        100.0
        * gas_balance["relative_residual"]
    )

    lines = [
        "% Automatically generated by scripts/generate_report_figures.py.",
        "% Do not edit manually.",
        "",
        r"\begin{table}[htbp]",
        r"    \centering",
        (
            r"    \caption{Global molar balance for the "
            r"synthetic gas-pressure-venting verification problem.}"
        ),
        r"    \label{tab:synthetic_gas_balance}",
        "",
        r"    \begin{tabular}{lr}",
        r"        \toprule",
        r"        Quantity & Amount \\",
        r"        \midrule",
        (
            "        Initial gas inventory"
            f" & {initial_mmol:.6f}"
            r" \si{\milli\mole} \\"
        ),
        (
            "        Generated gas"
            f" & {generated_mmol:.6f}"
            r" \si{\milli\mole} \\"
        ),
        (
            "        Vented gas"
            f" & {vented_mmol:.6f}"
            r" \si{\milli\mole} \\"
        ),
        (
            "        Final gas inventory"
            f" & {final_mmol:.6f}"
            r" \si{\milli\mole} \\"
        ),
        r"        \midrule",
        (
            "        Molar-balance residual"
            f" & {residual_umol:.6f}"
            r" \si{\micro\mole} \\"
        ),
        r"        \bottomrule",
        r"    \end{tabular}",
        "",
        (
            r"    \vspace{0.5em}"
        ),
        "",
        (
            r"    Relative residual: "
            f"{relative_residual_percent:.6f}"
            r"\%."
        ),
        r"\end{table}",
        "",
    ]

    output_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

def build_reference_problem() -> tuple[
    ThermalProblem,
    object,
]:
    """Build the Hu 2020 reference thermal-runaway problem."""
    cell = cell_21700_generic()

    reaction_network = hu2020_reaction_network(cell)

    chemistry_backend = ReactionNetworkBackend(
        reaction_network=reaction_network,
        cell=cell,
    )

    problem = ThermalProblem(
        cell=cell,
        chemistry_backend=chemistry_backend,
        initial_temperature=480.0,
        initial_conversions=hu2020_initial_conversions(),
        ambient_temperature=298.15,
        convection_coefficient=10.0,
        duration=20.0,
        maximum_temperature=1200.0,
    )

    return problem, reaction_network


def calculate_heat_histories(
    reaction_network: object,
    time: np.ndarray,
    temperature: np.ndarray,
    conversions: np.ndarray,
) -> tuple[list[str], dict[str, np.ndarray], np.ndarray]:
    """Reconstruct reaction heat-generation histories from the solution."""
    reaction_names = [
        reaction.name
        for reaction in reaction_network.reactions
    ]

    heat_by_reaction = {
        name: np.zeros_like(time, dtype=float)
        for name in reaction_names
    }

    for time_index, current_temperature in enumerate(temperature):
        current_conversions = [
            float(conversions[reaction_index, time_index])
            for reaction_index in range(conversions.shape[0])
        ]

        current_heat = reaction_network.heat_generation_by_reaction(
            temperature=float(current_temperature),
            conversions=current_conversions,
        )

        for name in reaction_names:
            heat_by_reaction[name][time_index] = current_heat[name]

    total_heat = np.sum(
        np.vstack(
            [
                heat_by_reaction[name]
                for name in reaction_names
            ]
        ),
        axis=0,
    )

    return reaction_names, heat_by_reaction, total_heat

def calculate_reaction_summary(
    time: np.ndarray,
    conversions: np.ndarray,
    reaction_names: list[str],
    heat_by_reaction: dict[str, np.ndarray],
) -> list[dict[str, float | str]]:
    """Calculate integrated energy and conversion data for each reaction."""
    summary: list[dict[str, float | str]] = []

    integrated_energies = [
        float(
            np.trapezoid(
                heat_by_reaction[reaction_name],
                time,
            )
            / 1000.0
        )
        for reaction_name in reaction_names
    ]

    total_energy = sum(integrated_energies)

    for index, (
        reaction_name,
        label,
        energy,
    ) in enumerate(
        zip(
            reaction_names,
            REACTION_LABELS,
            integrated_energies,
            strict=True,
        )
    ):
        if total_energy > 0.0:
            energy_fraction = energy / total_energy
        else:
            energy_fraction = 0.0

        summary.append(
            {
                "reaction_name": reaction_name,
                "label": label,
                "initial_conversion": float(
                    conversions[index, 0]
                ),
                "final_conversion": float(
                    conversions[index, -1]
                ),
                "integrated_energy_kJ_per_kg": energy,
                "energy_fraction": energy_fraction,
            }
        )

    return summary

def calculate_energy_balance(
    problem: ThermalProblem,
    time: np.ndarray,
    temperature: np.ndarray,
    total_heat: np.ndarray,
) -> tuple[dict[str, float], np.ndarray]:
    """Calculate the thermal energy balance of the simulation."""
    thermal_model = LumpedThermalModel(
        cell=problem.cell,
        convection_coefficient=problem.convection_coefficient,
        ambient_temperature=problem.ambient_temperature,
    )

    convective_loss_rate = np.asarray(
        [
            thermal_model.heat_loss(
                temperature=float(current_temperature)
            )
            / problem.cell.mass
            for current_temperature in temperature
        ],
        dtype=float,
    )

    reaction_energy = float(
        np.trapezoid(
            total_heat,
            time,
        )
        / 1000.0
    )

    convective_loss = float(
        np.trapezoid(
            convective_loss_rate,
            time,
        )
        / 1000.0
    )

    sensible_energy = float(
        (
            problem.cell.thermal_capacity
            * (
                float(temperature[-1])
                - float(temperature[0])
            )
            / problem.cell.mass
        )
        / 1000.0
    )

    residual = (
        reaction_energy
        - convective_loss
        - sensible_energy
    )

    if reaction_energy != 0.0:
        relative_residual = (
            residual / reaction_energy
        )
    else:
        relative_residual = 0.0

    energy_balance = {
        "reaction_energy_kJ_per_kg": reaction_energy,
        "sensible_energy_kJ_per_kg": sensible_energy,
        "convective_loss_kJ_per_kg": convective_loss,
        "residual_kJ_per_kg": residual,
        "relative_residual": relative_residual,
    }

    return energy_balance, convective_loss_rate

def save_reaction_summary_csv(
    summary: list[dict[str, float | str]],
) -> None:
    """Save the reaction-wise energy summary."""
    output_path = DATA_DIR / "hu2020_reaction_summary.csv"

    with output_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as file:
        writer = csv.writer(file)

        writer.writerow(
            [
                "reaction",
                "initial_conversion",
                "final_conversion",
                "integrated_energy_kJ_per_kg",
                "energy_fraction",
            ]
        )

        for row in summary:
            writer.writerow(
                [
                    row["label"],
                    row["initial_conversion"],
                    row["final_conversion"],
                    row["integrated_energy_kJ_per_kg"],
                    row["energy_fraction"],
                ]
            )

def save_reaction_summary_latex(
    summary: list[dict[str, float | str]],
) -> None:
    """Write the reaction-wise energy summary as a LaTeX table."""
    output_path = (
        GENERATED_DIR
        / "hu2020_reaction_summary.tex"
    )

    latex_labels = {
        "SEI": "SEI",
        "Anode-electrolyte": "Anode--electrolyte",
        "Cathode": "Cathode",
        "Electrolyte": "Electrolyte",
    }

    lines = [
        "% Automatically generated by scripts/generate_report_figures.py.",
        "% Do not edit manually.",
        "",
        r"\begin{table}[htbp]",
        r"    \centering",
        r"    \caption{Reaction-wise energy release in the Hu et al.\ reference simulation.}",
        r"    \label{tab:hu_reaction_summary}",
        "",
        r"    \begin{tabular}{lrrrr}",
        r"        \toprule",
        r"        Reaction",
        r"        & $\alpha_0$",
        r"        & $\alpha_f$",
        r"        & Energy (\si{\kilo\joule\per\kilogram})",
        r"        & Fraction (\%) \\",
        r"        \midrule",
    ]

    for row in summary:
        label = latex_labels[str(row["label"])]

        initial_conversion = float(
            row["initial_conversion"]
        )
        final_conversion = float(
            row["final_conversion"]
        )
        energy = float(
            row["integrated_energy_kJ_per_kg"]
        )
        fraction_percent = (
            100.0
            * float(row["energy_fraction"])
        )

        lines.append(
            "        "
            f"{label}"
            f" & {initial_conversion:.3f}"
            f" & {final_conversion:.3f}"
            f" & {energy:.2f}"
            f" & {fraction_percent:.1f}"
            r" \\"
        )

    total_energy = sum(
        float(
            row["integrated_energy_kJ_per_kg"]
        )
        for row in summary
    )

    lines.extend(
        [
            r"        \midrule",
            (
                r"        \textbf{Total}"
                r" & --"
                r" & --"
                f" & {total_energy:.2f}"
                r" & 100.0 \\"
            ),
            r"        \bottomrule",
            r"    \end{tabular}",
            r"\end{table}",
            "",
        ]
    )

    output_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def save_temperature_figure(
    time: np.ndarray,
    temperature: np.ndarray,
) -> None:
    """Save the reference temperature-history figure."""
    figure, axis = plt.subplots(figsize=(6.4, 4.0))

    axis.plot(
        time,
        temperature,
        linewidth=1.8,
    )

    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Temperature (K)")
    axis.grid(True, alpha=0.25)

    figure.tight_layout()

    figure.savefig(
        FIGURE_DIR / "hu2020_temperature.pdf",
        bbox_inches="tight",
    )

    plt.close(figure)


def save_conversion_figure(
    time: np.ndarray,
    conversions: np.ndarray,
) -> None:
    """Save the reaction-conversion histories."""
    figure, axis = plt.subplots(figsize=(6.4, 4.2))

    for reaction_index, label in enumerate(REACTION_LABELS):
        axis.plot(
            time,
            conversions[reaction_index],
            label=label,
            linewidth=1.6,
        )

    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Reaction conversion")
    axis.set_ylim(-0.02, 1.02)
    axis.grid(True, alpha=0.25)
    axis.legend()

    figure.tight_layout()

    figure.savefig(
        FIGURE_DIR / "hu2020_conversion.pdf",
        bbox_inches="tight",
    )

    plt.close(figure)


def save_heat_generation_figure(
    time: np.ndarray,
    reaction_names: list[str],
    heat_by_reaction: dict[str, np.ndarray],
    total_heat: np.ndarray,
) -> None:
    """Save the reaction heat-generation histories."""
    figure, axis = plt.subplots(figsize=(6.4, 4.2))

    for reaction_index, reaction_name in enumerate(reaction_names):
        axis.plot(
            time,
            heat_by_reaction[reaction_name],
            label=REACTION_LABELS[reaction_index],
            linewidth=1.4,
        )

    axis.plot(
        time,
        total_heat,
        label="Total",
        linewidth=2.0,
        linestyle="--",
    )

    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Heat generation (W/kg)")
    axis.grid(True, alpha=0.25)
    axis.legend()

    figure.tight_layout()

    figure.savefig(
        FIGURE_DIR / "hu2020_heat_generation.pdf",
        bbox_inches="tight",
    )

    plt.close(figure)


def save_csv(
    time: np.ndarray,
    temperature: np.ndarray,
    conversions: np.ndarray,
    reaction_names: list[str],
    heat_by_reaction: dict[str, np.ndarray],
    total_heat: np.ndarray,
) -> None:
    """Save numerical data used to construct the report figures."""
    output_path = DATA_DIR / "hu2020_reference.csv"

    with output_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as file:
        writer = csv.writer(file)

        writer.writerow(
            [
                "time_s",
                "temperature_K",
                "alpha_sei",
                "alpha_anode_electrolyte",
                "alpha_cathode",
                "alpha_electrolyte",
                "qdot_sei_W_per_kg",
                "qdot_anode_electrolyte_W_per_kg",
                "qdot_cathode_W_per_kg",
                "qdot_electrolyte_W_per_kg",
                "qdot_total_W_per_kg",
            ]
        )

        for time_index in range(len(time)):
            writer.writerow(
                [
                    time[time_index],
                    temperature[time_index],
                    conversions[0, time_index],
                    conversions[1, time_index],
                    conversions[2, time_index],
                    conversions[3, time_index],
                    heat_by_reaction[reaction_names[0]][time_index],
                    heat_by_reaction[reaction_names[1]][time_index],
                    heat_by_reaction[reaction_names[2]][time_index],
                    heat_by_reaction[reaction_names[3]][time_index],
                    total_heat[time_index],
                ]
            )


def save_latex_results(
    final_temperature: float,
    integrated_energy_kj_per_kg: float,
    reaction_summary: list[dict[str, float | str]],
    energy_balance: dict[str, float],
) -> None:
    """Write numerical result macros consumed by the LaTeX report."""
    output_path = GENERATED_DIR / "hu2020_results.tex"

    contents = [
        "% Automatically generated by scripts/generate_report_figures.py.",
        "% Do not edit manually.",
        "",
        (
            "\\newcommand{\\HuFinalTemperature}"
            f"{{{final_temperature:.1f}}}"
        ),
        (
            "\\newcommand{\\HuIntegratedEnergy}"
            f"{{{integrated_energy_kj_per_kg:.1f}}}"
        ),
    ]

    macro_names = {
        "SEI": "Sei",
        "Anode-electrolyte": "Anode",
        "Cathode": "Cathode",
        "Electrolyte": "Electrolyte",
    }

    for row in reaction_summary:
        label = str(row["label"])
        macro_name = macro_names[label]

        energy = float(
            row["integrated_energy_kJ_per_kg"]
        )
        fraction = (
            100.0
            * float(row["energy_fraction"])
        )
        final_conversion = float(
            row["final_conversion"]
        )

        contents.extend(
            [
                (
                    f"\\newcommand{{\\Hu{macro_name}Energy}}"
                    f"{{{energy:.1f}}}"
                ),
                (
                    f"\\newcommand{{\\Hu{macro_name}Fraction}}"
                    f"{{{fraction:.1f}}}"
                ),
                (
                    f"\\newcommand{{\\Hu{macro_name}FinalConversion}}"
                    f"{{{final_conversion:.3f}}}"
                ),
            ]
        )

    contents.extend(
        [
            (
                "\\newcommand{\\HuSensibleEnergy}"
                f"{{{energy_balance['sensible_energy_kJ_per_kg']:.1f}}}"
            ),
            (
                "\\newcommand{\\HuConvectiveLoss}"
                f"{{{energy_balance['convective_loss_kJ_per_kg']:.1f}}}"
            ),
            (
                "\\newcommand{\\HuEnergyResidual}"
                f"{{{energy_balance['residual_kJ_per_kg']:.4f}}}"
            ),
            (
                "\\newcommand{\\HuRelativeEnergyResidual}"
                f"{{{100.0 * energy_balance['relative_residual']:.4f}}}"
            ),
        ]
    )

    contents.append("")

    output_path.write_text(
        "\n".join(contents),
        encoding="utf-8",
    )

def save_thermal_power_balance_figure(
    time: np.ndarray,
    energy_histories: dict[str, np.ndarray],
) -> None:
    """Save the instantaneous thermal-power balance."""
    figure, axis = plt.subplots(
        figsize=(6.4, 4.2)
    )

    axis.plot(
        time,
        energy_histories[
            "reaction_rate_W_per_kg"
        ]
        / 1000.0,
        label="Reaction heat generation",
        linewidth=1.8,
    )

    axis.plot(
        time,
        energy_histories[
            "sensible_heating_rate_W_per_kg"
        ]
        / 1000.0,
        label="Sensible heating",
        linewidth=1.6,
        linestyle="--",
    )

    axis.plot(
        time,
        energy_histories[
            "convective_loss_rate_W_per_kg"
        ]
        / 1000.0,
        label="Convective heat loss",
        linewidth=1.4,
        linestyle=":",
    )

    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Specific thermal power (kW/kg)")
    axis.grid(True, alpha=0.25)
    axis.legend()

    figure.tight_layout()

    figure.savefig(
        FIGURE_DIR
        / "hu2020_thermal_power_balance.pdf",
        bbox_inches="tight",
    )

    plt.close(figure)

def save_cumulative_energy_balance_figure(
    time: np.ndarray,
    energy_histories: dict[str, np.ndarray],
) -> None:
    """Save the cumulative thermal-energy balance."""
    figure, axis = plt.subplots(
        figsize=(6.4, 4.2)
    )

    axis.plot(
        time,
        energy_histories[
            "reaction_energy_kJ_per_kg"
        ],
        label="Reaction energy",
        linewidth=2.0,
    )

    axis.plot(
        time,
        energy_histories[
            "accounted_energy_kJ_per_kg"
        ],
        label="Sensible energy + convective loss",
        linewidth=1.7,
        linestyle="--",
    )

    axis.plot(
        time,
        energy_histories[
            "convective_loss_kJ_per_kg"
        ],
        label="Convective loss",
        linewidth=1.4,
        linestyle=":",
    )

    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Specific energy (kJ/kg)")
    axis.grid(True, alpha=0.25)
    axis.legend()

    figure.tight_layout()

    figure.savefig(
        FIGURE_DIR
        / "hu2020_energy_balance_cumulative.pdf",
        bbox_inches="tight",
    )

    plt.close(figure)

def save_energy_balance_csv(
    time: np.ndarray,
    energy_histories: dict[str, np.ndarray],
) -> None:
    """Save instantaneous and cumulative thermal-energy data."""
    output_path = (
        DATA_DIR
        / "hu2020_energy_balance.csv"
    )

    with output_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as file:
        writer = csv.writer(file)

        writer.writerow(
            [
                "time_s",
                "reaction_rate_W_per_kg",
                "convective_loss_rate_W_per_kg",
                "sensible_heating_rate_W_per_kg",
                "reaction_energy_kJ_per_kg",
                "convective_loss_kJ_per_kg",
                "sensible_energy_kJ_per_kg",
                "accounted_energy_kJ_per_kg",
                "residual_kJ_per_kg",
            ]
        )

        for index in range(len(time)):
            writer.writerow(
                [
                    time[index],
                    energy_histories[
                        "reaction_rate_W_per_kg"
                    ][index],
                    energy_histories[
                        "convective_loss_rate_W_per_kg"
                    ][index],
                    energy_histories[
                        "sensible_heating_rate_W_per_kg"
                    ][index],
                    energy_histories[
                        "reaction_energy_kJ_per_kg"
                    ][index],
                    energy_histories[
                        "convective_loss_kJ_per_kg"
                    ][index],
                    energy_histories[
                        "sensible_energy_kJ_per_kg"
                    ][index],
                    energy_histories[
                        "accounted_energy_kJ_per_kg"
                    ][index],
                    energy_histories[
                        "residual_kJ_per_kg"
                    ][index],
                ]
            )

def calculate_energy_histories(
    problem: ThermalProblem,
    time: np.ndarray,
    temperature: np.ndarray,
    total_heat: np.ndarray,
    convective_loss_rate: np.ndarray,
) -> dict[str, np.ndarray]:
    """Calculate instantaneous and cumulative thermal-energy histories."""
    sensible_heating_rate = (
        total_heat - convective_loss_rate
    )

    cumulative_reaction_energy = (
        cumulative_trapezoid(
            total_heat,
            time,
            initial=0.0,
        )
        / 1000.0
    )

    cumulative_convective_loss = (
        cumulative_trapezoid(
            convective_loss_rate,
            time,
            initial=0.0,
        )
        / 1000.0
    )

    specific_heat_capacity = (
        problem.cell.thermal_capacity
        / problem.cell.mass
    )

    sensible_energy = (
        specific_heat_capacity
        * (
            temperature
            - temperature[0]
        )
        / 1000.0
    )

    accounted_energy = (
        sensible_energy
        + cumulative_convective_loss
    )

    residual = (
        cumulative_reaction_energy
        - accounted_energy
    )

    return {
        "reaction_rate_W_per_kg": total_heat,
        "convective_loss_rate_W_per_kg": convective_loss_rate,
        "sensible_heating_rate_W_per_kg": sensible_heating_rate,
        "reaction_energy_kJ_per_kg": cumulative_reaction_energy,
        "convective_loss_kJ_per_kg": cumulative_convective_loss,
        "sensible_energy_kJ_per_kg": sensible_energy,
        "accounted_energy_kJ_per_kg": accounted_energy,
        "residual_kJ_per_kg": residual,
    }

def save_energy_balance_latex(
    energy_balance: dict[str, float],
) -> None:
    """Write the thermal energy-balance table for the report."""
    output_path = (
        GENERATED_DIR
        / "hu2020_energy_balance.tex"
    )

    reaction_energy = energy_balance[
        "reaction_energy_kJ_per_kg"
    ]
    sensible_energy = energy_balance[
        "sensible_energy_kJ_per_kg"
    ]
    convective_loss = energy_balance[
        "convective_loss_kJ_per_kg"
    ]
    residual = energy_balance[
        "residual_kJ_per_kg"
    ]
    relative_residual_percent = (
        100.0
        * energy_balance["relative_residual"]
    )

    lines = [
        "% Automatically generated by scripts/generate_report_figures.py.",
        "% Do not edit manually.",
        "",
        r"\begin{table}[htbp]",
        r"    \centering",
        (
            r"    \caption{Thermal energy-balance closure "
            r"for the Hu et al.\ reference simulation.}"
        ),
        r"    \label{tab:hu_energy_balance}",
        "",
        r"    \begin{tabular}{lr}",
        r"        \toprule",
        (
            r"        Quantity"
            r" & Energy (\si{\kilo\joule\per\kilogram}) \\"
        ),
        r"        \midrule",
        (
            "        Reaction energy"
            f" & {reaction_energy:.3f}"
            r" \\"
        ),
        (
            "        Sensible-energy increase"
            f" & {sensible_energy:.3f}"
            r" \\"
        ),
        (
            "        Convective heat loss"
            f" & {convective_loss:.3f}"
            r" \\"
        ),
        r"        \midrule",
        (
            "        Energy-balance residual"
            f" & {residual:.6f}"
            r" \\"
        ),
        r"        \bottomrule",
        r"    \end{tabular}",
        "",
        (
            r"    \vspace{0.5em}"
        ),
        "",
        (
            r"    Relative residual: "
            f"{relative_residual_percent:.4f}"
            r"\%."
        ),
        r"\end{table}",
        "",
    ]

    output_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

def save_synthetic_vent_composition_figure(
    histories: dict[str, np.ndarray],
    vent_open_time: float,
) -> None:
    """Save the internal gas-composition history."""
    figure, axis = plt.subplots(
        figsize=(6.4, 4.2)
    )

    axis.plot(
        histories["time_s"],
        histories["n2_mole_fraction"],
        label="N2",
        linewidth=1.8,
    )

    axis.plot(
        histories["time_s"],
        histories["co2_mole_fraction"],
        label="CO2",
        linewidth=1.8,
    )

    axis.axvline(
        vent_open_time,
        label="Vent opening",
        linestyle="--",
    )

    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Gas mole fraction (-)")
    axis.set_ylim(0.0, 1.0)
    axis.grid(True, alpha=0.25)
    axis.legend()

    figure.tight_layout()

    figure.savefig(
        FIGURE_DIR
        / "synthetic_vent_composition.pdf",
        bbox_inches="tight",
    )

    plt.close(figure)

def main() -> None:
    """Run the Hu 2020 reference simulation and generate report outputs."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    GENERATED_DIR.mkdir(parents=True, exist_ok=True)

    problem, reaction_network = build_reference_problem()

    results = ScipySolver().solve(problem)

    time = np.asarray(
        results.time,
        dtype=float,
    )

    temperature = np.asarray(
        results.temperature,
        dtype=float,
    )

    conversions = np.vstack(
        [
            np.asarray(
                results.get_variable(f"conversion_{reaction_index}"),
                dtype=float,
            )
            for reaction_index in range(
                len(reaction_network.reactions)
            )
        ]
    )

    reaction_names, heat_by_reaction, total_heat = (
        calculate_heat_histories(
            reaction_network=reaction_network,
            time=time,
            temperature=temperature,
            conversions=conversions,
        )
    )

    reaction_summary = calculate_reaction_summary(
        time=time,
        conversions=conversions,
        reaction_names=reaction_names,
        heat_by_reaction=heat_by_reaction,
    )

    energy_balance, convective_loss_rate = calculate_energy_balance(
        problem=problem,
        time=time,
        temperature=temperature,
        total_heat=total_heat,
    )

    energy_histories = calculate_energy_histories(
        problem=problem,
        time=time,
        temperature=temperature,
        total_heat=total_heat,
        convective_loss_rate=convective_loss_rate,
    )

    save_energy_balance_latex(
        energy_balance=energy_balance,
    )

    integrated_energy_kj_per_kg = (
        np.trapezoid(
            total_heat,
            time,
        )
        / 1000.0
    )

    final_temperature = float(
        temperature[-1]
    )

    save_temperature_figure(
        time=time,
        temperature=temperature,
    )

    save_conversion_figure(
        time=time,
        conversions=conversions,
    )

    save_heat_generation_figure(
        time=time,
        reaction_names=reaction_names,
        heat_by_reaction=heat_by_reaction,
        total_heat=total_heat,
    )

    save_csv(
        time=time,
        temperature=temperature,
        conversions=conversions,
        reaction_names=reaction_names,
        heat_by_reaction=heat_by_reaction,
        total_heat=total_heat,
    )

    save_reaction_summary_csv(
        summary=reaction_summary,
    )

    save_reaction_summary_latex(
        summary=reaction_summary,
    )

    save_latex_results(
        final_temperature=final_temperature,
        integrated_energy_kj_per_kg=integrated_energy_kj_per_kg,
        reaction_summary=reaction_summary,
        energy_balance=energy_balance,
    )

    save_thermal_power_balance_figure(
        time=time,
        energy_histories=energy_histories,
    )

    save_cumulative_energy_balance_figure(
        time=time,
        energy_histories=energy_histories,
    )

    save_energy_balance_csv(
        time=time,
        energy_histories=energy_histories,
    )

#    save_synthetic_vent_composition_figure(
#        histories=synthetic_histories,
#        vent_open_time=vent_open_time,
#    )













    print("Generated report data and figures:")
    print(
        "  report/figures/hu2020_temperature.pdf"
    )
    print(
        "  report/figures/hu2020_conversion.pdf"
    )
    print(
        "  report/figures/hu2020_heat_generation.pdf"
    )
    print(
        "  report/data/hu2020_reference.csv"
    )
    print(
        "  report/generated/hu2020_results.tex"
    )
    print(
        "  report/data/hu2020_reaction_summary.csv"
    )
    print(
        "  report/generated/hu2020_reaction_summary.tex"
    )
    print("  report/generated/hu2020_energy_balance.tex")
    print("  report/figures/hu2020_thermal_power_balance.pdf")
    print("  report/figures/hu2020_energy_balance_cumulative.pdf")
    print("  report/data/hu2020_energy_balance.csv")
    print("  report/figures/synthetic_vent_composition.pdf")
    print()
    print(
        f"Final temperature: "
        f"{final_temperature:.3f} K"
    )
    print(
        f"Integrated reaction energy: "
        f"{integrated_energy_kj_per_kg:.3f} kJ/kg"
    )

    print()
    print("Reaction energy summary:")
    print(
        f"{'Reaction':<22}"
        f"{'alpha_0':>10}"
        f"{'alpha_f':>10}"
        f"{'Energy (kJ/kg)':>18}"
        f"{'Fraction (%)':>15}"
    )

    for row in reaction_summary:
        print(
            f"{row['label']!s:<22}"
            f"{float(row['initial_conversion']):>10.3f}"
            f"{float(row['final_conversion']):>10.3f}"
            f"{float(row['integrated_energy_kJ_per_kg']):>18.3f}"
            f"{100.0 * float(row['energy_fraction']):>15.2f}"
        )

    print()
    print("Thermal energy balance:")
    print(
        "Reaction energy:       "
        f"{energy_balance['reaction_energy_kJ_per_kg']:.6f} kJ/kg"
    )
    print(
        "Sensible energy:       "
        f"{energy_balance['sensible_energy_kJ_per_kg']:.6f} kJ/kg"
    )
    print(
        "Convective heat loss:  "
        f"{energy_balance['convective_loss_kJ_per_kg']:.6f} kJ/kg"
    )
    print(f"Energy residual:       {energy_balance['residual_kJ_per_kg']:.6f} kJ/kg")
    print(f"Relative residual:     {100.0 * energy_balance['relative_residual']:.6f} %")

    closed_vent_problem = build_synthetic_vent_problem(
        with_vent=False,
    )

    vented_problem = build_synthetic_vent_problem(
        with_vent=True,
    )

    closed_vent_results = ScipySolver().solve(closed_vent_problem)

    vented_results = ScipySolver().solve(vented_problem)

    synthetic_histories = calculate_synthetic_vent_histories(
        problem=vented_problem,
        results=vented_results,
    )

    synthetic_summary = calculate_synthetic_vent_summary(
        histories=synthetic_histories,
        closed_results=closed_vent_results,
    )

    synthetic_gas_balance = calculate_synthetic_gas_balance(
        problem=vented_problem,
        histories=synthetic_histories,
    )

    synthetic_flow_summary = calculate_synthetic_flow_regime_summary(
        problem=vented_problem,
        histories=synthetic_histories,
    )

    
    vent_open_time = synthetic_summary["vent_open_time_s"]

    vent_open_pressure = vented_problem.vent_open_pressure

    if vent_open_pressure is None:
        raise RuntimeError("Synthetic vent opening pressure is undefined.")

    save_synthetic_vent_temperature_figure(
        histories=synthetic_histories,
        vent_open_time=vent_open_time,
    )

    save_synthetic_vent_pressure_figure(
        histories=synthetic_histories,
        vent_open_time=vent_open_time,
        vent_open_pressure=vent_open_pressure,
    )

    save_synthetic_vent_gas_inventory_figure(
        histories=synthetic_histories,
        vent_open_time=vent_open_time,
    )

    save_synthetic_vent_flow_figure(
        histories=synthetic_histories,
        vent_open_time=vent_open_time,
        choked_end_time=synthetic_flow_summary["choked_end_time_s"],
    )

    save_synthetic_vent_csv(
        histories=synthetic_histories,
        closed_results=closed_vent_results,
    )

    save_synthetic_vent_results(
        summary=synthetic_summary,
        synthetic_gas_balance=synthetic_gas_balance,
        synthetic_flow_summary=synthetic_flow_summary,
    )

    save_synthetic_gas_balance_latex(
        gas_balance=synthetic_gas_balance,
    )
    
    save_synthetic_vent_composition_figure(
        histories=synthetic_histories,
        vent_open_time=vent_open_time,
    )

    print()
    print("Synthetic venting summary:")
    print(f"Vent opening time:      {synthetic_summary['vent_open_time_s']:.6f} s")
    print(
        f"Temperature at opening: {synthetic_summary['vent_open_temperature_K']:.3f} K"
    )
    print(
        "Peak pressure:          "
        f"{synthetic_summary['peak_pressure_Pa'] / 1000.0:.3f} kPa"
    )
    print(
        "Final vented pressure:  "
        f"{synthetic_summary['final_pressure_Pa'] / 1000.0:.3f} kPa"
    )
    print(
        "Final closed pressure:  "
        f"{synthetic_summary['closed_final_pressure_Pa'] / 1000.0:.3f} kPa"
    )
    print(
        "Peak vent mass flow:    "
        f"{synthetic_summary['peak_mass_flow_kg_per_s'] * 1.0e6:.3f} mg/s"
    )
    print(
        "Maximum gas inventory:  "
        f"{synthetic_summary['maximum_total_moles'] * 1000.0:.6f} mmol"
    )
    print(
        "Final gas inventory:    "
        f"{synthetic_summary['final_total_moles'] * 1000.0:.6f} mmol"
    )
    print()
    print("Synthetic gas molar balance:")
    print(
        "Initial gas:           "
        f"{synthetic_gas_balance['initial_moles'] * 1000.0:.6f} mmol"
    )
    print(
        "Generated gas:         "
        f"{synthetic_gas_balance['generated_moles'] * 1000.0:.6f} mmol"
    )
    print(
        "Vented gas:            "
        f"{synthetic_gas_balance['vented_moles'] * 1000.0:.6f} mmol"
    )
    print(
        "Final gas:             "
        f"{synthetic_gas_balance['final_moles'] * 1000.0:.6f} mmol"
    )
    print(
        "Molar balance residual:"
        f" {synthetic_gas_balance['residual_moles'] * 1.0e6:.6f} umol"
    )
    print(
        "Relative residual:     "
        f"{100.0 * synthetic_gas_balance['relative_residual']:.6f} %"
    )
    print("  report/figures/synthetic_vent_temperature.pdf")
    print("  report/figures/synthetic_vent_pressure.pdf")
    print("  report/figures/synthetic_vent_gas_inventory.pdf")
    print("  report/figures/synthetic_vent_flow.pdf")
    print("  report/data/synthetic_vent_reference.csv")
    print("  report/generated/synthetic_vent_results.tex")
    print("  report/generated/synthetic_gas_balance.tex")

    print()
    print("Synthetic vent flow regime:")
    print(
        "Critical pressure ratio: "
        f"{synthetic_flow_summary['critical_pressure_ratio']:.6f}"
    )
    print(
        "Critical upstream pressure: "
        f"{synthetic_flow_summary['critical_upstream_pressure_Pa'] / 1000.0:.3f} kPa"
    )
    print(
        "Choked-flow start:       "
        f"{synthetic_flow_summary['choked_start_time_s']:.6f} s"
    )
    print(
        f"Choked-flow end:         {synthetic_flow_summary['choked_end_time_s']:.6f} s"
    )
    print(
        "Choked-flow duration:    "
        f"{synthetic_flow_summary['choked_duration_s'] * 1.0e6:.3f} us"
    )
    print(
        f"Unchoked transition:     {synthetic_flow_summary['transition_time_s']:.6f} s"
    )
    print(
        "Transition pressure:     "
        f"{synthetic_flow_summary['transition_pressure_Pa'] / 1000.0:.3f} kPa"
    )
    print(
        "Peak CO2 mole fraction:  "
        f"{synthetic_flow_summary['peak_co2_mole_fraction']:.6f}"
    )
    print(
        "Final CO2 mole fraction: "
        f"{synthetic_flow_summary['final_co2_mole_fraction']:.6f}"
    )




if __name__ == "__main__":
    main()