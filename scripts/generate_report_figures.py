"""Generate reproducible numerical figures for the LiionTR report."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from liiontr.chemistry.reaction_backend import ReactionNetworkBackend
from liiontr.library.cells import cell_21700_generic
from liiontr.library.hu2020 import (
    hu2020_initial_conversions,
    hu2020_reaction_network,
)
from liiontr.problems.thermal import ThermalProblem
from liiontr.solver.scipy_solver import ScipySolver


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
) -> None:
    """Write numerical result macros consumed by the LaTeX report."""
    output_path = GENERATED_DIR / "hu2020_results.tex"

    contents = (
        "% Automatically generated by "
        "scripts/generate_report_figures.py.\n"
        "% Do not edit manually.\n\n"
        f"\\newcommand{{\\HuFinalTemperature}}"
        f"{{{final_temperature:.1f}}}\n"
        f"\\newcommand{{\\HuIntegratedEnergy}}"
        f"{{{integrated_energy_kj_per_kg:.1f}}}\n"
    )

    output_path.write_text(
        contents,
        encoding="utf-8",
    )


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

    save_latex_results(
        final_temperature=final_temperature,
        integrated_energy_kj_per_kg=integrated_energy_kj_per_kg,
    )

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
    print()
    print(
        f"Final temperature: "
        f"{final_temperature:.3f} K"
    )
    print(
        f"Integrated reaction energy: "
        f"{integrated_energy_kj_per_kg:.3f} kJ/kg"
    )


if __name__ == "__main__":
    main()