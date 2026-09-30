"""Sweep the elemental release time constant."""

from __future__ import annotations

import csv
import os
import re
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt

TAU_VALUES = [
    0.1,
    0.3,
    1.0,
    3.0,
    10.0,
]

PROJECT_ROOT = Path(__file__).resolve().parents[1]

DIAGNOSTIC_SCRIPT = (
    PROJECT_ROOT
    / "scripts"
    / "check_hu2020_howard2025_cantera_vent.py"
)

OUTPUT_DIRECTORY = (
    PROJECT_ROOT
    / "results"
    / "element_release_sweep"
)


def extract_value(
    pattern: str,
    output: str,
) -> float:
    """Extract one floating-point value from diagnostic output."""
    match = re.search(
        pattern,
        output,
    )

    if match is None:
        raise RuntimeError(
            f"Could not extract value using pattern: {pattern}"
        )

    return float(
        match.group(1)
    )


def run_case(
    time_constant: float,
) -> dict[str, float]:
    """Run one elemental-release diagnostic case."""
    environment = os.environ.copy()

    environment[
        "LIIONTR_RELEASE_TIME_CONSTANT"
    ] = str(time_constant)

    completed_process = subprocess.run(
        [
            sys.executable,
            str(DIAGNOSTIC_SCRIPT),
        ],
        cwd=PROJECT_ROOT,
        env=environment,
        check=True,
        capture_output=True,
        text=True,
    )

    output = completed_process.stdout

    maximum_temperature = extract_value(
        r"Maximum temperature:\s+"
        r"([0-9.eE+-]+)",
        output,
    )

    maximum_pressure = extract_value(
        r"Maximum pressure:\s+"
        r"([0-9.eE+-]+)",
        output,
    )

    maximum_pressure_time = extract_value(
        r"Maximum pressure time:\s+"
        r"([0-9.eE+-]+)",
        output,
    )

    maximum_release = extract_value(
        r"Maximum release:\s+"
        r"([0-9.eE+-]+)",
        output,
    )

    maximum_pending_carbon = extract_value(
        r"Maximum pending C:\s+"
        r"([0-9.eE+-]+)",
        output,
    )

    vent_open_time = extract_value(
        r"Vent opening time:\s+"
        r"([0-9.eE+-]+)",
        output,
    )

    final_pending_carbon = extract_value(
        r"Final pending elemental inventory"
        r"[\s\S]*?"
        r"\nC\s+([0-9.eE+-]+)",
        output,
    )

    return {
        "tau_s": time_constant,
        "maximum_pressure_bar": (
            maximum_pressure
        ),
        "maximum_pressure_time_s": (
            maximum_pressure_time
        ),
        "maximum_temperature_K": (
            maximum_temperature
        ),
        "vent_open_time_s": (
            vent_open_time
        ),
        "maximum_release_C_mol_s": (
            maximum_release
        ),
        "maximum_pending_C_mol": (
            maximum_pending_carbon
        ),
        "final_pending_C_mol": (
            final_pending_carbon
        ),
    }


def write_csv(
    rows: list[dict[str, float]],
) -> Path:
    """Write sweep results to CSV."""
    output_path = (
        OUTPUT_DIRECTORY
        / "element_release_sweep.csv"
    )

    field_names = list(
        rows[0].keys()
    )

    with output_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as file_handle:
        writer = csv.DictWriter(
            file_handle,
            fieldnames=field_names,
        )

        writer.writeheader()
        writer.writerows(
            rows
        )

    return output_path


def plot_variable(
    rows: list[dict[str, float]],
    variable_name: str,
    ylabel: str,
    filename: str,
) -> None:
    """Plot one sweep variable against release time constant."""
    tau_values = [
        row["tau_s"]
        for row in rows
    ]

    values = [
        row[variable_name]
        for row in rows
    ]

    figure, axis = plt.subplots()

    axis.plot(
        tau_values,
        values,
        marker="o",
    )

    axis.set_xscale(
        "log"
    )

    axis.set_xlabel(
        "Release time constant [s]"
    )

    axis.set_ylabel(
        ylabel
    )

    axis.grid(
        True,
        which="both",
    )

    figure.tight_layout()

    figure.savefig(
        OUTPUT_DIRECTORY
        / filename,
        dpi=200,
    )

    plt.close(
        figure
    )


def main() -> None:
    """Run the complete release-time sensitivity sweep."""
    OUTPUT_DIRECTORY.mkdir(
        parents=True,
        exist_ok=True,
    )

    rows: list[dict[str, float]] = []

    print(
        "Element release time sensitivity"
    )

    print(
        "================================"
    )

    for time_constant in TAU_VALUES:
        print(
            f"Running tau = "
            f"{time_constant:g} s"
        )

        row = run_case(
            time_constant
        )

        rows.append(
            row
        )

        print(
            f"  Pmax = "
            f"{row['maximum_pressure_bar']:.6f} bar"
        )

        print(
            f"  Tmax = "
            f"{row['maximum_temperature_K']:.6f} K"
        )

        print(
            f"  tvent = "
            f"{row['vent_open_time_s']:.6f} s"
        )

    csv_path = write_csv(
        rows
    )

    plot_variable(
        rows=rows,
        variable_name=(
            "maximum_pressure_bar"
        ),
        ylabel="Maximum pressure [bar]",
        filename="maximum_pressure.png",
    )

    plot_variable(
        rows=rows,
        variable_name=(
            "maximum_temperature_K"
        ),
        ylabel="Maximum temperature [K]",
        filename="maximum_temperature.png",
    )

    plot_variable(
        rows=rows,
        variable_name=(
            "vent_open_time_s"
        ),
        ylabel="Vent opening time [s]",
        filename="vent_open_time.png",
    )

    plot_variable(
        rows=rows,
        variable_name=(
            "maximum_release_C_mol_s"
        ),
        ylabel=(
            "Maximum carbon release rate "
            "[mol-atoms/s]"
        ),
        filename="maximum_carbon_release.png",
    )

    print()

    print(
        f"CSV written to: "
        f"{csv_path}"
    )

    print(
        f"Figures written to: "
        f"{OUTPUT_DIRECTORY}"
    )


if __name__ == "__main__":
    main()