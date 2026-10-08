"""Calibration tools for reduced-order internal release models."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares


@dataclass(slots=True, frozen=True)
class FirstOrderReleaseFit:
    """Result of a first-order internal-release fit."""

    time_constant: float
    rmse: float
    maximum_absolute_error: float


def fit_first_order_release_time(
    time: np.ndarray,
    released_fraction: np.ndarray,
) -> FirstOrderReleaseFit:
    """Fit a first-order release time constant.

    The fitted reduced-order model is

        F(t) = 1 - exp(-t / tau)

    where F is the cumulative released fraction and tau is the
    characteristic internal release time.

    Parameters
    ----------
    time : numpy.ndarray
        Time values in seconds.
    released_fraction : numpy.ndarray
        Cumulative released fraction, between zero and one.

    Returns
    -------
    FirstOrderReleaseFit
        Fitted time constant and residual diagnostics.
    """
    time_values = np.asarray(
        time,
        dtype=float,
    )

    fraction_values = np.asarray(
        released_fraction,
        dtype=float,
    )

    if time_values.ndim != 1:
        raise ValueError(
            "Time values must be one-dimensional."
        )

    if fraction_values.ndim != 1:
        raise ValueError(
            "Released fractions must be one-dimensional."
        )

    if time_values.size != fraction_values.size:
        raise ValueError(
            "Time and released-fraction arrays must "
            "have the same length."
        )

    if time_values.size < 3:
        raise ValueError(
            "At least three samples are required."
        )

    if not np.all(
        np.isfinite(time_values)
    ):
        raise ValueError(
            "Time values must be finite."
        )

    if not np.all(
        np.isfinite(fraction_values)
    ):
        raise ValueError(
            "Released fractions must be finite."
        )

    if np.any(
        np.diff(time_values) <= 0.0
    ):
        raise ValueError(
            "Time values must be strictly increasing."
        )

    if np.any(
        fraction_values < 0.0
    ) or np.any(
        fraction_values > 1.0
    ):
        raise ValueError(
            "Released fractions must be between zero and one."
        )

    if abs(
        float(fraction_values[0])
    ) > 1.0e-6:
        raise ValueError(
            "Released fraction must begin at zero."
        )

    if np.any(
        np.diff(fraction_values)
        < -1.0e-10
    ):
        raise ValueError(
            "Released fraction must be nondecreasing."
        )

    elapsed_time = (
        time_values
        - time_values[0]
    )

    characteristic_index = int(
        np.argmin(
            np.abs(
                fraction_values
                - (
                    1.0
                    - np.exp(-1.0)
                )
            )
        )
    )

    initial_time_constant = max(
        float(
            elapsed_time[
                characteristic_index
            ]
        ),
        float(
            elapsed_time[-1]
        )
        / 10.0,
        1.0e-12,
    )

    def residual(
        parameters: np.ndarray,
    ) -> np.ndarray:
        time_constant = float(
            parameters[0]
        )

        predicted_fraction = (
            1.0
            - np.exp(
                -elapsed_time
                / time_constant
            )
        )

        return (
            predicted_fraction
            - fraction_values
        )

    fit = least_squares(
        residual,
        x0=np.asarray(
            [
                initial_time_constant,
            ],
            dtype=float,
        ),
        bounds=(
            np.asarray(
                [
                    1.0e-12,
                ]
            ),
            np.asarray(
                [
                    np.inf,
                ]
            ),
        ),
    )

    time_constant = float(
        fit.x[0]
    )

    predicted_fraction = (
        1.0
        - np.exp(
            -elapsed_time
            / time_constant
        )
    )

    errors = (
        predicted_fraction
        - fraction_values
    )

    rmse = float(
        np.sqrt(
            np.mean(
                errors**2
            )
        )
    )

    maximum_absolute_error = float(
        np.max(
            np.abs(
                errors
            )
        )
    )

    return FirstOrderReleaseFit(
        time_constant=time_constant,
        rmse=rmse,
        maximum_absolute_error=(
            maximum_absolute_error
        ),
    )