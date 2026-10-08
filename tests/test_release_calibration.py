"""Tests for internal-release calibration."""

import numpy as np
import pytest

from liiontr.chemistry.release_calibration import (
    fit_first_order_release_time,
)


def test_first_order_release_fit_recovers_time_constant():
    """Recover tau from an exact first-order release curve."""
    expected_time_constant = 2.5

    time = np.linspace(
        0.0,
        15.0,
        151,
    )

    released_fraction = (
        1.0
        - np.exp(
            -time
            / expected_time_constant
        )
    )

    fit = fit_first_order_release_time(
        time=time,
        released_fraction=(
            released_fraction
        ),
    )

    assert fit.time_constant == pytest.approx(
        expected_time_constant,
        rel=1.0e-6,
    )

    assert fit.rmse < 1.0e-8

    assert (
        fit.maximum_absolute_error
        < 1.0e-8
    )


def test_release_fit_rejects_nonmonotonic_fraction():
    """Reject a decreasing cumulative release fraction."""
    with pytest.raises(
        ValueError,
        match="nondecreasing",
    ):
        fit_first_order_release_time(
            time=np.asarray(
                [
                    0.0,
                    1.0,
                    2.0,
                ]
            ),
            released_fraction=np.asarray(
                [
                    0.0,
                    0.8,
                    0.7,
                ]
            ),
        )