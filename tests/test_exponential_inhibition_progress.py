import math

import pytest

from liiontr.kinetics import (
    ExponentialInhibitionProgress,
    PowerLawProgress,
)
from liiontr.reactions import ReactionContext


def test_exponential_inhibition_reduces_progress():
    model = ExponentialInhibitionProgress(
        progress_model=PowerLawProgress(
            order=1.0,
        ),
        variable_name="sei_thickness_ratio",
    )

    context = ReactionContext(
        variables={
            "sei_thickness_ratio": 1.0,
        }
    )

    factor = model.factor(
        conversion=0.25,
        context=context,
    )

    expected = 0.75 * math.exp(-1.0)

    assert factor == pytest.approx(expected)


def test_zero_inhibition_preserves_progress():
    model = ExponentialInhibitionProgress(
        progress_model=PowerLawProgress(
            order=1.0,
        ),
        variable_name="sei_thickness_ratio",
    )

    context = ReactionContext(
        variables={
            "sei_thickness_ratio": 0.0,
        }
    )

    factor = model.factor(
        conversion=0.25,
        context=context,
    )

    assert factor == pytest.approx(0.75)


def test_exponential_inhibition_requires_context():
    model = ExponentialInhibitionProgress(
        progress_model=PowerLawProgress(
            order=1.0,
        ),
        variable_name="sei_thickness_ratio",
    )

    with pytest.raises(
        ValueError,
        match="Reaction context is required",
    ):
        model.factor(
            conversion=0.25,
        )


def test_exponential_inhibition_uses_context_variable():
    model = ExponentialInhibitionProgress(
        progress_model=PowerLawProgress(
            order=1.0,
        ),
        variable_name="sei_thickness_ratio",
    )

    context = ReactionContext(
        variables={
            "sei_thickness_ratio": 2.0,
        }
    )

    factor = model.factor(
        conversion=0.0,
        context=context,
    )

    assert factor == pytest.approx(math.exp(-2.0))
