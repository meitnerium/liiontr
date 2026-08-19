import pytest

from liiontr.kinetics import (
    PowerLawProgress,
    ThresholdProgress,
)
from liiontr.reactions import ReactionContext


def test_threshold_progress_is_inactive_above_threshold():
    model = ThresholdProgress(
        progress_model=PowerLawProgress(
            order=1.0,
        ),
        reaction_name="SEI decomposition",
        remaining_below=0.10,
    )

    context = ReactionContext(
        conversions={
            "SEI decomposition": 0.80,
        }
    )

    factor = model.factor(
        conversion=0.25,
        context=context,
    )

    assert factor == pytest.approx(0.0)


def test_threshold_progress_activates_below_threshold():
    model = ThresholdProgress(
        progress_model=PowerLawProgress(
            order=1.0,
        ),
        reaction_name="SEI decomposition",
        remaining_below=0.10,
    )

    context = ReactionContext(
        conversions={
            "SEI decomposition": 0.95,
        }
    )

    factor = model.factor(
        conversion=0.25,
        context=context,
    )

    assert factor == pytest.approx(0.75)


def test_threshold_progress_requires_context():
    model = ThresholdProgress(
        progress_model=PowerLawProgress(
            order=1.0,
        ),
        reaction_name="SEI decomposition",
        remaining_below=0.10,
    )

    with pytest.raises(
        ValueError,
        match="Reaction context is required",
    ):
        model.factor(
            conversion=0.25,
        )


@pytest.mark.parametrize(
    "remaining_below",
    [
        -0.1,
        1.1,
    ],
)
def test_threshold_must_be_between_zero_and_one(
    remaining_below: float,
):
    with pytest.raises(
        ValueError,
        match="Remaining fraction threshold must be between 0 and 1",
    ):
        ThresholdProgress(
            progress_model=PowerLawProgress(
                order=1.0,
            ),
            reaction_name="SEI decomposition",
            remaining_below=remaining_below,
        )
