from __future__ import annotations

import pytest

from liiontr.kinetics import Arrhenius, ProgressModel
from liiontr.reactions import Reaction, ReactionContext


class SEIDependentProgress(ProgressModel):
    """
    Test progress model depending on SEI state.
    """

    def factor(
        self,
        conversion: float,
        context: ReactionContext | None = None,
    ) -> float:
        if context is None:
            raise ValueError("Reaction context is required.")

        sei_remaining = context.remaining_fraction("SEI decomposition")

        return (1.0 - conversion) * sei_remaining


def test_reaction_passes_context_to_progress_model():
    reaction = Reaction(
        name="Anode-electrolyte",
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=500000.0,
        progress_model=SEIDependentProgress(),
    )

    context = ReactionContext(
        conversions={
            "SEI decomposition": 0.80,
        }
    )

    rate = reaction.progress_rate(
        temperature=500.0,
        conversion=0.25,
        context=context,
    )

    expected = reaction.rate(500.0) * 0.75 * 0.20

    assert rate == pytest.approx(expected)
