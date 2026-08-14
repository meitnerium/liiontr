import pytest

from liiontr.kinetics import Arrhenius, ProgressModel
from liiontr.reactions import (
    Reaction,
    ReactionContext,
    ReactionNetwork,
)


class SEIDependentProgress(ProgressModel):
    """
    Test progress model depending on SEI conversion.
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


def test_network_passes_context_to_progress_models():
    sei = Reaction(
        name="SEI decomposition",
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=200000.0,
    )

    anode = Reaction(
        name="Anode-electrolyte",
        kinetics=Arrhenius(
            activation_energy=90000.0,
            pre_exponential_factor=2.0e5,
        ),
        enthalpy=500000.0,
        progress_model=SEIDependentProgress(),
    )

    network = ReactionNetwork(
        reactions=[
            sei,
            anode,
        ]
    )

    conversions = [
        0.80,
        0.25,
    ]

    rates = network.progress_rates(
        temperature=500.0,
        conversions=conversions,
    )

    expected_anode_rate = anode.rate(500.0) * 0.75 * 0.20

    assert rates[1] == pytest.approx(expected_anode_rate)


def test_network_passes_context_to_heat_generation():
    sei = Reaction(
        name="SEI decomposition",
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=200000.0,
    )

    anode = Reaction(
        name="Anode-electrolyte",
        kinetics=Arrhenius(
            activation_energy=90000.0,
            pre_exponential_factor=2.0e5,
        ),
        enthalpy=500000.0,
        mass_fraction=0.10,
        progress_model=SEIDependentProgress(),
    )

    network = ReactionNetwork(
        reactions=[
            sei,
            anode,
        ]
    )

    conversions = [
        0.80,
        0.25,
    ]

    total_heat = network.heat_generation(
        temperature=500.0,
        conversions=conversions,
    )

    sei_heat = sei.heat_generation(
        temperature=500.0,
        conversion=0.80,
    )

    expected_anode_heat = (
        anode.enthalpy * anode.mass_fraction * anode.rate(500.0) * 0.75 * 0.20
    )

    assert total_heat == pytest.approx(sei_heat + expected_anode_heat)
