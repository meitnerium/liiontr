import pytest

from liiontr.kinetics import (
    Arrhenius,
    AutocatalyticProgress,
)
from liiontr.reactions import (
    MultiChannelReaction,
    ReactionChannel,
)


def test_multichannel_reaction_progress_rate():
    channel_1 = ReactionChannel(
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=200000.0,
    )

    channel_2 = ReactionChannel(
        kinetics=Arrhenius(
            activation_energy=90000.0,
            pre_exponential_factor=2.0e5,
        ),
        enthalpy=300000.0,
    )

    reaction = MultiChannelReaction(
        name="Synthetic multichannel reaction",
        channels=[
            channel_1,
            channel_2,
        ],
        mass_fraction=0.20,
        progress_model=AutocatalyticProgress(
            autocatalytic_order=1.0,
            remaining_order=1.0,
        ),
    )

    temperature = 500.0
    conversion = 0.25

    factor = conversion * (1.0 - conversion)

    expected = (
        channel_1.kinetics.rate(temperature) + channel_2.kinetics.rate(temperature)
    ) * factor

    rate = reaction.progress_rate(
        temperature=temperature,
        conversion=conversion,
    )

    assert rate == pytest.approx(expected)


def test_multichannel_reaction_heat_generation():
    channel_1 = ReactionChannel(
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        enthalpy=200000.0,
    )

    channel_2 = ReactionChannel(
        kinetics=Arrhenius(
            activation_energy=90000.0,
            pre_exponential_factor=2.0e5,
        ),
        enthalpy=300000.0,
    )

    reaction = MultiChannelReaction(
        name="Synthetic multichannel reaction",
        channels=[
            channel_1,
            channel_2,
        ],
        mass_fraction=0.20,
        progress_model=AutocatalyticProgress(
            autocatalytic_order=1.0,
            remaining_order=1.0,
        ),
    )

    temperature = 500.0
    conversion = 0.25

    factor = conversion * (1.0 - conversion)

    expected = (
        reaction.mass_fraction
        * factor
        * (
            channel_1.enthalpy * channel_1.kinetics.rate(temperature)
            + channel_2.enthalpy * channel_2.kinetics.rate(temperature)
        )
    )

    heat_generation = reaction.heat_generation(
        temperature=temperature,
        conversion=conversion,
    )

    assert heat_generation == pytest.approx(expected)


def test_multichannel_reaction_requires_channel():
    with pytest.raises(
        ValueError,
        match="At least one reaction channel is required",
    ):
        MultiChannelReaction(
            name="Invalid reaction",
            channels=[],
            mass_fraction=0.20,
            progress_model=AutocatalyticProgress(
                autocatalytic_order=1.0,
                remaining_order=1.0,
            ),
        )
