from liiontr.kinetics import (
    Arrhenius,
    AutocatalyticProgress,
)
from liiontr.reactions import (
    MultiChannelReaction,
    ReactionChannel,
    ReactionNetwork,
)


def test_reaction_network_accepts_multichannel_reaction():
    reaction = MultiChannelReaction(
        name="Cathode decomposition",
        channels=[
            ReactionChannel(
                kinetics=Arrhenius(
                    activation_energy=80000.0,
                    pre_exponential_factor=1.0e5,
                ),
                enthalpy=200000.0,
            ),
            ReactionChannel(
                kinetics=Arrhenius(
                    activation_energy=90000.0,
                    pre_exponential_factor=2.0e5,
                ),
                enthalpy=300000.0,
            ),
        ],
        mass_fraction=0.20,
        progress_model=AutocatalyticProgress(
            autocatalytic_order=1.0,
            remaining_order=1.0,
        ),
    )

    network = ReactionNetwork(
        reactions=[
            reaction,
        ]
    )

    rates = network.progress_rates(
        temperature=500.0,
        conversions=[
            0.04,
        ],
    )

    heat_generation = network.heat_generation(
        temperature=500.0,
        conversions=[
            0.04,
        ],
    )

    assert len(rates) == 1
    assert rates[0] > 0.0
    assert heat_generation > 0.0
