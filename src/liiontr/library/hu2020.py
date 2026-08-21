from __future__ import annotations

from typing import TYPE_CHECKING

from liiontr.kinetics import (
    Arrhenius,
    AutocatalyticProgress,
    ExponentialInhibitionProgress,
    PowerLawProgress,
    TemperatureThresholdKinetics,
    ThresholdProgress,
)
from liiontr.reactions import (
    MultiChannelReaction,
    ReactionChannel,
    ReactionNetwork,
    RemainingFractionRatioVariable,
)

from .reactions import VolumetricReactionParameters

if TYPE_CHECKING:
    from liiontr.cells.cell import Cell


HU2020_REFERENCE = "Hu et al., ACS Omega (2020), DOI: 10.1021/acsomega.0c01862"


def hu2020_sei_decomposition() -> VolumetricReactionParameters:
    """
    Return the SEI decomposition parameters reported by Hu et al. (2020).

    Original reported parameters
    ----------------------------
    A:
        1.667e15 1/s

    Ea:
        1.3508e5 J/mol

    H:
        257 J/g

    Wc:
        6.104e5 g/m3

    Initial remaining fraction:
        0.15

    Reaction order:
        1.0

    LiionTR unit conversions
    ------------------------
    H:
        257 J/g -> 257000 J/kg

    Wc:
        6.104e5 g/m3 -> 610.4 kg/m3

    LiionTR uses conversion alpha, where:

        alpha = 1 - remaining_fraction

    Therefore:

        alpha_initial = 1 - 0.15 = 0.85
    """

    return VolumetricReactionParameters(
        name="SEI decomposition",
        activation_energy=1.3508e5,
        pre_exponential_factor=1.667e15,
        enthalpy=257000.0,
        specific_content=610.4,
        initial_remaining_fraction=0.15,
        reaction_order=1.0,
        reference=HU2020_REFERENCE,
    )


def hu2020_anode_electrolyte() -> VolumetricReactionParameters:
    """
    Return the anode-electrolyte reaction parameters
    reported by Hu et al. (2020).

    Original reported parameters
    ----------------------------
    A:
        2.5e13 1/s

    Ea:
        1.3508e5 J/mol

    H:
        1714 J/g

    Wc:
        6.104e5 g/m3

    Initial remaining fraction:
        0.75

    Reaction order:
        1.0

    LiionTR unit conversions
    ------------------------
    H:
        1714 J/g -> 1.714e6 J/kg

    Wc:
        6.104e5 g/m3 -> 610.4 kg/m3

    Initial conversion:

        alpha_initial = 1 - 0.75 = 0.25

    The SEI-dependent threshold and exponential inhibition
    are represented separately by hu2020_anode_progress_model().
    """

    return VolumetricReactionParameters(
        name="Anode-electrolyte",
        activation_energy=1.3508e5,
        pre_exponential_factor=2.5e13,
        enthalpy=1.714e6,
        specific_content=610.4,
        initial_remaining_fraction=0.75,
        reaction_order=1.0,
        reference=HU2020_REFERENCE,
    )


def hu2020_anode_progress_model() -> ThresholdProgress:
    """
    Return the coupled progress model for the Hu 2020
    anode-electrolyte reaction.

    The base reaction follows:

        f(alpha) = (1 - alpha)

    and is inhibited by:

        exp(-sei_thickness_ratio)

    The reaction is active only when the remaining SEI
    fraction is below 0.10.
    """

    return ThresholdProgress(
        progress_model=ExponentialInhibitionProgress(
            progress_model=PowerLawProgress(
                order=1.0,
            ),
            variable_name="sei_thickness_ratio",
        ),
        reaction_name="SEI decomposition",
        remaining_below=0.10,
    )


def hu2020_initial_conversions() -> list[float]:
    """
    Return initial reaction conversions for the Hu 2020 network.

    Reaction ordering:

        0: SEI decomposition
        1: Anode-electrolyte
        2: Cathode decomposition
    """

    sei_parameters = hu2020_sei_decomposition()
    anode_parameters = hu2020_anode_electrolyte()

    return [
        sei_parameters.initial_conversion,
        anode_parameters.initial_conversion,
        hu2020_cathode_initial_conversion(),
    ]


def hu2020_reaction_network(
    cell: Cell,
) -> ReactionNetwork:
    """
    Build the Hu 2020 SEI, anode, and cathode reaction network.

    Reaction ordering:

        0: SEI decomposition
        1: Anode-electrolyte
        2: Cathode decomposition

    The SEI thickness ratio is treated as an algebraic variable:

        sei_thickness_ratio
            = current_SEI_remaining
            / reference_SEI_remaining
    """

    sei_parameters = hu2020_sei_decomposition()

    sei_reaction = sei_parameters.build(
        cell=cell,
    )

    anode_parameters = hu2020_anode_electrolyte()

    anode_reaction = anode_parameters.build(
        cell=cell,
        progress_model=hu2020_anode_progress_model(),
    )

    cathode_reaction = hu2020_cathode_decomposition(
        cell=cell,
    )

    return ReactionNetwork(
        reactions=[
            sei_reaction,
            anode_reaction,
            cathode_reaction,
        ],
        context_variables={
            "sei_thickness_ratio": RemainingFractionRatioVariable(
                reaction_name="SEI decomposition",
                reference_remaining_fraction=(
                    sei_parameters.initial_remaining_fraction
                ),
            )
        },
    )


def hu2020_cathode_initial_conversion() -> float:
    """
    Return the initial cathode reaction conversion reported
    by Hu et al. (2020).
    """

    return 0.04


def hu2020_cathode_decomposition(
    cell: Cell,
) -> MultiChannelReaction:
    """
    Build the Hu 2020 cathode decomposition reaction.

    The cathode is represented by two parallel kinetic channels
    sharing a single reaction conversion alpha.

    Both channels are active only at temperatures greater than
    or equal to 393.15 K.

    Reported cathode content:

        Wp = 1.221e6 g/m3
           = 1221 kg/m3

    Channel 1
    ---------
    A:
        1.75e9 1/s

    Ea:
        1.1495e5 J/mol

    H:
        77 J/g
        = 77000 J/kg

    Channel 2
    ---------
    A:
        1.077e12 1/s

    Ea:
        1.5888e5 J/mol

    H:
        84 J/g
        = 84000 J/kg

    Progress law
    ------------
    f(alpha) = alpha * (1 - alpha)
    """

    cell_density = cell.material.density(298.15)

    mass_fraction = 1221.0 / cell_density

    if mass_fraction > 1.0:
        raise ValueError("Cathode volumetric content exceeds total cell mass.")

    channel_1 = ReactionChannel(
        kinetics=TemperatureThresholdKinetics(
            kinetics=Arrhenius(
                activation_energy=1.1495e5,
                pre_exponential_factor=1.75e9,
            ),
            minimum_temperature=393.15,
        ),
        enthalpy=77000.0,
    )

    channel_2 = ReactionChannel(
        kinetics=TemperatureThresholdKinetics(
            kinetics=Arrhenius(
                activation_energy=1.5888e5,
                pre_exponential_factor=1.077e12,
            ),
            minimum_temperature=393.15,
        ),
        enthalpy=84000.0,
    )

    return MultiChannelReaction(
        name="Cathode decomposition",
        channels=[
            channel_1,
            channel_2,
        ],
        mass_fraction=mass_fraction,
        progress_model=AutocatalyticProgress(
            autocatalytic_order=1.0,
            remaining_order=1.0,
        ),
    )
