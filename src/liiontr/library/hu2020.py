from __future__ import annotations

from typing import TYPE_CHECKING

from liiontr.kinetics import (
    ExponentialInhibitionProgress,
    PowerLawProgress,
    ThresholdProgress,
)
from liiontr.reactions import (
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
    """

    sei_parameters = hu2020_sei_decomposition()
    anode_parameters = hu2020_anode_electrolyte()

    return [
        sei_parameters.initial_conversion,
        anode_parameters.initial_conversion,
    ]


def hu2020_reaction_network(
    cell: Cell,
) -> ReactionNetwork:
    """
    Build the Hu 2020 SEI + anode reaction network.

    The SEI thickness ratio is treated as an algebraic variable:

        sei_thickness_ratio
            = current_SEI_remaining
            / reference_SEI_remaining

    with a reference remaining SEI fraction of 0.15.
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

    return ReactionNetwork(
        reactions=[
            sei_reaction,
            anode_reaction,
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
