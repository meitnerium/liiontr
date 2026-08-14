from __future__ import annotations

from .reactions import VolumetricReactionParameters


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
