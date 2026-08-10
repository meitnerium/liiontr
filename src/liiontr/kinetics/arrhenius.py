from __future__ import annotations

import math

from .reaction import Reaction


R = 8.314462618


def reaction_rate(
    reaction: Reaction,
    temperature: float,
) -> float:
    """
    Arrhenius reaction rate.
    """

    return reaction.pre_exponential_factor * math.exp(
        -reaction.activation_energy / (R * temperature)
    )


def heat_generation(
    reaction: Reaction,
    temperature: float,
) -> float:
    """
    Heat generation rate [W/kg].
    """

    return (
        reaction.enthalpy
        * reaction.mass_fraction
        * reaction_rate(
            reaction,
            temperature,
        )
    )
