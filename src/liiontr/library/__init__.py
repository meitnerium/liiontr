"""Reference cells and literature-based thermal runaway models."""

from .cells import cell_18650_lco_hatchard2001, cell_21700_generic
from .hatchard2001 import (
    HATCHARD2001_REFERENCE,
    hatchard2001_anode_electrolyte,
    hatchard2001_boundary_heat_transfer_coefficient,
    hatchard2001_cathode_reaction,
    hatchard2001_initial_conversions,
    hatchard2001_lumped_heat_transfer_coefficient,
    hatchard2001_reaction_network,
    hatchard2001_sei_decomposition,
)
from .hu2020 import (
    HU2020_REFERENCE,
    hu2020_anode_electrolyte,
    hu2020_anode_progress_model,
    hu2020_cathode_decomposition,
    hu2020_cathode_initial_conversion,
    hu2020_electrolyte_decomposition,
    hu2020_electrolyte_reaction,
    hu2020_initial_conversions,
    hu2020_reaction_network,
    hu2020_sei_decomposition,
)
from .reactions import (
    ReactionParameters,
    VolumetricConversionReactionParameters,
    VolumetricReactionParameters,
)

__all__ = [
    "HATCHARD2001_REFERENCE",
    "HU2020_REFERENCE",
    "ReactionParameters",
    "VolumetricConversionReactionParameters",
    "VolumetricReactionParameters",
    "cell_18650_lco_hatchard2001",
    "cell_21700_generic",
    "hatchard2001_anode_electrolyte",
    "hatchard2001_boundary_heat_transfer_coefficient",
    "hatchard2001_cathode_reaction",
    "hatchard2001_initial_conversions",
    "hatchard2001_lumped_heat_transfer_coefficient",
    "hatchard2001_reaction_network",
    "hatchard2001_sei_decomposition",
    "hu2020_anode_electrolyte",
    "hu2020_anode_progress_model",
    "hu2020_cathode_decomposition",
    "hu2020_cathode_initial_conversion",
    "hu2020_electrolyte_decomposition",
    "hu2020_electrolyte_reaction",
    "hu2020_initial_conversions",
    "hu2020_reaction_network",
    "hu2020_sei_decomposition",
]
