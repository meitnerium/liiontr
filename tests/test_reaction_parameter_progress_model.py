from liiontr.kinetics import (
    ExponentialInhibitionProgress,
    PowerLawProgress,
)
from liiontr.library import (
    VolumetricReactionParameters,
    cell_21700_generic,
)


def test_volumetric_parameters_accept_progress_model():
    cell = cell_21700_generic()

    parameters = VolumetricReactionParameters(
        name="Synthetic reaction",
        activation_energy=80000.0,
        pre_exponential_factor=1.0e5,
        enthalpy=500000.0,
        specific_content=500.0,
        initial_remaining_fraction=0.75,
    )

    progress_model = ExponentialInhibitionProgress(
        progress_model=PowerLawProgress(
            order=1.0,
        ),
        variable_name="sei_thickness_ratio",
    )

    reaction = parameters.build(
        cell=cell,
        progress_model=progress_model,
    )

    assert reaction.progress_model is progress_model
