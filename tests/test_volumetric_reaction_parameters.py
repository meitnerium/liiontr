import pytest

from liiontr.library import (
    VolumetricReactionParameters,
    cell_21700_generic,
)


def test_volumetric_parameters_convert_to_mass_fraction():
    cell = cell_21700_generic()

    parameters = VolumetricReactionParameters(
        name="Synthetic reaction",
        activation_energy=80000.0,
        pre_exponential_factor=1.0e5,
        enthalpy=500000.0,
        specific_content=500.0,
        reaction_order=1.0,
        initial_remaining_fraction=0.25,
        reference="Synthetic test data",
    )

    reaction = parameters.build(
        cell=cell,
    )

    expected_mass_fraction = parameters.specific_content / cell.material.density(298.15)

    assert reaction.mass_fraction == pytest.approx(expected_mass_fraction)


def test_initial_conversion_is_derived_from_remaining_fraction():
    parameters = VolumetricReactionParameters(
        name="Synthetic reaction",
        activation_energy=80000.0,
        pre_exponential_factor=1.0e5,
        enthalpy=500000.0,
        specific_content=500.0,
        initial_remaining_fraction=0.25,
    )

    assert parameters.initial_conversion == pytest.approx(0.75)


@pytest.mark.parametrize(
    "initial_remaining_fraction",
    [
        -0.1,
        1.1,
    ],
)
def test_initial_remaining_fraction_must_be_bounded(
    initial_remaining_fraction: float,
):
    with pytest.raises(
        ValueError,
        match="Initial remaining fraction must be between 0 and 1",
    ):
        VolumetricReactionParameters(
            name="Synthetic reaction",
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
            enthalpy=500000.0,
            specific_content=500.0,
            initial_remaining_fraction=initial_remaining_fraction,
        )


def test_volumetric_content_cannot_exceed_cell_mass():
    cell = cell_21700_generic()

    parameters = VolumetricReactionParameters(
        name="Synthetic reaction",
        activation_energy=80000.0,
        pre_exponential_factor=1.0e5,
        enthalpy=500000.0,
        specific_content=3000.0,
    )

    with pytest.raises(
        ValueError,
        match="Volumetric reactant content exceeds total cell mass",
    ):
        parameters.build(
            cell=cell,
        )
