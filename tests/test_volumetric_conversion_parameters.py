import pytest

from liiontr.library import (
    VolumetricConversionReactionParameters,
    cell_21700_generic,
)


def test_volumetric_conversion_parameters_store_initial_conversion():
    parameters = VolumetricConversionReactionParameters(
        name="Synthetic autocatalytic reaction",
        activation_energy=80000.0,
        pre_exponential_factor=1.0e5,
        enthalpy=500000.0,
        specific_content=500.0,
        initial_conversion=0.04,
    )

    assert parameters.initial_conversion == pytest.approx(0.04)


def test_volumetric_conversion_parameters_build_reaction():
    cell = cell_21700_generic()

    parameters = VolumetricConversionReactionParameters(
        name="Synthetic autocatalytic reaction",
        activation_energy=80000.0,
        pre_exponential_factor=1.0e5,
        enthalpy=500000.0,
        specific_content=500.0,
        initial_conversion=0.04,
    )

    reaction = parameters.build(
        cell=cell,
    )

    expected_mass_fraction = 500.0 / cell.material.density(298.15)

    assert reaction.mass_fraction == pytest.approx(expected_mass_fraction)


@pytest.mark.parametrize(
    "initial_conversion",
    [
        -0.1,
        1.1,
    ],
)
def test_initial_conversion_must_be_between_zero_and_one(
    initial_conversion: float,
):
    with pytest.raises(
        ValueError,
        match="Initial conversion must be between 0 and 1",
    ):
        VolumetricConversionReactionParameters(
            name="Synthetic autocatalytic reaction",
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
            enthalpy=500000.0,
            specific_content=500.0,
            initial_conversion=initial_conversion,
        )
