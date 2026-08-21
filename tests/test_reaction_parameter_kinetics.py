from liiontr.kinetics import (
    Arrhenius,
    TemperatureThresholdKinetics,
)
from liiontr.library import (
    VolumetricConversionReactionParameters,
    VolumetricReactionParameters,
    cell_21700_generic,
)


def test_volumetric_conversion_parameters_accept_custom_kinetics():
    cell = cell_21700_generic()

    parameters = VolumetricConversionReactionParameters(
        name="Synthetic cathode reaction",
        activation_energy=80000.0,
        pre_exponential_factor=1.0e5,
        enthalpy=500000.0,
        specific_content=500.0,
        initial_conversion=0.04,
    )

    kinetics = TemperatureThresholdKinetics(
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        minimum_temperature=393.15,
    )

    reaction = parameters.build(
        cell=cell,
        kinetics=kinetics,
    )

    assert reaction.kinetics is kinetics


def test_volumetric_parameters_accept_custom_kinetics():
    cell = cell_21700_generic()

    parameters = VolumetricReactionParameters(
        name="Synthetic electrolyte reaction",
        activation_energy=80000.0,
        pre_exponential_factor=1.0e5,
        enthalpy=500000.0,
        specific_content=500.0,
        initial_remaining_fraction=1.0,
    )

    kinetics = TemperatureThresholdKinetics(
        kinetics=Arrhenius(
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
        ),
        minimum_temperature=473.15,
    )

    reaction = parameters.build(
        cell=cell,
        kinetics=kinetics,
    )

    assert reaction.kinetics is kinetics
