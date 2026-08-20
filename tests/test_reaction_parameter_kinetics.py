from liiontr.kinetics import (
    Arrhenius,
    TemperatureThresholdKinetics,
)
from liiontr.library import (
    VolumetricConversionReactionParameters,
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
