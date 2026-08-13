import pytest

from liiontr.library import ReactionParameters


def test_reaction_parameters_store_kinetic_data():
    parameters = ReactionParameters(
        name="Synthetic reaction",
        activation_energy=80000.0,
        pre_exponential_factor=1.0e5,
        enthalpy=500000.0,
        mass_fraction=0.10,
        reaction_order=2.0,
        reference="Synthetic test data",
    )

    assert parameters.name == "Synthetic reaction"
    assert parameters.activation_energy == 80000.0
    assert parameters.pre_exponential_factor == 1.0e5
    assert parameters.enthalpy == 500000.0
    assert parameters.mass_fraction == 0.10
    assert parameters.reaction_order == 2.0
    assert parameters.reference == "Synthetic test data"


def test_reaction_parameters_build_reaction():
    parameters = ReactionParameters(
        name="Synthetic reaction",
        activation_energy=80000.0,
        pre_exponential_factor=1.0e5,
        enthalpy=500000.0,
        mass_fraction=0.10,
        reaction_order=2.0,
    )

    reaction = parameters.build()

    assert reaction.name == parameters.name
    assert reaction.enthalpy == parameters.enthalpy
    assert reaction.mass_fraction == parameters.mass_fraction

    assert reaction.kinetics.activation_energy == parameters.activation_energy

    assert reaction.kinetics.pre_exponential_factor == parameters.pre_exponential_factor


@pytest.mark.parametrize(
    "activation_energy",
    [
        0.0,
        -1.0,
    ],
)
def test_activation_energy_must_be_positive(
    activation_energy: float,
):
    with pytest.raises(
        ValueError,
        match="Activation energy must be greater than zero",
    ):
        ReactionParameters(
            name="Synthetic reaction",
            activation_energy=activation_energy,
            pre_exponential_factor=1.0e5,
            enthalpy=500000.0,
        )


@pytest.mark.parametrize(
    "pre_exponential_factor",
    [
        0.0,
        -1.0,
    ],
)
def test_pre_exponential_factor_must_be_positive(
    pre_exponential_factor: float,
):
    with pytest.raises(
        ValueError,
        match="Pre-exponential factor must be greater than zero",
    ):
        ReactionParameters(
            name="Synthetic reaction",
            activation_energy=80000.0,
            pre_exponential_factor=pre_exponential_factor,
            enthalpy=500000.0,
        )


@pytest.mark.parametrize(
    "mass_fraction",
    [
        -0.1,
        1.1,
    ],
)
def test_mass_fraction_must_be_between_zero_and_one(
    mass_fraction: float,
):
    with pytest.raises(
        ValueError,
        match="Mass fraction must be between 0 and 1",
    ):
        ReactionParameters(
            name="Synthetic reaction",
            activation_energy=80000.0,
            pre_exponential_factor=1.0e5,
            enthalpy=500000.0,
            mass_fraction=mass_fraction,
        )
