import pytest

from liiontr.gases import ReactionGasYield


def test_reaction_gas_yield_stores_species_yields():
    gas_yield = ReactionGasYield(
        reaction_name="Electrolyte decomposition",
        species_yields={
            "CO2": 2.0,
            "H2": 0.5,
        },
    )

    assert gas_yield.reaction_name == ("Electrolyte decomposition")

    assert gas_yield.yield_of("CO2") == pytest.approx(2.0)

    assert gas_yield.yield_of("H2") == pytest.approx(0.5)


def test_missing_species_has_zero_yield():
    gas_yield = ReactionGasYield(
        reaction_name="Reaction",
        species_yields={
            "CO2": 1.0,
        },
    )

    assert gas_yield.yield_of("H2") == pytest.approx(0.0)


def test_reaction_gas_yield_computes_generation_rates():
    gas_yield = ReactionGasYield(
        reaction_name="Reaction",
        species_yields={
            "CO2": 2.0,
            "H2": 0.5,
        },
    )

    rates = gas_yield.generation_rates(
        cell_mass=0.060,
        reaction_mass_fraction=0.20,
        progress_rate=0.50,
    )

    reacted_mass_rate = 0.060 * 0.20 * 0.50

    assert rates["CO2"] == pytest.approx(2.0 * reacted_mass_rate)

    assert rates["H2"] == pytest.approx(0.5 * reacted_mass_rate)


def test_negative_species_yield_is_rejected():
    with pytest.raises(
        ValueError,
        match="Gas yield",
    ):
        ReactionGasYield(
            reaction_name="Reaction",
            species_yields={
                "CO2": -1.0,
            },
        )


def test_empty_reaction_name_is_rejected():
    with pytest.raises(
        ValueError,
        match="Reaction name",
    ):
        ReactionGasYield(
            reaction_name="",
            species_yields={
                "CO2": 1.0,
            },
        )


def test_negative_progress_rate_is_rejected():
    gas_yield = ReactionGasYield(
        reaction_name="Reaction",
        species_yields={
            "CO2": 1.0,
        },
    )

    with pytest.raises(
        ValueError,
        match="Progress rate",
    ):
        gas_yield.generation_rates(
            cell_mass=0.060,
            reaction_mass_fraction=0.20,
            progress_rate=-1.0,
        )
