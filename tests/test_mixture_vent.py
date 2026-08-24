import pytest

from liiontr.gases import (
    CompressibleVentFlowModel,
    GasInventory,
    GasSpecies,
    MixtureVentFlowModel,
)


def test_mixture_vent_distributes_flow_by_mole_fraction():
    inventory = GasInventory(
        species=[
            GasSpecies(
                name="CO2",
                molar_mass=44.0e-3,
            ),
            GasSpecies(
                name="H2",
                molar_mass=2.0e-3,
            ),
        ],
        moles={
            "CO2": 2.0e-3,
            "H2": 3.0e-3,
        },
    )

    flow_model = CompressibleVentFlowModel(
        vent_area=1.0e-6,
        discharge_coefficient=0.8,
        heat_capacity_ratio=1.30,
    )

    vent = MixtureVentFlowModel(
        flow_model=flow_model,
        downstream_pressure=101325.0,
    )

    rates = vent.species_molar_flow_rates(
        inventory=inventory,
        upstream_pressure=2.0e6,
        temperature=600.0,
    )

    total_rate = sum(
        rates.values()
    )

    assert rates["CO2"] / total_rate == pytest.approx(
        2.0 / 5.0
    )

    assert rates["H2"] / total_rate == pytest.approx(
        3.0 / 5.0
    )


def test_species_rates_sum_to_total_molar_flow():
    inventory = GasInventory(
        species=[
            GasSpecies(
                name="CO2",
                molar_mass=44.0e-3,
            ),
            GasSpecies(
                name="H2",
                molar_mass=2.0e-3,
            ),
        ],
        moles={
            "CO2": 2.0e-3,
            "H2": 3.0e-3,
        },
    )

    vent = MixtureVentFlowModel(
        flow_model=CompressibleVentFlowModel(
            vent_area=1.0e-6,
            discharge_coefficient=0.8,
            heat_capacity_ratio=1.30,
        ),
    )

    rates = vent.species_molar_flow_rates(
        inventory=inventory,
        upstream_pressure=2.0e6,
        temperature=600.0,
    )

    total_rate = vent.total_molar_flow_rate(
        inventory=inventory,
        upstream_pressure=2.0e6,
        temperature=600.0,
    )

    assert sum(
        rates.values()
    ) == pytest.approx(
        total_rate
    )


def test_species_rates_preserve_mass_flow():
    inventory = GasInventory(
        species=[
            GasSpecies(
                name="CO2",
                molar_mass=44.0e-3,
            ),
            GasSpecies(
                name="H2",
                molar_mass=2.0e-3,
            ),
        ],
        moles={
            "CO2": 2.0e-3,
            "H2": 3.0e-3,
        },
    )

    flow_model = CompressibleVentFlowModel(
        vent_area=1.0e-6,
        discharge_coefficient=0.8,
        heat_capacity_ratio=1.30,
    )

    vent = MixtureVentFlowModel(
        flow_model=flow_model,
    )

    rates = vent.species_molar_flow_rates(
        inventory=inventory,
        upstream_pressure=2.0e6,
        temperature=600.0,
    )

    species_by_name = {
        species.name: species
        for species in inventory.species
    }

    mass_flow_from_species = sum(
        rate
        * species_by_name[name].molar_mass
        for name, rate in rates.items()
    )

    expected_mass_flow = flow_model.mass_flow_rate(
        upstream_pressure=2.0e6,
        downstream_pressure=101325.0,
        temperature=600.0,
        molar_mass=inventory.mean_molar_mass,
    )

    assert mass_flow_from_species == pytest.approx(
        expected_mass_flow
    )


def test_empty_inventory_has_zero_vent_flow():
    inventory = GasInventory(
        species=[
            GasSpecies(
                name="CO2",
                molar_mass=44.0e-3,
            ),
        ]
    )

    vent = MixtureVentFlowModel(
        flow_model=CompressibleVentFlowModel(
            vent_area=1.0e-6,
        ),
    )

    rates = vent.species_molar_flow_rates(
        inventory=inventory,
        upstream_pressure=2.0e6,
        temperature=600.0,
    )

    assert rates == {
        "CO2": 0.0,
    }

    assert vent.total_molar_flow_rate(
        inventory=inventory,
        upstream_pressure=2.0e6,
        temperature=600.0,
    ) == pytest.approx(
        0.0
    )


def test_no_flow_at_ambient_pressure():
    inventory = GasInventory(
        species=[
            GasSpecies(
                name="N2",
                molar_mass=28.0134e-3,
            ),
        ],
        moles={
            "N2": 1.0e-3,
        },
    )

    vent = MixtureVentFlowModel(
        flow_model=CompressibleVentFlowModel(
            vent_area=1.0e-6,
        ),
        downstream_pressure=101325.0,
    )

    rates = vent.species_molar_flow_rates(
        inventory=inventory,
        upstream_pressure=101325.0,
        temperature=400.0,
    )

    assert rates["N2"] == pytest.approx(
        0.0
    )