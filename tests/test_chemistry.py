from liiontr.chemistry.nmc import NMC811


def test_nmc811():
    chemistry = NMC811()

    assert chemistry.name == "NMC811"
    assert chemistry.nominal_voltage > 3.0
    assert chemistry.specific_capacity > 100.0
