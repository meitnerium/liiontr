"""Nickel-manganese-cobalt lithium-ion battery chemistries."""

from dataclasses import dataclass

from .chemistry import Chemistry


@dataclass(slots=True)
class NMC811(Chemistry):
    """Represent a generic NMC811 lithium-ion battery chemistry.

    Notes
    -----
    The current implementation defines nominal bulk properties for an
    NMC811 chemistry:

    - nominal voltage: 3.65 V
    - specific capacity: 200 Ah/kg
    - density: 4800 kg/m³

    These values currently act as generic chemistry parameters rather
    than a cell-specific calibrated parameter set.
    """

    def __init__(self) -> None:
        """Initialize the default NMC811 chemistry."""
        super().__init__(
            name="NMC811",
            nominal_voltage=3.65,
            specific_capacity=200.0,
            density=4800.0,
        )
