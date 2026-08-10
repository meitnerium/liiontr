from dataclasses import dataclass

from .chemistry import Chemistry


@dataclass(slots=True)
class NMC811(Chemistry):
    def __init__(self):
        super().__init__(
            name="NMC811",
            nominal_voltage=3.65,
            specific_capacity=200.0,
            density=4800.0,
        )
