from dataclasses import dataclass


@dataclass(slots=True, frozen=True)
class ConstantProperty:
    """
    Temperature-independent material property.
    """

    value: float
    unit: str

    def __call__(self, temperature: float) -> float:
        return self.value
