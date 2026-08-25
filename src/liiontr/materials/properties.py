"""Thermophysical property models for LiionTR materials."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True, frozen=True)
class ConstantProperty:
    """Represent a temperature-independent material property.

    Parameters
    ----------
    value : float
        Numerical value of the property.
    unit : str
        Physical unit associated with the property.
    """

    value: float
    unit: str

    def __call__(self, temperature: float) -> float:
        """Return the constant property value.

        Parameters
        ----------
        temperature : float
            Temperature in K. The value is accepted for API
            compatibility but does not affect a constant property.

        Returns
        -------
        float
            Temperature-independent property value.
        """
        return self.value
