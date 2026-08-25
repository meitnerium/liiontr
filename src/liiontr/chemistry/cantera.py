"""Cantera-based thermochemical backend for LiionTR."""

from __future__ import annotations

from .backend import ChemistryBackend


class CanteraBackend(ChemistryBackend):
    """Provide the interface for Cantera-based thermochemistry.

    Notes
    -----
    The backend currently stores the Cantera mechanism identifier but
    does not yet implement heat-generation calculations. Future
    versions are expected to use Cantera for gas-phase thermodynamics
    and chemical-equilibrium calculations.
    """

    def __init__(
        self,
        mechanism: str,
    ) -> None:
        """Initialize the Cantera backend.

        Parameters
        ----------
        mechanism : str
            Cantera mechanism file or mechanism identifier to use for
            thermochemical calculations.
        """
        self.mechanism = mechanism

    def heat_generation(
        self,
        temperature: float,
        conversions: list[float] | None = None,
    ) -> float:
        """Return the total cell heat generation rate in W.

        Parameters
        ----------
        temperature : float
            Cell temperature in K.
        conversions : list[float] or None, optional
            Reaction conversion values. This argument is reserved for
            compatibility with the chemistry backend interface.

        Returns
        -------
        float
            Total cell heat generation rate in W.

        Raises
        ------
        NotImplementedError
            Always raised because Cantera heat generation is not yet
            implemented.
        """
        raise NotImplementedError
