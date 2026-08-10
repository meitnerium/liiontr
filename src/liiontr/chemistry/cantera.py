from __future__ import annotations

from .backend import ChemistryBackend


class CanteraBackend(ChemistryBackend):
    """
    Cantera based thermal chemistry.
    """

    def __init__(
        self,
        mechanism: str,
    ) -> None:
        self.mechanism = mechanism

    def heat_generation(
        self,
        temperature: float,
    ) -> float:
        raise NotImplementedError
