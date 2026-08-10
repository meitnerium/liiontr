from __future__ import annotations

from abc import ABC

import numpy as np


class Field(ABC):
    """
    Base class for every simulation field.
    """

    def __init__(self, values):
        self.values = np.asarray(values, dtype=float)

    @property
    def size(self):
        return self.values.size

    def copy(self):
        return self.__class__(self.values.copy())
