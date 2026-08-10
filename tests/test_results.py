import numpy as np

from liiontr.core.results import Results


def test_results_storage():
    results = Results()

    temperature = np.array([300.0, 310.0, 350.0])

    results.add_variable(
        "temperature",
        temperature,
    )

    assert "temperature" in results.variables

    assert results.get("temperature")[1] == 310.0
