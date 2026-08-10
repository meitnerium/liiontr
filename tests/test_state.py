from liiontr.core.state import State
from liiontr.core.variable import Variable


def test_state():
    state = State()

    temperature = Variable(
        name="temperature",
        unit="K",
    )

    state.add(temperature)

    assert len(state) == 1

    assert "temperature" in state

    assert state.get("temperature").unit == "K"
