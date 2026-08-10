from liiontr.core.domain import Domain


class DummyDomain(Domain):
    pass


def test_domain():
    domain = DummyDomain(
        name="18650",
        dimension=1,
    )

    assert domain.name == "18650"

    assert domain.dimension == 1

    assert domain.describe() == "18650 (1D)"
