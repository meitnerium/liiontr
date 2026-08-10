from liiontr.core.problem import Problem


class DummyProblem(Problem):
    pass


def test_problem_creation():
    problem = DummyProblem()

    assert problem.parameters is None

    assert problem.domain is None

    assert problem.physics == []
