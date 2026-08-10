from liiontr.library import cell_21700_generic
from liiontr.problems import ThermalProblem
from liiontr.solver.scipy_solver import ScipySolver


cell = cell_21700_generic()


problem = ThermalProblem(
    cell=cell,
    initial_temperature=350.0,
    duration=3600.0,
)


solver = ScipySolver()


solution = solver.solve(problem)


print(solution.time)
print(solution.temperature)
