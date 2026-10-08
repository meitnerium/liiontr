Quick Start
===========

LiionTR simulations are assembled from a cell, physical models, a problem
definition, and a numerical solver.

Minimal Thermal Simulation
--------------------------

The bundled ``examples/first_simulation.py`` demonstrates the smallest thermal
problem:

.. code-block:: python

   from liiontr.library import cell_21700_generic
   from liiontr.problems import ThermalProblem
   from liiontr.solver import ScipySolver

   cell = cell_21700_generic()

   problem = ThermalProblem(
       cell=cell,
       initial_temperature=350.0,
       duration=3600.0,
   )

   results = ScipySolver().solve(problem)

   print(results.time)
   print(results.temperature)

Run it from the repository root:

.. code-block:: bash

   python examples/first_simulation.py

Hu et al. (2020) Reaction Network
---------------------------------

A literature-based thermal-runaway calculation combines the generic 21700 cell
with the Hu et al. reaction network:

.. code-block:: python

   from liiontr.chemistry import ReactionNetworkBackend
   from liiontr.library import cell_21700_generic
   from liiontr.library.hu2020 import (
       hu2020_initial_conversions,
       hu2020_reaction_network,
   )
   from liiontr.problems import ThermalProblem
   from liiontr.solver import ScipySolver

   cell = cell_21700_generic()
   network = hu2020_reaction_network(cell=cell)

   backend = ReactionNetworkBackend(
       reaction_network=network,
       cell=cell,
   )

   problem = ThermalProblem(
       cell=cell,
       chemistry_backend=backend,
       initial_temperature=480.0,
       initial_conversions=hu2020_initial_conversions(),
       ambient_temperature=298.15,
       convection_coefficient=10.0,
       duration=20.0,
       maximum_temperature=1200.0,
   )

   results = ScipySolver().solve(problem)

   print(results.time[-1])
   print(results.temperature[-1])

The complete example, including reaction-energy post-processing, is available
in ``examples/hu2020_simulation.py``.

Gas Pressure and Venting
------------------------

``examples/vented_thermal_runaway.py`` demonstrates the empirical gas path:
reaction progress generates species, the gas inventory determines pressure,
and an irreversible vent-opening event switches to compressible discharge.

.. code-block:: bash

   python examples/vented_thermal_runaway.py

Elemental Thermochemistry
-------------------------

The elemental thermochemistry path is more advanced and currently assembled
through the lower-level model interfaces. Its conceptual workflow is

.. code-block:: text

   reaction or reaction-energy release
              |
              v
      elemental generation
              |
              v
      elemental inventory
              |
      optional release delay
              |
              v
        Cantera equilibrium
              |
       +------+------+ 
       |             |
       v             v
    pressure     gas energy
       |             |
       +------> venting

See :doc:`theory/thermochemistry` and the Cantera-related tests for the current
interfaces.

Results
-------

Simulation outputs are stored in :class:`liiontr.core.results.Results`.
Depending on the configured problem, they may include temperature, reaction
conversions, gas inventories, pressure, vent state, elemental inventories,
thermal energy, and cumulative vented quantities.
