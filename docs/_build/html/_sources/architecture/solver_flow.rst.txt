Solver and State Flow
=====================

Overview
--------

The numerical solver integrates the complete LiionTR simulation state.

The current solver uses SciPy ``solve_ivp`` with a BDF integration
method suitable for stiff thermal runaway kinetics.


State Vector
------------

For a thermal problem containing :math:`N` reactions and :math:`M` gas
species, the solver state is conceptually

.. math::

   \mathbf{y}
   =
   \left[
       T,
       \alpha_1,
       \ldots,
       \alpha_N,
       n_1,
       \ldots,
       n_M
   \right].

The state contains:

* one temperature variable;
* one conversion variable per reaction;
* one mole amount per gas species.


Temperature State
-----------------

The first state variable is the lumped cell temperature,

.. math::

   y_0 = T.

Its derivative is supplied by the thermal model,

.. math::

   \frac{dT}{dt}
   =
   \frac{
       \dot{Q}_{\mathrm{gen}}
       -
       \dot{Q}_{\mathrm{loss}}
   }{
       C_{\mathrm{th}}
   }.


Reaction States
---------------

Each thermal runaway reaction contributes one integrated conversion,

.. math::

   \alpha_i.

The corresponding derivative is

.. math::

   \frac{d\alpha_i}{dt}
   =
   k_i(T)
   f_i(\boldsymbol{\alpha}).


Gas States
----------

When gas generation is enabled, each gas species contributes one mole
state,

.. math::

   n_s.

Before vent opening,

.. math::

   \frac{dn_s}{dt}
   =
   \dot{n}_{s,\mathrm{gen}}.

After vent opening,

.. math::

   \frac{dn_s}{dt}
   =
   \dot{n}_{s,\mathrm{gen}}
   -
   \dot{n}_{s,\mathrm{vent}}.


Right-Hand-Side Evaluation
--------------------------

For each solver evaluation, LiionTR conceptually performs the following
sequence:

.. code-block:: text

   current state y
        |
        v
   extract T, alpha, gas moles
        |
        v
   build reaction context
        |
        v
   evaluate reaction rates
        |
        +---------> reaction progress rates
        |
        +---------> heat generation
        |
        v
   evaluate thermal model
        |
        +---------> dT/dt
        |
        v
   evaluate gas generation
        |
        +---------> dn_gen/dt
        |
        v
   evaluate pressure
        |
        v
   if vent open:
        evaluate vent flow
        |
        +---------> dn_vent/dt
        |
        v
   assemble dy/dt


Closed Phase
------------

The simulation begins with the vent closed.

The gas balance is

.. math::

   \frac{dn_s}{dt}
   =
   \dot{n}_{s,\mathrm{gen}}.

During this phase, internal pressure is continuously reconstructed from
the current temperature and gas inventory.


Vent-Opening Event
------------------

The solver monitors the residual

.. math::

   g_{\mathrm{vent}}
   =
   P_{\mathrm{vent}}
   -
   P.

The event occurs when

.. math::

   g_{\mathrm{vent}} = 0.

This corresponds to

.. math::

   P = P_{\mathrm{vent}}.


Open Phase
----------

After the vent-opening event, LiionTR starts a second integration phase.

The gas balance becomes

.. math::

   \frac{dn_s}{dt}
   =
   \dot{n}_{s,\mathrm{gen}}
   -
   \dot{n}_{s,\mathrm{vent}}.

The vent remains open for the remainder of the simulation.


Irreversible State Transition
-----------------------------

The solver architecture currently uses an irreversible transition

.. code-block:: text

   CLOSED
      |
      | P = P_vent
      v
    OPEN

There is no transition from ``OPEN`` back to ``CLOSED``.

This prevents numerical oscillation around the opening threshold and
represents a mechanical vent that does not reseal during the event.


Maximum-Temperature Event
-------------------------

A simulation may define a maximum temperature.

The solver monitors a temperature event residual that reaches zero when

.. math::

   T = T_{\mathrm{max}}.

This event may terminate integration.


Maximum-Pressure Event
----------------------

A simulation may also define a maximum pressure.

The corresponding event is triggered when

.. math::

   P = P_{\mathrm{max}}.

This threshold is distinct from the vent-opening pressure.


Stiff Integration
-----------------

Thermal runaway kinetics can be highly stiff because Arrhenius reaction
rates may change by many orders of magnitude over a relatively small
temperature interval.

LiionTR therefore uses a BDF integration method by default.

BDF methods are implicit multistep schemes that are commonly used for
stiff systems of ordinary differential equations.


Solver Tolerances
-----------------

The solver exposes relative and absolute integration tolerances.

The relative tolerance controls errors proportional to the magnitude of
the state variables.

The absolute tolerance provides a lower error scale for states close to
zero.

Appropriate tolerances are particularly important when the state vector
contains quantities with different physical magnitudes, such as:

* temperature in hundreds of kelvin;
* conversion between zero and one;
* gas inventories that may initially be very small.


Result Reconstruction
---------------------

After numerical integration, the solver reconstructs simulation outputs
from the integrated state.

These may include:

* time;
* temperature;
* reaction conversions;
* gas-species mole amounts;
* internal pressure;
* vent state.

The two integration phases are combined into one continuous result when
vent opening occurs.


Physical Model Ownership
------------------------

The solver does not own the physical laws.

Instead, it orchestrates calls to:

* reaction models;
* chemistry backend;
* thermal model;
* gas-generation model;
* pressure model;
* vent-flow model.

This separation keeps numerical integration independent from individual
physical model implementations.


Benefits of the Current Structure
---------------------------------

The current solver architecture supports:

* independent unit testing of physical models;
* event-driven state transitions;
* interchangeable physics components;
* future alternative numerical solvers;
* clear ownership of integrated state;
* straightforward addition of new state variables.


Future Extensions
-----------------

Potential solver extensions include:

* analytical or sparse Jacobians;
* solver-specific state scaling;
* adaptive event handling for rupture;
* phase-change state variables;
* dynamic vent geometry;
* spatial thermal states;
* propagation between cells;
* thermochemical equilibrium states;
* coupling to PyBaMM electrochemical variables;
* restart and checkpoint support;
* sensitivity analysis and parameter estimation.


Implementation
--------------

The SciPy solver is implemented by

:class:`liiontr.solver.scipy_solver.ScipySolver`.

The thermal problem definition is

:class:`liiontr.problems.thermal.ThermalProblem`.

Simulation results are represented by

:class:`liiontr.core.results.Results`.