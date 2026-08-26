Architecture Overview
=====================

Design Goals
------------

LiionTR is designed as a modular scientific framework for lithium-ion
battery thermal runaway modeling.

The architecture separates:

* cell geometry and effective material properties;
* reaction kinetics and progress laws;
* thermal runaway reaction networks;
* heat generation;
* gas generation;
* gas inventory and pressure;
* vent-flow models;
* numerical integration;
* optional thermochemical and electrochemical backends.

This separation allows individual physical models to be replaced or
extended without rewriting the complete simulation framework.


High-Level Data Flow
--------------------

The current reduced-order architecture can be represented as

.. code-block:: text

   Cell + Material
          |
          v
   Reaction Network
      |         |
      |         +----------------> Heat Generation
      |                                  |
      |                                  v
      |                            Thermal Model
      |                                  |
      |                                  v
      +----------------------------> Temperature
      |
      v
   Gas Generation
      |
      v
   Gas Inventory
      |
      v
   Pressure Model
      |
      v
   Vent Flow
      |
      +----------------------------> Gas Inventory

                    ScipySolver
             integrates all state variables


Cells
-----

A battery cell combines three main concepts:

* geometry;
* effective material properties;
* chemistry.

The cell provides quantities required by other models, including:

* volume;
* surface area;
* mass;
* thermal capacity.

The cell itself does not perform thermal runaway calculations.


Geometry
--------

Geometry models provide geometric quantities independently from material
or chemistry models.

For cylindrical cells, LiionTR currently evaluates quantities such as

.. math::

   V
   =
   \pi r^2 h

and the external surface area required by the lumped thermal model.

This separation makes it possible to add prismatic or pouch-cell
geometries without modifying the thermal or reaction subsystems.


Materials
---------

Material models provide effective thermophysical properties such as:

* density;
* specific heat capacity;
* thermal conductivity.

The current implementation supports constant effective properties.

These properties are intentionally represented independently from the
cell geometry.


Reaction Kinetics
-----------------

Kinetic models determine the intrinsic temperature dependence of
thermal runaway reactions.

A common example is the Arrhenius relation

.. math::

   k(T)
   =
   A
   \exp\left(
       -\frac{E_a}{RT}
   \right).

The kinetic model does not contain the complete physical reaction
definition.


Progress Models
---------------

Progress models describe how reaction state modifies the kinetic rate.

For example,

.. math::

   f(\alpha)
   =
   (1-\alpha)^n.

Progress models may also depend on the state of other reactions through
a reaction context.

Separating kinetics from progress makes the same kinetic law reusable
with multiple physical reaction formulations.


Reaction Models
---------------

A reaction combines:

* a kinetic model;
* a progress model;
* reaction enthalpy;
* reacting-material mass fraction.

The reaction provides:

* a conversion rate;
* a mass-specific heat-generation rate.

The reaction object does not integrate its own state.


Reaction Network
----------------

A reaction network groups multiple reaction models and evaluates them
using a shared reaction state.

The network provides:

* reaction progress rates;
* heat generation by reaction;
* total heat generation;
* context information for coupled reactions.

The numerical solver owns the actual conversion vector.


Chemistry Backend
-----------------

The chemistry backend separates the physical problem from the mechanism
used to calculate heat generation.

The current reaction-network backend maps a
:class:`liiontr.reactions.network.ReactionNetwork` to the generic
chemistry interface.

This design leaves room for additional chemistry implementations.


Thermal Model
-------------

The thermal model converts net heat flow into a temperature derivative.

For the current lumped model,

.. math::

   \frac{dT}{dt}
   =
   \frac{
       \dot{Q}_{\mathrm{gen}}
       -
       hA(T-T_{\infty})
   }{
       C_{\mathrm{th}}
   }.

The thermal model does not calculate reaction rates.


Gas Generation
--------------

The current gas-generation subsystem maps reaction progress to
prescribed gas-species generation rates.

It uses:

* reaction progress rates;
* cell mass;
* reaction mass fraction;
* empirical species yields.

Gas generation is therefore coupled to the reaction network but remains
independent from pressure and vent-flow calculations.


Gas Inventory
-------------

The gas inventory stores the instantaneous mole amount of each species.

It provides derived quantities such as:

* total moles;
* total gas mass;
* mole fractions;
* mean molar mass.

The inventory forms the connection between gas generation, pressure,
and venting.


Pressure Model
--------------

The current pressure model applies the ideal-gas equation

.. math::

   P
   =
   \frac{nRT}{V}.

It uses the instantaneous gas inventory and temperature but does not
modify those states directly.


Vent Flow
---------

The vent-flow subsystem calculates gas discharge when the internal
pressure exceeds the configured vent-opening threshold.

The model consumes:

* internal pressure;
* downstream pressure;
* temperature;
* gas-mixture properties.

It returns vent molar-flow rates that become sink terms in the gas
species balances.


Thermal Problem
---------------

A thermal problem acts as the simulation configuration container.

It assembles the models and conditions required by the solver,
including:

* cell;
* thermal model;
* chemistry backend;
* initial temperature;
* reaction conversions;
* gas models;
* pressure model;
* vent model;
* event thresholds;
* simulation time span.

The problem definition does not itself perform numerical integration.


Numerical Solver
----------------

The numerical solver is responsible for advancing the complete physical
state.

The current implementation uses SciPy ``solve_ivp`` with a BDF method.

The solver evaluates the configured physical models to construct the
right-hand side of the governing ODE system.


Separation of State and Physics
-------------------------------

An important architectural principle in LiionTR is that physical model
objects do not mutate their own integration state.

Instead, the solver stores state variables such as

.. math::

   T,
   \alpha_1,
   \ldots,
   \alpha_N,
   n_1,
   \ldots,
   n_M.

Models receive the current state and return rates or derived
properties.

This design has several advantages:

* deterministic right-hand-side evaluation;
* compatibility with standard ODE solvers;
* easier testing;
* reduced hidden coupling;
* simpler restart and event handling;
* easier future implementation of alternative solvers.


Optional Scientific Backends
----------------------------

LiionTR is intended to support external scientific engines without
making them mandatory dependencies.

Cantera
~~~~~~~

Cantera is intended for future thermochemical calculations such as:

* equilibrium gas composition;
* mixture thermodynamic properties;
* heat capacities;
* chemical equilibrium;
* possibly energy-consistent gas chemistry.

PyBaMM
~~~~~~

PyBaMM may eventually provide electrochemical states before or during
abuse simulations, including:

* state of charge;
* electrode concentrations;
* voltage;
* current;
* electrochemical heat generation;
* degradation state.

LiionTR remains responsible for thermal runaway and safety-oriented
physics.


Intended Long-Term Architecture
-------------------------------

A future coupled workflow may conceptually follow

.. code-block:: text

   PyBaMM
      |
      | electrochemical state
      v
   LiionTR Reaction Models
      |
      +-----------> Heat
      |
      +-----------> Element Sources
                        |
                        v
                     Cantera
                        |
                        +----> Gas Composition
                        +----> Thermodynamic Properties
                        +----> Pressure Support
                        |
                        v
                     Venting

The external packages therefore complement LiionTR rather than replace
its orchestration layer.


Implementation Map
------------------

The major source packages are:

.. list-table::
   :header-rows: 1

   * - Package
     - Responsibility

   * - ``liiontr.cells``
     - Battery-cell definitions

   * - ``liiontr.geometry``
     - Cell geometry

   * - ``liiontr.materials``
     - Effective material properties

   * - ``liiontr.kinetics``
     - Kinetic and reaction-progress models

   * - ``liiontr.reactions``
     - Reactions and reaction networks

   * - ``liiontr.chemistry``
     - Chemistry backend interfaces

   * - ``liiontr.thermal``
     - Thermal models

   * - ``liiontr.gases``
     - Gas generation, pressure, inventory, and venting

   * - ``liiontr.problems``
     - Physical problem definitions

   * - ``liiontr.solver``
     - Numerical integration

   * - ``liiontr.library``
     - Reference cells and literature models

   * - ``liiontr.core``
     - Generic framework infrastructure