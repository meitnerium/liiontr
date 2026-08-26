Gas Pressure Model
==================

Overview
--------

LiionTR currently calculates the internal cell pressure using an
ideal-gas model.

The pressure model is coupled to the gas inventory and uses the gas
temperature, free volume, and total gas amount to determine the
absolute internal pressure.

The current formulation is deliberately simple and provides a
consistent pressure state for vent-opening and pressure-termination
events.


Ideal-Gas Equation
------------------

The pressure model is based on the ideal-gas relation

.. math::

   P V
   =
   n R T,

where

.. list-table::
   :header-rows: 1
   :widths: 20 45 20

   * - Symbol
     - Description
     - SI unit

   * - :math:`P`
     - Absolute gas pressure
     - Pa

   * - :math:`V`
     - Gas free volume
     - m\ :sup:`3`

   * - :math:`n`
     - Total amount of gas
     - mol

   * - :math:`R`
     - Universal gas constant
     - J mol\ :sup:`-1` K\ :sup:`-1`

   * - :math:`T`
     - Absolute gas temperature
     - K


Free Volume
-----------

The pressure model uses a prescribed free volume

.. math::

   V_{\mathrm{free}}.

This volume represents the region available to the gas phase inside the
cell or modeled containment volume.

The current implementation assumes that this free volume remains
constant during the simulation.

Thus,

.. math::

   V
   =
   V_{\mathrm{free}}.


Initial Gas Inventory
---------------------

The initial gas amount is determined from the specified initial
pressure, temperature, and free volume.

Using the ideal-gas equation,

.. math::

   n_0
   =
   \frac{
       P_0 V_{\mathrm{free}}
   }{
       R T_0
   },

where :math:`P_0` is the initial absolute pressure and :math:`T_0` is
the initial gas temperature.

This initial inventory may represent air, inert gas, vapor, or another
effective gas phase depending on the model configuration.


Pressure from Generated Gas
---------------------------

For a closed system, the total gas amount may be expressed as

.. math::

   n_{\mathrm{tot}}
   =
   n_0
   +
   n_{\mathrm{generated}}.

The pressure is then

.. math::

   P
   =
   \frac{
       n_{\mathrm{tot}} R T
   }{
       V_{\mathrm{free}}
   }.

Substituting the total gas amount gives

.. math::

   P
   =
   \frac{
       \left(
           n_0 + n_{\mathrm{generated}}
       \right)
       R T
   }{
       V_{\mathrm{free}}
   }.

This formulation allows both gas generation and temperature increase
to raise the internal pressure.


Pressure from Total Gas Amount
------------------------------

During venting, the gas inventory is no longer equal to the initial
inventory plus cumulative generation because gas is simultaneously
removed from the system.

The solver therefore uses the instantaneous total gas amount,

.. math::

   n_{\mathrm{tot}}
   =
   \sum_s n_s,

and evaluates

.. math::

   P
   =
   \frac{
       n_{\mathrm{tot}} R T
   }{
       V_{\mathrm{free}}
   }.

This formulation is used during the coupled gas-generation and venting
calculation.


Temperature Coupling
--------------------

The pressure model currently assumes that the gas temperature is equal
to the cell temperature,

.. math::

   T_{\mathrm{gas}}
   =
   T_{\mathrm{cell}}.

Consequently, internal pressure may rise for two independent reasons:

* the total gas amount increases;
* the gas temperature increases.

For constant volume and fixed gas amount,

.. math::

   P \propto T.

For constant temperature and fixed volume,

.. math::

   P \propto n.


Closed-System Pressure Evolution
--------------------------------

Before vent opening,

.. math::

   \frac{dn_{\mathrm{tot}}}{dt}
   =
   \sum_s
   \dot{n}_{s,\mathrm{gen}}.

The pressure therefore evolves according to both the thermal and gas
generation states.

Conceptually,

.. math::

   P(t)
   =
   \frac{
       n_{\mathrm{tot}}(t)
       R
       T(t)
   }{
       V_{\mathrm{free}}
   }.


Pressure During Venting
-----------------------

After the vent opens, gas generation and discharge occur
simultaneously.

For each species,

.. math::

   \frac{dn_s}{dt}
   =
   \dot{n}_{s,\mathrm{gen}}
   -
   \dot{n}_{s,\mathrm{vent}}.

The total remaining gas amount is therefore time dependent and may
either increase or decrease depending on the relative magnitudes of
generation and discharge.

The instantaneous pressure continues to follow

.. math::

   P
   =
   \frac{
       n_{\mathrm{tot}} R T
   }{
       V_{\mathrm{free}}
   }.


Vent-Opening Criterion
----------------------

The internal pressure may be compared against a prescribed vent-opening
pressure,

.. math::

   P_{\mathrm{vent}}.

The opening event occurs when

.. math::

   P
   =
   P_{\mathrm{vent}}.

Before this event, the vent is treated as closed.

After the event, LiionTR switches permanently to the open-vent
governing equations.


Maximum-Pressure Event
----------------------

A thermal problem may also define a maximum allowable pressure,

.. math::

   P_{\mathrm{max}}.

The solver can terminate integration when

.. math::

   P
   =
   P_{\mathrm{max}}.

This event is independent of the vent-opening threshold.

It can therefore be used to represent a numerical safety limit,
containment failure criterion, or other pressure-based termination
condition.


Absolute Pressure
-----------------

All pressures used by the current gas model are absolute pressures.

This distinction is important for the ideal-gas equation and
compressible-flow calculations.

Atmospheric pressure may be represented approximately as

.. math::

   P_{\mathrm{atm}}
   \approx
   101325
   \ \mathrm{Pa}.

Gauge pressure is not used internally by the present implementation.


Units
-----

LiionTR uses SI units internally.

.. list-table::
   :header-rows: 1
   :widths: 45 30

   * - Quantity
     - Unit

   * - Pressure
     - Pa

   * - Free volume
     - m\ :sup:`3`

   * - Gas amount
     - mol

   * - Temperature
     - K

   * - Universal gas constant
     - J mol\ :sup:`-1` K\ :sup:`-1`


Current Assumptions
-------------------

The present pressure model assumes that:

* the gas behaves ideally;
* the gas phase is spatially uniform;
* the gas temperature equals the cell temperature;
* the free volume remains constant;
* all gas species share a common temperature and pressure;
* liquid-vapor equilibrium is neglected;
* condensation is neglected;
* gas dissolution in condensed phases is neglected;
* compressibility effects are represented only in the vent-flow model,
  not in the pressure equation itself.


Limitations
-----------

The ideal-gas approximation is suitable for reduced-order modeling but
may become inaccurate at elevated pressure, for strongly non-ideal gas
mixtures, or when condensable species are important.

A real battery undergoing thermal runaway may also experience changes
in available free volume because of swelling, separator collapse,
electrode deformation, liquid displacement, or vent-component motion.

These effects are not currently represented.


Future Extensions
-----------------

Future pressure-model extensions may include:

* variable free volume;
* real-gas equations of state;
* gas-liquid phase equilibrium;
* vaporization and condensation;
* pressure-dependent thermodynamic properties;
* coupling to mechanical deformation;
* rupture and structural failure models;
* thermochemical properties supplied directly by Cantera or another
  equilibrium backend.


Implementation
--------------

The ideal-gas pressure model is implemented by

:class:`liiontr.gases.ideal.IdealGasPressureModel`.

The gas inventory supplying the total mole amount is represented by

:class:`liiontr.gases.inventory.GasInventory`.

Pressure-based events are handled by

:class:`liiontr.solver.scipy_solver.ScipySolver`.

The pressure model is configured through

:class:`liiontr.problems.thermal.ThermalProblem`.