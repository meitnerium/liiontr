Thermal Model
=============

Overview
--------

LiionTR currently uses a zero-dimensional lumped-capacitance thermal
model for individual battery cells.

The model assumes that the temperature is spatially uniform throughout
the cell. Internal heat generation is balanced by heat transfer from the
cell surface to the surrounding environment.

The thermal model is intentionally separated from the reaction models.
Reaction kinetics determine the internal heat-generation rate, while the
thermal model converts the net heat flow into a cell temperature rate.


Energy Balance
--------------

The cell energy balance is

.. math::

   C_{\mathrm{th}} \frac{dT}{dt}
   =
   \dot{Q}_{\mathrm{gen}}
   -
   \dot{Q}_{\mathrm{loss}},

where

.. math::

   \dot{Q}_{\mathrm{loss}}
   =
   h A (T - T_{\infty}).

The resulting temperature evolution equation is therefore

.. math::

   \frac{dT}{dt}
   =
   \frac{
       \dot{Q}_{\mathrm{gen}}
       -
       h A (T - T_{\infty})
   }{
       C_{\mathrm{th}}
   }.


Definitions
-----------

The quantities appearing in the thermal balance are:

.. list-table::
   :header-rows: 1
   :widths: 20 45 20

   * - Symbol
     - Description
     - SI unit

   * - :math:`T`
     - Uniform cell temperature
     - K

   * - :math:`T_{\infty}`
     - Ambient temperature
     - K

   * - :math:`\dot{Q}_{\mathrm{gen}}`
     - Total internal heat-generation rate
     - W

   * - :math:`\dot{Q}_{\mathrm{loss}}`
     - Heat-loss rate to the surroundings
     - W

   * - :math:`h`
     - Convective heat-transfer coefficient
     - W m\ :sup:`-2` K\ :sup:`-1`

   * - :math:`A`
     - External cell surface area
     - m\ :sup:`2`

   * - :math:`C_{\mathrm{th}}`
     - Cell thermal capacity
     - J K\ :sup:`-1`


Thermal Capacity
----------------

LiionTR evaluates the lumped thermal capacity as

.. math::

   C_{\mathrm{th}}
   =
   m c_p,

where :math:`m` is the cell mass and :math:`c_p` is the effective
specific heat capacity of the cell material.

For the current cell implementation,

.. math::

   m
   =
   \rho V,

so that

.. math::

   C_{\mathrm{th}}
   =
   \rho V c_p.

The current implementation evaluates density and heat capacity at
298.15 K when calculating cell mass and thermal capacity.

This treatment therefore assumes constant effective thermophysical
properties during the simulation.


Heat Generation
---------------

The thermal model receives the total internal heat-generation rate

.. math::

   \dot{Q}_{\mathrm{gen}}

from the chemistry or reaction backend.

For a reaction-network model, the total heat-generation rate may be
written conceptually as

.. math::

   \dot{Q}_{\mathrm{gen}}
   =
   m_{\mathrm{cell}}
   \sum_i
   \dot{q}_i,

where :math:`\dot{q}_i` is the mass-specific heat-generation rate of
reaction :math:`i`.

The thermal model itself does not determine reaction kinetics or
reaction enthalpies. This separation allows the same thermal model to be
used with different thermal-runaway reaction mechanisms.


Convective Heat Transfer
------------------------

Heat loss to the surroundings is currently represented using Newton's
law of cooling,

.. math::

   \dot{Q}_{\mathrm{conv}}
   =
   h A (T - T_{\infty}).

When

.. math::

   T > T_{\infty},

the term is positive and removes heat from the cell.

When

.. math::

   T < T_{\infty},

the same formulation produces a negative heat-loss term and therefore
represents heating of the cell by the environment.


Lumped-Capacitance Assumption
-----------------------------

The lumped model assumes that internal temperature gradients are
negligible and that the entire cell can be represented using a single
temperature state.

A conventional criterion for assessing this approximation is the Biot
number,

.. math::

   \mathrm{Bi}
   =
   \frac{h L_c}{k},

where :math:`L_c` is a characteristic length and :math:`k` is the
effective thermal conductivity of the cell.

The lumped-capacitance approximation is generally most appropriate when

.. math::

   \mathrm{Bi} \ll 1.

LiionTR does not currently enforce a Biot-number validity criterion.
Users are responsible for determining whether the lumped approximation
is appropriate for a particular cell geometry, heat-transfer condition,
and thermal-runaway regime.


Current Assumptions
-------------------

The present thermal model makes the following assumptions:

* the cell temperature is spatially uniform;
* effective material properties are constant;
* heat transfer to the environment is purely convective;
* the convection coefficient is uniform over the complete cell surface;
* the ambient temperature is constant;
* radiative heat transfer is neglected;
* conductive coupling to neighbouring cells or structures is neglected;
* phase changes are not represented explicitly;
* internal thermal gradients are not represented.


Limitations
-----------

The lumped approximation is useful for model development, parameter
studies, and coupling thermal runaway chemistry to gas and pressure
models, but it cannot reproduce spatially resolved temperature
distributions.

During severe thermal runaway, large internal temperature gradients may
develop because reaction rates can increase strongly with temperature.
For such conditions, a spatially resolved thermal model may eventually
be required.

Future extensions may therefore include:

* temperature-dependent thermophysical properties;
* radiative heat transfer;
* conductive heat transfer between cells;
* spatially resolved one-dimensional or multidimensional thermal models;
* thermal propagation between neighbouring cells;
* coupling to phase-change and venting enthalpy effects.


Implementation
--------------

The lumped thermal balance is implemented by

:class:`liiontr.thermal.lumped.LumpedThermalModel`.

The corresponding abstract interface is

:class:`liiontr.thermal.model.ThermalModel`.

The thermal model returns the temperature derivative in K/s to the
numerical solver.