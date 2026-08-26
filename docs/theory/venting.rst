Vent Flow Model
===============

Overview
--------

LiionTR models gas discharge from a battery cell using a compressible
ideal-gas vent-flow model.

The flow rate depends on the upstream cell pressure, downstream ambient
pressure, gas temperature, vent area, gas molar mass, discharge
coefficient, and heat-capacity ratio.

Two flow regimes are represented:

* choked flow;
* unchoked compressible flow.

The model does not permit reverse flow into the cell.


Pressure Definitions
--------------------

The vent model distinguishes between the upstream pressure

.. math::

   P_u

and downstream pressure

.. math::

   P_d.

In the current battery application,

.. math::

   P_u
   =
   P_{\mathrm{cell}}

and the downstream pressure typically represents atmospheric pressure,

.. math::

   P_d
   =
   P_{\mathrm{ambient}}.

Both pressures are absolute pressures expressed in Pa.


Pressure Ratio
--------------

The downstream-to-upstream pressure ratio is

.. math::

   r_p
   =
   \frac{P_d}{P_u}.

This ratio determines whether the flow is choked or unchoked.


Specific Gas Constant
---------------------

The vent-flow equations require the specific gas constant of the gas
mixture,

.. math::

   R_s
   =
   \frac{R}{M},

where

.. list-table::
   :header-rows: 1
   :widths: 20 45 20

   * - Symbol
     - Description
     - SI unit

   * - :math:`R_s`
     - Specific gas constant
     - J kg\ :sup:`-1` K\ :sup:`-1`

   * - :math:`R`
     - Universal gas constant
     - J mol\ :sup:`-1` K\ :sup:`-1`

   * - :math:`M`
     - Gas molar mass
     - kg mol\ :sup:`-1`

For a gas mixture, LiionTR currently uses the mole-weighted mean molar
mass of the gas inventory.


Heat-Capacity Ratio
-------------------

The compressible-flow model uses the heat-capacity ratio

.. math::

   \gamma
   =
   \frac{c_p}{c_v}.

The current model treats :math:`\gamma` as a prescribed constant.

Its value affects both the critical pressure ratio and the predicted
mass flow rate.


Critical Pressure Ratio
-----------------------

For ideal-gas isentropic flow, the critical downstream-to-upstream
pressure ratio is

.. math::

   r_{\mathrm{crit}}
   =
   \left(
       \frac{2}{\gamma + 1}
   \right)^{
       \frac{\gamma}{\gamma - 1}
   }.

The flow is treated as choked when

.. math::

   \frac{P_d}{P_u}
   \leq
   r_{\mathrm{crit}}.

When this condition is satisfied, the gas velocity reaches the local
speed of sound at the restrictive section of the vent.


Choked Flow
-----------

For choked flow, the mass discharge rate is calculated as

.. math::

   \dot{m}
   =
   C_d A_v P_u
   \sqrt{
       \frac{\gamma}
            {R_s T}
   }
   \left(
       \frac{2}
            {\gamma + 1}
   \right)^{
       \frac{\gamma + 1}
            {2(\gamma - 1)}
   },

where

.. list-table::
   :header-rows: 1
   :widths: 20 45 20

   * - Symbol
     - Description
     - SI unit

   * - :math:`\dot{m}`
     - Gas mass flow rate
     - kg/s

   * - :math:`C_d`
     - Discharge coefficient
     - dimensionless

   * - :math:`A_v`
     - Effective vent area
     - m\ :sup:`2`

   * - :math:`P_u`
     - Upstream absolute pressure
     - Pa

   * - :math:`T`
     - Upstream gas temperature
     - K

   * - :math:`R_s`
     - Specific gas constant
     - J kg\ :sup:`-1` K\ :sup:`-1`

   * - :math:`\gamma`
     - Heat-capacity ratio
     - dimensionless


Physical Meaning of Choking
---------------------------

Once the flow becomes choked, decreasing the downstream pressure
further does not increase the mass flow rate predicted by the
isentropic nozzle relation.

The discharge rate is then controlled primarily by the upstream
thermodynamic state and the effective vent geometry.

This behavior is particularly relevant during rapid battery venting
when the cell pressure is substantially greater than ambient pressure.


Unchoked Flow
-------------

When

.. math::

   \frac{P_d}{P_u}
   >
   r_{\mathrm{crit}},

the flow is treated as unchoked.

The mass flow rate is calculated using

.. math::

   \dot{m}
   =
   C_d A_v P_u
   \sqrt{
       \frac{
           2\gamma
       }{
           R_s T
           (\gamma - 1)
       }
       \left[
           \left(
               \frac{P_d}{P_u}
           \right)^{
               \frac{2}{\gamma}
           }
           -
           \left(
               \frac{P_d}{P_u}
           \right)^{
               \frac{\gamma+1}{\gamma}
           }
       \right]
   }.

As the upstream and downstream pressures approach one another, the
predicted mass flow rate approaches zero.


No Reverse Flow
---------------

The current model does not simulate gas entering the battery through
the vent.

Therefore, when

.. math::

   P_u
   \leq
   P_d,

the discharge rate is taken to be zero,

.. math::

   \dot{m}
   =
   0.

This assumption is appropriate for the current thermal-runaway venting
model, which focuses on outward discharge following internal
pressurization.


Mass-to-Molar Flow Conversion
-----------------------------

The compressible-flow equations produce a mass flow rate.

LiionTR converts this value to total molar flow using

.. math::

   \dot{n}_{\mathrm{vent}}
   =
   \frac{\dot{m}}
        {M},

where :math:`M` is the gas or mixture molar mass.

The resulting molar flow rate has units of mol/s.


Mixture Vent Flow
-----------------

The current mixture model assumes that the gas inside the cell free
volume is perfectly mixed.

The total vent molar flow rate is calculated using the mean molar mass

.. math::

   \overline{M}
   =
   \sum_s x_s M_s.

Thus,

.. math::

   \dot{n}_{\mathrm{vent}}
   =
   \frac{
       \dot{m}_{\mathrm{vent}}
   }{
       \overline{M}
   }.


Species Vent Rates
------------------

The total molar flow is distributed among gas species according to the
current mole fractions.

For species :math:`s`,

.. math::

   \dot{n}_{s,\mathrm{vent}}
   =
   x_s
   \dot{n}_{\mathrm{vent}}.

Consequently,

.. math::

   \sum_s
   \dot{n}_{s,\mathrm{vent}}
   =
   \dot{n}_{\mathrm{vent}}.

This treatment assumes that no species is preferentially retained or
discharged through the vent.


Species Balance During Venting
------------------------------

After vent opening, each species follows

.. math::

   \frac{dn_s}{dt}
   =
   \dot{n}_{s,\mathrm{gen}}
   -
   \dot{n}_{s,\mathrm{vent}}.

The pressure is recalculated from the remaining gas inventory after the
generation and discharge terms are applied.


Vent Opening
------------

The vent remains closed until the internal pressure reaches the
configured vent-opening pressure,

.. math::

   P_{\mathrm{cell}}
   =
   P_{\mathrm{vent}}.

The numerical solver treats this condition as an event.

Before this event,

.. math::

   \dot{n}_{s,\mathrm{vent}}
   =
   0.

After the event, the venting terms are enabled.


Irreversible Vent State
-----------------------

LiionTR currently treats vent opening as irreversible.

The simulation therefore consists conceptually of two phases.

Closed phase:

.. math::

   \frac{dn_s}{dt}
   =
   \dot{n}_{s,\mathrm{gen}}.

Open phase:

.. math::

   \frac{dn_s}{dt}
   =
   \dot{n}_{s,\mathrm{gen}}
   -
   \dot{n}_{s,\mathrm{vent}}.

Even if the cell pressure subsequently decreases below the original
vent-opening pressure, the vent remains open.


Coupling to the Numerical Solver
--------------------------------

The SciPy solver integrates the closed system until the vent-opening
event is detected.

A second integration phase is then started using the open-vent
governing equations.

Conceptually,

.. math::

   \text{closed RHS}
   \rightarrow
   P = P_{\mathrm{vent}}
   \rightarrow
   \text{open RHS}.

The two integration segments are combined into a single simulation
result.


Coupling to Gas Generation
--------------------------

Gas generation and vent discharge occur simultaneously after vent
opening.

A thermal runaway reaction can therefore continue increasing the gas
inventory even while material is being discharged.

Whether the pressure rises or falls depends on the competition between

.. math::

   \dot{n}_{\mathrm{generation}}

and

.. math::

   \dot{n}_{\mathrm{vent}},

as well as on the simultaneous evolution of temperature.


Coupling to Temperature
-----------------------

Temperature affects vent flow through the compressible-flow equations.

For a fixed pressure state, the mass-flow expressions contain a factor
proportional to

.. math::

   \frac{1}{\sqrt{T}}.

However, temperature also affects internal pressure and reaction-driven
gas generation.

The complete coupled response therefore cannot generally be inferred
from this direct temperature dependence alone.


Units
-----

LiionTR uses SI units internally.

.. list-table::
   :header-rows: 1
   :widths: 45 30

   * - Quantity
     - Unit

   * - Absolute pressure
     - Pa

   * - Temperature
     - K

   * - Vent area
     - m\ :sup:`2`

   * - Mass flow rate
     - kg/s

   * - Molar flow rate
     - mol/s

   * - Molar mass
     - kg/mol

   * - Specific gas constant
     - J kg\ :sup:`-1` K\ :sup:`-1`

   * - Heat-capacity ratio
     - dimensionless

   * - Discharge coefficient
     - dimensionless


Current Assumptions
-------------------

The present vent model assumes that:

* the gas behaves ideally;
* the gas mixture is spatially uniform;
* the upstream gas has a single temperature;
* flow through the vent is quasi-steady;
* the vent is represented by a fixed effective area;
* the discharge coefficient is constant;
* the heat-capacity ratio is constant;
* the mean gas molar mass adequately represents mixture flow;
* the gas composition discharged through the vent equals the internal
  gas composition;
* reverse flow is neglected;
* vent opening is irreversible;
* vent geometry does not evolve after opening;
* multiphase discharge is neglected.


Limitations
-----------

Real lithium-ion battery venting can be substantially more complex than
single-phase ideal-gas nozzle flow.

Battery vent streams may contain:

* permanent gases;
* electrolyte vapor;
* liquid droplets;
* aerosols;
* condensed particles;
* electrode or separator decomposition products.

The vent itself may also open progressively, deform, rupture, or change
effective area during the event.

These effects are not represented by the current reduced-order model.


Future Extensions
-----------------

Potential improvements include:

* temperature-dependent mixture heat capacities;
* dynamically calculated :math:`\gamma`;
* gas properties supplied by Cantera;
* variable vent area;
* explicit vent-opening dynamics;
* two-phase or multiphase discharge;
* flashing of liquid electrolyte;
* non-equilibrium gas composition;
* momentum and jet models outside the cell;
* coupling to fire and combustion models.


Implementation
--------------

Compressible gas discharge is implemented by

:class:`liiontr.gases.vent.CompressibleVentFlowModel`.

Gas-mixture discharge is implemented by

:class:`liiontr.gases.vent.MixtureVentFlowModel`.

The gas composition is supplied by

:class:`liiontr.gases.inventory.GasInventory`.

Pressure is evaluated by

:class:`liiontr.gases.ideal.IdealGasPressureModel`.

Vent opening and the transition between closed and open integration
phases are handled by

:class:`liiontr.solver.scipy_solver.ScipySolver`.