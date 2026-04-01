Overview
========

.. toctree::
  :maxdepth: 1
  :titlesonly:
  :hidden:

  tutorial.rst
  theory.rst
  numerical_implementation.rst
  verification/index.rst

The :class:`cedar.HeatTransfer` model solves the 3D transient heat transfer equation.

The theory and implementation are based heavily on [Mazumder_2016]_.

.. math::
    \rho(\vec{r}) c_p(T, \vec{r}) \frac{\partial T}{\partial t} - \nabla (k(T, \vec{r}) \nabla T) = \dot{q}'''(t)

Where

| :math:`\rho` is mass density :math:`[\frac{kg}{m^3}]`
| :math:`\vec{r}` is spatial position vector :math:`[m]`
| :math:`c_p` is specific heat capacity :math:`[\frac{J}{kg K}]`
| :math:`T` is temperature :math:`[K]`
| :math:`t` is time :math:`[s]`
| :math:`k` is thermal conductivity :math:`[\frac{W}{m K}]`
| :math:`\dot{q}'''` is volumetric internal heat source :math:`[\frac{W}{m^3}]`

.. rubric:: Fields
  :heading-level: 2

.. list-table:: 
   :header-rows: 1
   :align: center

   * - Field
     - Description
     - Units
     - Type
     - Requires IC
   * - T
     - Temperature
     - :math:`K`
     - Cell
     - X
   * - Qgen
     - Heat Generation
     - :math:`W`
     - Cell
     - 
   * - volQgen
     - Volumetric Heat Generation
     - :math:`\frac{W}{m^3}`
     - Cell
     - 
   * - k
     - Thermal Conductivity
     - :math:`\frac{W}{m-K}`
     - Cell
     - 
   * - cp
     - Specific Heat Capacity
     - :math:`\frac{J}{kg-K}`
     - Cell
     - 
   * - T_wall
     - Boundary temperature
     - :math:`K`
     - Boundary Face
     - 
   * - J
     - Boundary heat flux
     - :math:`\frac{W}{m^2}`
     - Boundary Face
     - 
   * - Qdot
     - Boundary heat transfer
     - :math:`W`
     - Boundary Face
     - 


.. rubric:: Boundary Conditions
  :heading-level: 2

.. list-table:: 
   :header-rows: 1
   :align: center

   * - Condition
     - Description
     - Transfer Targets
   * - :class:`cedar.AdiabaticBC`
     - Enforces zero normal heat flux across the boundary.
     - N/A
   * - :class:`cedar.KnownTemperatureBC`
     - Enforces a prescribed boundary temperature.
     - T
   * - :class:`cedar.KnownFluxBC`
     - Enforces a prescribed heat flux across the boundary.
     - N/A
   * - :class:`cedar.ConvectiveBC`
     - Enforces a heat flux based on convection to a fluid.
     - T_flow, htc

.. rubric:: Sources
  :heading-level: 2

.. list-table:: 
   :header-rows: 1
   :align: center

   * - Source
     - Description
     - Transfer Targets
   * - :class:`cedar.HeatTransferTotalSource`
     - Total heating value is distributed using optional shape function.
     - N/A
   * - :class:`cedar.HeatTransferVolumetricSource`
     - Volumetric heating is applied to cells.
     - N/A

.. rubric:: References
  :heading-level: 2

.. [Mazumder_2016] S. Mazumder, Numerical Methods for Partial Differential
       Equations: Finite Difference and Finite Volume Methods. 2016.

.. |space| unicode:: 32