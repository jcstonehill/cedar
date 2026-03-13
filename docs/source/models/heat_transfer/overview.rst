Overview
========

.. toctree::
   :maxdepth: 1
   :hidden:

   tutorial.rst
   theory.rst
   numerical_implementation.rst
   verification/index.rst

The :class:`cedar.HeatTransfer` model solves the heat transfer equation on a 3D mesh. 

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

Fields
------

.. list-table:: 
   :header-rows: 1
   :align: center

   * - Field
     - Description
     - Units
     - Requires IC
   * - T
     - Temperature
     - :math:`K`
     - Yes
   * - Qgen
     - Heat Generation
     - :math:`W`
     - No
   * - volQgen
     - Volumetric Heat Generation
     - :math:`\frac{W}{m^3}`
     - No
   * - k
     - Thermal Conductivity
     - :math:`\frac{W}{m-K}`
     - No
   * - cp
     - Specific Heat Capacity
     - :math:`\frac{J}{kg-K}`
     - No
   * - J
     - Boundary heat flux
     - :math:`\frac{W}{m^2}`
     - No
   * - Qdot
     - Boundary heat flux
     - :math:`W`
     - No


Boundary Conditions
-------------------

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

Sources
-------

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

References
----------

.. [Mazumder_2016] S. Mazumder, Numerical Methods for Partial Differential
       Equations: Finite Difference and Finite Volume Methods. 2016.

.. |space| unicode:: 32