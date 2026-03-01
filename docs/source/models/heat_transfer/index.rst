Heat Transfer
=============

Documentation
-------------

.. toctree::
   :maxdepth: 1

   tutorial
   theory
   numerical_implementation
   ../../benchmarks/index

Overview
--------

The :class:`cedar.HeatTransfer` model solves the heat conduction equation on a 3D mesh. 

The theory and implementation are based heavily on [Mazumder_2016]_.

.. math::
    \rho(\vec{r}) c_p(T, \vec{r}) \frac{\partial T}{\partial t} - \nabla (k(T, \vec{r}) \nabla T) = \dot{q}'''

Where

| :math:`\rho` is mass density :math:`[\frac{kg}{m^3}]`
| :math:`c_p` is specific heat capacity :math:`[\frac{J}{kg K}]`
| :math:`T` is temperature :math:`[K]`
| :math:`k` is thermal conductivity :math:`[\frac{W}{m K}]`
| :math:`\dot{q}'''` is volumetric internal heat source :math:`[\frac{W}{m^3}]`

References
----------

.. [Mazumder_2016] S. Mazumder, Numerical Methods for Partial Differential
       Equations: Finite Difference and Finite Volume Methods. 2016.

.. |space| unicode:: 32