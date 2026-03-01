Models
======

.. toctree::
   :hidden:

   heat_transfer/index.rst
   flow/index.rst

Below is a directory to technical documentation for all models currently included in Cedar.

Heat Transfer
-------------

.. list-table::

   * - .. figure:: img/heat_transfer.png
          :target: heat_transfer/index.html
          :scale: 25 %
          :align: center

          :doc:`View Documentation <heat_transfer/index>`

     - | Solves the transient heat conduction equation on a 3D mesh.
       | 
       | :math:`\rho(\vec{r}) c_p(T, \vec{r}) \frac{\partial T}{\partial t} - \nabla (k(T, \vec{r}) \nabla T) = \dot{q}'''`
       
Flow
-------------

.. list-table::

   * - .. figure:: img/heat_transfer.png
          :target: heat_transfer/index.html
          :scale: 25 %
          :align: center

          :doc:`View Documentation <flow/index>`

     - | Solves the compressible 1D area-averaged Euler flow equations.
       | 
       | :math:`\frac{\partial \rho A}{\partial t} + \frac{\partial \rho u A}{\partial x} = 0`
       | 
       | :math:`\frac{\partial \rho u A}{\partial t} + \frac{\partial (\rho u^2 + P) A}{\partial x} = P \frac{\partial A}{\partial x} - F A + \rho g_x A`
       | 
       | :math:`\frac{\partial \rho E A}{\partial t} + \frac{\partial u (\rho E + P) A}{\partial x} = \rho u g_x A + \dot{q}''' A`

Neutronics
-------------

.. list-table::

   * - .. figure:: img/neutronics.png
          :target: neutronics/index.html
          :scale: 25 %
          :align: center

          :doc:`View Documentation <neutronics/index>`

     - | Solves the neutron transport equation using OpenMC and the point reactor kinetics equations.
       | 
       | :math:`\frac{\partial \rho A}{\partial t} + \frac{\partial \rho u A}{\partial x} = 0`
       | 
       | :math:`\frac{\partial \rho u A}{\partial t} + \frac{\partial (\rho u^2 + P) A}{\partial x} = P \frac{\partial A}{\partial x} - F A + \rho g_x A`
       | 
       | :math:`\frac{\partial \rho E A}{\partial t} + \frac{\partial u (\rho E + P) A}{\partial x} = \rho u g_x A + \dot{q}''' A`

