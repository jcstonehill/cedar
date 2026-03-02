Materials
==========

All material property objects are created such that they will provide continuous
and positive values for all temperature dependent functions from 1e-12 K to
10,000 K. The properties are extended to such extreme bounds to aid in problem
convergence. For example, during iterations the temperature might move outside
of the correlation validity bounds, even though the final solution is inside of
the validity bounds.

.. note::
   When a thermophysical property is requested at temperature outside the validity
   bounds for that property, the property value at the bound is returned instead.

   For example, if the validity bounds for a property are from 300 K to 1000 K, 
   but the property value at 1100 K is requested, the property value at 1000 K
   is returned.

Once the solution has converged, all material property objects are checked to
see if the solution temperature violates any of the correlations used to achieve
the solution. If so, a warning is provided to the user.

Below are all of the material property objects currently available in Cedar.

Properties
----------

All material objects include the following properties:

.. list-table::
   :header-rows: 1

   * - Property
     - Method
     - Units
     - Notes
   * - Room Temperature Density
     - `Material.rho_rt()`
     - :math:`\frac{kg}{m^3}`
     - 
   * - Thermal Conductivity
     - `Material.k(T)`
     - :math:`\frac{W}{m-K}`
     - Temperature-Dependent
   * - Specific Heat Capacity
     - `Material.cp(T)`
     - :math:`\frac{J}{kg-K}`
     - Temperature-Dependent

Available Materials
-------------------

.. list-table::
   :header-rows: 1

   * - Name
     - Formula
     - API
     - Ref
   * - :doc:`be`
     - Be
     - :class:`cedar.materials.Be`
     - [1]_
   * - :doc:`beo`
     - BeO
     - :class:`cedar.materials.BeO`
     - [1]_
   * - :doc:`constantmaterial`
     - Constant
     - :class:`cedar.materials.ConstantMaterial`
     - 
   * - :doc:`g348`
     - G348
     - :class:`cedar.materials.G348`
     - [2]_   [3]_
   * - :doc:`u10mo`
     - U10Mo
     - :class:`cedar.materials.U10Mo`
     - [4]_
   * - :doc:`uc_zrc_c`
     - (U, Zr)C-C
     - :class:`cedar.materials.uc_zrc_c`
     - [5]_
   * - :doc:`un`
     - UN
     - :class:`cedar.materials.un`
     - [1]_
   * - :doc:`uo2`
     - UO2
     - :class:`cedar.materials.uo2`
     - [6]_   [7]_
   * - :doc:`yh188`
     - YH188
     - :class:`cedar.materials.yh188`
     - [8]_
   * - :doc:`zrc_c`
     - ZrC-C
     - :class:`cedar.materials.zrc_c`
     - [1]_
   * - :doc:`zrc`
     - ZrC
     - :class:`cedar.materials.zrc`
     - [1]_
   
References
----------

.. [1] SNP-HDBK-0008 SNP Material Property Handbook
      https://ntrs.nasa.gov/citations/20240004217

.. [2] D. M. McEligot, et. al., "Thermal Properties of
      G-348 Graphite,” Idaho National Laboratory, Report
      INL/EXT-16-38241, Idaho Falls, Idaho, May. 2016.

.. [3] A. T. D. Butland, et. al., "The Specific Heat of
      Graphite: An Evaluation of Measurements," Journal of Nuclear
      Materials 49 (1973/74) 45-46, Jun. 1973.

.. [4] D. E. Burkes, et. al., "Thermophysical Properties of U-10Mo
      Alloy," Idaho National Laboratory, Report INL/EXT-10-19373, Idaho Falls,
      ID, Nov. 2010.

.. [5] L. L. Lyon, “Performance of (U, Zr)C-Graphite (Composite) and of (U, Zr)C
      (Carbide) Fuel Elements in the Nuclear Furance 1 Rest Reactor,” Los Alamos
      National Laboratory, Report LA-5398-MS Vol 1, Los Alamos, NM, Sept. 1973.

.. [6] J. V. Miller, “Estimating Thermal Conductivity of CERMET Fuel Materials
      for Nuclear Reactor Application,” Lewis Research Center, Report NASA TND-3898,
      Cleveland, OH, Apr. 1967. https://ntrs.nasa.gov/citations/19670013537

.. [7] S. G. Popov, et. al., "Thermophysical Properties of MOX and UO2 Fuel
      Including the Effects of Radiation," Oak Ridge National Laboratory, Report
      ORNL/TM-2000/351, Oak Ridge, TN, Nov. 2000.

.. [8] X. Hu, et. al., “Handbook on the Material Properties of Yttrium Hydride
      for High-Temperature Moderator Applications,” Oak Ridge National Laboratory,
      Report ORNL/TM-2021/2052, Oak Ridge, TN, Jun. 2021.

.. |space| unicode:: 32