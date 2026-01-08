Benchmarks
==========

.. toctree::
   :hidden:

   thmcase1.rst
   adiabatictransient0d.rst
   dirichletsteady1d.rst
   nafemstransient1d.rst

.. list-table::
   :header-rows: 1

   * - Model
     - Benchmark
     - Reference
     - Error
   * - Flow
     - :doc:`thmcase1`
     - MOOSE THM
     - 0.1%
   * - Thermal
     - :doc:`adiabatictransient0d`
     - Analytical
     - 0.0%
   * - Thermal
     - :doc:`dirichletsteady1d`
     - Analytical
     - 2.2%
   * - Thermal
     - :doc:`nafemstransient1d`
     - NAFEMS
     - 0.4%