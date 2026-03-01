Benchmarks
==========

.. toctree::
  :hidden:

  thmcase1.rst
  adiabatictransient0d.rst
  dirichletsteady1d.rst
  nafemstransient1d.rst

Flow
----

.. list-table::
   :header-rows: 1

   * - Benchmark
     - Reference
     - Error
   * - :doc:`thmcase1`
     - MOOSE THM
     - 0.1%

Heat Transfer
-------------

.. list-table::
   :header-rows: 1

   * - Benchmark
     - Reference
     - Error
   * - :doc:`adiabatictransient0d`
     - Analytical
     - 0.0%
   * - :doc:`dirichletsteady1d`
     - Analytical
     - 2.2%
   * - :doc:`nafemstransient1d`
     - NAFEMS
     - 0.4%