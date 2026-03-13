Verification
============

.. toctree::
  :hidden:

  adiabatictransient0d.rst
  dirichletsteady1d.rst
  nafemstransient1d.rst

Below are the verification cases for the :class:`cedar.HeatTransfer` model.

.. list-table::
   :header-rows: 1
   :align: center

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