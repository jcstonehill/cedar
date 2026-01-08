API
===

cedar
-----
.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: class.rst

   cedar.FlowStateVar
   cedar.Log
   cedar.Mesh1D
   cedar.Mesh3D
   cedar.Mesh1DVar
   cedar.Mesh3DVar
   cedar.Problem
   cedar.ScalarVar

cedar.adapters
--------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: class.rst

   cedar.adapters.Direct
   cedar.adapters.NearestValue
   cedar.adapters.Summation

cedar.base
----------
.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: class.rst

   cedar.base.Adapter
   cedar.base.Benchmark
   cedar.base.Fluid
   cedar.base.Material
   cedar.base.Mesh
   cedar.base.Model
   cedar.base.Var

cedar.benchmarks
----------------
.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: childclass.rst

   cedar.benchmarks.flow.THMCase1
   cedar.benchmarks.thermal.AdiabaticTransient0D
   cedar.benchmarks.thermal.DirichletSteady1D
   cedar.benchmarks.thermal.NAFEMSTransient1D

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: function.rst

   cedar.benchmarks.plot
   cedar.benchmarks.run

cedar.fluids
---------------
.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: childclass.rst

   cedar.fluids.Hydrogen
   cedar.fluids.IdealGas
   cedar.fluids.Parahydrogen

cedar.materials
---------------
.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: childclass.rst

   cedar.materials.ConstantMaterial
   cedar.materials.UC_ZrC_C
   cedar.materials.ZrC
   cedar.materials.ZrC_C

cedar.models
------------
.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: childclass.rst

   cedar.models.Flow
   cedar.models.ScalarFromFunc
   cedar.models.Thermal