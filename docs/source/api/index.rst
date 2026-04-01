API
===

cedar
-----
.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: class.rst

   cedar.Benchmark
   cedar.Field
   cedar.FieldView
   cedar.Log
   cedar.Material
   cedar.Mesh
   cedar.Mesh1D
   cedar.Mesh3D
   cedar.Model
   cedar.Problem
   cedar.Term
   cedar.Transfer

cedar.models
------------
.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: childclass.rst

   cedar.HeatTransfer
   cedar.HeatTransferBC
   cedar.AdiabaticBC
   cedar.KnownTemperatureBC
   cedar.KnownFluxBC
   cedar.ConvectiveBC
   cedar.HeatTransferSource
   cedar.HeatTransferTotalSource
   cedar.HeatTransferVolumetricSource
   cedar.Flow
   cedar.FlowInletBC
   cedar.FlowOutletBC
   cedar.FlowSource
   cedar.FlowQdotSource
   cedar.FlowHeatFluxSource
   
cedar.materials
---------------
.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: childclass.rst

   cedar.materials.Be
   cedar.materials.BeO
   cedar.materials.ConstantMaterial
   cedar.materials.G348
   cedar.materials.U10Mo
   cedar.materials.UC_ZrC_C
   cedar.materials.UN
   cedar.materials.UO2
   cedar.materials.ZrC_C
   cedar.materials.ZrC
   
cedar.functions
---------------
.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: function.rst

   cedar.functions.format_computation_time
   cedar.functions.gmsh_to_vtkhdf
   cedar.functions.MAPE
   cedar.functions.print_vtkhdf
   cedar.functions.residual
   cedar.functions.sort3
   cedar.functions.tetra_vol
   cedar.functions.triangle_area
   cedar.functions.triangle_center
   cedar.functions.triangle_normal