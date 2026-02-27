Materials
==========

.. toctree::
   :maxdepth: 1

   be.rst
   beo.rst
   constant.rst
   g348.rst
   u10mo.rst
   uc_zrc_c.rst
   un.rst
   uo2.rst
   yh188.rst
   zrc_c.rst
   zrc.rst

All material property objects are created such that they will provide continuous
and positive values for all temperature dependent functions from 1e-12 K to
10,000 K. The properties are extended to such extreme bounds to aid in
convergence. For example, during iterations the temperature might move outside
of the correlation validity bounds, even though the final solution is inside of
the validity bounds.

Once the solution has converged, all material property objects are checked to
see if the solution temperature violates any of the correlations used to achieve
the solution. If so, a warning is provided to the user.

Below are all of the material property objects currently available in Cedar.