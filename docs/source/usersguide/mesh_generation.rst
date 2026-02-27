Mesh generation
===============

Definitions
-----------

**Boundary**: A collection of faces that are connected to only one cell. These
boundaries must have a boundary condition applied to close the problem.

**Cell**: A small, discrete control volume over which the heat conduction equation
is integrated and solved. Typically, this is a tetrahedron.

**Face**: A surface connecting two neighboring cells. Heat fluxes are computed and
exchanged between cells.

**Mesh**: A collection of cells and faces that represent the computational domain of
the problem.

**Region**: A collection of cells. Typically, regions are created to distinguish
between materials or separated areas of the computational domain.