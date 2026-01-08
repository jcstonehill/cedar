Thermal
=======

The ``Thermal`` model solves the heat conduction equation on a 3D mesh. 


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

Theory
------

Heat conduction is the process by which thermal energy moves through a material
due to microscopic interactions between its atoms. When one area is hotter
than another, faster-moving particles in the hot area transfer energy to
neighboring, slower-moving particles through collisions or vibrations. This
results in a net flow of heat from high temperatures to low temperatures,
without any bulk motion of the material.

At the macroscopic level, this behavior is described by Fourier's law.

.. math::
    \rho(\vec{r}) c_p(T, \vec{r}) \frac{\partial T}{\partial t} - \nabla (k(T, \vec{r}) \nabla T) = q'''

Where

| :math:`\rho` is mass density :math:`[\frac{kg}{m^3}]`
| :math:`c_p` is specific heat capacity :math:`[\frac{J}{kg K}]`
| :math:`T` is temperature :math:`[K]`
| :math:`k` is thermal conductivity :math:`[\frac{W}{m K}]`
| :math:`q'''` is volumetric internal heat source :math:`[\frac{W}{m^3}]`

The mass density, specific heat capacity, and thermal conductivity are all
properties of the material in question. The internal heat source can come from
many physical phenomena including nuclear fission, induction heating, and
chemical reactions.

The possible boundary conditions for the heat transfer equation are as follows:

1. Adiabatic
    No heat is transferred through the boundary.

    .. math::
        \nabla T = 0

2. Known Value (Dirichlet)
    The temperature at the boundary is known. It can be varying with time or
    space.

    .. math::
        T = f(t, \vec{r})

3. Known Flux (Neumann)
    The heat flux at the bounday is known. It can be varying with time or space.

    .. math::
        k(T, \vec{r}) \nabla T = f(t, \vec{r})

4. Convective (Robin)
    The heat flux is a function of the temperature and some bulk fluid
    temperature.

    .. math::
        k(T, \vec{r}) \nabla T = h (T - T_b)

    Where

    | :math:`h` is the heat transfer coefficient :math:`[\frac{W}{m^2 K}]`
    | :math:`T_b` bulk fluid temperature :math:`[K]`

5. Radiative
    The heat flux is equal to the net heat loss through thermal radiation from
    the boundary.

    .. math::
        k(T, \vec{r}) \nabla T = \epsilon \sigma (T^4 - T^4_env)

    Where

    | :math:`\epsilon` is surface emissivity
    | :math:`\sigma` is Stefan-Boltzmann constant :math:`[\frac{W}{m^2 K^4}]`
    | :math:`T_env` is the effective environment temperature :math:`[\frac{W}{m^2 K}]`

Numerical Implementation
------------------------

1. Finite Volume Method
    The Finite Volume Method (FVM) is applied to the heat transfer equation.
    
    Starting from the Partial Differential Equation (PDE):

    .. math::
        \rho(\vec{r}) c_p(T, \vec{r}) \frac{\partial T}{\partial t} - \nabla (k(T, \vec{r}) \nabla T) = q'''

    Rewrite in terms of flux.

    .. math::
        \rho(\vec{r}) c_p(T, \vec{r}) \frac{\partial T}{\partial t} - \nabla J = q_i'''

    Where

    | :math:`J` is the flux :math:`[\frac{W}{m^2}]`

    Now, we separate the computational domain into discrete volumes, known as
    cells. The governing equation is rewritten to be applied to a single cell,
    :math:`i`, by integrating all terms over the cell volume.

    .. math::
        \left[\rho(\vec{r}) c_p(T, \vec{r}) \frac{\partial T}{\partial t}\right]_i V_i - \sum^{m}_{j=1} A_j J_j = q_i

    Where

    | :math:`i` denotes the cell of interest
    | :math:`j` denotes a face of a cell
    | :math:`m` is the number of faces of cell :math:`i`
    | :math:`V` is cell volume :math:`[m^3]`
    | :math:`A` is area of a face :math:`[m^2]`
    | :math:`q` is internal heat source :math:`[W]`

2. Face Calculations
    Flux is defined by the gradient of temperature and the thermal conductivity.

    .. math::
        J = -k(T, \vec{r}) \nabla T

    This means that for every face we must approximate the thermal conductivity
    :math:`k` and the temperature gradient at a face using cell values. Once we
    do this, we can easily rearrange to a matrix form. The details depend on the
    condition at the face.

    Internal Face
        If a face is shared by two cells, then it is an internal face.

        The thermal conductivity at the face is found by interpolation between
        the two cell values of thermal conductivity.
        
        For internal faces, the thermal conductivity at the face is found by
        interpolation between the two cell values of thermal conductivity.

        .. math::
            w_f = \frac{1}{1+\frac{d_1}{d_2}}

        .. math::
            k_f = w_f k_1 + (1-w_f) k_2

        Where

        | :math:`1` denotes the first cell that uses the face :math:`f`
        | :math:`2` denotes the second cell that uses the face :math:`f`
        | :math:`d` is the distance from the cell center to the face center :math:`f`

        .. note::
            The order of the first and second cell is completely arbitrary and is determined by the mesh generation software and the order they're read in by the mesh reader.

        Next, we calculate the characteristic transport distance. This comes
        from dotting the vector pointing from one cell center to another with
        the normal of the face.

        .. math::
            L = (C_2 - C_1) \cdot \hat{n_f}

        Where

        | :math:`C` is the cell center
        | :math:`n` is the normal vector of the face

        Now we can define the temperature gradient.

        .. math::
            \nabla T = \frac{T_j - T_i}{L}

        Which gives us the flux through an internal face.

        .. math::
            J = -k(T, \vec{r}) \nabla T

        .. math::
            J = -k_f \frac{T_j - T_i}{L}
        
        .. important::
            :math:`k_f` is a function of the k for each cell, which are a function of temperature. So, this step makes the problem inherently nonlinear.

    Boundary Face: Adiabatic
        In the case of an adiabatic boundary, the heat flux is zero.

        .. math::
            J = 0 = -k_f(T, \vec{r}) \nabla T

        For all boundary faces the thermal conductivity is assumed to be the
        same at the face as it is at the cell center.

        .. math::
            k_f = k_i

        Since the thermal conductivity is non-zero, the temperature gradient
        must be zero.

        .. math::
            \nabla T = 0

        .. note::
            This is the default treatment of boundary faces unless a different boundary condition is set.

        .. note::
            Since :math:`J` of adiabatic faces are not a function of any cell temperatures, they do not contribute to the matrix form.

    Boundary Face: Known Value (Dirichlet) 
        The flux is computed using the temperature gradient calculated from a
        known temperature at the face and the computed temperature at the cell
        center.

        .. math::
            T_f = f(t, \vec{r})

        .. note::
            Even though the temperature is a known value, it can still change with time or position.

        Like all other boundary faces, the thermal conductivity is assumed to be
        the same at the face as it is at the cell center.

        .. math::
            k_f = k_i

        The temperature gradient is calculated in a similar fashion to internal
        faces with some key changes. The characteristic length, :math:`L` is
        calculated using the vector pointing from the face center to the cell
        center instead of from one cell center to a neighbors cell center. Also,
        the temperature value comes from the boundary condition instead of a
        neighboring cell.

        .. math::
            L = (C_f - C_1) \cdot \hat{n_f}

            \nabla T = \frac{T_f - T_i}{L}

    Boundary Face: Known Flux (Neumann)
        The face flux is directly known.

        .. math::
            J = f(t, \vec{r})

        Like all other boundary faces, the thermal conductivity is assumed to be
        the same at the face as it is at the cell center.

        .. math::
            k_f = k_i

        The temperature gradient is calculated in the exact same way as a face
        with a Dirichlet boundary condition.

        .. math::
            L = (C_f - C_1) \cdot \hat{n_f}

            \nabla T = \frac{T_f - T_i}{L}
    
3. Temporal Discretization
    If the problem is solving to steady state, the time derivative is zero.

    .. math::
        \frac{\partial T}{\partial t} = 0

    If the problem is doing a transient solve, the Backward Euler method (also
    known as Implicit Euler Method) is applied to the time derivatives.

    .. math::
        \rho_i {c^{n+1}_p}_i(T^{n+1}_i) V_i \frac{T^{n+1}_i - T^{n}_i}{\Delta t}  - \sum^{m}_{j=1} A_j J^{n+1}_j = q^{n+1}_i

    Where

    | :math:`n` denotes the timestep

    The only reference to the previous timestep's temperature value in cell
    :math:`i` comes from the time derivative. This value must be moved to the
    :math:`b` matrix since it is a known value.