Numerical Implementation
========================

.. warning::
    The :class:`cedar.HeatTransfer` model does not account for tangential cell flux at this time. It will be incorporated as future work.

Governing Equation
------------------

The Finite Volume Method is applied to the heat transfer equation.

Starting from the Partial Differential Equation:

.. math::
    \rho(\vec{r}) c_p(T, \vec{r}) \frac{\partial T}{\partial t} - \nabla (k(T, \vec{r}) \nabla T) = \dot{q}'''(t)

| :math:`\rho` is mass density :math:`[\frac{kg}{m^3}]`
| :math:`\vec{r}` is spatial position vector :math:`[m]`
| :math:`c_p` is specific heat capacity :math:`[\frac{J}{kg K}]`
| :math:`T` is temperature :math:`[K]`
| :math:`t` is time :math:`[s]`
| :math:`k` is thermal conductivity :math:`[\frac{W}{m K}]`
| :math:`\dot{q}'''` is volumetric internal heat source :math:`[\frac{W}{m^3}]`

Use divergence theorem to rewrite in terms of flux.

.. math::
    \rho(\vec{r}) c_p(T, \vec{r}) \frac{\partial T}{\partial t} + \nabla J = \dot{q}_i'''(t)

Where

| :math:`J` is the flux out of a cell :math:`[\frac{W}{m^2}]`

Now, we separate the computational domain into discrete volumes, known as
cells. The governing equation is rewritten to be applied to a single cell,
:math:`i`, by integrating all terms over the cell volume.

.. math::
    \boxed{\left[\rho(\vec{r}) c_p(T, \vec{r}) \frac{\partial T}{\partial t}\right]_i V_i + \sum^{m}_{j=1} A_j J_j = \dot{q}_i(t)}

Where

| :math:`i` denotes the cell of interest
| :math:`j` denotes a face of a cell
| :math:`m` is the number of faces of cell :math:`i`
| :math:`V` is cell volume :math:`[m^3]`
| :math:`A` is area of a face :math:`[m^2]`
| :math:`\dot{q}` is internal heat source :math:`[W]`

Examining the boxed equation above, further manipulationg is required. We need
expressions for the face flux and the temporal terms.

Face Flux
---------

For every face condition, we must derive an expression for flux to be used in
the governing equation. Once we do this, it is trivial to rearrange to a matrix
form. The flux derivations for each face condition are shown below.

1. Internal Face
    If a face is shared by two cells, then it is an internal face.

    The thermal conductivity at an internal face is derived with consideration
    of flux balance at a face. This is especially important for internal faces
    where the two connecting cells have dissimilar material properties.

    The derivation starts by considering the approximations of heat flux as the
    face is approached from each connecting cell. To ensure flux balance, the
    face flux calculated from the neighboring cell centered values must be the
    same as the flux found when approaching the face from each side.

    .. math::
        J_f = k_f \frac{T_i - T_j}{\delta_i + \delta_j} = k_i \frac{T_i - T_f}{\delta_i} = k_j \frac{T_f - T_j}{\delta_j}

    .. math::
        \delta = d \cdot \hat{n_f}

    Where

    | :math:`i` denotes the first cell that uses the face :math:`f`
    | :math:`j` denotes the second cell that uses the face :math:`f`
    | :math:`d` is the distance from the cell center to the face center :math:`f`
    | :math:`\hat{n}` is the normal vector of face :math:`f`

    Rearrange to express temperature change in terms of heat flux.

    .. math::

        T_i - T_f = J_f \frac{\delta_i}{k_i}

    .. math::

        T_f - T_j = J_f \frac{\delta_j}{k_j}

    Add together.

    .. math::

        T_i - T_j = J_f \left( \frac{\delta_i}{k_i} + \frac{\delta_j}{k_j}\right)

    Now, substitute this expression for :math:`T_i-T_j` into the expression of
    heat flux using cell centered values.

    .. math::
    
        J_f = k_f \frac{J_f \left( \frac{\delta_i}{k_i} + \frac{\delta_j}{k_j}\right)}{\delta_i + \delta_j}

    Assume heat flux is non zero on internal faces.

    .. math::
    
        1 = k_f \frac{\frac{\delta_i}{k_i} + \frac{\delta_j}{k_j}}{\delta_i + \delta_j}

    Rearrange.

    .. math::

        k_f = \frac{\delta_i + \delta_j}{\frac{\delta_i}{k_i}+\frac{\delta_j}{k_j}}

    Now, define a face weighting which is a function of geometry only.

    .. math::
        w_f = \frac{\delta_j}{\delta_i+\delta_j}

    Where

    | :math:`w` is geometric weighting

    Rewrite :math:`k_f` expression in terms of :math:`w`.

    .. math::

        k_f = \frac{k_i k_j}{\frac{k_j(1-w_f)}{k_i} + k_j w_f}

    Now that we have an expression for the thermal conductivity at the internal
    face, we can derive the expression for face flux.

    .. math::
        \nabla T = -\frac{T_i - T_j}{\delta_i + \delta_j}

    .. math::
        J_f = -k_f(T, \vec{r}) \nabla T

    .. math::
        \boxed{J_f = k_f \frac{T_i - T_j}{\delta_i + \delta_j}}
    
    .. important::
        :math:`k_f` is a function of the k for each cell, which is a function of temperature. So, this step makes the problem inherently nonlinear.

2. Boundary Face: Adiabatic
    In the case of an adiabatic boundary, the heat flux is zero.

    .. math::
        J_f = 0 = -k_f(T, \vec{r}) \nabla T

    .. note::
        This is the default treatment of boundary faces unless a different boundary condition is set.

    .. note::
        Since :math:`J_f` of adiabatic faces are not a function of any cell temperatures, they do not contribute to the linear A-b matrices.

3. Boundary Face: Known Value (Dirichlet) 
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
    faces with some key changes. The characteristic length is now
    :math:`\delta_i` instead of :math:`\delta_i + \delta_j`. In other words, the
    characteristic length is effectively the distance from the cell center to
    the face center, rather than from the cell center to the neighboring cell
    center.
    
    The temperature value comes from the boundary condition
    instead of a neighboring cell.

    .. math::

        \nabla T = -\frac{T_i - T_f}{\delta_i}

    .. math::
        J_f = -k_f \nabla T

    .. math::
        J_f = k_f \frac{T_i - T_f}{\delta_i}

4. Boundary Face: Known Flux (Neumann)
    The face flux is directly prescribed.

    .. math::
        J_f = f(t, \vec{r})

5. Boundary Face: Convective (Robin)
    The convection boundary condition enforces a flux that is a function of the
    boundary face temperature and some bulk flow temperature associated with the
    boundary.

    .. math::
        J_f = htc_f (T_f - T_\infty)

    Where

    | :math:`\infty` denotes the bulk flow temperature.
    | :math:`htc` is the heat transfer coefficient :math:`\frac{W}{m^2 K}`

    However, since this expression of flux is in terms of boundary face
    temperature, we must use a conservation of heat flux to write this in terms
    of cell-centered temperatures.

    .. math::
        J_f = k_f \frac{T_i - T_f}{\delta_i} = htc_f \cdot (T_f - T_\infty)

    We rearrange the convective boundary flux expression to obtain boundary face
    temperature as a function of flux.

    .. math::
        T_f = \frac{J_f}{htc_f} + T_\infty

    Plug this into the expression of heat flux in the cell.

    .. math::
        J_f = k_f \frac{T_i - T_f}{\delta_i}

    .. math::
        J_f = \frac{k_f}{\delta_i} (T_i - \frac{J_f}{htc_f} - T_\infty)

    Now, solve for heat flux.

    .. math::
        \boxed{J_f = \frac{htc_f \cdot k_f}{\delta_i \cdot htc_f + k_f}(T_i - T_\infty)}

Temporal Terms
--------------

If the problem is solving to steady state, the time derivative is zero.

.. math::
    \frac{\partial T}{\partial t} = 0

If the problem is doing a transient solve, the Backward Euler method (also
known as Implicit Euler Method) is applied to the time derivatives.

.. math::
    \rho_i {c^{n+1}_p}_i(T^{n+1}_i) V_i \frac{T^{n+1}_i - T^{n}_i}{\Delta t}  - \sum^{m}_{j=1} A_j J^{n+1}_j = \dot{q}^{n+1}_i

Where

| :math:`n` denotes the timestep

The only reference to the previous timestep's temperature value in cell
:math:`i` comes from the time derivative. This value must be moved to the
:math:`b` matrix since it is a known value.