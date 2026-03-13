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
        J_f = k_f \frac{T_1 - T_2}{\delta_1 + \delta_2} = k_1 \frac{T_1 - T_f}{\delta_1} = k_2 \frac{T_f - T_2}{\delta_2}

    .. math::
        \delta = d \cdot \hat{n_f}

    Where

    | :math:`1` denotes the first cell that uses the face :math:`f`
    | :math:`2` denotes the second cell that uses the face :math:`f`
    | :math:`d` is the distance from the cell center to the face center :math:`f`
    | :math:`\hat{n}` is the normal vector of face :math:`f`

    Rearrange to express temperature change in terms of heat flux.

    .. math::

        T_1 - T_f = J_f \frac{\delta_1}{k_1}

    .. math::

        T_f - T_2 = J_f \frac{\delta_2}{k_2}

    Add together.

    .. math::

        T_1 - T_2 = J_f \left( \frac{\delta_1}{k_1} + \frac{\delta_2}{k_2}\right)

    Now, substitute this expression for :math:`T_1-T_2` into the expression of
    heat flux using cell centered values.

    .. math::
    
        J_f = k_f \frac{J_f \left( \frac{\delta_1}{k_1} + \frac{\delta_2}{k_2}\right)}{\delta_1 + \delta_2}

    Assume heat flux is non zero on internal faces.

    .. math::
    
        1 = k_f \frac{\frac{\delta_1}{k_1} + \frac{\delta_2}{k_2}}{\delta_1 + \delta_2}

    Rearrange.

    .. math::

        k_f = \frac{d_1 + d_2}{\frac{d_1}{k_1}+\frac{d_2}{k_2}}

    Now, define a face weighting which is a function of geometry only.

    .. math::
        w_f = \frac{d_2}{d_1+d_2}

    Where

    | :math:`w` is geometric weighting

    Rewrite :math:`k_f` expression in terms of :math:`w`.

    .. math::

        \boxed{k_f = \frac{k_1 k_2}{\frac{k_2(1-w_f)}{k_1} + k_2 w_f}}

    Now we can define the temperature gradient.

    .. math::
        \boxed{\nabla T = \frac{T_j - T_i}{\delta_1 + \delta_2}}

    Which gives us the flux through an internal face.

    .. math::
        J = -k(T, \vec{r}) \nabla T

    .. math::
        J = -k_f \frac{T_j - T_i}{\delta_1 + \delta_2}
    
    .. important::
        :math:`k_f` is a function of the k for each cell, which is a function of temperature. So, this step makes the problem inherently nonlinear.

2. Boundary Face: Adiabatic
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
        Since :math:`J` of adiabatic faces are not a function of any cell temperatures, they do not contribute to the linear A-b matrices.

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
    :math:`\delta_i` instead of :math:`\delta_1 + \delta_2`. In other words, the
    characteristic length is effectively the distance from the cell center to
    the face center, rather than from the cell center to the neighboring cell
    center.
    
    The temperature value comes from the boundary condition
    instead of a neighboring cell.

    .. math::

        \nabla T = \frac{T_f - T_i}{\delta_i}

4. Boundary Face: Known Flux (Neumann)
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

        \nabla T = \frac{T_f - T_i}{\delta_i}
    
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