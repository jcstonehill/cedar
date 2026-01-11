Flow
=======

The ``Flow`` model solves the 1D compressible Euler equations for internal flow.


Theory
------

The governing equations for 1D compressible Euler flow are below [Berry_2018]_.

Conservation of Mass
^^^^^^^^^^^^^^^^^^^^
.. math::
    \frac{\partial \rho A}{\partial t} + \frac{\partial \rho u A}{\partial x} = 0

Where

| :math:`\rho` is mass density :math:`[\frac{kg}{m^3}]`
| :math:`A` is cross sectional flow area :math:`[m^2]`
| :math:`t` time :math:`[s]`
| :math:`u` time :math:`[\frac{m}{s}]`
| :math:`x` position :math:`[m]`

Conservation of Momentum
^^^^^^^^^^^^^^^^^^^^^^^^

.. math::
    \frac{\partial \rho u A}{\partial t} + \frac{\partial (\rho u^2 + P) A}{\partial x} = P \frac{\partial A}{\partial x} - F A + \rho g_x A

.. math::
    F = \frac{f \rho u |u|}{2 D_h}

.. math::
    D_h = \frac{4 A}{P_{wall}}

Where

| :math:`P` is pressure :math:`[Pa]`
| :math:`F` is viscous drag force density :math:`[\frac{N}/{m^3}]`
| :math:`g_x` is gravitaitonal acceleration :math:`[\frac{m}{s^2}]`
| :math:`f` is friction factor
| :math:`D_h` is hydraulic diameter :math:`[m]`
| :math:`P_{wall}` is wetted perimeter :math:`[m]`

Conservation of Energy
^^^^^^^^^^^^^^^^^^^^^^

.. math::
    \frac{\partial \rho E A}{\partial t} + \frac{\partial u (\rho E + P) A}{\partial x} = \rho u g_x A + \dot{q}''' A

Where

| :math:`E` is specific total energy :math:`[\frac{J}{kg}]`
| :math:`\dot{q}'''` is volumectric heat addition :math:`[\frac{W}{m^3}]`

Numerical Implementation
------------------------

Key assumptions are applied to the Partial Differential Equations (PDEs) above:

- All flow models are quasi-static. That is, the physics resolve much faster
  than other physics in the ``Problem``. Therefore, the time derivatives are zero.
- The cross sectional flow area is uniform throughout.
- Gravity is negligible.
- The values at a cell face can be approximated by the cell value upwind of the cell face.

The Finite Volume Method is applied to the three PDEs.

Conservation of Mass
^^^^^^^^^^^^^^^^^^^^

The conservation of mass equation is solved for cell velocity, while assuming
that all other values are fixed.

Starting from the conservation equation.

.. math::
    \frac{\partial \rho A}{\partial t} + \frac{\partial \rho u A}{\partial x} = 0

Time derivatives are zero with quasi-static assumption.

.. math::
    \frac{\partial \rho u A}{\partial x} = 0

Apply Finite Volume Method to a cell, :math:`i`.

.. math::
    \sum^{2}_{j=1} J_j = 0

.. math::
    J = \rho u A

Where

| :math:`j` denotes a face of a cell
| :math:`J` is the mass flux :math:`[\frac{kg}{s}]`

Apply an upwind scheme to determine the flux at a face in terms of cell values.

.. math::
    J_{i-1} + J_i = 0

.. note::
    :math:`J_{i-1}` will have a negative value here because its vector is dotted with the normal vector of the cell face.

Substitute in definition of :math:`J`.

.. math::
    - \rho_{i-1} u_{i-1} A + \rho_i u_i A = 0

Assume that the cross sectional area is uniform everywhere.

.. math::
    [- \rho_{i-1}] u_{i-1} + [\rho_i] u_i = 0

This is the final form of the equation. It is transformed to a sparse matrix,
the residual is evaluated, and then it is solved with ``scipy.spsolve``.

Only one boundary condition is required. The mass flow rate at the inlet is a
a Dirichlet boundary condition. From there, the inlet velocity can be obtained.

.. math::
    u_{inlet} = \frac{\dot{m}_{inlet}}{\rho_{inlet} A}

Conservation of Momentum
^^^^^^^^^^^^^^^^^^^^^^^^

The conservation of momentum equation is solved for cell pressure, while
assuming that all other values are fixed.

Starting from the conservation equation.

.. math::
    \frac{\partial \rho u A}{\partial t} + \frac{\partial (\rho u^2 + P) A}{\partial x} = P \frac{\partial A}{\partial x} - F A + \rho g_x A

Time derivatives are zero with quasi-static assumption.

.. math::
    \frac{\partial (\rho u^2 + P) A}{\partial x} = P \frac{\partial A}{\partial x} - F A + \rho g_x A

The cross-sectional flow area does not change with position.

.. math::
    \frac{\partial (\rho u^2 + P) A}{\partial x} = - F A + \rho g_x A

Neglect gravity.

.. math::
    \frac{\partial (\rho u^2 + P) A}{\partial x} = - F A

Apply Finite Volume Method to a cell, :math:`i`.

.. math::
    \sum^{2}_{j=1} J_j = - F \Delta x A

.. math::
    J = (\rho u^2 + P) A

Where

| :math:`\Delta x` is the length of a cell :math:`[m]`

Apply an upwind scheme to determine the flux at a face in terms of cell values.

.. math::
    J_i + J_{i+1}  = - F \Delta x A

.. note::
    :math:`J_i` will have a negative value here because its vector is dotted with the normal vector of the cell face.

Substitute in definition of :math:`J`.

.. math::
    - (\rho_i u_i^2 + P_i) A + (\rho_{i+1} u_{i+1}^2 + P_{i+1}) A = - F_i \Delta x A

Assume that the cross sectional area is uniform everywhere.

.. math::
    - \rho_i u_i^2 - P_i + \rho_{i+1} u_{i+1}^2 + P_{i+1}  = - F_i \Delta x

Rearrange.

.. math::
    [1] P_i + [-1] P_{i+1} = F_i \Delta x + \rho_{i+1} u_{i+1}^2 - \rho_i u_i^2

.. math::
    F_i = \frac{f_i \rho_i u_i^2}{2 D_h}

This is the final form of the equation. It is transformed to a sparse matrix,
the residual is evaluated, and then it is solved with ``scipy.spsolve``.

Only one boundary condition is required. The stagnation pressure at the outlet is a
a Dirichlet boundary condition. From there, the outlet pressure can be obtained.

.. math::
    P_{outlet} = {P_0}_{outlet} - \frac{1}{2} \rho_{outlet} u^2_{outlet}

Conservation of Energy
^^^^^^^^^^^^^^^^^^^^^^

The conservation of energy equation is solved for total specific energy, while
assuming that all other values are fixed.

Starting from the conservation equation.

.. math::
    \frac{\partial \rho E A}{\partial t} + \frac{\partial u (\rho E + P) A}{\partial x} = \rho u g_x A + \dot{q}''' A

Time derivatives are zero with quasi-static assumption.

.. math::
    \frac{\partial u (\rho E + P) A}{\partial x} = \rho u g_x A + \dot{q}''' A

Neglect gravity.

.. math::
    \frac{\partial u (\rho E + P) A}{\partial x} = \dot{q}''' A

Apply Finite Volume Method to a cell, :math:`i`.

.. math::
    \sum^{2}_{j=1} J_j = \dot{q}_i

Where

| :math:`\dot{q}` is internal heat source :math:`[W]`

.. math::
    J = u (\rho E + P) A

Apply an upwind scheme to determine the flux at a face in terms of cell values.

.. math::
    J_{i-1} + J_i = \dot{q}_i

.. note::
    :math:`J_{i-1}` will have a negative value here because its vector is dotted with the normal vector of the cell face.

Substitute in definition of :math:`J`.

.. math::
    - u_{i-1} (\rho_{i-1} E_{i-1} + P_{i-1}) A + u_i (\rho_i E_i + P_i) A = \dot{q}_i

Assume that the cross sectional area is uniform everywhere.

.. math::
    - u_{i-1} (\rho_{i-1} E_{i-1} + P_{i-1}) + u_i (\rho_i E_i + P_i) = \frac{\dot{q}_i}{A}

Rearrange.

.. math::
    [-u_{i-1} \rho_{i-1}] E_{i-1} + [u_i \rho_i] E_i = \frac{\dot{q}_i}{A} + u_{i-1} P_{i-1} - u_i P_i

This is the final form of the equation. It is transformed to a sparse matrix,
the residual is evaluated, and then it is solved with ``scipy.spsolve``.

Only one boundary condition is required. The stagnation temperature at the inlet
is a a Dirichlet boundary condition. The specific internal energy is found from
fluid properties as a function of temperature and pressure. The total internal
energy at the inlet is then calculated.

.. math::
    T_{inlet} = {T_0}_{inlet} - \frac{1}{2 {c_p}_{inlet}}  {u^2}_{inlet}

Where

| :math:`c_p` is specific heat capacity :math:`[\frac{J}{kg \cdot K}]`

.. math::
    e_{inlet} = f(T_{inlet}, P_{inlet})

.. math::
    E_{inlet} = e_{inlet} + \frac{1}{2} {u^2}_{inlet}

Friction Factor
^^^^^^^^^^^^^^^
A friction factor correlation is required to close the momentum equation. The
Churchill correlation is implemented for this [Churchill_1977]_.

.. math::
    f = 8\left[\left(\frac{8}{Re}\right)^{12}+\frac{1}{(A+B)^{1.5}}\right]^{\frac{1}{12}}

.. math::
    A = \left[2.457 \cdot ln\left( \frac{1}{(\frac{7}{Re})^{0.9} + 0.27 \frac{\epsilon}{D_h}} \right) \right]^{16}

.. math::
    B = (\frac{37530}{Re})^{16}

.. math::
    Re = \frac{\rho u D_h}{\mu}

Where

| :math:`Re` is Reynolds number.
| :math:`\epsilon` is wall roughness :math:`[m]`
| :math:`\mu` is dynamic viscosity :math:`[Pa \cdot s]`

Wall Temperature
^^^^^^^^^^^^^^^^
The ``Flow`` model uses a known :math:`\dot{q}_i` transferred to the fluid as
an input. The wall temperature is calculated from Newton's law of cooling as an
output during a post-processing step after the ``Fluid`` model is iterated.

Starting from Newton's Law of Cooling.

.. math::
    \dot{q}_i = h_i P_{wall} \Delta x ({T_{wall}}_i - T_i)

Where

| :math:`h` is heat transfer coefficient :math:`[\frac{W}{m^2 \cdot K}]`
| :math:`T_{wall}` is wall temperature :math:`[K]`

Rearrange.

.. math::
    {T_{wall}}_i = \frac{\dot{q}_i}{h_i P_{wall} \Delta x} + T_i

The heat transfer coefficient is calculated from the definition of Nusselt
number.

.. math::
    {Nu}_i = \frac{h_i D_h}{k_i}

Rearrange.

.. math::
    h_i = \frac {{Nu}_i k_i}{D_h}

A Nusselt number correlation is required to close the wall temperature
calculation. The Westinghouse-Modified McCarthy-Wolf correlation is implemented
for this [Thomas_1962]_.

.. math::
    {Nu}_i = 0.025 Re^{0.8}_i Pr^{0.4}_i \left( \frac{{T_{wall}}_i}{T_i}\right)^{-0.55} \left[1 + 0.3 \left( \frac{x_i}{D_h}\right)^{-0.7} \right]

.. math::

    Pr_i = \frac{\mu_i {c_p}_i}{k_i}

Where

| :math:`Pr` is Prandtl number
| :math:`x` is distance from inlet :math:`[m]`

References
----------

.. [Berry_2018] R. Berry et al., RELAP-7 Theory Manual, Apr. 2018.
       doi:10.2172/1498781 

.. [Churchill_1977] S. W. Churchill, “Friction Factor Equations
    Spans All Fluid-Flow Regimes,” Chemical Engineering Journal,
    84, 91-92 (1977).

.. [Thomas_1962] G. R. Thomas, "An Interim Study of Single Phase Heat Transfer
          Correlations Using Hydrogen," WANL-TNR-056, Apr. 1962.

.. |space| unicode:: 32