Heat Transfer
=============

.. toctree::
   :maxdepth: 1

   heat_transfer/implementation.rst     

The :class:`Heat Transfer` model solves the heat conduction equation on a 3D mesh. 

Heat conduction is the process by which thermal energy moves through a material
due to microscopic interactions between its atoms. When one area is hotter
than another, faster-moving particles in the hot area transfer energy to
neighboring, slower-moving particles through collisions or vibrations. This
results in a net flow of heat from high temperatures to low temperatures,
without any bulk motion of the material.

At the macroscopic level, this behavior is described by Fourier's law.

.. math::
    \rho(\vec{r}) c_p(T, \vec{r}) \frac{\partial T}{\partial t} - \nabla (k(T, \vec{r}) \nabla T) = \dot{q}'''

Where

| :math:`\rho` is mass density :math:`[\frac{kg}{m^3}]`
| :math:`c_p` is specific heat capacity :math:`[\frac{J}{kg K}]`
| :math:`T` is temperature :math:`[K]`
| :math:`k` is thermal conductivity :math:`[\frac{W}{m K}]`
| :math:`\dot{q}'''` is volumetric internal heat source :math:`[\frac{W}{m^3}]`

The mass density, specific heat capacity, and thermal conductivity are all
properties of the material in question. The internal heat source can come from
many physical phenomena including nuclear fission, induction heating, and
chemical reactions.

The possible boundary conditions for the heat transfer equation are as follows:

1. Adiabatic
    No heat is transferred through the boundary.

    .. math::
        \nabla T = 0

2. Temperature (Dirichlet)
    The temperature at the boundary is known. It can be varying with time or
    space.

    .. math::
        T = f(t, \vec{r})

3. Flux (Neumann)
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