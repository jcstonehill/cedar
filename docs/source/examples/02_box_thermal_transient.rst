02 - Box Thermal Transient
==========================

This example is the exact same as 01 - Box Thermal except that the
``Problem.solve()`` method is called with the ``dt`` and ``t_end`` parameters
specified. When these parameters are specified, ``Problem.solve()`` initiates a
transient solve instead of a steady state solve.

A transient solve and a steady state solve with Cedar are very similar behind
the scenes, but there are some key differences.


For Steady State:

- ``Var.initial`` corresponds to the initial guess of ``Var.val``.
- Internally, all Models are solved with any time derivatives set to zero.
- ``Problem`` has one "timestep" (to steady state), and ``Var.val`` corresponds to the steady state value after solve is complete.
- For relevant types of ``Var``, only a single ``.vtk`` file will be created in the Outputs folder.

For Transient:

- ``Var.initial`` corresponds to the value of ``Var.val`` at the start of the timestep AND is also the initial guess of ``Var.val``.
- Internally, all Models are solved with relevant time derivatives included.
- ``Problem`` has many timesteps, where ``Var.val`` is the value after some ``dt``, starting from ``Var.initial``. Once a timestep is solved, the ``Problem`` marches to a new timestep by setting ``Var.initial`` = ``Var.val`` and repeating the process until ``t_end``.
- For relevant types of ``Var``, a ``.vtk`` file will be created in the Outputs folder for each timestep.

Specifying a Transient Solve
----------------------------

Specifying a transient solve is as simple as defining the ``dt`` and ``t_end`` parameters in the ``Problem.solve()`` method.

So...

.. code-block:: python

    problem.solve()

becomes...

.. code-block:: python

    problem.solve(dt = 100, t_end = 10000)

No other changes are necessary.

Entire Example Problem File
---------------------------

.. code-block:: python

    import cedar

    # Instantiate Problem Object
    problem = cedar.Problem("box_thermal")

    # Build Mesh for Thermal Model
    mesh = cedar.Mesh3D("box_thermal.msh")

    # Define Material Properties
    ZrC = cedar.materials.ZrC()

    # Instantiate Model
    thermal = cedar.models.Thermal("thermal", mesh, {"volume" : ZrC})

    # Set Heat Source
    thermal.set_Qdot(5e4)

    # Set Boundary Conditions (Default is Adiabatic)
    thermal.set_bc("left", "dirichlet", 300)
    thermal.set_bc("right", "dirichlet", 500)

    # Set Initial Conditions
    thermal.vars.T.set_initial(500)

    # Add Problem to Model and Solve
    problem.add_model(thermal)
    problem.solve(dt = 100, t_end = 10000)

Outputs
-------

Below is a temperature contour plot animation. You can see how the temperature
distribution evolves through time.

.. figure:: img/02_box_thermal_transient_T_anim.gif
   :scale: 50 %
   :align: center

   Temperature contour plot animation.