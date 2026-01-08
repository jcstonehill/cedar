import sys
import os, shutil
import numpy as np
import time

import cedar


class Problem:
    """
    Top-level problem definition and solver.

    A ``Problem`` coordinates one or more models, their meshes, and the
    coupling of variables between models. It is responsible for:

    - managing output directories and logging,
    - checking and setting up models,
    - applying inter-model variable coupling,
    - driving steady-state or transient solution loops,
    - exporting results.

    The ``Problem`` class does **not** solve physics directly; instead, it
    controls the iteration and data exchange between individual
    :class:`cedar.base.Model` instances.

    Parameters
    ----------
    name : str
        Name of the problem.
    restart_from_i : int, optional
        Restart index. If provided, the problem is treated as a restart.
        Restart capability is not yet implemented.
    create_outputs : bool, optional
        If ``True``, output directories and log files are created and results
        are exported.

    Attributes
    ----------
    name : str
        Name of the problem.
    path : str
        Filesystem path of the main executing script.
    models : list of cedar.base.Model
        Models participating in the problem.
    meshes : list of cedar.base.Mesh
        Meshes used by the models.
    connections : list of tuple
        Variable couplings defined as ``(src, tgt, adapter)``.
    is_restart : bool
        Indicates whether this problem is a restart.
    create_outputs : bool
        Controls whether outputs and logs are generated.

    Notes
    -----
    - Variable coupling is applied once per nonlinear iteration.
    - Convergence is based on the maximum residual returned by all models.
    - Restart capability is planned but not currently implemented.
    """

    def __init__(self, name, restart_from_i = None, create_outputs = True):
        """
        Initialize a problem instance.

        Parameters
        ----------
        name : str
            Name of the problem.
        restart_from_i : int, optional
            Restart index.
        create_outputs : bool, optional
            Enable or disable output generation.
        """

        self.is_restart = restart_from_i is not None
        self.i = restart_from_i
        self.create_outputs = create_outputs
        
        self.name = name
        main_file_path = os.path.abspath(str(sys.modules["__main__"].__file__))
        self.path = os.path.dirname(main_file_path)

        self.models: list[cedar.base.Model] = []
        self.meshes: list[cedar.base.Mesh] = []
        self.connections: list = []

    def add_model(self, model: cedar.base.Model):
        """
        Add a model to the problem.

        Parameters
        ----------
        model : cedar.base.Model
            Model to be added.
        """
         
        self.models.append(model)

    def couple(self, src, tgt, adapter = cedar.adapters.Direct):
        """
        Couple a source variable to a target variable.

        Parameters
        ----------
        src : cedar.base.Var
            Source variable.
        tgt : cedar.base.Var
            Target variable.
        adapter : cedar.base.Adapter, optional
            Adapter class used to transfer values from src to tgt.

        Notes
        -----
        Coupling is applied during every nonlinear iteration before model
        updates.
        """
        self.connections.append((src, tgt, adapter))

    def export(self, i = None, t = None):
        """
        Export model results.

        Parameters
        ----------
        i : int, optional
            Timestep index.
        t : float, optional
            Physical time.

        Notes
        -----
        Export is skipped if ``create_outputs`` is ``False``.
        """

        if not self.create_outputs:
            return
        
        model: cedar.base.Model
        for model in self.models:
            model.export(self.path + "/outputs", i, t)

    def restart(self, i):
        pass

    def solve(self, dt = None, t0 = 0, t_end = None):
        """
        Solve the problem.

        Parameters
        ----------
        dt : float, optional
            Time step size. If provided, a transient solve is performed.
            Otherwise, the problem is solved to steady state.
        t0 : float, optional
            Initial time.
        t_end : float, optional
            End time for transient simulations.

        Notes
        -----
        - Steady-state problems solve with no time derivatives.
        - Transient problems advance in time until ``t_end`` is reached.
        """

        is_transient = dt is not None

        if self.create_outputs:
            self._create_outputs_dir()
            cedar.Log.create(self.path)
        
        cedar.Log.message(f"Starting Problem: \"{self.name}\"")
        cedar.Log.add_level()
        self._log_models()
        cedar.Log.remove_level()
        cedar.Log.line_break()
        cedar.Log.message("Checking models.")
        self._check()
        cedar.Log.message("Setting up.")
        self._setup(t0, is_transient)
        cedar.Log.message("Starting solver.")

        if not self.connections:
            cedar.Log.message("No variables are coupled. Skipping problem iterations.")
            cedar.Log.line_break()

        # model: cedar.base.Model
        # for model in self.models.values():
        #     var: cedar.base.Var
        #     for var in model.var_dict.values():
        #         var.set(var.initial)

        if is_transient:
            self._solve_transient(dt, t0, t_end)

        else:
            self._solve_steady_state()

    def _check(self):
        model: cedar.base.Model
        for model in self.models:
            model.check()

    def _create_outputs_dir(self):
        output_path = self.path + "/outputs"
        backup_path = self.path + "/outputs/_backup"

        if not self.is_restart:
            if not os.path.isdir(output_path):
                os.mkdir(output_path)

            if os.path.isdir(backup_path):
                shutil.rmtree(backup_path)

            os.mkdir(backup_path)

            for item in os.listdir(output_path):
                if item == "_backup":
                    continue

                source_item_path = os.path.join(output_path, item)
                destination_item_path = os.path.join(backup_path, item)

                shutil.move(source_item_path, destination_item_path)

        else:
            cedar.Log.error("Problem restart capability is not yet implemented.")

    def _log_models(self):
        model: cedar.base.Model
        for model in self.models:
            cedar.Log.message(f"{model.name} ({model.__class__.__name__})")

    def _setup(self, t0, is_transient):
        model: cedar.base.Model
        for model in self.models:
            model.setup()

            if self.create_outputs:
                model.create_output_files(self.path + "/outputs", t0, is_transient)

    def _solve_once(self, t0, dt, tol = 1e-6, f = 1):
        """
        Perform a single nonlinear solve loop.

        Parameters
        ----------
        t0 : float or None
            Current physical time.
        dt : float or None
            Time step size.
        tol : float, optional
            Convergence tolerance.
        f : float, optional
            Relaxation factor applied to coupled variables.

        Notes
        -----
        - Residuals are computed per model.
        - Convergence is determined by the maximum residual.
        """
        #if self.connections:
        #cedar.Log.message(f"{0:>4} |R|: N/A")
            # cedar.Log.add_level()

        f = 0.5

        for i in range(1000):
            residuals = []

            model: cedar.base.Model
            src: cedar.base.Var
            tgt: cedar.base.Var
            adapter: cedar.base.Adapter

            for model in self.models:
                for src, tgt, adapter in self.connections:
                    if model.name == tgt.model_name:
                        src.relax(f)
                        adapter.transfer(src, tgt)

                residuals.append(model.iterate(t0, dt))

            # for src, tgt, adapter in self.connections:
            #     residuals.append(np.sqrt(src.res2()))
            #print(residuals)
            res = np.max(residuals)

            cedar.Log.message(f"{i+1:>4} |R|: {res:.5e}")

            if res <= tol:
                break

        # if self.connections:
        #     cedar.Log.remove_level()

    def _solve_steady_state(self):
        """
        Solve the problem to steady state.

        Notes
        -----
        - Results are exported after convergence.
        - Total solve time is logged.
        """
        cedar.Log.line_break()
        cedar.Log.message("Solving to steady state.")
        cedar.Log.add_level()

        start_time = time.time()
        self._solve_once(None, None)
        end_time = time.time()

        self.export()

        cedar.Log.line_break()
        cedar.Log.message("SOLVE COMPLETE!")

        cedar.Log.message(cedar.helper.format_computation_time(end_time - start_time))
        cedar.Log.remove_level()
        cedar.Log.line_break()
        cedar.Log.line_break()

    def _solve_transient(self, dt, start = 0, end = None):
        """
        Solve the problem transiently.

        Parameters
        ----------
        dt : float
            Time step size.
        start : float, optional
            Start time.
        end : float, optional
            End time.

        Notes
        -----
        - Variables are advanced in time and marched after each step.
        - Results are exported at every timestep.
        """
        cedar.Log.message(f"Starting at t = {start} [s].")

        time_i = 0
        t0 = start
        end = t0 + dt if end is None else end

        cedar.Log.message(f"Ending at t = {end} [s].")
        cedar.Log.line_break()
        cedar.Log.line_break()

        while t0 < end:
            if t0 + dt > end:
                dt = end - t0

            cedar.Log.message("Time = " + str(t0+dt) + " [s], dt = " + str(dt) +" [s]")
            cedar.Log.add_level()

            start_time = time.time()
            self._solve_once(t0, dt)
            end_time = time.time()

            time_i += 1
            self.export(time_i, t0+dt)

            cedar.Log.line_break()
            cedar.Log.message("Timestep Converged!")

            cedar.Log.message(cedar.helper.format_computation_time(end_time - start_time))
            cedar.Log.remove_level()
            cedar.Log.line_break()
            cedar.Log.line_break()

            model: cedar.base.Model
            for model in self.models:

                var: cedar.base.Var
                for var in model.var_dict.values():
                    var.set_initial(var.val)

            t0 += dt

        cedar.Log.message("SOLVE COMPLETE!")
        cedar.Log.line_break()
        cedar.Log.line_break()

    # def _residual(self):
    #     res2 = 0

    #     tgt: cedar.base.Var
    #     for _, tgt, _ in self.connections:
    #         res2 += tgt.res2()

    #     return np.sqrt(res2)
