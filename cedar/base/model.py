from abc import ABC, abstractmethod
import numpy as np
import meshio
import os

import cedar
from cedar.base.mesh import Mesh
from cedar.base.var import Var


class Model(ABC):
    """
    Abstract base class for simulation models.

    A ``Model`` owns a mesh, a collection of variables, and the logic required
    to set up, solve, and export simulation results. Concrete models must
    implement validation, setup, and solution procedures.
    """

    model_type = "Model"

    @property
    def mesh(self):
        """
        Mesh associated with the model.
        """
        return self._mesh

    @property
    def params(self):
        """
        Model parameters container.
        """
        return self._params

    @property
    def vars(self):
        """
        Model variables container.
        """
        return self._vars

    def __init__(self, name: str, mesh: Mesh):
        """
        Initialize the model.

        Parameters
        ----------
        name : str
            Model name used for logging and output files.
        mesh : Mesh
            Computational mesh associated with the model.
        """
        self.name = name
        self._mesh = mesh

        self.var_dict: dict[Var] = None

        self._params = None
        self._vars = None

    def _check_var_name(self, name: str):
        """
        Ensure variable names are unique within the model.
        """
        if name in self.var_dict:
            self.log_error(
                f'Tried to add variable named "{name}" but a variable with that name already exists.'
            )

    def add_var(self, var: Var) -> Var:
        """
        Register a variable with the model.

        Parameters
        ----------
        var : Var
            Variable to be added.

        Returns
        -------
        Var
            The added variable.
        """
        if not hasattr(self, "var_dict"):
            self.var_dict = {}

        self._check_var_name(var.name)
        self.var_dict[var.name] = var

        return var

    def log_message(self, message: str):
        """
        Log an informational message with model context.
        """
        cedar.Log.message(f"{self.name} ({self.model_type}) :: {message}")

    def log_error(self, message: str):
        """
        Log an error message with model context.
        """
        cedar.Log.error(f"{self.name} ({self.model_type}) :: {message}")

    @abstractmethod
    def check(self):
        """
        Validate model configuration and inputs.
        """
        pass

    @abstractmethod
    def setup(self):
        """
        Perform model setup prior to solving.
        """
        pass

    @abstractmethod
    def iterate(self, t0=None, dt=None, res_reduc=1e-2) -> np.ndarray:
        """
        Perform one model iteration.

        Parameters
        ----------
        t0 : float, optional
            Initial time.
        dt : float, optional
            Time step size.
        res_reduc : float, optional
            Target residual reduction for convergence.

        Returns
        -------

        """
        pass

    def export(self, case_path, i=None, t=None):
        """
        Export model results to CSV and VTK formats.

        Parameters
        ----------
        case_path : str
            Output directory.
        i : int, optional
            Time step or iteration index.
        t : float, optional
            Simulation time.
        """
        base_filename = case_path + "/" + self.name

        self._export_to_csv(base_filename, self._scalars(), self._cell_arrays(), i, t)

        if self.mesh is not None:
            self._export_to_vtk(base_filename, self._scalars(), self._cell_arrays(), i, t)

    def create_output_files(self, case_path, t0, is_transient):
        """
        Create initial output files and headers.

        Parameters
        ----------
        case_path : str
            Output directory.
        t0 : float
            Initial simulation time.
        is_transient : bool
            Whether the simulation is transient.
        """
        scalars = self._scalars(use_initial=True)
        arrays = self._cell_arrays(use_initial=True)

        base_filename = case_path + "/" + self.name
        filename = base_filename + ".csv"

        with open(filename, "w") as file:
            file.write("i, t")
            for key in scalars.keys():
                file.write(", " + key)
            file.write("\n")

        if isinstance(self.mesh, cedar.Mesh1D):
            for key in arrays.keys():
                filename = base_filename + "_" + key + ".csv"
                with open(filename, "w") as file:
                    file.write("_, z_pos")
                    for cell_center in self.mesh.cell_centers:
                        file.write(", " + str(cell_center[2]))
                    file.write("\n")

                    if is_transient:
                        file.write(f"0, {t0}")
                        for val in arrays[key]:
                            file.write(", " + str(val))
                        file.write("\n")

        if is_transient:
            self._export_to_csv(base_filename, scalars, arrays, 0, t0)
            if self.mesh is not None:
                self._export_to_vtk(base_filename, scalars, arrays, 0, t0)

    def _scalars(self, use_initial=False):
        """
        Collect scalar-valued variables for output.
        """
        scalars = {}

        for var in self.var_dict.values():
            val = var.initial if use_initial else var.val

            if isinstance(var, cedar.ScalarVar):
                scalars[var.name] = val

            elif isinstance(var, cedar.FlowStateVar):
                scalars[var.name + "_T0"] = val[0]
                scalars[var.name + "_P0"] = val[1]
                scalars[var.name + "_mdot"] = val[2]

        return scalars

    def _cell_arrays(self, use_initial=False):
        """
        Collect cell-based array variables for output.
        """
        arrays = {}

        for var in self.var_dict.values():
            val = var.initial if use_initial else var.val

            if (
                isinstance(var, (cedar.Mesh1DVar, cedar.Mesh3DVar))
                and var.boundary is None
            ):
                arrays[var.name] = val

        return arrays

    def _export_to_csv(self, base_filename, scalars, arrays, i=None, t=None):
        """
        Append results to CSV output files.
        """
        i = 0 if i is None else i
        t = "Steady" if t is None else t

        filename = base_filename + ".csv"
        with open(filename, "a") as file:
            file.write(f"{i}, {t}")
            for val in scalars.values():
                file.write(", " + str(val))
            file.write("\n")

        if isinstance(self.mesh, cedar.Mesh1D):
            for key, array in arrays.items():
                filename = base_filename + "_" + key + ".csv"
                with open(filename, "a") as file:
                    file.write(f"{i}, {t}")
                    for val in array:
                        file.write(", " + str(val))

    def _export_to_vtk(self, base_filename, scalars, arrays, i=None, t=None):
        """
        Export results to VTK format for visualization.
        """
        suffix = "" if i is None else "_" + f"{i:05d}"
        filename = base_filename + suffix + ".vtk"

        cell_type = "tetra" if isinstance(self.mesh, cedar.Mesh3D) else "line"
        cells = [(cell_type, self.mesh.cell_point_ids)]

        data = {}
        for key, val in scalars.items():
            data[key] = [val * np.ones(self.mesh.N_cells)]
        for key, val in arrays.items():
            data[key] = [val]

        mesh0 = meshio.Mesh(
            points=self.mesh.points,
            cells=cells,
            cell_data=data,
        )
        meshio.write(filename, mesh0)