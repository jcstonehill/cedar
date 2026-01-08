from __future__ import annotations
from abc import ABC, abstractmethod
import numpy as np

import cedar


class Var(ABC):
    """
    Abstract base class for model variables.

    A ``Var`` represents a named physical quantity associated with a model.
    Variables may be scalar-valued or defined over a mesh, region, or boundary.
    Direct mutation of stored values is prohibited to ensure consistency and
    enable controlled updates (e.g. relaxation).

    Notes
    -----
    - Values are always stored internally as NumPy arrays.
    - Accessors return copies to prevent accidental in-place modification.
    - Concrete subclasses define how values are interpreted and mapped to space.
    """

    @property
    def val(self) -> np.ndarray:
        """
        Current value of the variable.

        Returns
        -------
        ndarray or None
            Copy of the current value array.
        """
        return np.copy(self._val) if self._val is not None else None

    @val.setter
    def val(self, val):
        """
        Prevent direct assignment to ``val``.
        """
        self.log_error(
            'Tried to set "Var.val" attribute directly. Use var.set() instead.'
        )

    @property
    def initial(self) -> np.ndarray:
        """
        Initial value of the variable.

        Returns
        -------
        ndarray or None
            Copy of the initial value array.
        """
        return np.copy(self._initial) if self._initial is not None else None

    @initial.setter
    def initial(self, initial):
        """
        Prevent direct assignment to ``initial``.
        """
        self.log_error(
            'Tried to set "Var.initial" attribute directly. Use var.set_initial() instead.'
        )

    @property
    def mesh(self) -> cedar.base.Mesh:
        """
        Mesh associated with the variable, if any.
        """
        return self._mesh

    @property
    def region(self) -> str:
        """
        Region name if the variable is defined on a mesh region.
        """
        return self._region

    @property
    def boundary(self) -> str:
        """
        Boundary name if the variable is defined on a mesh boundary.
        """
        return self._boundary

    @property
    def N(self) -> int:
        """
        Number of degrees of freedom for the variable.
        """
        return self._N

    def __init__(
        self,
        model_name: str,
        name: str,
        initial,
        units: str = None,
        mesh: cedar.base.Mesh = None,
        region: str = None,
        boundary: str = None,
    ):
        """
        Initialize a variable.

        Parameters
        ----------
        model_name : str
            Name of the owning model (used for logging).
        name : str
            Variable name.
        initial : ndarray
            Initial value(s) of the variable.
        units : str, optional
            Physical units of the variable.
        mesh : Mesh, optional
            Mesh on which the variable is defined.
        region : str, optional
            Mesh region name.
        boundary : str, optional
            Mesh boundary name.

        Notes
        -----
        A variable may be associated with *either* a region or a boundary,
        but not both.
        """
        if region is not None and boundary is not None:
            self.log_error(
                "Both region and boundary can not be defined. "
                "A variable can only be a region variable or a boundary variable."
            )

        if mesh is None:
            self._N = 1
        else:
            if region is not None:
                self._N = mesh.region_N[region]
            elif boundary is not None:
                self._N = mesh.boundary_N[boundary]
            else:
                self._N = mesh.N_cells

        self.model_name = model_name
        self.name = name
        self.units = units

        self._mesh = mesh
        self._region = region
        self._boundary = boundary

        self.prev_val = self.interpret_val(initial)
        self._val = self.interpret_val(initial)
        self._initial = self.interpret_val(initial)

    def log_message(self, message: str):
        """
        Log an informational message with variable context.
        """
        cedar.Log.message(f"{self.model_name} :: {self.name} :: {message}")

    def log_error(self, message: str):
        """
        Log an error message with variable context.
        """
        cedar.Log.error(f"{self.model_name} :: {self.name} :: {message}")

    def relax(self, f: float):
        """
        Apply under-relaxation to the variable.

        Parameters
        ----------
        f : float
            Relaxation factor in the range [0, 1].
        """
        prev_val = np.copy(self.prev_val)
        self.set(f * self.val + (1 - f) * self.prev_val)
        self.prev_val = prev_val

    def set(self, val: np.ndarray):
        """
        Update the variable value.

        Parameters
        ----------
        val : ndarray
            New value(s) of the variable.
        """
        self.prev_val = self._val
        self._val = self.interpret_val(val)

        if self.initial is None:
            self.set_initial(val)

    def set_initial(self, initial: np.ndarray):
        """
        Set the initial value of the variable.

        Parameters
        ----------
        initial : ndarray
            Initial value(s).
        """
        self._initial = self.interpret_val(initial)

    @abstractmethod
    def interpret_val(self, val: np.ndarray) -> np.ndarray:
        """
        Interpret and validate raw input values.

        Returns
        -------
        ndarray
            Validated and correctly shaped value array.
        """
        pass

    @abstractmethod
    def points(self) -> np.ndarray:
        """
        Return spatial points associated with the variable.

        Returns
        -------
        ndarray
            Coordinates corresponding to each degree of freedom.
        """
        pass

    # @abstractmethod
    # def res2(self) -> np.ndarray:
    #     """
    #     Return squared residuals for the variable.

    #     Returns
    #     -------
    #     ndarray
    #         Squared residual values.
    #     """
    #     pass