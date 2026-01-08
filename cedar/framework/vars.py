import numbers
import numpy as np
import meshio

import cedar


class ScalarVar(cedar.base.Var):
    """
    Scalar variable with no spatial dependence.

    This class represents a single scalar quantity that is not associated
    with a mesh. Internally, values are stored as NumPy arrays of dtype
    ``float64`` to ensure numerical consistency with other variable types.

    Notes
    -----
    - Residuals are computed using a normalized squared difference.
    """

    def interpret_val(self, val: np.ndarray) -> np.ndarray:
        """
        Interpret and store a scalar value.

        Parameters
        ----------
        val : ndarray
            Input value.

        Returns
        -------
        ndarray
            Scalar value converted to a ``float64`` NumPy array.
        """
        return np.array(val, dtype=np.float64)

    def points(self) -> list[tuple]:
        """
        Return the spatial point associated with this variable.

        Returns
        -------
        list of tuple
            A single dummy point at the origin.
        """
        return [(0, 0, 0)]

    def res2(self) -> np.ndarray:
        """
        Compute the normalized squared residual.

        Returns
        -------
        ndarray
            Normalized squared residual between the current and previous
            values.
        """
        return (self.val - self.prev_val) ** 2 / self.val ** 2
    
class FlowStateVar(cedar.base.Var):
    """
    Lumped flow-state variable.

    This variable represents a three-component flow state consisting of
    stagnation temperature, stagnation pressure, and mass flow rate.

    The value vector is ordered as:

    - ``T0`` : stagnation temperature [K]
    - ``P0`` : stagnation pressure [Pa]
    - ``m_dot`` : mass flow rate [kg/s]

    Notes
    -----
    - This variable is not associated with a mesh.
    - The internal state vector always has length 3.
    """

    def __init__(self, model_name: str, name: str):
        """
        Initialize a flow-state variable.

        Parameters
        ----------
        model_name : str
            Name of the model this variable belongs to.
        name : str
            Name of the variable.
        """
        self._N = 3

        self.model_name = model_name
        self.name = name
        self.units = ("K", "Pa", "kg/s")

        self.prev_val = np.array((300.0, 101325.0, 1.0))
        self._val = np.array((300.0, 101325.0, 1.0))
        self._initial = np.array((300.0, 101325.0, 1.0))

    def interpret_val(self, val: np.ndarray) -> np.ndarray:
        """
        Interpret and validate the flow-state vector.

        Parameters
        ----------
        val : ndarray
            Flow-state values ``(T0, P0, m_dot)``.

        Returns
        -------
        ndarray
            Flow-state vector as a ``float64`` NumPy array.

        Raises
        ------
        Exception
            If the input does not contain exactly three values.
        """
        if len(val) != 3:
            self.log_error(
                "Tried to set value using incorrect size. "
                "FlowStateVar.val must be of size 3. (T0, P0, m_dot)."
            )

        return np.array(val, dtype=np.float64)

    def points(self):
        """
        Return the spatial point associated with this variable.

        Returns
        -------
        list of tuple
            A single dummy point at the origin.
        """
        return [(0, 0, 0)]

    def res2(self) -> np.ndarray:
        """
        Compute the normalized squared residual.

        Returns
        -------
        ndarray
            Normalized squared residual between the current and previous
            flow-state vectors.
        """
        return np.sum((self.val - self.prev_val) ** 2) / np.sum(self.val ** 2)

class Mesh1DVar(cedar.base.Var):
    """
    Variable defined on a one-dimensional mesh.

    This variable may be associated with:
    - all cells,
    - a specific region, or
    - a specific boundary.

    Scalar inputs are automatically expanded to match the number of mesh
    entities.
    """

    def interpret_val(self, val: np.ndarray) -> np.ndarray:
        """
        Interpret values defined on a 1D mesh.

        Parameters
        ----------
        val : array_like
            Input value(s).

        Returns
        -------
        numpy.ndarray
            Array of values with length equal to the number of mesh cells
            or faces.
        """
        val = np.array(val, dtype=np.float64)

        if val.size == 1:
            return val * np.ones(self.N)
        else:
            return val

    def points(self):
        """
        Return spatial locations associated with this variable.

        Returns
        -------
        numpy.ndarray
            Coordinates of mesh points corresponding to the variable
            definition.
        """
        if self.boundary is not None:
            return self.mesh.face_centers[self.mesh.boundaries[self.boundary]]

        elif self.region is not None:
            return self.mesh.cell_centers[self.mesh.regions[self.region]]

        else:
            return self.mesh.cell_centers

    def res2(self) -> np.ndarray:
        """
        Compute the normalized squared residual.

        Returns
        -------
        numpy.ndarray
            Normalized squared residual between current and previous values.
        """
        return np.sum((self.val - self.prev_val) ** 2) / np.sum(self.val ** 2)

class Mesh3DVar(cedar.base.Var):
    """
    Variable defined on a three-dimensional mesh.

    This variable may be associated with:
    - all cells,
    - a specific region, or
    - a specific boundary.

    Scalar inputs are automatically expanded to match the number of mesh
    entities.
    """

    def interpret_val(self, val: np.ndarray) -> np.ndarray:
        """
        Interpret values defined on a 3D mesh.

        Parameters
        ----------
        val : ndarray
            Input value(s).

        Returns
        -------
        ndarray
            Array of values with length equal to the number of mesh cells
            or faces.
        """
        val = np.array(val, dtype=np.float64)

        if val.size == 1:
            return val * np.ones(self.N)
        else:
            return val

    def points(self):
        """
        Return spatial locations associated with this variable.

        Returns
        -------
        ndarray
            Coordinates of mesh points corresponding to the variable
            definition.
        """
        if self.boundary is not None:
            return self.mesh.face_centers[self.mesh.boundaries[self.boundary]]

        elif self.region is not None:
            return self.mesh.cell_centers[self.mesh.regions[self.region]]

        else:
            return self.mesh.cell_centers

    def res2(self) -> np.ndarray:
        """
        Compute the normalized squared residual.

        Returns
        -------
        ndarray
            Normalized squared residual between current and previous values.
        """
        return np.sum((self.val - self.prev_val) ** 2) / np.sum(self.val ** 2)