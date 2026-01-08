from abc import ABC, abstractmethod
import cedar
from cedar.base.var import Var


class Adapter(ABC):
    """
    Abstract base class for data adapters between models.

    An Adapter defines a standardized interface for transferring data
    from a source variable to a target variable. It also provides
    logging utilities for tracking operations and errors during transfers.
    """

    name = "Adapter"

    @classmethod
    @abstractmethod
    def transfer(cls, src: Var, tgt: Var):
        """
        Transfer data from a source variable to a target variable.

        This is an abstract method that must be implemented by subclasses.

        Parameters
        ----------
        src : Var
            The source variable providing the data.
        tgt : Var
            The target variable receiving the data.

        Raises
        ------
        NotImplementedError
            If the subclass does not implement this method.
        """
        pass

    @classmethod
    def log_message(cls, src: Var, tgt: Var, message: str):
        """
        Log a standard informational message for the adapter.

        Parameters
        ----------
        src : Var
            The source variable involved in the operation.
        tgt : Var
            The target variable involved in the operation.
        message : str
            Custom message describing the operation or status.
        """
        cedar.Log.message(
            f"{cls.__name__} Adapter :: "
            f"{src.model_name}.{src.name} -> {tgt.model_name}.{tgt.name} :: {message}"
        )

    @classmethod
    def log_error(cls, src: Var, tgt: Var, message: str):
        """
        Log an error message for the adapter.

        Parameters
        ----------
        src : Var
            The source variable involved in the operation.
        tgt : Var
            The target variable involved in the operation.
        message : str
            Custom error message describing the issue.
        """
        cedar.Log.error(
            f"{cls.__name__} Adapter :: "
            f"{src.model_name}.{src.name} -> {tgt.model_name}.{tgt.name} :: {message}"
        )
