import numpy as np

import cedar


class Direct(cedar.base.Adapter):
    """
    Direct variable adapter.

    This adapter performs a direct transfer of a variable value from a
    source variable to a target variable. It assumes that both variables
    represent the same physical quantity and therefore must be of the same
    concrete :class:`cedar.base.Var` subclass.
    """

    @classmethod
    def transfer(cls, src: cedar.base.Var, tgt: cedar.base.Var):
        """
        Transfer the value of a source variable directly to a target variable.

        Parameters
        ----------
        src : cedar.base.Var
            Source variable providing the value to be transferred.
        tgt : cedar.base.Var
            Target variable that will receive the value.

        Raises
        ------
        RuntimeError
            If the source and target variables are not of the same type.
        """

        # Enforce strict type equivalence between source and target variables
        if type(src) != type(tgt):
            cls.log_error(
                src,
                tgt,
                (
                    "Direct requires the src and tgt be the same variable type. "
                    f"type({src.model_name}.{src.name}) = {src.__class__.__name__} "
                    f"and type({tgt.model_name}.{tgt.name}) = {tgt.__class__.__name__}."
                ),
            )

        # Perform direct value assignment
        tgt.set(src.val)