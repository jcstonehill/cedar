import numpy as np

import cedar


class NearestValue(cedar.base.Adapter):
    """
    Nearest-neighbor value transfer adapter.

    This adapter transfers values from a lower-dimensional or coarser variable
    to a higher-dimensional target variable by selecting the *nearest available
    source value* for each target location.

    Notes
    -----
    - No interpolation is performed.
    - Mapping is based solely on positional proximity.
    - This adapter is useful when source resolution is insufficient for
      interpolation but nearest-value mapping is acceptable.
    """

    @classmethod
    def transfer(cls, src: cedar.base.Var, tgt: cedar.base.Var):
        """
        Transfer values using a nearest-value strategy.

        Parameters
        ----------
        src : cedar.base.Var
            Source variable. Must be either a :class:`cedar.ScalarVar` or
            :class:`cedar.Mesh1DVar`.
        tgt : cedar.base.Var
            Target variable. Must be compatible with the selected transfer
            strategy.

        Raises
        ------
        RuntimeError
            If the source/target variable combination is not supported.
        """

        # Scalar → any mesh variable: broadcast scalar value
        if isinstance(src, cedar.ScalarVar):
            data_for_tgt = src.val * np.ones(tgt.N)

        # Mesh1D → Mesh3D: nearest-value mapping along the z-axis
        elif isinstance(src, cedar.Mesh1DVar):
            if isinstance(tgt, cedar.Mesh3DVar):
                data_for_tgt = np.zeros(tgt.N)

                # Target spatial coordinates (x, y, z) for each cell/node
                points = tgt.points()

                # Source 1D mesh used to locate nearest z-index
                mesh: cedar.Mesh1D = src.mesh

                for j in range(tgt.N):
                    # Find nearest 1D mesh index for the target z-position
                    i = mesh.get_i_for_z_pos(points[j][2])
                    data_for_tgt[j] = src.val[i]

            else:
                cls.log_error(
                    src,
                    tgt,
                    "When source is a Mesh1DVar, NearestValue can only transfer to a Mesh3DVar.",
                )

        else:
            cls.log_error(
                src,
                tgt,
                "Source must be a ScalarVar or a Mesh1DVar.",
            )

        # Assign computed values to the target variable
        tgt.set(data_for_tgt)