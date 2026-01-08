import numpy as np

import cedar


class Summation(cedar.base.Adapter):
    """
    Summation-based variable adapter.

    This adapter transfers data by *summing* source values and mapping the
    result to a lower-dimensional target variable.

    Notes
    -----
    - No normalization or averaging is performed.
    - Spatial aggregation for ``Mesh3D → Mesh1D`` is one-directional along
      the z-axis.
    - This adapter is typically used for conserved quantities (e.g. mass,
      energy) where summation is physically meaningful.
    - ``Mesh3DVar → ScalarVar``:
      Sum of all values in the 3D mesh.
    - ``Mesh3DVar → Mesh1DVar``:
      Values are summed into 1D bins based on the z-location of each 3D point
      relative to the target mesh face centers.
    - ``Mesh1DVar → ScalarVar``:
      Sum of all values in the 1D mesh.
    """

    @classmethod
    def transfer(cls, src: cedar.base.Var, tgt: cedar.base.Var):
        """
        Transfer values using summation-based aggregation.

        Parameters
        ----------
        src : cedar.base.Var
            Source variable. Must be a :class:`cedar.Mesh3DVar` or
            :class:`cedar.Mesh1DVar`.
        tgt : cedar.base.Var
            Target variable. Must be compatible with the selected summation
            strategy.

        Raises
        ------
        RuntimeError
            If the source/target variable combination is not supported.
        """

        # Mesh3D source: reduce dimensionality via summation
        if isinstance(src, cedar.Mesh3DVar):
            if isinstance(tgt, cedar.ScalarVar):
                # Sum all values in the 3D mesh
                data_for_tgt = np.sum(src.val)

            elif isinstance(tgt, cedar.Mesh1DVar):
                # Initialize accumulation array for 1D bins
                data_for_tgt = np.zeros(tgt.N)

                # z-coordinates of interior face centers in the target mesh
                face_z = [
                    center[2]
                    for center in tgt.mesh.face_centers[tgt.mesh.face_is_interior]
                ]

                # Accumulate 3D values into the appropriate 1D bin
                for point, val_from_src in zip(src.points(), src.val):
                    for i, z in enumerate(face_z):
                        if point[2] <= z:
                            data_for_tgt[i] += val_from_src
                            break
                    else:
                        # Values above the highest face center go to the last bin
                        data_for_tgt[-1] += val_from_src

            else:
                cls.log_error(
                    src,
                    tgt,
                    "When source is a Mesh3DVar, Summation can only transfer "
                    "to a Mesh1DVar or a ScalarVar.",
                )

        # Mesh1D source: total summation
        elif isinstance(src, cedar.Mesh1DVar):
            if isinstance(tgt, cedar.ScalarVar):
                # Sum all values in the 1D mesh
                data_for_tgt = np.sum(src.val)

            else:
                cls.log_error(
                    src,
                    tgt,
                    "When source is a Mesh1DVar, Summation can only transfer "
                    "to a ScalarVar.",
                )

        else:
            cls.log_error(
                src,
                tgt,
                "Source must be a Mesh3DVar or a Mesh1DVar",
            )

        # Assign aggregated result to the target variable
        tgt.set(data_for_tgt)