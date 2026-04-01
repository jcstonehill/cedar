import h5py as h5
import numpy as np
import os
import time

import cedar


class Problem:
    def __init__(self):
        self.models: list[cedar.Model] = []

        self.meshes: dict[str, cedar.Mesh] = {}

        self.output: str = "output.vtkhdf"

        self.tol = 1e-6
        self.max_iter = 1000

        self.t_start: float = 0.0
        self.dt: float = 1.0
        self.t_end: float = -1.0

        self._is_built = False
        self._is_initialized = False

    def add_model(self, model: cedar.Model):
        self.models.append(model)

        if model.mesh.id not in self.meshes:
            self.meshes[model.mesh.id] = model.mesh

    def _build(self):
        if self._is_built:
            return

        cedar.Log.message("Building problem...")
        cedar.Log.add_level()

        # Get all unique meshes in problem.
        for model in self.models:
            if model.mesh.id not in self.meshes:
                self.meshes[model.mesh.id] = model.mesh

        cedar.Log.message("Building meshes...")
        cedar.Log.add_level()
        mesh_ids = self.meshes.keys()
        ordered_mesh_ids = sorted(mesh_ids)

        # Build all meshes
        for id in ordered_mesh_ids:
            cedar.Log.message(
                f"Mesh ID: {id:02d} ({self.meshes[id].__class__.__name__})", end=""
            )
            start = time.time()
            self.meshes[id]._build()
            end = time.time()
            cedar.Log.message(cedar.functions.format_computation_time(end - start))

        # All meshes are built, so now we can build the models
        for model in self.models:
            model._build()

            for field in model.all_fields.values():
                field._build(model.mesh)

            for boundary, bc in model.bcs.items():
                N = model.mesh.boundary_N[boundary]
                face_i = model.mesh.boundary_i[boundary]
                pts = model.mesh.face_centers[face_i]
                face_areas = model.mesh.face_areas[face_i]

                bc._build(pts, N, face_areas)

            for region, source in model.sources.items():
                N = model.mesh.region_N[region]
                cell_i = model.mesh.region_i[region]
                pts = model.mesh.cell_centers[cell_i]
                cell_vols = model.mesh.cell_vols[cell_i]

                source._build(pts, N, cell_vols)

        cedar.Log.remove_level()
        cedar.Log.message("Done")
        cedar.Log.remove_level()
        cedar.Log.line_break()

        self._is_built = True

    def _initialize(self):
        if self._is_initialized:
            return

        cedar.Log.message("Initializing...")
        cedar.Log.add_level()

        self.transfers: list[cedar.Transfer] = []

        # Initialize model.
        for model in self.models:
            model._initialize()

        # Initialize fields, BCs, and sources.
        for model in self.models:
            for field in model.all_fields.values():
                field._initialize()

            for source in model.sources.values():
                source._initialize()

            for bc in model.bcs.values():
                bc._initialize()

        # Now that everything is initialized, pull out transfers from BCs and
        # sources.
        for model in self.models:
            for source in model.sources.values():
                self.transfers.extend(source.transfers)

            for bc in model.bcs.values():
                self.transfers.extend(bc.transfers)

        cedar.Log.message("Done")
        cedar.Log.remove_level()
        cedar.Log.line_break()

        self._is_initialized = True

    def _post_commands(self):
        try:
            with open("post_commands.txt", "r") as f:
                for line in f.readlines():
                    os.system(line)

        except:
            pass

    def solve(self):
        self._build()
        self._initialize()

        is_steady = self.t_end < 0

        if self.output:
            self._create_output()
            self._append_metadata_to_output()

        if is_steady:
            cedar.Log.message("Solving to steady state...")
            cedar.Log.line_break()
            cedar.Log.add_level()

            start = time.time()
            self._solve_once(0)
            end = time.time()

            cedar.Log.line_break()
            cedar.Log.message(
                f"SOLVE COMPLETED     {cedar.functions.format_computation_time(end - start)}"
            )
            cedar.Log.line_break()

            cedar.Log.remove_level()

            if self.output:
                self._append_step_to_output(0)
                self._append_N_steps_to_output(1)

        else:
            cedar.Log.message("Solving transient...")
            cedar.Log.add_level()
            cedar.Log.message(f"{"t_start":<7}   = {self.t_start:>7.2f} [s]")
            cedar.Log.message(f"{"dt":<7}   = {self.dt:>7.2f} [s]")
            cedar.Log.message(f"{"t_end":<7}   = {self.t_end:>7.2f} [s]")
            cedar.Log.line_break()

            i = 0
            t = self.t_start

            if self.output:
                self._append_step_to_output(self.t_start)

            start = time.time()
            while t < self.t_end:
                t += self.dt
                i += 1
                cedar.Log.message(f"Time Step {i},   t = {t} [s],   dt = {self.dt} [s]")
                cedar.Log.add_level()

                self._step()
                self._solve_once(t, self.dt)

                cedar.Log.remove_level()
                cedar.Log.line_break()

                if self.output:
                    self._append_step_to_output(t)

            end = time.time()

            cedar.Log.message(
                f"SOLVE COMPLETED     {cedar.functions.format_computation_time(end - start)}"
            )
            cedar.Log.line_break()

            cedar.Log.remove_level()

            if self.output:
                self._append_N_steps_to_output(i)

        self._post_commands()

    def _step(self):
        for model in self.models:

            # Cell fields are the only ones that have initial conditions, so
            # they're the only ones that need to step
            for field in model.cell_fields.values():
                field._step()

    def _solve_once(self, t: float, dt: float = None) -> bool:
        for i in range(self.max_iter):
            residuals = []

            for model in self.models:
                for source in model.sources.values():
                    source._update(t)

                for bc in model.bcs.values():
                    bc._update(t)

                residuals.append(model._iterate(1e-2, dt))

            for transfer in self.transfers:
                residuals.append(transfer.residual())

            res = np.max(residuals)

            cedar.Log.message(f"{i+1:>4} |R|: {res:.5e}")

            if res <= self.tol:
                return i + 1

        cedar.Log.error(
            f"Problem was not solved using the maximum number of iterations: {self.max_iter}"
        )

    def _append_metadata_to_output(self):
        with h5.File(self.output, "a") as f:
            assembly = f["VTKHDF"]["Assembly"]

            block: h5.Group
            for block in assembly.values():
                steps = block.create_group("Steps")
                steps.create_dataset(
                    "Values", shape=(0,), maxshape=(None,), dtype="f"
                )  # time values

                steps.create_dataset(
                    "PartOffsets", shape=(0,), maxshape=(None,), dtype="i8"
                )
                steps.create_dataset(
                    "NumberOfParts", shape=(0,), maxshape=(None,), dtype="i8"
                )
                steps.create_dataset(
                    "ConnectivityIdOffsets", shape=(0,), maxshape=(None,), dtype="i8"
                )
                steps.create_dataset(
                    "CellOffsets", shape=(0,), maxshape=(None,), dtype="i8"
                )
                steps.create_dataset(
                    "PointOffsets", shape=(0,), maxshape=(None,), dtype="i8"
                )

                steps.create_group("CellDataOffsets")

            for model in self.models:
                # Cell Fields
                for field_name in model.cell_fields:
                    if field_name in model.output_fields:
                        for block in model.mesh.regions:
                            assembly[block]["CellData"].create_dataset(
                                field_name, (0,), maxshape=(None,), dtype="f"
                            )
                            assembly[block]["Steps"]["CellDataOffsets"].create_dataset(
                                field_name, (0,), maxshape=(None,), dtype="f"
                            )

                # Boundary Face Fields
                for field_name in model.boundary_face_fields:
                    if field_name in model.output_fields:
                        for block in model.mesh.boundaries:
                            assembly[block]["CellData"].create_dataset(
                                field_name, (0,), maxshape=(None,), dtype="f"
                            )
                            assembly[block]["Steps"]["CellDataOffsets"].create_dataset(
                                field_name, (0,), maxshape=(None,), dtype="f"
                            )

    def _append_N_steps_to_output(self, N_steps: int):
        with h5.File(self.output, "a") as f:
            assembly = f["VTKHDF"]["Assembly"]

            for block in assembly.values():
                block["Steps"].attrs["NSteps"] = N_steps

        cedar.Log.message(f"Output File: {self.output}")
        cedar.Log.line_break()

    def _append_step_to_output(self, t: float):
        with h5.File(self.output, "a") as f:
            assembly = f["VTKHDF"]["Assembly"]

            for block in assembly.values():
                steps = block["Steps"]

                self._extend_dataset(steps["Values"], (t,))

                self._extend_dataset(steps["PartOffsets"], (0,))
                self._extend_dataset(steps["NumberOfParts"], (1,))
                self._extend_dataset(steps["PointOffsets"], (0,))
                self._extend_dataset(steps["CellOffsets"], (0,))
                self._extend_dataset(steps["ConnectivityIdOffsets"], (0,))

            for model in self.models:
                # Cell Fields
                for field_name, field in model.cell_fields.items():
                    if field_name in model.output_fields:
                        for block in model.mesh.regions:
                            data = field.get(block)

                            offsets = assembly[block]["Steps"]["CellDataOffsets"]

                            self._extend_dataset(
                                offsets[field_name],
                                (assembly[block]["CellData"][field_name].shape[0],),
                            )

                            cell_data = assembly[block]["CellData"]
                            self._extend_dataset(cell_data[field_name], data)

                # Boundary Face Fields
                for field_name, field in model.boundary_face_fields.items():
                    if field_name in model.output_fields:
                        for block in model.mesh.boundaries:
                            data = field.get(block)

                            offsets = assembly[block]["Steps"]["CellDataOffsets"]

                            self._extend_dataset(
                                offsets[field_name],
                                (assembly[block]["CellData"][field_name].shape[0],),
                            )

                            cell_data = assembly[block]["CellData"]
                            self._extend_dataset(cell_data[field_name], data)

    def _create_output(self):
        VERSION = (2, 5)

        with h5.File(self.output, "w") as f:
            root = f.create_group("VTKHDF", track_order=True)
            root.attrs["Version"] = VERSION
            root.attrs["Type"] = "MultiBlockDataSet"

            assembly = root.create_group("Assembly", track_order=True)

            block_index = 0

            for mesh in self.meshes.values():
                vtkhdf_dict = mesh._vtkhdf_dict()

                N_pts_key = None
                pts_key = None

                key: str
                for key, value in vtkhdf_dict.items():
                    if key.endswith("_NumberOfPoints"):
                        N_pts_key = key
                        root.create_dataset(key, data=value, dtype="i8")

                    if key.endswith("_Points"):
                        pts_key = key
                        root.create_dataset(key, data=value, dtype="f")

                for block_name, block_dict in vtkhdf_dict["blocks"].items():
                    block = root.create_group(block_name, track_order=True)
                    block.attrs["Version"] = VERSION
                    block.attrs["Type"] = "UnstructuredGrid"
                    block.attrs["Index"] = block_index
                    block_index += 1

                    block["NumberOfPoints"] = h5.SoftLink("/VTKHDF/" + N_pts_key)
                    block["Points"] = h5.SoftLink("/VTKHDF/" + pts_key)

                    block.create_group("CellData", track_order=True)

                    for key, (value, dtype) in block_dict.items():
                        block.create_dataset(key, data=value, dtype=dtype)

                    assembly[block_name] = h5.SoftLink("/VTKHDF/" + block_name)

                    # steps = block.create_group("Steps")
                    # steps.create_dataset(
                    #     "Values", shape=(0,), maxshape=(None,), dtype="f"
                    # )  # time values

                    # steps.create_dataset("PartOffsets", shape=(0,), maxshape=(None,), dtype="i8")
                    # steps.create_dataset("NumberOfParts", shape=(0,), maxshape=(None,), dtype="i8")
                    # steps.create_dataset(
                    #     "ConnectivityIdOffsets", shape=(0,), maxshape=(None,), dtype="i8"
                    # )
                    # steps.create_dataset("CellOffsets", shape=(0,), maxshape=(None,), dtype="i8")
                    # steps.create_dataset("PointOffsets", shape=(0,), maxshape=(None,), dtype="i8")

            #         steps.create_group("CellDataOffsets")

            # for model in self.models:
            #     for field_name in model.output_fields:
            #         field = model.fields[field_name]

            #         for block, data in field.values.items():
            #             root[block]["CellData"].create_dataset(field_name, (0,), maxshape=(None,), dtype="f")
            #             root[block]["Steps"]["CellDataOffsets"].create_dataset(field_name, (0,), maxshape=(None,), dtype="f")

    def _extend_dataset(self, dataset: h5.Dataset, array: np.ndarray):
        """Resize a chunked dataset, adding `array` at the end"""
        original_size = dataset.shape[0]
        dataset.resize(original_size + len(array), axis=0)
        dataset[original_size:] = array
