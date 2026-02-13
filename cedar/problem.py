import h5py as h5
import numpy as np
import time

import cedar


class Problem:
    def __init__(
        self,
        models: list[cedar.Model] = [],
    ):
        self.models: list[cedar.Model] = models

        self.meshes: dict[str, cedar.Mesh] = {}

        self.output: str = "output.vtkhdf"

        self.tol = 1e-6
        self.max_iter = 1000

        self.t_start: float = 0.0
        self.dt: float = 1.0
        self.t_end: float = -1.0

        self._is_built = False
        self._is_initialized = False

    def build(self):
        if self._is_built:
            return
        
        cedar.Log.message("Building problem...")
        cedar.Log.add_level()

        cedar.Log.message("Retrieving unique meshes")
        # Get all unique meshes
        for model in self.models:
            if model.mesh.id not in self.meshes:
                self.meshes[model.mesh.id] = model.mesh

        cedar.Log.message("Building meshes...")
        cedar.Log.add_level()
        mesh_ids = self.meshes.keys()
        ordered_mesh_ids = sorted(mesh_ids)

        # Build all meshes
        for id in ordered_mesh_ids:
            cedar.Log.message(f"Mesh ID: {id:02d} ({self.meshes[id].__class__.__name__})", end = "")
            start = time.time()
            self.meshes[id].build()
            end = time.time()
            cedar.Log.message(cedar.functions.format_computation_time(end-start))
            
        cedar.Log.remove_level()
        cedar.Log.message("Done")
        cedar.Log.remove_level()
        cedar.Log.line_break()

        self._is_built = True

    def initialize(self):
        if self._is_initialized:
            return
        
        cedar.Log.message("Initializing...")
        cedar.Log.add_level()
        for model in self.models:
            model.initialize(self.t_start)

            model.source.mesh = model.mesh
            model.source.initialize(self.t_start)

            for field in model.fields.values():
                field.mesh = model.mesh
                field.initialize(self.t_start)

            for boundary, bc in model.bc.items():
                bc.boundary = boundary
                bc.mesh = model.mesh

                bc.initialize(self.t_start)

        cedar.Log.message("Done")
        cedar.Log.remove_level()
        cedar.Log.line_break()

        self._is_initialized = True

    # def check(self):
    #     for model in self.models:
    #         model.check()

    def solve(self):
        self.build()
        self.initialize()
        # self.check()

        is_steady = (self.t_end < 0)

        if is_steady:
            cedar.Log.message("Solving to steady state...")
            cedar.Log.line_break()
            cedar.Log.add_level()

            self._create_output()
            
            start = time.time()
            self._solve_once()
            end = time.time()

            cedar.Log.line_break()
            cedar.Log.message(f"SOLVE COMPLETED     {cedar.functions.format_computation_time(end - start)}")
            cedar.Log.line_break()

            cedar.Log.remove_level()

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

            self._create_output()
            self._append_step_to_output(self.t_start)

            start = time.time()
            while t < self.t_end:
                t += self.dt
                i += 1
                cedar.Log.message(f"Time Step {i},   t = {t} [s],   dt = {self.dt} [s]")
                cedar.Log.add_level()

                for model in self.models:
                    for bc in model.bc.values():
                        bc.step(t)

                    for field in model.fields.values():
                        field.step(t)

                    model.source.step(t)

                self._solve_once(self.dt)

                cedar.Log.remove_level()

                cedar.Log.line_break()

                self._append_step_to_output(t)

            end = time.time()

            cedar.Log.message(f"SOLVE COMPLETED     {cedar.functions.format_computation_time(end - start)}")
            cedar.Log.line_break()

            cedar.Log.remove_level()

            self._append_N_steps_to_output(i)
            
    def _solve_once(self, dt = None) -> bool:
        for i in range(self.max_iter):
            residuals = []

            for model in self.models:
                residuals.append(model.iterate(1e-2, dt))

            res = np.max(residuals)

            cedar.Log.message(f"{i+1:>4} |R|: {res:.5e}")

            if res <= self.tol:
                return i+1
            
        cedar.Log.error(f"Problem was not solved using the maximum number of iterations: {self.max_iter}")

    def _append_N_steps_to_output(self, N_steps: int):
        if self.output is None:
            return
        
        with h5.File(self.output, "a") as f:
            assembly = f["VTKHDF"]["Assembly"]
            
            for block in assembly.values():
                block["Steps"].attrs["NSteps"] = N_steps

        cedar.Log.message(f"Output File: {self.output}")
        cedar.Log.line_break()

    def _append_step_to_output(self, t: float):
        if self.output is None:
            return
        
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
                for field_name in model.output_fields:
                    field = model.fields[field_name]
                    
                    for block_name, data in field.values.items():
                        offsets = assembly[block_name]["Steps"]["CellDataOffsets"]
                        self._extend_dataset(offsets[field_name], (assembly[block_name]["CellData"][field_name].shape[0],))

                        cell_data = assembly[block_name]["CellData"]
                        self._extend_dataset(cell_data[field_name], data)

    def _create_output(self):
        if self.output is None:
            return
        
        VERSION = (2,5)

        with h5.File(self.output, "w") as f:
            root = f.create_group("VTKHDF", track_order = True)
            root.attrs["Version"] = VERSION
            root.attrs["Type"] = "MultiBlockDataSet"

            assembly = root.create_group("Assembly", track_order = True)

            block_index = 0

            for mesh in self.meshes.values():
                vtkhdf_dict = mesh.vtkhdf_dict()

                N_pts_key = None
                pts_key = None

                key: str
                for key, value in vtkhdf_dict.items():
                    if key.endswith("_NumberOfPoints"):
                        N_pts_key = key
                        root.create_dataset(key, data = value, dtype = "i8")

                    if key.endswith("_Points"):
                        pts_key = key
                        root.create_dataset(key, data = value, dtype = "f")

                for block_name, block_dict in vtkhdf_dict["blocks"].items():
                    block = root.create_group(block_name, track_order = True)
                    block.attrs["Version"] = VERSION
                    block.attrs["Type"] = "UnstructuredGrid"
                    block.attrs["Index"] = block_index
                    block_index += 1

                    block["NumberOfPoints"] = h5.SoftLink("/VTKHDF/" + N_pts_key)
                    block["Points"] = h5.SoftLink("/VTKHDF/" + pts_key)
                    
                    block.create_group("CellData", track_order = True)

                    for key, (value, dtype) in block_dict.items():
                        block.create_dataset(key, data = value, dtype = dtype)

                    assembly[block_name] = h5.SoftLink("/VTKHDF/" + block_name)

                    steps = block.create_group("Steps")
                    steps.create_dataset(
                        "Values", shape=(0,), maxshape=(None,), dtype="f"
                    )  # time values

                    steps.create_dataset("PartOffsets", shape=(0,), maxshape=(None,), dtype="i8")
                    steps.create_dataset("NumberOfParts", shape=(0,), maxshape=(None,), dtype="i8")
                    steps.create_dataset(
                        "ConnectivityIdOffsets", shape=(0,), maxshape=(None,), dtype="i8"
                    )
                    steps.create_dataset("CellOffsets", shape=(0,), maxshape=(None,), dtype="i8")
                    steps.create_dataset("PointOffsets", shape=(0,), maxshape=(None,), dtype="i8")

                    steps.create_group("CellDataOffsets")

            for model in self.models:
                for field_name in model.output_fields:
                    field = model.fields[field_name]

                    for block, data in field.values.items():
                        root[block]["CellData"].create_dataset(field_name, (0,), maxshape=(None,), dtype="f")
                        root[block]["Steps"]["CellDataOffsets"].create_dataset(field_name, (0,), maxshape=(None,), dtype="f")

    def _extend_dataset(self, dataset: h5.Dataset, array: np.ndarray):
        """Resize a chunked dataset, adding `array` at the end"""
        original_size = dataset.shape[0]
        dataset.resize(original_size + len(array), axis=0)
        dataset[original_size:] = array   