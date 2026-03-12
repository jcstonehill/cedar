from __future__ import annotations
from abc import ABC, abstractmethod
import numpy as np

import cedar


class Transfer:
    def __init__(self, src: cedar.Field | cedar.FieldView, mapping: np.ndarray):
        self.src: cedar.Field | cedar.FieldView = src
        self.mapping: np.ndarray = mapping

        self.prev = None
        self.current = None

    def get(self) -> np.ndarray:
        self.prev = self.current
        self.current = self.mapping @ self.src.get()

        return self.current

    def residual(self) -> float:
        if self.prev is None:
            return 1e12

        return np.sum(np.abs((self.current - self.prev) / self.current))

    @classmethod
    def with_nearest_value_mapping(
        cls, src: cedar.Field | cedar.FieldView, tgt_pts: np.ndarray
    ) -> "Transfer":

        if isinstance(src, cedar.Field):
            src_pts = src.mesh.cell_centers

        elif isinstance(src, cedar.FieldView):
            if src.domain in src.mesh.regions:
                src_pts = src.mesh.cell_centers[src.mesh.region_i[src.domain]]

            else:
                src_pts = src.mesh.face_centers[src.mesh.boundary_i[src.domain]]

        # pairwise distances: (n_tgt, n_src)
        distances = np.linalg.norm(src_pts[None, :, :] - tgt_pts[:, None, :], axis=2)

        nearest = np.argmin(distances, axis=1)

        mapping = np.zeros((tgt_pts.shape[0], src_pts.shape[0]))
        mapping[np.arange(tgt_pts.shape[0]), nearest] = 1

        return Transfer(src, mapping)

    @classmethod
    def with_layered_area_avg_mapping(
        cls, src: cedar.FieldView, layer_mesh: cedar.Mesh1D
    ) -> Transfer:
        mapping = np.zeros((layer_mesh.N_cells, src.mesh.boundary_N[src.domain]))
        face_i = src.mesh.boundary_i[src.domain]

        for i in range(src.mesh.boundary_N[src.domain]):
            face_pts = src.mesh.pts[src.mesh.face_pts_i[face_i[i]]]
            area_total = src.mesh.face_areas[face_i[i]]

            for j in range(layer_mesh.N_cells):
                area_overlap = cedar.functions.face_area_1d_overlap(
                    face_pts, layer_mesh.pts[j], layer_mesh.pts[j + 1]
                )
                mapping[j][i] += area_overlap / area_total

        for i in range(layer_mesh.N_cells):
            mapping[i][:] = mapping[i][:] / np.sum(mapping[i][:])

        return Transfer(src, mapping)

    @classmethod
    def with_layered_summation_mapping(
        cls, src: cedar.Field | cedar.FieldView, layer_mesh: cedar.Mesh1D
    ) -> "Transfer":
        mapping = np.zeros((layer_mesh.N_cells, src.mesh.boundary_N[src.domain]))
        face_i = src.mesh.boundary_i[src.domain]

        for i in range(src.mesh.boundary_N[src.domain]):
            face_pts = src.mesh.pts[src.mesh.face_pts_i[face_i[i]]]
            area_total = src.mesh.face_areas[face_i[i]]

            for j in range(layer_mesh.N_cells):
                area_overlap = cedar.functions.face_area_1d_overlap(
                    face_pts, layer_mesh.pts[j], layer_mesh.pts[j + 1]
                )
                mapping[j][i] += area_overlap / area_total

        return Transfer(src, mapping)

    @classmethod
    def with_summation_mapping(cls, src: cedar.Field | cedar.FieldView) -> "Transfer":
        # TODO This shouldn't only work with boundaries, we should make it work with regions or fields too.
        mapping = np.ones(src.mesh.boundary_N[src.domain], dtype=np.float64)

        return Transfer(src, mapping)
