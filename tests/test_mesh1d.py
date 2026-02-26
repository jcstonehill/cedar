import numpy as np
import pytest

import cedar   # adjust import as needed


@pytest.fixture
def simple_mesh() -> cedar.Mesh1D:
    mesh = cedar.Mesh1D(
        N_cells=4,
        L = 4.0,
        region="domain",
        start_boundary="left",
        end_boundary="right",
        start=np.array([0.0, 0.0, 0.0]),
        basis = "x"
    )
    mesh.build()
    return mesh

def test_basic_attributes(simple_mesh):
    mesh: cedar.Mesh1D = simple_mesh

    assert mesh.dim == 1
    assert mesh.N_cells == 4
    assert mesh.ds == pytest.approx(1.0)


def test_number_of_points(simple_mesh):
    mesh: cedar.Mesh1D = simple_mesh

    assert mesh.N_pts == mesh.N_cells + 1
    assert mesh.pts.shape == (5, 3)


def test_points_are_evenly_spaced(simple_mesh):
    mesh: cedar.Mesh1D = simple_mesh

    expected_x = np.linspace(0.0, 4.0, 5)
    np.testing.assert_allclose(mesh.pts[:, 0], expected_x)
    np.testing.assert_allclose(mesh.pts[:, 1:], 0.0)


def test_cell_point_connectivity(simple_mesh):
    mesh: cedar.Mesh1D = simple_mesh

    expected = np.array([
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 4],
    ])
    np.testing.assert_array_equal(mesh.cell_pts_i, expected)
    np.testing.assert_array_equal(mesh.cell_faces_i, expected)


def test_faces(simple_mesh):
    mesh: cedar.Mesh1D = simple_mesh

    assert mesh.N_faces == mesh.N_pts
    np.testing.assert_array_equal(mesh.face_pts_i, np.arange(mesh.N_pts))

    assert mesh.face_is_interior[0] == False
    assert mesh.face_is_interior[-1] == False
    assert np.all(mesh.face_is_interior[1:-1])


def test_cell_centers(simple_mesh):
    mesh: cedar.Mesh1D = simple_mesh

    expected_centers = np.array([
        [0.5, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [2.5, 0.0, 0.0],
        [3.5, 0.0, 0.0],
    ])
    np.testing.assert_allclose(mesh.cell_centers, expected_centers)


def test_face_centers(simple_mesh):
    mesh: cedar.Mesh1D = simple_mesh

    np.testing.assert_allclose(mesh.face_centers, mesh.pts)


def test_region_data(simple_mesh):
    mesh: cedar.Mesh1D = simple_mesh

    assert "domain" in mesh.region_i
    np.testing.assert_array_equal(mesh.region_i["domain"], np.arange(4))
    assert mesh.region_N["domain"] == 4


def test_boundary_data(simple_mesh):
    mesh: cedar.Mesh1D = simple_mesh

    assert "left" in mesh.boundary_i
    assert "right" in mesh.boundary_i

    np.testing.assert_array_equal(mesh.boundary_i["left"], np.array([0]))
    np.testing.assert_array_equal(mesh.boundary_i["right"], np.array([-1]))

    assert mesh.boundary_N["left"] == 1
    assert mesh.boundary_N["right"] == 1


def test_vtkhdf_dict_structure(simple_mesh):
    mesh: cedar.Mesh1D = simple_mesh
    vtk = mesh.vtkhdf_dict()

    assert f"{mesh.id}_NumberOfPoints" in vtk
    assert f"{mesh.id}_Points" in vtk
    assert "blocks" in vtk

    blocks = vtk["blocks"]
    assert "domain" in blocks
    assert "left" in blocks
    assert "right" in blocks


def test_vtkhdf_domain_block(simple_mesh):
    mesh: cedar.Mesh1D = simple_mesh
    block = mesh.vtkhdf_dict()["blocks"]["domain"]

    assert block["NumberOfCells"][0] == (mesh.N_cells,)
    assert np.all(block["Types"][0] == 3)  # VTK_LINE
    assert block["NumberOfConnectivityIds"][0] == (2 * mesh.N_cells,)


def test_vtkhdf_boundary_blocks(simple_mesh):
    mesh: cedar.Mesh1D = simple_mesh
    vtk = mesh.vtkhdf_dict()

    left = vtk["blocks"]["left"]
    right = vtk["blocks"]["right"]

    assert left["Types"][0] == (1,)   # VTK_VERTEX
    assert right["Types"][0] == (1,)

    assert left["Connectivity"][0] == (0,)
    assert right["Connectivity"][0] == (mesh.N_cells,)