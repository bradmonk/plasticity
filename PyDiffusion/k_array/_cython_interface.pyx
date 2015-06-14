"""Module for code intended to be used from MATLAB.

The goal is to hide all Python OOP and just use data (i.e.
floats, ints, matrices) that both languages (and C) shared
happily.
"""

import itertools
import numpy as np
cimport cython
from cython cimport view

from particle_diffusion_on_mesh import Mesh
from particle_diffusion_on_mesh import Point


@cython.boundscheck(False)
@cython.wraparound(False)
def advance_one_step(double[:, ::1] _xyz_loc, long[:, ::1] _face_indices,
                     double[:, ::1] _k, double[::1] _initial_point,
                     long initial_face_index, double[:, ::1] _all_vertices,
                     long[:, ::1] _triangles, double[:, ::1] _face_local_bases,
                     long[:, ::1] _neighbor_faces):
    """Custom method for advancing simulation by one-step.

    This is a bare-bones method which accepts and returns only simple types
    (matrices) to avoid surfacing Python OOP in it's interface.

    The main values being advanced are the points and faces corresponding to
    those points. In addition, the remaining (7) arguments are used to
    construct a `Mesh` object.

    Args:
        _xyz_loc: An Mx3 array (2D) of point locations.
        _face_indices: An Mx1 vector (2D so MATLAB is happy) of face indices
                       corresponding to each point in `_xyz_loc`.
        k: Float, the "global" diffusion rate in the mesh.
        _initial_point: A 3-vector (1D) of floats containing the initial
                        starting point on the mesh.
        initial_face_index: Integer. The face containing `_initial_point`.
        _all_vertices: A Vx3 matrix of floats containing all the vertices in
                       the mesh.
        _triangles: A Tx3 matrix of integers containing the vertex indices of
                    each triple of vertices in a mesh triangle exterior face.
        _face_local_bases: A Tx6 matrix of floats containing a pre-computed
                           orthogonal basis of the plane going through each
                           triangle. (The rows of `_face_local_bases` correspond
                           to the rows of `_triangles`).
        _neighbor_faces: A Tx3 matrix of integers containing indices of the
                         neighbors of each face in `_triangles`. A triangle has
                         3 sides hence 3 neighbors. (The rows of
                         `_neighbor_faces` correspond to the rows of
                         `_triangles`).

    Returns:
        Returns two arrays, the updated version of _xyz_loc and _face_indices.
    """
    face_indices = np.asarray(_face_indices)
    # NOTE: If passed in from C, columns and rows will be swapped.
    if face_indices.shape[0] == 1:
        face_indices = face_indices.T

    xyz_loc = np.asarray(_xyz_loc)
    # NOTE: If passed in from C, columns and rows will be swapped.
    if xyz_loc.shape[0] == 3:
        xyz_loc = xyz_loc.T

    k = np.asarray(_k)
    # NOTE: If passed in from C, columns and rows will be swapped.
    if k.shape[0] == 1:
        k = k.T

    initial_point = np.asarray(_initial_point)

    all_vertices = np.asarray(_all_vertices)
    # NOTE: If passed in from C, columns and rows will be swapped.
    if all_vertices.shape[0] == 3:
        all_vertices = all_vertices.T

    triangles = np.asarray(_triangles)
    # NOTE: If passed in from C, columns and rows will be swapped.
    if triangles.shape[0] == 3:
        triangles = triangles.T

    face_local_bases = np.asarray(_face_local_bases)
    # NOTE: If passed in from C, columns and rows will be swapped.
    if face_local_bases.shape[0] == 6:
        face_local_bases = face_local_bases.T

    neighbor_faces = np.asarray(_neighbor_faces)
    # NOTE: If passed in from C, columns and rows will be swapped.
    if neighbor_faces.shape[0] == 3:
        neighbor_faces = neighbor_faces.T

    mesh_wrapper = Mesh(k, initial_point, initial_face_index, all_vertices,
                        triangles, face_local_bases, neighbor_faces)

    # NOTE: We assume but don't check that xyz_loc and face_indices
    #       have the same number of rows.
    points = []
    for pt, face_index in itertools.izip(xyz_loc, face_indices):
        face_index, = face_index  # Unpacking a row of an Mx1 matrix.
        points.append(Point(mesh_wrapper, point=pt, face_index=face_index))

    # After all points are set, we are OK to move them.
    for pt in points:
        pt.move()

    cdef int i
    for i, pt in enumerate(points):
        # Update the xyz locations in place.
        xyz_loc[i, 0] = pt.point[0]
        xyz_loc[i, 1] = pt.point[1]
        xyz_loc[i, 2] = pt.point[2]
        # Update the face indices in place.
        face_indices[i, 0] = pt.face.face_index


cdef public void advance_one_step_c(
        size_t num_points, size_t num_vertices, size_t num_triangles,
        double* xyz_loc, long* face_indices, double* k, double* initial_point,
        long initial_face_index, double* all_vertices,
        long* triangles, double* face_local_bases, long* neighbor_faces):
    # NOTE: This is very crucial. We swap rows and columns since the data is
    #       in Fortran order (from MATLAB) but the views are in C order.
    cdef view.array py_k = <double[:1, :num_points]> k
    cdef view.array py_face_indices = <long[:1, :num_points]> face_indices
    cdef view.array py_xyz_loc = <double[:3, :num_points]> xyz_loc
    cdef view.array py_initial_point = <double[:3]> initial_point
    cdef view.array py_all_vertices = <double[:3, :num_vertices]> all_vertices
    cdef view.array py_triangles = <long[:3, :num_triangles]> triangles
    cdef view.array py_face_local_bases = \
        <double[:6, :num_triangles]> face_local_bases
    cdef view.array py_neighbor_faces = \
        <long[:3, :num_triangles]> neighbor_faces
    advance_one_step(py_xyz_loc, py_face_indices, py_k, py_initial_point,
                     initial_face_index, py_all_vertices, py_triangles,
                     py_face_local_bases, py_neighbor_faces)
