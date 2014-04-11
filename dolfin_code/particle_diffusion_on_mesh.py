import dolfin
import numpy as np


def convert_point_to_array(point_object):
  return np.array([point_object.x(), point_object.y(), point_object.z()])


def get_neighbor(facet, edge, facets_list):
  """Gets neighbor for a facet across a given edge."""
  neighbor_faces = []
  for index in edge.entities(2):
    # Skip current facet.
    if index == facet.index():
      continue
    # Skip non-exterior facets.
    if not facets_list[index].exterior():
      continue
    neighbor_faces.append(index)

  if len(neighbor_faces) != 1:
    raise ValueError('Wrong number of edge neighbors.')

  return neighbor_faces[0]


class FaceWrapper(object):

  def __init__(self, facet, parent_mesh_wrapper):
    self.facet = facet
    self.check_facet_type()

    self.parent_mesh_wrapper = parent_mesh_wrapper

    self.set_vertices()
    self.set_neighbors()
    self.compute_gram_schmidt_directions()

  def check_facet_type(self):
    if self.facet.dim() != 2:
      raise ValueError('Expected triangular facet.')

  def set_vertices(self):
    # This will fail if not exactly 3 vertices.
    self.a_index, self.b_index, self.c_index = self.facet.entities(0)

    vertex_list = self.parent_mesh_wrapper.vertex_list
    self.a = convert_point_to_array(vertex_list[self.a_index].point())
    self.b = convert_point_to_array(vertex_list[self.b_index].point())
    self.c = convert_point_to_array(vertex_list[self.c_index].point())

  def find_missing_vertex(self, index_list):
    return_values = []
    if self.a_index not in index_list:
      return_values.append('a')
    if self.b_index not in index_list:
      return_values.append('b')
    if self.c_index not in index_list:
      return_values.append('c')

    if len(return_values) != 1:
      raise ValueError('Edge did not contain exactly two of three vertices.')

    return return_values[0]

  def set_neighbors(self):
    """Sets neighbor facet indices.

    NOTE: This behavior is subject to change.

    Determines the face across an edge from each vertex a, b, c and
    labels accordingly.

    This ignores the neighbors which go through a vertex.
    """
    facet_edge_indices = self.facet.entities(1)

    edge_vertex_pairing = []
    # NOTE: We don't check that there are 3 edges, but could / should.
    for edge_index in facet_edge_indices:
      curr_edge = self.parent_mesh_wrapper.edge_list[edge_index]
      index_list = curr_edge.entities(0)
      vertex_str = self.find_missing_vertex(index_list)
      edge_vertex_pairing.append((vertex_str, curr_edge))

    # Sort pairs and check they are correct.
    edge_vertex_pairing.sort(key=lambda pair: pair[0])
    vertex_strs = tuple(pair[0] for pair in edge_vertex_pairing)
    if vertex_strs != ('a', 'b', 'c'):
      raise ValueError('Edges are not opposite the expected vertices.')

    # Get the edges, sorted now by the sides they are opposite to.
    sorted_edges = [pair[1] for pair in edge_vertex_pairing]

    # Find indices of faces opposite our vertices.
    facets_list = self.parent_mesh_wrapper.facets_list
    self.face_opposite_a = get_neighbor(self.facet, sorted_edges[0],
                                        facets_list)
    self.face_opposite_b = get_neighbor(self.facet, sorted_edges[1],
                                        facets_list)
    self.face_opposite_c = get_neighbor(self.facet, sorted_edges[2],
                                        facets_list)

  def compute_gram_schmidt_directions(self):
    v1 = self.b - self.a
    v2 = self.c - self.a
    n = convert_point_to_array(self.facet.normal())

    vec_as_cols = np.array([v1, v2, n]).T
    Q, _ = np.linalg.qr(vec_as_cols)
    # NOTE: We don't do this, but we could / should check that
    #       w3 = Q[2, :] is either n or -n (same direction, unit length).
    #       For example, we could do this by checking:
    #           np.allclose(w3.dot(n), 1) or np.allclose(w3.dot(n), -1)
    self.w1 = Q[:, 0]
    self.w2 = Q[:, 1]

  def compute_vertex_angle(self, delta):
    component1 = self.w1.dot(delta)
    component2 = self.w2.dot(delta)
    return np.arctan2(component2, component1)

  def compute_vertex_angles(self, point):
    a_theta = self.compute_vertex_angle(point - self.a)
    b_theta = self.compute_vertex_angle(point - self.b)
    c_theta = self.compute_vertex_angle(point - self.c)
    return a_theta, b_theta, c_theta


class MeshWrapper(object):

  def __init__(self, mesh):
    self._faces_added = False

    # Add the mesh and make sure it is the right kind of mesh.
    self.mesh = mesh
    self.mesh.init()  # This will do nothing if already called.
    self.check_mesh_type()

    # Compute extra data not stored on the object.
    self.vertex_list = list(dolfin.vertices(self.mesh))
    self.edge_list = list(dolfin.edges(self.mesh))
    self.facets_list = list(dolfin.facets(self.mesh))
    self.add_faces()

  def check_mesh_type(self):
    if self.mesh.geometry().dim() != 3:
      raise ValueError('Expecting 3D mesh.')
    if self.mesh.cells().shape[1] != 4:
      raise ValueError('Expecting tetrahedral mesh.')

  def add_faces(self):
    """Adds faces from mesh to stored list on object.

    Only adds exterior faces.
    """
    if self._faces_added:
      return

    self.faces = {}
    # Use facets since they have facet.exterior() set.
    for facet in self.facets_list:
      if not facet.exterior():
        continue
      self.faces[facet.index()] = FaceWrapper(facet, self)

    self._faces_added = True


class Point(object):

  def __init__(self, x, y, z, mesh_wrapper):
    pass


def sample_code():
  resolution = 96
  mesh_full_filename = 'mesh_res_%d_full.xml' % resolution
  mesh_3d = dolfin.Mesh(mesh_full_filename)
  return MeshWrapper(mesh_3d)
