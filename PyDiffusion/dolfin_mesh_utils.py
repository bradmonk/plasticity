import dolfin
import numpy as np

from particle_diffusion_on_mesh import Mesh


def get_face(point, mesh, cell_list):
  face_intersection = dolfin.cpp.mesh.intersect(mesh, point)
  # Require a unique intersection. This is not actually necessary as some
  # points may lie on an edge / vertex or may not be on the mesh at all.
  if face_intersection.intersected_cells().size != 1:
    raise ValueError('Point does not intersect mesh in a unique cell.')
  cell_index = face_intersection.intersected_cells()[0]

  matched_cell = cell_list[cell_index]
  exterior_facets = [facet for facet in dolfin.facets(matched_cell)
                     if facet.exterior()]
  if len(exterior_facets) != 1:
    raise ValueError('Number of facets on cell marked exterior != 1.')
    print 'exterior_facets:', exterior_facets

  return exterior_facets[0].index()


def check_facet_type(facet):
  if facet.dim() != 2:
    raise ValueError('Expected triangular facet.')


def check_mesh_type(mesh):
  if mesh.geometry().dim() != 3:
    raise ValueError('Expecting 3D mesh.')
  if mesh.cells().shape[1] != 4:
    raise ValueError('Expecting tetrahedral mesh.')


def convert_point_to_array(point_object):
  return np.array([point_object.x(), point_object.y(), point_object.z()])


def get_face_vertices(facet, vertex_list):
  # This will fail if not exactly 3 vertices.
  a_index, b_index, c_index = facet.entities(0)

  a = convert_point_to_array(vertex_list[a_index].point())
  b = convert_point_to_array(vertex_list[b_index].point())
  c = convert_point_to_array(vertex_list[c_index].point())

  return a, b, c, a_index, b_index, c_index


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


def find_missing_vertex(a_index, b_index, c_index, index_list):
  return_values = []
  if a_index not in index_list:
    return_values.append('a')
  if b_index not in index_list:
    return_values.append('b')
  if c_index not in index_list:
    return_values.append('c')

  if len(return_values) != 1:
    raise ValueError('Edge did not contain exactly two of three vertices.')

  return return_values[0]


def get_neighbor_faces(facet, edge_list, facets_list,
                       a_index, b_index, c_index):
  """Sets neighbor facet indices.

  NOTE: This behavior is subject to change.

  Determines the face across an edge from each vertex a, b, c and
  labels accordingly.

  This ignores the neighbors which go through a vertex.
  """
  facet_edge_indices = facet.entities(1)

  edge_vertex_pairing = []
  for edge_index in facet_edge_indices:
    curr_edge = edge_list[edge_index]
    index_list = curr_edge.entities(0)
    vertex_str = find_missing_vertex(a_index, b_index, c_index, index_list)
    edge_vertex_pairing.append((vertex_str, curr_edge))

  # Sort pairs and check they are correct.
  edge_vertex_pairing.sort(key=lambda pair: pair[0])
  vertex_strs = tuple(pair[0] for pair in edge_vertex_pairing)
  if vertex_strs != ('a', 'b', 'c'):
    raise ValueError('Edges are not opposite the expected vertices.')

  # Get the edges, sorted now by the sides they are opposite to.
  sorted_edges = [pair[1] for pair in edge_vertex_pairing]

  # Find indices of faces opposite our vertices.
  face_opposite_a = get_neighbor(facet, sorted_edges[0], facets_list)
  face_opposite_b = get_neighbor(facet, sorted_edges[1], facets_list)
  face_opposite_c = get_neighbor(facet, sorted_edges[2], facets_list)
  return face_opposite_a, face_opposite_b, face_opposite_c


def get_face_properties(facet, vertex_list, edge_list, facets_list):
  a, b, c, a_index, b_index, c_index = get_face_vertices(facet, vertex_list)

  face_opposite_a, face_opposite_b, face_opposite_c = get_neighbor_faces(
      facet, edge_list, facets_list, a_index, b_index, c_index)

  return [
      (a, a_index, face_opposite_a),
      (b, b_index, face_opposite_b),
      (c, c_index, face_opposite_c),
  ]


def compute_gram_schmidt_directions(facet, directed_side):
  w1 = directed_side / np.linalg.norm(directed_side)
  n = convert_point_to_array(facet.normal())
  # n acts as z == e3, w1 acts as x == e1, hence y == e2 (== e3 x e1):
  w2 = np.cross(n, w1)
  return w1, w2


def overwrite_numpy_value(key, value, curr_dict):
  if key in curr_dict:
    if not np.allclose(value, curr_dict[key]):
      raise ValueError('Values do not match.')
  else:
    curr_dict[key] = value


class RawFaceData(object):

  def __init__(self, facet, vertex_list, edge_list, facets_list):
    check_facet_type(facet)
    self.face_index = facet.index()

    [
        (self.a, self.a_index, self.face_opposite_a),
        (self.b, self.b_index, self.face_opposite_b),
        (self.c, self.c_index, self.face_opposite_c),
    ] = get_face_properties(facet, vertex_list, edge_list, facets_list)

    self.w1, self.w2 = compute_gram_schmidt_directions(facet,
                                                       self.b - self.a)

  def add_vertices(self, vertices_dict):
    overwrite_numpy_value(self.a_index, self.a, vertices_dict)
    overwrite_numpy_value(self.b_index, self.b, vertices_dict)
    overwrite_numpy_value(self.c_index, self.c, vertices_dict)

  def add_triangle(self, triangles):
    triangles.append((self.a_index, self.b_index, self.c_index))


def get_vertex_array(vertices_dict):
  # We re-number the vertices as well.
  renumber_vertex_dict = {key: i
                          for i, key in enumerate(vertices_dict.keys())}

  # Use mapping in reverse so we can create a (N x 3) array.
  renumbered_vertices = {}
  for key, value in vertices_dict.iteritems():
    new_key = renumber_vertex_dict[key]
    renumbered_vertices[new_key] = value

  num_vertices = len(vertices_dict)
  all_vertices = np.vstack([renumbered_vertices[i]
                            for i in xrange(num_vertices)])

  return renumber_vertex_dict, all_vertices


def get_face_data(vertex_list, edge_list, facets_list, initial_face_index):
  """Gets faces from dolfin.Mesh and turns them into Face objects.

  Only adds exterior faces.
  """
  vertices_dict = {}
  triangles = []
  face_local_bases = []
  neighbor_faces = []

  renumber_face_dict = {}
  new_face_index = 0
  for facet in facets_list:
    if not facet.exterior():
      continue
    raw_face_data = RawFaceData(facet, vertex_list, edge_list, facets_list)
    renumber_face_dict[raw_face_data.face_index] = new_face_index
    new_face_index += 1  # Increment to prepare for next one.

    raw_face_data.add_vertices(vertices_dict)
    raw_face_data.add_triangle(triangles)

    # Add the local orthogonal basis for the face as a 1D vector in R^6.
    face_local_bases.append(np.hstack([raw_face_data.w1, raw_face_data.w2]))
    # Add neighbor faces (before re-numbering).
    neighbor_faces.append([raw_face_data.face_opposite_a,
                           raw_face_data.face_opposite_b,
                           raw_face_data.face_opposite_c])

  # Re-number the vertices
  renumber_vertex_dict, all_vertices = get_vertex_array(vertices_dict)

  # Convert `triangles` to a numpy array and re-number based on vertex
  # re-numbering.
  triangles = np.array(triangles, dtype=np.int32)
  def renumber_triangles(vertex_index):
    return renumber_vertex_dict[vertex_index]
  renumber_triangles = np.vectorize(renumber_triangles)
  triangles = renumber_triangles(triangles)

  # Turn `face_local_bases` into a numpy array.
  face_local_bases = np.vstack(face_local_bases)

  # Convert `neighbor_faces` to a numpy array and re-number based on face
  # re-numbering.
  neighbor_faces = np.array(neighbor_faces, dtype=np.int32)
  def renumber_faces(face_index):
    return renumber_face_dict[face_index]
  renumber_faces = np.vectorize(renumber_faces)
  neighbor_faces = renumber_faces(neighbor_faces)

  # Update the index based on re-numbering.
  initial_face_index = renumber_face_dict[initial_face_index]

  return (all_vertices, triangles, face_local_bases,
          neighbor_faces, initial_face_index)


def from_mesh(cls, mesh, initial_point, k):
  print 'Creating Mesh from dolfin.Mesh data.'

  # Make sure it is the right kind of mesh.
  print 'Initializing mesh attributes (edges, faces, etc.)'
  mesh.init()  # This will do nothing if already called.
  check_mesh_type(mesh)

  # Compute extra data not stored on the object.
  print 'Reading vertex list from the mesh object'
  vertex_list = list(dolfin.vertices(mesh))
  print 'Reading edge list from the mesh object'
  edge_list = list(dolfin.edges(mesh))
  print 'Reading facets list from the mesh object'
  # Use facets since they have facet.exterior() set.
  facets_list = list(dolfin.facets(mesh))
  # Get values specific to motion on the mesh.
  print 'Reading cell list from the mesh object'
  cell_list = list(dolfin.cells(mesh))
  initial_face_index = get_face(dolfin.Point(*initial_point),
                                mesh, cell_list)

  print 'Parsing exterior faces and creating Face objects'
  (all_vertices, triangles, face_local_bases, neighbor_faces,
   initial_face_index) = get_face_data(vertex_list, edge_list,
                                       facets_list, initial_face_index)

  return cls(k, initial_point, initial_face_index,
             all_vertices, triangles, face_local_bases, neighbor_faces)
