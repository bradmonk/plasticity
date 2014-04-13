import dolfin
import numpy as np
import random


def convert_point_to_array(point_object):
  return np.array([point_object.x(), point_object.y(), point_object.z()])


def get_neighbor(face, edge, faces_list):
  """Gets neighbor for a face across a given edge."""
  neighbor_faces = []
  for index in edge.entities(2):
    # Skip current face.
    if index == face.index():
      continue
    neighbor_faces.append(index)

  if len(neighbor_faces) != 1:
    raise ValueError('Wrong number of edge neighbors.')

  return neighbor_faces[0]


def find_intersection(center0, direction0, center1, direction1):
  """Finds the intersection of two lines in R^3.

  Solves
      c0 + d0 t = c1 + d1 s
  and returns (t, s). If the lines do not cross then returns (None, None).

  Args:
    center0: 1D NumPy are with 3 elements.
    direction0: 1D NumPy are with 3 elements.
    center1: 1D NumPy are with 3 elements.
    direction1: 1D NumPy are with 3 elements.

  Returns:
    Two real values s and t.
  """
  # c0 + d0 t = c1 + d1 s
  # (-d0) t + (d1) s = c0 - c1
  # [-d0, d1] [t,s]^T = delta
  A = np.array([-direction0, direction1]).T
  delta = center0 - center1
  (t, s), residue, rank, _ = np.linalg.lstsq(A, delta)
  if not (rank == 2 and np.allclose(residue, 0)):
    raise ValueError('Solution not unique')

  return t, s


class FaceWrapper(object):

  def __init__(self, face, parent_mesh_wrapper):
    self.face = face
    self.check_face_type()

    self.parent_mesh_wrapper = parent_mesh_wrapper

    self.set_vertices()
    self.set_neighbors()
    self.compute_gram_schmidt_directions()

    self.points = {}

  def __str__(self):
    return 'Face(%d)' % self.face.index()

  def __repr__(self):
    return str(self)

  def check_face_type(self):
    if self.face.dim() != 2:
      raise ValueError('Expected triangular face.')

  def add_point(self, point):
    self.points[point.point_index] = point

  def remove_point(self, point):
    self.points.pop(point.point_index)

  def set_vertices(self):
    # This will fail if not exactly 3 vertices.
    self.a_index, self.b_index, self.c_index = self.face.entities(0)

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
    """Sets neighbor face indices.

    NOTE: This behavior is subject to change.

    Determines the face across an edge from each vertex a, b, c and
    labels accordingly.

    This ignores the neighbors which go through a vertex.
    """
    face_edge_indices = self.face.entities(1)

    edge_vertex_pairing = []
    # NOTE: We don't check that there are 3 edges, but could / should.
    for edge_index in face_edge_indices:
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
    faces_list = self.parent_mesh_wrapper.faces_list
    self.face_opposite_a = get_neighbor(self.face, sorted_edges[0],
                                        faces_list)
    self.face_opposite_b = get_neighbor(self.face, sorted_edges[1],
                                        faces_list)
    self.face_opposite_c = get_neighbor(self.face, sorted_edges[2],
                                        faces_list)

  def compute_gram_schmidt_directions(self):
    v1 = self.b - self.a
    v2 = self.c - self.a
    n = convert_point_to_array(self.face.cell_normal())

    vec_as_cols = np.array([v1, v2, n]).T
    Q, _ = np.linalg.qr(vec_as_cols)
    # NOTE: We don't do this, but we could / should check that
    #       w3 = Q[2, :] is either n or -n (same direction, unit length).
    #       For example, we could do this by checking:
    #           np.allclose(w3.dot(n), 1) or np.allclose(w3.dot(n), -1)
    self.w1 = Q[:, 0]
    self.w2 = Q[:, 1]

  @staticmethod
  def check_intersection(t, s):
    """Check that intersection is on the right line.

    t and s come from solving
        p + t * d = a + (b - a) * s
    where p is the point we are moving from, d is the normal direction
    we are moving and a and b are vertices of a side of a triangle.

    Check that (t >= 0), i.e. we are actually moving in the direction
    given by `d`. Also checks that (0 <= s <= 1), i.e. that the intersection
    is between `a` and `b`.
    """
    # Allow wiggle room nearby 0.
    if not (np.allclose(t, 0) or 0 <= t):
      return False

    # Allow wiggle room nearby 0.
    if not (np.allclose(s, 0) or 0 <= s):
      return False

    # Allow wiggle room nearby 1.
    if not (np.allclose(s, 1) or s <= 1):
      return False

    return True

  def move_toward_side(self, list_of_moves,
                       particle_center, particle_direction,
                       side_first_vertex, side_second_vertex,
                       move_length, next_face_index):
    center_side_line = side_first_vertex
    direction_side_line = side_second_vertex - side_first_vertex

    t, s = find_intersection(particle_center, particle_direction,
                             center_side_line, direction_side_line)

    if not self.check_intersection(t, s):
      # Do nothing; i.e. don't augment `list_of_moves`.
      return

    # If we haven't returned, this side is a valid choice.
    actual_move_length = move_length
    remaining_length = 0
    next_face = self.face.index()
    # If the move takes us past the intersection, we need to change to
    # a different face and can't use the full length of the move.
    if move_length > t:
      actual_move_length = t
      remaining_length = move_length - t
      next_face = next_face_index

    next_point = particle_center + actual_move_length * particle_direction

    direction_new = None  # To be computed.
    list_of_moves.append((next_face, next_point,
                          remaining_length, direction_new))

  def move(self, point, L, particle_direction):
    # In case the point lies on a vertex, we may have more than
    # one possible move.
    move_choices = []

    # Side opposite A is B-->C.
    self.move_toward_side(move_choices, point, particle_direction,
                          self.b, self.c, L, self.face_opposite_a)
    # Side opposite B is C-->A.
    self.move_toward_side(move_choices, point, particle_direction,
                          self.c, self.a, L, self.face_opposite_b)
    # Side opposite C is A-->B.
    self.move_toward_side(move_choices, point, particle_direction,
                          self.a, self.b, L, self.face_opposite_c)

    # We expect to travel in the direction of either 1 or 2 sides (if
    # we are moving in the direction of a vertex).
    if len(move_choices) not in (1, 2):
      raise ValueError('Unexpected number of possible moves on face.')
    # The goal in the case of a tie is to correctly determine which face
    # to go into *through* the vertex. We can do this by taking each triangle
    # and determining the angle contributed at the vertex.

    # NOTE: In the case of multiple moves, we expect that `next_point` and
    #       `L_new` should be identical, but don't check this below.
    next_face, next_point, L_new, direction_new = random.choice(move_choices)
    return next_face, next_point, L_new, direction_new


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
    # In a triangular 3D mesh, the cells are faces.
    self.faces_list = list(dolfin.cells(self.mesh))
    self.add_faces()

  def __str__(self):
    return 'Mesh(num_faces=%d)' % len(self.faces)

  def __repr__(self):
    return str(self)

  def check_mesh_type(self):
    if self.mesh.geometry().dim() != 3:
      raise ValueError('Expecting 3D mesh.')
    if self.mesh.cells().shape[1] != 3:
      raise ValueError('Expecting triangular / boundary mesh.')

  def add_faces(self):
    """Adds faces from mesh to stored list on object.

    Only adds exterior faces.
    """
    if self._faces_added:
      return

    self.faces = []
    for face in self.faces_list:
      self.faces.append(FaceWrapper(face, self))

    self._faces_added = True


def get_face(mesh_wrapper, point):
  face_intersection = dolfin.cpp.mesh.intersect(mesh_wrapper.mesh, point)
  # Require a unique intersection. This is not actually necessary as some
  # points may lie on an edge or a vertex or may not be on the mesh at all.
  if face_intersection.intersected_cells().size != 1:
    raise ValueError('Point does not intersect the mesh on a unique face.')
  cell_index = face_intersection.intersected_cells()[0]
  return mesh_wrapper.faces[cell_index]


def get_random_components(k):
  x_rand = k * np.random.randn()
  y_rand = k * np.random.randn()
  L = np.linalg.norm([x_rand, y_rand])
  theta = np.arctan2(y_rand, x_rand)
  return L, theta


class Point(object):

  def __init__(self, point_index, x, y, z, k, mesh_wrapper):
    self.point_index = point_index
    # NOTE: We don't need to store this since `self.face.parent_mesh_wrapper`
    #       will also hold this value.
    self.mesh_wrapper = mesh_wrapper

    # Start with Null face.
    self.face = None

    dolfin_point = dolfin.Point(x, y, z)
    self.change_face(get_face(mesh_wrapper, dolfin_point))
    self.point = np.array([x, y, z])

    self.k = k

  def __str__(self):
    return 'Point(%d, face=%d)' % (self.point_index,
                                   self.face.face.index())

  def __repr__(self):
    return str(self)

  def change_face(self, face):
    """Updates the face on current object and adds point to face.

    Args:
      face: A FaceWrapper object.
    """
    if self.face is not None:
      self.face.remove_point(self)

    self.face = face
    self.face.add_point(self)

  def move(self):
    L, theta = get_random_components(self.k)
    direction = np.cos(theta) * self.face.w1 + np.sin(theta) * self.face.w2
    next_face, next_point, L_new, direction_new = self.face.move(
        self.point, L, direction)
    # Currently we ignore `L_new`, simply change the face and move on.
    # We intend to use `L_new` to move in the other face, but currently
    # hold off since it requires computing a new theta.
    self.point = next_point
    if next_face != self.face.face.index():
      self.change_face(self.mesh_wrapper.faces[next_face])


def sample_code():
  resolution = 96
  mesh_full_filename = 'mesh_res_%d_boundary.xml' % resolution
  mesh_3d = dolfin.Mesh(mesh_full_filename)
  mesh_wrapper = MeshWrapper(mesh_3d)

  # NOTE: This is temporary. These are parameters of the mesh (when it was
  #       created in full_dendrite_mesh.py) and we should package them in a
  #       different way.
  SCALE_FACTOR = 50.0
  STARTING_X = SCALE_FACTOR * 0.0
  STARTING_Y = SCALE_FACTOR * 0.0
  STARTING_Z = SCALE_FACTOR * 1.0
  STARTING_K = SCALE_FACTOR * 0.01

  points = [Point(i, STARTING_X, STARTING_Y, STARTING_Z,
                  STARTING_K, mesh_wrapper)
            for i in xrange(10)]
  return points
