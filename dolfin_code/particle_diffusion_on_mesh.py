import matplotlib
matplotlib.use('TKAgg')
import dolfin
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import random


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


def get_face_vertices(facet, vertex_list):
  # This will fail if not exactly 3 vertices.
  a_index, b_index, c_index = facet.entities(0)

  a = convert_point_to_array(vertex_list[a_index].point())
  b = convert_point_to_array(vertex_list[b_index].point())
  c = convert_point_to_array(vertex_list[c_index].point())

  return a, b, c, a_index, b_index, c_index


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
      (a, face_opposite_a),
      (b, face_opposite_b),
      (c, face_opposite_c),
  ]


def compute_gram_schmidt_directions(facet, directed_side):
  w1 = directed_side / np.linalg.norm(directed_side)
  n = convert_point_to_array(facet.normal())
  # n acts as z == e3, w1 acts as x == e1, hence y == e2 (== e3 x e1):
  w2 = np.cross(n, w1)
  return w1, w2


def check_facet_type(facet):
  if facet.dim() != 2:
    raise ValueError('Expected triangular facet.')


class FaceWrapper(object):

  def __init__(self, facet_index, a, b, c,
               face_opposite_a, face_opposite_b, face_opposite_c,
               w1, w2):
    self.facet_index = facet_index

    self.a = a
    self.b = b
    self.c = c

    self.face_opposite_a = face_opposite_a
    self.face_opposite_b = face_opposite_b
    self.face_opposite_c = face_opposite_c

    self.w1 = w1
    self.w2 = w2

    self.points = {}

  @classmethod
  def from_mesh_data(cls, facet, vertex_list, edge_list, facets_list):
    check_facet_type(facet)

    [
        (a, face_opposite_a),
        (b, face_opposite_b),
        (c, face_opposite_c),
    ] = get_face_properties(facet, vertex_list, edge_list, facets_list)

    w1, w2 = compute_gram_schmidt_directions(facet, b - a)
    return cls(facet.index(), a, b, c,
               face_opposite_a, face_opposite_b, face_opposite_c,
               w1, w2)

  def __str__(self):
    return 'Face(%d)' % self.facet_index

  def __repr__(self):
    return str(self)

  def add_point(self, point):
    self.points[point.point_index] = point

  def remove_point(self, point):
    self.points.pop(point.point_index)

  def add_mesh_wrapper(self, mesh_wrapper):
    self.parent_mesh_wrapper = mesh_wrapper

  def flatten_data(self):
    return np.hstack([self.a, self.b, self.c, self.w1, self.w2,
                      self.facet_index, self.face_opposite_a,
                      self.face_opposite_b, self.face_opposite_c])

  @classmethod
  def from_flattened_data(cls, flattened_data):
    a = flattened_data[:3]
    b = flattened_data[3:6]
    c = flattened_data[6:9]

    w1 = flattened_data[9:12]
    w2 = flattened_data[12:15]

    facet_index = int(flattened_data[15])

    face_opposite_a = int(flattened_data[16])
    face_opposite_b = int(flattened_data[17])
    face_opposite_c = int(flattened_data[18])

    return cls(facet_index, a, b, c,
               face_opposite_a, face_opposite_b, face_opposite_c,
               w1, w2)

  def compute_angle(self, direction):
    """Computes angle of `direction` in current orthogonal coordinates.

    Assumes `direction` lies in plane spanned by w1 and w2. If this
    is not true, we'll still get a number, but it won't make sense.

    NOTE: We could check the validity of this by computing
              ||direction|| cos(theta) w1 + ||direction|| sin(theta) w2
          and comparing it to `direction`.
    """
    scaled_cosine = self.w1.dot(direction)  # ||direction|| cos(theta)
    scaled_sine = self.w2.dot(direction)  # ||direction|| sin(theta)
    return np.arctan2(scaled_sine, scaled_cosine)

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

  def compute_new_direction_theta(self, fixed_direction,
                                  direction_to_change, new_face_index):
    """Computes the angle of direction of motion in new face.

    Assumes `fixed_direction` is the vector along an edge of the triangle
    and `direction_to_change` should maintain the same angle between that
    vector on the current and next face. Also assumes `new_face_index`
    describes the other face with `fixed_direction` as an edge.
    """
    angle_change = (self.compute_angle(direction_to_change) -
                    self.compute_angle(fixed_direction))

    new_face = self.parent_mesh_wrapper.faces[new_face_index]
    return angle_change + new_face.compute_angle(fixed_direction)

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
    next_face = self.facet_index
    # We don't need these attributes if the move stays on the same face.
    remaining_length = theta_new = None
    # If the move takes us past the intersection, we need to change to
    # a different face and can't use the full length of the move.
    if move_length > t:
      actual_move_length = t
      remaining_length = move_length - t
      next_face = next_face_index
      theta_new = self.compute_new_direction_theta(direction_side_line,
                                                   particle_direction,
                                                   next_face)

    next_point = particle_center + actual_move_length * particle_direction
    list_of_moves.append((next_face, next_point,
                          remaining_length, theta_new))

  def choose_move(self, point, move_choices):
    num_choices = len(move_choices)
    if num_choices == 1:
      return move_choices[0]

    # We expect to travel in the direction of either 1 or 2 sides (if
    # we are moving in the direction of a vertex).
    if num_choices != 2:
      raise ValueError('Unexpected number of possible moves on face.')

    first_choice, second_choice = move_choices
    no_move_first = np.allclose(first_choice[1], point)
    no_move_second = np.allclose(second_choice[1], point)

    if no_move_first and not no_move_second:
      return second_choice
    elif not no_move_first and no_move_second:
      return first_choice
    else:
      # NOTE: In the case of multiple moves, we expect that `next_point` and
      #       `L_new` should be identical for all tuples in `move_choices`,
      #       but don't check this here.
      # BEGIN: Temporary statements to showcase issues with `random.choice`.
      print 'Used random.choice'
      # END: Temporary statements to showcase issues with `random.choice`.
      # The goal in the case of a tie is to correctly determine which face
      # to go into *through* the vertex. We can do this by taking each
      # triangle and determining the angle contributed at the vertex.
      return random.choice(move_choices)

  def move(self, point, L, theta):
    particle_direction = np.cos(theta) * self.w1 + np.sin(theta) * self.w2
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

    return self.choose_move(point, move_choices)


def check_mesh_type(mesh):
  if mesh.geometry().dim() != 3:
    raise ValueError('Expecting 3D mesh.')
  if mesh.cells().shape[1] != 4:
    raise ValueError('Expecting tetrahedral mesh.')


def get_face(point, mesh, cell_list, face_dictionary):
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
    raise ValueError('Multiple facets on cell marked exterior.')

  exterior_facet_index = exterior_facets[0].index()
  return face_dictionary[exterior_facet_index]


def get_faces(vertex_list, edge_list, facets_list):
  """Gets faces from mesh and turns them into FaceWrapper objects.

  Only adds exterior faces.
  """
  faces = {}
  # Use facets since they have facet.exterior() set.
  for facet in facets_list:
    if not facet.exterior():
      continue
    faces[facet.index()] = FaceWrapper.from_mesh_data(facet, vertex_list,
                                                      edge_list, facets_list)

  return faces


class MeshWrapper(object):

  def __init__(self, k, initial_point, initial_face, faces):
    self.k = k
    self.initial_point = initial_point
    self.initial_face = initial_face

    self.faces = faces
    for face in faces.values():
      face.add_mesh_wrapper(self)

  @classmethod
  def from_mesh(cls, mesh, initial_point, k, return_facets=False):
    # Make sure it is the right kind of mesh.
    mesh.init()  # This will do nothing if already called.
    check_mesh_type(mesh)

    # Compute extra data not stored on the object.
    vertex_list = list(dolfin.vertices(mesh))
    edge_list = list(dolfin.edges(mesh))
    facets_list = list(dolfin.facets(mesh))
    faces = get_faces(vertex_list, edge_list, facets_list)

    # Set values specific to motion on the mesh.
    cell_list = list(dolfin.cells(mesh))
    initial_face = get_face(dolfin.Point(*initial_point), mesh,
                            cell_list, faces)

    if return_facets:
      return cls(k, initial_point, initial_face, faces), facets_list
    else:
      return cls(k, initial_point, initial_face, faces)

  def __str__(self):
    return 'Mesh(num_faces=%d)' % len(self.faces)

  def __repr__(self):
    return str(self)


def get_random_components(k):
  x_rand = k * np.random.randn()
  y_rand = k * np.random.randn()
  L = np.linalg.norm([x_rand, y_rand])
  theta = np.arctan2(y_rand, x_rand)
  return L, theta


class Point(object):

  def __init__(self, point_index, mesh_wrapper):
    self.point_index = point_index
    # NOTE: We don't need to store this since `self.face.parent_mesh_wrapper`
    #       will also hold this value.
    self.mesh_wrapper = mesh_wrapper

    # Start with Null `Face` object.
    self.face = None

    self.change_face(mesh_wrapper.initial_face)
    self.point = np.array(mesh_wrapper.initial_point)

    self.k = mesh_wrapper.k

    self.move_counter = 0
    self.values = [(self.face.facet_index, self.point)]

  def __str__(self):
    return 'Point(%d, face=%d)' % (self.point_index,
                                   self.face.facet_index)

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

  def _move(self, L, theta):
    next_face, next_point, L_new, theta_new = self.face.move(
        self.point, L, theta)
    self.point = next_point
    if next_face != self.face.facet_index:
      self.change_face(self.mesh_wrapper.faces[next_face])
      # Continue to move until we stay on the same face.
      self._move(L_new, theta_new)

  def move(self):
    L, theta = get_random_components(self.k)
    self._move(L, theta)

    self.move_counter += 1
    self.values.append((self.face.facet_index, self.point))


# Mostly borrowed from particle_diffusion.py (without using dolfin).
def plot_simulation(num_points, mesh_wrapper,
                    num_frames=200, print_frequency=None,
                    interval=30, filename=None):
  points = [Point(i, mesh_wrapper) for i in xrange(num_points)]

  # Attaching 3D axis to the figure
  fig = plt.figure()
  ax = p3.Axes3D(fig)

  # Create lines with a single point as a scatter.
  all_points = [ax.plot([pt.point[0]], [pt.point[1]], [pt.point[2]],
                        c='b', marker='o')[0]
                for pt in points]

  def update_plot(step_num):
    if print_frequency is not None and step_num % print_frequency == 0:
      print 'Step Number:', step_num

    for pt, point_container in zip(points, all_points):
      pt.move()
      # point_container.set_color(color)
      point_container.set_data([pt.point[0]], [pt.point[1]])
      point_container.set_3d_properties([pt.point[2]])

    return all_points

  # Setting the axes properties
  ax.set_xlim3d([-17.0, 17.0])  # Just slightly bigger than 15
  ax.set_xlabel('X')

  ax.set_ylim3d([-17.0, 17.0])
  ax.set_ylabel('Y')

  ax.set_zlim3d([0.0, 50.0])
  ax.set_zlabel('Z')

  anim = animation.FuncAnimation(fig, update_plot,
                                 repeat=False, frames=num_frames,
                                 interval=interval, blit=False)
                                 # interval=interval, blit=True)

  # Could also use plt.draw() for interactive work, not for 3D though.
  plt.show(block=False)


def sample_code():
  resolution = 96
  mesh_full_filename = 'mesh_res_%d_full.xml' % resolution
  mesh_3d = dolfin.Mesh(mesh_full_filename)

  # NOTE: This is temporary. These are parameters of the mesh (when it was
  #       created in full_dendrite_mesh.py) and we should package them in a
  #       different way.
  SCALE_FACTOR = 50.0
  STARTING_X = SCALE_FACTOR * 0.0
  STARTING_Y = SCALE_FACTOR * 0.0
  STARTING_Z = SCALE_FACTOR * 1.0
  STARTING_K = SCALE_FACTOR * 0.01

  initial_point = np.array((STARTING_X, STARTING_Y, STARTING_Z))
  mesh_wrapper = MeshWrapper.from_mesh(mesh_3d, initial_point, STARTING_K)

  points = [Point(i, mesh_wrapper) for i in xrange(10)]

  for i in xrange(5):
    points[0].move()

  return points


def error_off_plane(face_id, point, facets_list):
  facet = facets_list[face_id]
  n = convert_point_to_array(facet.normal())
  midpoint = convert_point_to_array(facet.midpoint())
  return np.dot(point - midpoint, n)


def test_accurary_on_face(mesh_wrapper=None, facets_list=None, num_steps=1000):
  SCALE_FACTOR = 50.0
  STARTING_X = SCALE_FACTOR * 0.0
  STARTING_Y = SCALE_FACTOR * 0.0
  STARTING_Z = SCALE_FACTOR * 1.0
  STARTING_K = SCALE_FACTOR * 0.01

  if mesh_wrapper is None or facets_list is None:
    resolution = 96
    mesh_full_filename = 'mesh_res_%d_full.xml' % resolution
    mesh_3d = dolfin.Mesh(mesh_full_filename)
    initial_point = np.array((STARTING_X, STARTING_Y, STARTING_Z))
    mesh_wrapper, facets_list = MeshWrapper.from_mesh(
        mesh_3d, initial_point, STARTING_K, return_facets=True)

  point = Point(0, mesh_wrapper)

  for i in xrange(num_steps):
    point.move()

  errors = [error_off_plane(face_id, pt, facets_list)
            for face_id, pt in point.values]
  print 'Max Error after %d steps' % num_steps
  print np.max(np.abs(errors))
