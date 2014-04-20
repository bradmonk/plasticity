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


def check_facet_type(facet):
  if facet.dim() != 2:
    raise ValueError('Expected triangular facet.')


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


class Face(object):

  def __init__(self, face_index, mesh_wrapper):
    self.face_index = face_index
    self.mesh_wrapper = mesh_wrapper

    self.points = {}

  @property
  def a(self):
    vertex_index = self.mesh_wrapper.triangles[self.face_index, 0]
    return self.mesh_wrapper.all_vertices[vertex_index, :]

  @property
  def b(self):
    vertex_index = self.mesh_wrapper.triangles[self.face_index, 1]
    return self.mesh_wrapper.all_vertices[vertex_index, :]

  @property
  def c(self):
    vertex_index = self.mesh_wrapper.triangles[self.face_index, 2]
    return self.mesh_wrapper.all_vertices[vertex_index, :]

  @property
  def face_opposite_a(self):
    return self.mesh_wrapper.neighbor_faces[self.face_index, 0]

  @property
  def face_opposite_b(self):
    return self.mesh_wrapper.neighbor_faces[self.face_index, 1]

  @property
  def face_opposite_c(self):
    return self.mesh_wrapper.neighbor_faces[self.face_index, 2]

  @property
  def w1(self):
    return self.mesh_wrapper.face_local_bases[self.face_index, :3]

  @property
  def w2(self):
    return self.mesh_wrapper.face_local_bases[self.face_index, 3:]

  def __str__(self):
    return 'Face(%d)' % self.face_index

  def __repr__(self):
    return str(self)

  def add_point(self, point):
    self.points[point.point_index] = point

  def remove_point(self, point):
    self.points.pop(point.point_index)

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

    new_face = self.mesh_wrapper.faces[new_face_index]
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
    next_face = self.face_index
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
    raise ValueError('Multiple facets on cell marked exterior.')

  exterior_facet_index = exterior_facets[0].index()
  return exterior_facet_index


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
  """Gets faces from mesh and turns them into Face objects.

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


class Mesh(object):

  def __init__(self, k, initial_point, initial_face_index,
               all_vertices, triangles, face_local_bases, neighbor_faces):
    self.k = k
    self.initial_point = initial_point
    self.initial_face_index = initial_face_index

    self.all_vertices = all_vertices
    self.triangles = triangles
    self.face_local_bases = face_local_bases
    self.neighbor_faces = neighbor_faces

    self.faces = {}
    for face_index in xrange(triangles.shape[0]):
      self.faces[face_index] = Face(face_index, self)

    # Set initial_face object based on newly created faces.
    self.initial_face = self.faces[self.initial_face_index]
    # Counter to keep track of points on the mesh.
    self.current_point_index = -1

  @classmethod
  def from_mesh(cls, mesh, initial_point, k):
    print 'Creating Mesh from mesh data.'

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

  def next_index(self):
    self.current_point_index += 1
    return self.current_point_index

  def serialize_mesh(self, filename):
    print 'Saving mesh to', filename
    if self.current_point_index != -1:
      print 'Points on mesh will not be serialized.'

    np.savez(filename, k=self.k, initial_point=self.initial_point,
             initial_face_index=self.initial_face_index,
             all_vertices=self.all_vertices, triangles=self.triangles,
             face_local_bases=self.face_local_bases,
             neighbor_faces=self.neighbor_faces)

  @classmethod
  def from_file(cls, filename):
    print 'Loading mesh data from NPZ file', filename
    npzfile = np.load(filename)

    k = npzfile['k'].item()
    initial_point = npzfile['initial_point']
    initial_face_index = npzfile['initial_face_index'].item()

    all_vertices = npzfile['all_vertices']
    triangles = npzfile['triangles']
    face_local_bases = npzfile['face_local_bases']
    neighbor_faces = npzfile['neighbor_faces']

    return cls(k, initial_point, initial_face_index,
               all_vertices, triangles, face_local_bases, neighbor_faces)

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

  def __init__(self, mesh_wrapper):
    self.point_index = mesh_wrapper.next_index()
    # NOTE: We don't need to store this since `self.face.mesh_wrapper`
    #       will also hold this value.
    self.mesh_wrapper = mesh_wrapper

    # Start with Null `Face` object.
    self.face = None

    self.change_face(mesh_wrapper.initial_face)
    self.point = np.array(mesh_wrapper.initial_point)

    self.k = mesh_wrapper.k

    self.move_counter = 0
    self.values = [(self.face.face_index, self.point)]

  def __str__(self):
    return 'Point(%d, face=%d)' % (self.point_index,
                                   self.face.face_index)

  def __repr__(self):
    return str(self)

  @property
  def x(self):
    return self.point[0]

  @property
  def y(self):
    return self.point[1]

  @property
  def z(self):
    return self.point[2]

  def change_face(self, face):
    """Updates the face on current object and adds point to face.

    Args:
      face: A Face object.
    """
    if self.face is not None:
      self.face.remove_point(self)

    self.face = face
    self.face.add_point(self)

  def _move(self, L, theta):
    next_face, next_point, L_new, theta_new = self.face.move(
        self.point, L, theta)
    self.point = next_point
    if next_face != self.face.face_index:
      self.change_face(self.mesh_wrapper.faces[next_face])
      # Continue to move until we stay on the same face.
      self._move(L_new, theta_new)

  def move(self):
    L, theta = get_random_components(self.k)
    self._move(L, theta)

    self.move_counter += 1
    self.values.append((self.face.face_index, self.point))


class PlotBoundary(object):
  """Container object to define a plot area for an animation.

  Example:
    >>> # 17.0 is just slightly bigger than 15
    >>> plot_boundary = PlotBoundary(-17.0, 17.0, -17.0, 17.0, 0.0, 50.0)
  """

  def __init__(self, min_x, max_x, min_y, max_y, min_z, max_z):
    self.min_x = min_x
    self.max_x = max_x
    self.min_y = min_y
    self.max_y = max_y
    self.min_z = min_z
    self.max_z = max_z


def default_color_function(point):
  return 'b'


# Mostly borrowed from particle_diffusion.py (without using dolfin).
def plot_simulation(num_points, mesh_wrapper, plot_boundary,
                    color_function=default_color_function, show_mesh=False,
                    num_frames=200, print_frequency=None,
                    interval=30, filename=None):
  points = [Point(mesh_wrapper) for _ in xrange(num_points)]

  # Attaching 3D axis to the figure
  fig = plt.figure()
  ax = p3.Axes3D(fig)
  # Set default view.
  ax.view_init(elev=12, azim=-37)

  # Create lines with a single point as a scatter.
  all_points = [ax.plot([pt.x], [pt.y], [pt.z], c='b', marker='o')[0]
                for pt in points]
  if show_mesh:
    # Add permanent (unchanged) features to plot.
    x = mesh_wrapper.all_vertices[:, 0]
    y = mesh_wrapper.all_vertices[:, 1]
    z = mesh_wrapper.all_vertices[:, 2]
    ax.plot_trisurf(x, y, z, triangles=mesh_wrapper.triangles,
                    color='w', linewidth=0.05)

  def update_plot(step_num):
    if print_frequency is not None and step_num % print_frequency == 0:
      print 'Step Number:', step_num

    for pt, point_container in zip(points, all_points):
      pt.move()
      point_container.set_color(color_function(pt))
      point_container.set_data([pt.x], [pt.y])
      point_container.set_3d_properties([pt.z])

    return all_points

  # Setting the axes properties
  ax.set_xlim3d([plot_boundary.min_x, plot_boundary.max_x])
  ax.set_xlabel('X')

  ax.set_ylim3d([plot_boundary.min_y, plot_boundary.max_y])
  ax.set_ylabel('Y')

  ax.set_zlim3d([plot_boundary.min_z, plot_boundary.max_z])
  ax.set_zlabel('Z')

  create_gif = filename is not None
  anim = animation.FuncAnimation(fig, update_plot,
                                 repeat=False, frames=num_frames,
                                 interval=interval, blit=create_gif)

  if create_gif:
    anim.save(filename, writer='imagemagick_file')
  else:
    # Could also use plt.draw() for interactive work, not for 3D though.
    plt.show(block=False)
