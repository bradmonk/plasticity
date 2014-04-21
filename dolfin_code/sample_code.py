import numpy as np

from particle_diffusion_on_mesh import Mesh
from particle_diffusion_on_mesh import PlotBoundary
from particle_diffusion_on_mesh import Point
from particle_diffusion_on_mesh import plot_simulation


def save_serialized_mesh():
  print 'Importing dolfin, takes a bit of time...'
  import dolfin
  print 'Done importing dolfin.'

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
  mesh_wrapper = Mesh.from_mesh(mesh_3d, initial_point, STARTING_K)

  serialized_mesh_filename = 'serialized_mesh_res_%d.npz' % resolution
  mesh_wrapper.serialize_mesh(serialized_mesh_filename)


def sample_code():
  resolution = 96
  serialized_mesh_filename = 'serialized_mesh_res_%d.npz' % resolution
  mesh_wrapper = Mesh.from_file(serialized_mesh_filename)

  points = [Point(mesh_wrapper) for _ in xrange(10)]

  for i in xrange(5):
    points[0].move()

  return points


def error_off_plane(face, point):
  # Built so that cross(x, y) == z
  n = np.cross(face.w1, face.w2)
  return np.dot(point - face.a, n)


def test_accurary_on_face(mesh_wrapper=None, num_steps=1000):
  if mesh_wrapper is None:
    resolution = 96

    serialized_mesh_filename = 'serialized_mesh_res_%d.npz' % resolution
    mesh_wrapper = Mesh.from_file(serialized_mesh_filename)

  point = Point(mesh_wrapper)

  for i in xrange(num_steps):
    point.move()

  errors = [error_off_plane(mesh_wrapper.faces[face_id], pt)
            for face_id, pt in point.values]
  print 'Max Error after %d steps' % num_steps
  print np.max(np.abs(errors))


def plot_custom():
  X_CENTER1 = 0.0
  X_CENTER2 = 100.0
  X_CENTER3 = 200.0
  X_CENTER4 = 300.0
  MAX_RADIUS = 15.0
  def in_box(value, center):
    delta = value - center
    if np.allclose(delta, MAX_RADIUS) or np.allclose(delta, - MAX_RADIUS):
      return True

    return -MAX_RADIUS <= delta <= MAX_RADIUS

  def in_x_y_box(point, x_center):
    return in_box(point.x, x_center) and in_box(point.y, 0.0)

  def custom_color_function(point):
    if np.allclose(point.z, 50.0):
      return 'b'

    if np.allclose(point.z, 0) or point.z >= 0:
      if in_x_y_box(point, X_CENTER1) or in_x_y_box(point, X_CENTER3):
        return 'r'
      elif in_x_y_box(point, X_CENTER2) or in_x_y_box(point, X_CENTER4):
        return 'y'

    return 'g'

  resolution = 96
  serialized_mesh_filename = 'serialized_mesh_res_%d.npz' % resolution
  mesh_wrapper = Mesh.from_file(serialized_mesh_filename)
  mesh_wrapper.k = 6.0

  x = mesh_wrapper.all_vertices[:, 0]
  y = mesh_wrapper.all_vertices[:, 1]
  z = mesh_wrapper.all_vertices[:, 2]
  # Consider putting this into `plot_simulation`.
  plot_boundary = PlotBoundary(np.min(x), np.max(x),
                               np.min(y), np.max(y),
                               np.min(z), np.max(z))

  plot_simulation(10, mesh_wrapper, plot_boundary,
                  color_function=custom_color_function,
                  num_frames=500, print_frequency=20, show_mesh=True,
                  filename='10points_500steps_bigger_k.gif')
