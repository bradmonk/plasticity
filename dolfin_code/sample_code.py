import dolfin
import numpy as np

from particle_diffusion_on_mesh import Mesh
from particle_diffusion_on_mesh import PlotBoundary
from particle_diffusion_on_mesh import Point
from particle_diffusion_on_mesh import convert_point_to_array
from particle_diffusion_on_mesh import plot_simulation


def save_serialized_mesh():
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
  n = np.cross(face.c - face.a, face.b - face.a)
  n = n / np.linalg.norm(n)
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
  def custom_color_function(point):
    if np.allclose(point.z, 50.0):
      return 'b'
    else:
      return 'r'

  resolution = 96
  serialized_mesh_filename = 'serialized_mesh_res_%d.npz' % resolution
  mesh_wrapper = Mesh.from_file(serialized_mesh_filename)

  plot_boundary = PlotBoundary(-17.0, 17.0, -17.0, 17.0, 0.0, 50.0)

  plot_simulation(200, mesh_wrapper, plot_boundary,
                  color_function=custom_color_function,
                  print_frequency=20, show_mesh=True,
                  filename='200points_200steps_colored_points.gif')
