from mpl_toolkits.mplot3d import Axes3D  # Needs to be imported for 3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

from particle_diffusion_on_mesh import Mesh


resolution = 96
serialized_mesh_filename = 'serialized_mesh_res_%d.npz' % resolution
mesh_wrapper = Mesh.from_file(serialized_mesh_filename)

points_dict = {}
triangles_list = []
for face in mesh_wrapper.faces.itervalues():
  # Don't actually check if they already exist.
  points_dict[face.a_index] = face.a
  points_dict[face.b_index] = face.b
  points_dict[face.c_index] = face.c

  triangles_list.append((face.a_index, face.b_index, face.c_index))

# Renumber the points as 0,1,2,...
renumber_dict = {key: i for i, key in enumerate(points_dict.keys())}

old_points_dict = points_dict
points_dict = {}
for key, value in old_points_dict.iteritems():
  new_key = renumber_dict[key]
  points_dict[new_key] = value

old_triangles_list = triangles_list
triangles_list = []
for a_index, b_index, c_index in old_triangles_list:
  new_a_index = renumber_dict[a_index]
  new_b_index = renumber_dict[b_index]
  new_c_index = renumber_dict[c_index]

  triangles_list.append((new_a_index, new_b_index, new_c_index))

num_points = len(points_dict)
if sorted(points_dict.keys()) != range(num_points):
  raise ValueError('Renumbering did not work.')

stacked_points = np.vstack([points_dict[i] for i in xrange(num_points)])
x = stacked_points[:, 0]
y = stacked_points[:, 1]
z = stacked_points[:, 2]

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot_trisurf(x, y, z, triangles=triangles_list, color='w',
                linewidth=0.05)

plt.show()

# http://stackoverflow.com/questions/7965743/
# python-matplotlib-setting-aspect-ratio
