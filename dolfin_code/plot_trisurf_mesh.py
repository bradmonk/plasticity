from mpl_toolkits.mplot3d import Axes3D  # Needs to be imported for 3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

from particle_diffusion_on_mesh import Mesh


resolution = 96
serialized_mesh_filename = 'serialized_mesh_res_%d.npz' % resolution
mesh_wrapper = Mesh.from_file(serialized_mesh_filename)

x = mesh_wrapper.all_vertices[:, 0]
y = mesh_wrapper.all_vertices[:, 1]
z = mesh_wrapper.all_vertices[:, 2]

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot_trisurf(x, y, z, triangles=mesh_wrapper.triangles,
                color=(0, 0, 0, 0), edgecolor='Gray', linewidth=0.05)

plt.show()

# http://stackoverflow.com/questions/7965743/
# python-matplotlib-setting-aspect-ratio
