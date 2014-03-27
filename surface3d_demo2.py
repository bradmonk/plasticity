# From
#     http://matplotlib.org/mpl_examples/mplot3d/surface3d_demo3.py
# linked to on
#    http://matplotlib.org/mpl_toolkits/mplot3d/tutorial.html
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

RADIUS = 4
CENTER_X = 0.0
CENTER_Y = 0.0
CENTER_Z = 4.0
Z_CUTOFF = 1.0


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

theta = np.linspace(0, 2 * np.pi, 100)
# Need Z_CUTOFF = RADIUS cos(phi) + CENTER_Z
phi_max = np.arccos((Z_CUTOFF - CENTER_Z) / RADIUS)
phi = np.linspace(0, phi_max, 100)

x = RADIUS * np.outer(np.cos(theta), np.sin(phi)) + CENTER_X
y = RADIUS * np.outer(np.sin(theta), np.sin(phi)) + CENTER_Y
z = RADIUS * np.outer(np.ones(np.size(theta)), np.cos(phi)) + CENTER_Z

ax.plot_surface(x, y, z, rstride=4, cstride=4, color='y')

# Add a cylinder
Z_BOTTOM_CYL = -1.0
RADIUS_CYL = np.sqrt(RADIUS**2 - (Z_CUTOFF - CENTER_Z)**2)
theta_cyl = np.linspace(0, 2 * np.pi, 100)
z_mesh_cyl = np.linspace(Z_BOTTOM_CYL, Z_CUTOFF, 100)
x_cyl = RADIUS_CYL * np.outer(np.cos(theta_cyl),
                              np.ones(np.size(z_mesh_cyl)))
y_cyl = RADIUS_CYL * np.outer(np.sin(theta_cyl),
                              np.ones(np.size(z_mesh_cyl)))
z_cyl = RADIUS_CYL * np.outer(np.ones(np.size(theta_cyl)), z_mesh_cyl)

ax.plot_surface(x_cyl, y_cyl, z_cyl, rstride=4, cstride=4, color='b')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()
