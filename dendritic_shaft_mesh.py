import dolfin
import numpy as np


if not dolfin.has_cgal():
  raise ImportError('Does not have CGAL.')
else:
  print 'Done importing'
  print '=' * 50


SCALE_FACTOR = 50  # dolfin.Mesh has issues below a certain scale
SPHERE_Z = SCALE_FACTOR * 1 - np.sqrt(0.05 * SCALE_FACTOR**2)
SPHERE_RADIUS = SCALE_FACTOR * 0.3
CYLINDER_BOTTOM = 0.0
CYLINDER_TOP = SPHERE_Z - np.sqrt(0.08 * SCALE_FACTOR**2)
CYLINDER_RADIUS = SCALE_FACTOR * 0.1
BOX_WIDTH = SCALE_FACTOR * 2.0
BOX_BOTTOM = SCALE_FACTOR * 1.0


sphere = dolfin.Sphere(dolfin.Point(0, 0, SPHERE_Z),
                       SPHERE_RADIUS)
cone = dolfin.Cone(dolfin.Point(0, 0, CYLINDER_BOTTOM),
                   dolfin.Point(0, 0, CYLINDER_TOP),
                   CYLINDER_RADIUS, CYLINDER_RADIUS)
box = dolfin.Box(-BOX_WIDTH / 2, -BOX_WIDTH / 2, BOX_BOTTOM,
                 BOX_WIDTH / 2, BOX_WIDTH / 2, BOX_BOTTOM + BOX_WIDTH)


geometry_3d = sphere + cone - box
dolfin.info(geometry_3d, True)
resolution = 32
print '=' * 50
print 'Creating mesh'
print '=' * 50
mesh_3d = dolfin.Mesh(geometry_3d, resolution)
dolfin.plot(mesh_3d, '3D mesh')
print '=' * 50
print 'Done plotting, entering interactive mode'
print '=' * 50
dolfin.interactive()
