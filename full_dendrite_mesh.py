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
# Make artificially lower, since it has to intersect the shaft
# beneath it.
CYLINDER_BOTTOM = SCALE_FACTOR * (-0.2)
CYLINDER_TOP = SPHERE_Z - np.sqrt(0.08 * SCALE_FACTOR**2)
CYLINDER_RADIUS = SCALE_FACTOR * 0.1
BOX_WIDTH = SCALE_FACTOR * 2.0
BOX_BOTTOM = SCALE_FACTOR * 1.0
DENDRITE_RADIUS = SCALE_FACTOR * 1.0
# Assumes we'll have 3 dendritic shafts, there was a problem with 4.
DENDRITE_LEFT = - DENDRITE_RADIUS
DENDRITE_RIGHT = 5 * DENDRITE_RADIUS
SHAFT_HORIZ_SHIFT = 2 * DENDRITE_RADIUS


def get_shaft(shaft_index):
  curr_shift = (shaft_index - 1) * SHAFT_HORIZ_SHIFT
  sphere = dolfin.Sphere(dolfin.Point(curr_shift, 0, SPHERE_Z),
                         SPHERE_RADIUS)
  cone = dolfin.Cone(dolfin.Point(curr_shift, 0, CYLINDER_BOTTOM),
                     dolfin.Point(curr_shift, 0, CYLINDER_TOP),
                     CYLINDER_RADIUS, CYLINDER_RADIUS)
  box = dolfin.Box(
      -BOX_WIDTH / 2 + curr_shift, -BOX_WIDTH / 2, BOX_BOTTOM,
      BOX_WIDTH / 2 + curr_shift, BOX_WIDTH / 2, BOX_BOTTOM + BOX_WIDTH)
  return (sphere - box) + cone


dendrite_cone = dolfin.Cone(
    dolfin.Point(DENDRITE_LEFT, 0, -DENDRITE_RADIUS),
    dolfin.Point(DENDRITE_RIGHT, 0, -DENDRITE_RADIUS),
    DENDRITE_RADIUS, DENDRITE_RADIUS)

geometry_first_shaft = get_shaft(1)
geometry_second_shaft = get_shaft(2)
geometry_third_shaft = get_shaft(3)

# Combine all geometries
geometry_3d = (geometry_first_shaft + geometry_second_shaft +
               geometry_third_shaft + dendrite_cone)
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
