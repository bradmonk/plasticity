import dolfin
import numpy as np


if not dolfin.has_cgal():
  raise ImportError('Does not have CGAL.')
else:
  print 'Done importing'
  print '=' * 50


SCALE_FACTOR = 50  # dolfin.Mesh has issues below a certain scale
# Top cone data for synapse
TOP_CONE_TOP_Z = SCALE_FACTOR * 1.0
TOP_CONE_TOP_RADIUS = SCALE_FACTOR * 0.2
# Cone boundary data for synapse
CONE_BOUNDARY_Z = SCALE_FACTOR * (1.0 - 0.16)
CONE_BOUNDARY_RADIUS = SCALE_FACTOR * 0.3
# Bottom cone data for synapse
BOTTOM_CONE_BOTTOM_Z = SCALE_FACTOR * (1.0 - 3 * 0.16)
BOTTOM_CONE_BOTTOM_RADIUS = SCALE_FACTOR * 0.1
# Cylinder data for dendritic shaft
CYLINDER_TOP = BOTTOM_CONE_BOTTOM_Z
# Make bottom artificially lower, since it has to intersect the shaft
# beneath it.
CYLINDER_BOTTOM = SCALE_FACTOR * (-0.1)
CYLINDER_RADIUS = SCALE_FACTOR * 0.1
# Cylinder data for dendrite (main line)
DENDRITE_RADIUS = SCALE_FACTOR * 1.0
# Assumes we'll have 4 dendritic shafts.
DENDRITE_LEFT = - DENDRITE_RADIUS
DENDRITE_RIGHT = 7 * DENDRITE_RADIUS
SHAFT_HORIZ_SHIFT = 2 * DENDRITE_RADIUS


def get_shaft(shaft_index):
  curr_shift = (shaft_index - 1) * SHAFT_HORIZ_SHIFT
  top_cone = dolfin.Cone(dolfin.Point(curr_shift, 0, CONE_BOUNDARY_Z),
                         dolfin.Point(curr_shift, 0, TOP_CONE_TOP_Z),
                         CONE_BOUNDARY_RADIUS, TOP_CONE_TOP_RADIUS)
  bottom_cone = dolfin.Cone(dolfin.Point(curr_shift, 0, BOTTOM_CONE_BOTTOM_Z),
                            dolfin.Point(curr_shift, 0, CONE_BOUNDARY_Z),
                            BOTTOM_CONE_BOTTOM_RADIUS, CONE_BOUNDARY_RADIUS)
  cylinder = dolfin.Cone(dolfin.Point(curr_shift, 0, CYLINDER_BOTTOM),
                         dolfin.Point(curr_shift, 0, CYLINDER_TOP),
                         CYLINDER_RADIUS, CYLINDER_RADIUS)
  return top_cone + bottom_cone + cylinder


def get_dendrite():
  return dolfin.Cone(
      dolfin.Point(DENDRITE_LEFT, 0, -DENDRITE_RADIUS),
      dolfin.Point(DENDRITE_RIGHT, 0, -DENDRITE_RADIUS),
      DENDRITE_RADIUS, DENDRITE_RADIUS)


def get_geometry():
  dendrite_cone = get_dendrite()

  geometry_first_shaft = get_shaft(1)
  geometry_second_shaft = get_shaft(2)
  geometry_third_shaft = get_shaft(3)
  geometry_fourth_shaft = get_shaft(4)

  # Combine all geometries
  return (geometry_first_shaft + geometry_second_shaft +
          geometry_third_shaft + geometry_fourth_shaft +
          dendrite_cone)


def main():
  geometry_3d = get_geometry()
  dolfin.info(geometry_3d, True)
  resolution = 96  # 32 * 3
  print '=' * 50
  print 'Creating mesh'
  print '=' * 50
  mesh_3d = dolfin.Mesh(geometry_3d, resolution)
  print '=' * 50
  print 'Plotting in interactive mode'
  print '=' * 50
  dolfin.plot(mesh_3d, '3D mesh', interactive=True)


if __name__ == '__main__':
  main()
