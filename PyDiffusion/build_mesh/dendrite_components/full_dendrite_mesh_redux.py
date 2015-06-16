import dolfin
import numpy as np
import os


if not dolfin.has_cgal():
  raise ImportError('Does not have CGAL.')
else:
  print 'Done importing'
  print '=' * 60


SCALE_FACTOR = 50    # dolfin.Mesh has issues below a certain scale


#---------------- SPINE OBJECT SHAPE A -----------------#
TOTAL_LENGTH = 1.5   # total length of spine

# Z-AXIS VALUES FOR HIGHEST POINT OF SPINE OBJECTS
Z_SYNTOP = SCALE_FACTOR * (TOTAL_LENGTH - 0.10)  # SYNAPSE_SURFACE
Z_PERSYN = SCALE_FACTOR * (TOTAL_LENGTH - 0.20)  # PERISYNAPTIC_AREA
Z_HEDTOP = SCALE_FACTOR * (TOTAL_LENGTH - 0.30)  # SPINE_HEAD_UPPER
Z_HEDMID = SCALE_FACTOR * (TOTAL_LENGTH - 0.50)  # SPINE_HEAD_MIDDLE
Z_HEDLOW = SCALE_FACTOR * (TOTAL_LENGTH - 0.70)  # SPINE_HEAD_LOWER
Z_SPNECK = SCALE_FACTOR * (TOTAL_LENGTH - 0.90)  # SPINE_NECK
Z_NKLINK = SCALE_FACTOR * (TOTAL_LENGTH - 1.40)  # NECK_LINKER
Z_SHLINK = SCALE_FACTOR * (TOTAL_LENGTH - 1.50) - .1  # SHAFT_LINKER
# SHAFT_LINKER must sit arbitrarily below dendritic shaft, hence "- .1"

# RADIUS AT TOP & BOTTOM OF SPINE OBJECTS
R_SYNTOP = SCALE_FACTOR * 0.20		# SURF SYNAPSE_SURFACE
R_tPERSYN = R_SYNTOP              # TOP PERISYNAPTIC_AREA
R_bPERSYN = SCALE_FACTOR * 0.25		# BOT PERISYNAPTIC_AREA
Z_tHEDTOP = R_bPERSYN				      # TOP SPINE_HEAD_UPPER
Z_bHEDTOP = SCALE_FACTOR * 0.30		# BOT SPINE_HEAD_UPPER
R_tHEDMID = Z_bHEDTOP				      # TOP SPINE_HEAD_MIDDLE
R_bHEDMID = SCALE_FACTOR * 0.26		# BOT SPINE_HEAD_MIDDLE
R_tHEDLOW = R_bHEDMID				      # TOP SPINE_HEAD_LOWER
R_bHEDLOW = SCALE_FACTOR * 0.09		# BOT SPINE_HEAD_LOWER
R_tSPNECK = R_bHEDLOW				      # TOP SPINE_NECK
R_bSPNECK = SCALE_FACTOR * 0.09		# BOT SPINE_NECK
R_tNKLINK = R_bSPNECK				      # TOP NECK_LINKER
R_bNKLINK = SCALE_FACTOR * 0.10		# BOT NECK_LINKER
R_tSHLINK = R_bNKLINK				      # TOP SHAFT_LINKER
R_bSHLINK = R_tSHLINK				      # BOT SHAFT_LINKER
# (LAST TWO MUST BE THE SAME VALUE)




#---------------- SPINE OBJECT SHAPE B -----------------#
TOTAL_LENGTH = 1.3   # total length of spine

# Z-AXIS VALUES FOR HIGHEST POINT OF SPINE OBJECTS
Z_SYNTOP2 = SCALE_FACTOR * (TOTAL_LENGTH - 0.20)  # SYNAPSE_SURFACE
Z_PERSYN2 = SCALE_FACTOR * (TOTAL_LENGTH - 0.25)  # PERISYNAPTIC_AREA
Z_HEDTOP2 = SCALE_FACTOR * (TOTAL_LENGTH - 0.30)  # SPINE_HEAD_UPPER
Z_HEDMID2 = SCALE_FACTOR * (TOTAL_LENGTH - 0.50)  # SPINE_HEAD_MIDDLE
Z_HEDLOW2 = SCALE_FACTOR * (TOTAL_LENGTH - 0.70)  # SPINE_HEAD_LOWER
Z_SPNECK2 = SCALE_FACTOR * (TOTAL_LENGTH - 0.95)  # SPINE_NECK
Z_NKLINK2 = SCALE_FACTOR * (TOTAL_LENGTH - 1.20)  # NECK_LINKER
Z_SHLINK2 = SCALE_FACTOR * (TOTAL_LENGTH - 1.30) - .1  #Z= 0.00  SHAFT_LINKER
# SHAFT_LINKER must sit arbitrarily below dendritic shaft, hence "- .1"


# RADIUS AT TOP & BOTTOM OF SPINE OBJECTS
R_SYNTOP2 = SCALE_FACTOR * 0.28		# SURF SYNAPSE_SURFACE
R_tPERSYN2 = R_SYNTOP2            # TOP PERISYNAPTIC_AREA
R_bPERSYN2 = SCALE_FACTOR * 0.40  # OT PERISYNAPTIC_AREA
Z_tHEDTOP2 = R_bPERSYN2           # TOP SPINE_HEAD_UPPER
Z_bHEDTOP2 = SCALE_FACTOR * 0.45	# BOT SPINE_HEAD_UPPER
R_tHEDMID2 = Z_bHEDTOP2				    # TOP SPINE_HEAD_MIDDLE
R_bHEDMID2 = SCALE_FACTOR * 0.40  # BOT SPINE_HEAD_MIDDLE
R_tHEDLOW2 = R_bHEDMID2				    # TOP SPINE_HEAD_LOWER
R_bHEDLOW2 = SCALE_FACTOR * 0.15	# BOT SPINE_HEAD_LOWER
R_tSPNECK2 = R_bHEDLOW2				    # TOP SPINE_NECK
R_bSPNECK2 = SCALE_FACTOR * 0.12	# BOT SPINE_NECK
R_tNKLINK2 = R_bSPNECK2				    # TOP NECK_LINKER
R_bNKLINK2 = SCALE_FACTOR * 0.15  # BOT NECK_LINKER
R_tSHLINK2 = R_bNKLINK2				    # TOP SHAFT_LINKER
R_bSHLINK2 = R_tSHLINK2				    # BOT SHAFT_LINKER
# (LAST TWO MUST BE THE SAME VALUE)








#---------------- DENDRITIC SHAFT OBJECT -----------------#
# [ DENDRITE_SHAFT | SHAFT_LINKER | NONE ]
DENDRITE_RADIUS = SCALE_FACTOR * 1.0
DENDRITE_LEFT = - DENDRITE_RADIUS  # Assumes we'll have 4 dendritic shafts.
DENDRITE_RIGHT = 7 * DENDRITE_RADIUS
SHAFT_HORIZ_SHIFT = 2 * DENDRITE_RADIUS


def get_shaft(shaft_index):
  curr_shift = (shaft_index - 1) * SHAFT_HORIZ_SHIFT

  # GENERIC EXAMPLE OF WHAT GETS SENT TO dolfin.Cone()
  # SPINE_OBJECT = dolfin.Cone(dolfin.Point(curr_shift, 0, Z_UPPER   ),
  #                            dolfin.Point(curr_shift, 0, Z_LOWER   ),
  #                                                        BOT_RAD ,
  #                                                        TOP_RAD )

  if ((shaft_index == 1) or (shaft_index == 3)):

    SPINE_PSD = dolfin.Cone(dolfin.Point(curr_shift, 0, Z_PERSYN   ),
                               dolfin.Point(curr_shift, 0, Z_SYNTOP   ),
                                                           R_bPERSYN ,
                                                           R_tPERSYN )


    SPINE_SYNAPSE = dolfin.Cone(dolfin.Point(curr_shift, 0, Z_HEDTOP   ),
                               dolfin.Point(curr_shift, 0, Z_PERSYN   ),
                                                           R_tHEDMID ,
                                                           R_bPERSYN )


    SPINE_HEAD_UPPER = dolfin.Cone(dolfin.Point(curr_shift, 0, Z_HEDMID   ),
                               dolfin.Point(curr_shift, 0, Z_HEDTOP   ),
                                                           R_bHEDMID ,
                                                           R_tHEDMID )

    SPINE_HEAD_LOWER = dolfin.Cone(dolfin.Point(curr_shift, 0, Z_HEDLOW   ),
                               dolfin.Point(curr_shift, 0, Z_HEDMID   ),
                                                           R_bHEDLOW ,
                                                           R_tHEDLOW )

    SPINE_NECK = dolfin.Cone(dolfin.Point(curr_shift, 0, Z_SPNECK   ),
                               dolfin.Point(curr_shift, 0, Z_HEDLOW   ),
                                                           R_bSPNECK ,
                                                           R_tSPNECK )

    NECK_LINKER = dolfin.Cone(dolfin.Point(curr_shift, 0, Z_NKLINK   ),
                               dolfin.Point(curr_shift, 0, Z_SPNECK   ),
                                                           R_bNKLINK ,
                                                           R_tNKLINK )

    SHAFT_LINKER = dolfin.Cone(dolfin.Point(curr_shift, 0, Z_SHLINK   ),
                               dolfin.Point(curr_shift, 0, Z_NKLINK   ),
                                                           R_bSHLINK ,
                                                           R_tSHLINK )

    SPINE_ALL = (SPINE_SYNAPSE + SPINE_HEAD_UPPER + SPINE_HEAD_LOWER +
                  SPINE_NECK + NECK_LINKER + SHAFT_LINKER + SPINE_PSD)

  else:

    SPINE_PSD = dolfin.Cone(dolfin.Point(curr_shift, 0, Z_PERSYN2   ),
                               dolfin.Point(curr_shift, 0, Z_SYNTOP2   ),
                                                           R_bPERSYN2 ,
                                                           R_tPERSYN2 )


    SPINE_SYNAPSE = dolfin.Cone(dolfin.Point(curr_shift, 0, Z_HEDTOP2   ),
                               dolfin.Point(curr_shift, 0, Z_PERSYN2   ),
                                                           R_tHEDMID2 ,
                                                           R_bPERSYN2 )


    SPINE_HEAD_UPPER = dolfin.Cone(dolfin.Point(curr_shift, 0, Z_HEDMID2   ),
                               dolfin.Point(curr_shift, 0, Z_HEDTOP2   ),
                                                           R_bHEDMID2 ,
                                                           R_tHEDMID2 )

    SPINE_HEAD_LOWER = dolfin.Cone(dolfin.Point(curr_shift, 0, Z_HEDLOW2   ),
                               dolfin.Point(curr_shift, 0, Z_HEDMID2   ),
                                                           R_bHEDLOW2 ,
                                                           R_tHEDLOW2 )

    SPINE_NECK = dolfin.Cone(dolfin.Point(curr_shift, 0, Z_SPNECK2   ),
                               dolfin.Point(curr_shift, 0, Z_HEDLOW2   ),
                                                           R_bSPNECK2 ,
                                                           R_tSPNECK2 )

    NECK_LINKER = dolfin.Cone(dolfin.Point(curr_shift, 0, Z_NKLINK2   ),
                               dolfin.Point(curr_shift, 0, Z_SPNECK2   ),
                                                           R_bNKLINK2 ,
                                                           R_tNKLINK2 )

    SHAFT_LINKER = dolfin.Cone(dolfin.Point(curr_shift, 0, Z_SHLINK2   ),
                               dolfin.Point(curr_shift, 0, Z_NKLINK2   ),
                                                           R_bSHLINK2 ,
                                                           R_tSHLINK2 )

    SPINE_ALL = (SPINE_SYNAPSE + SPINE_HEAD_UPPER + SPINE_HEAD_LOWER +
                  SPINE_NECK + NECK_LINKER + SHAFT_LINKER + SPINE_PSD)


  return SPINE_ALL


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
  resolution = 32  # 32 * 3
  # filename = 'data/mesh_res_%d_boundary.xml' % resolution
  filename = 'data/mesh_res_%d_boundary.xml' % resolution

  if os.path.exists(filename):
    print '=' * 60
    print 'Mesh file exists, loading from file'
    print '=' * 60
    mesh_3d = dolfin.Mesh(filename)
  else:
    geometry_3d = get_geometry()
    dolfin.info(geometry_3d, True)
    print '=' * 60
    print 'Creating mesh'
    print '=' * 60
    mesh_3d = dolfin.Mesh(geometry_3d, resolution)

    print '=' * 60
    print 'Converting to a boundary mesh'
    print '=' * 60
    mesh_3d = dolfin.BoundaryMesh(mesh_3d, 'exterior')
    # Since much smaller, it's much easier to refine the boundary mesh
    # as well. It may be worth doing this via:
    #     mesh_3d = dolfin.mesh.refine(mesh_3d)

    print '=' * 60
    print 'Saving to file:', filename
    print '=' * 60
    data_file = dolfin.File(filename)
    data_file << mesh_3d

  # Plot in either case.
  print '=' * 60
  print 'Plotting in interactive mode'
  print '=' * 60
  fh = dolfin.plot(mesh_3d, '3D mesh', interactive=True,
    window_width=900, window_height=700)


if __name__ == '__main__':
  main()
