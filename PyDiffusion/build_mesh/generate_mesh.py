#########################################################################
###                   GENERATE 3D TETRAHEDRAL MESH                    ###
#########################################################################

#####################      ENTER PARAMETERS       #######################

mesh_filename_base = 'dendritic_mesh'   # mesh filename without extension
mesh_dir = 'data/'                      # mesh directory/ relative to this file
resolution = 32

#########################################################################

##########################     PATH SETUP     ###########################
# REQUIRES:
# export PATH="$PATH:~/Library/Frameworks/Python.framework/Versions/2.7/bin"
# source "~/Applications/FEniCS.app/Contents/Resources/share/fenics/fenics.conf"
#
# DISABLE:
# # source "~/Library/Enthought/Canopy_64bit/User/bin/activate"
#########################################################################

########################     PYTHON IMPORTS     #########################
import dolfin
import numpy as np
import os
import scipy.io

if not dolfin.has_cgal():
  raise ImportError('Does not have CGAL.')
else:
  print('Done importing')
  print('================================')
#########################################################################


SCALE_FACTOR = 50    # dolfin.Mesh has issues below a certain scale
TOTAL_LENGTH = 1.2   # total length of spine

# Z-AXIS VALUES FOR HIGHEST POINT OF SPINE OBJECTS
Z_SYNTOP = SCALE_FACTOR * (TOTAL_LENGTH - 0.20)  # SYNAPSE_SURFACE
Z_PERSYN = SCALE_FACTOR * (TOTAL_LENGTH - 0.25)  # PERISYNAPTIC_AREA
Z_HEDTOP = SCALE_FACTOR * (TOTAL_LENGTH - 0.30)  # SPINE_HEAD_UPPER
Z_HEDMID = SCALE_FACTOR * (TOTAL_LENGTH - 0.50)  # SPINE_HEAD_MIDDLE
Z_HEDLOW = SCALE_FACTOR * (TOTAL_LENGTH - 0.70)  # SPINE_HEAD_LOWER
Z_SPNECK = SCALE_FACTOR * (TOTAL_LENGTH - 0.95)  # SPINE_NECK
Z_NKLINK = SCALE_FACTOR * (TOTAL_LENGTH - 1.20)  # NECK_LINKER
Z_SHLINK = SCALE_FACTOR * (-.1)  # SHAFT_LINKER
# SHAFT_LINKER must sit arbitrarily below dendritic shaft, hence "- .1"

# RADIUS AT TOP & BOTTOM OF SPINE OBJECTS
R_SYNTOP = SCALE_FACTOR * 0.22    # SURF SYNAPSE_SURFACE
R_tPERSYN = R_SYNTOP              # TOP PERISYNAPTIC_AREA
R_bPERSYN = SCALE_FACTOR * 0.30   # BOT PERISYNAPTIC_AREA
Z_tHEDTOP = R_bPERSYN             # TOP SPINE_HEAD_UPPER
Z_bHEDTOP = SCALE_FACTOR * 0.35   # BOT SPINE_HEAD_UPPER
R_tHEDMID = Z_bHEDTOP             # TOP SPINE_HEAD_MIDDLE
R_bHEDMID = SCALE_FACTOR * 0.32   # BOT SPINE_HEAD_MIDDLE
R_tHEDLOW = R_bHEDMID             # TOP SPINE_HEAD_LOWER
R_bHEDLOW = SCALE_FACTOR * 0.15   # BOT SPINE_HEAD_LOWER
R_tSPNECK = R_bHEDLOW             # TOP SPINE_NECK
R_bSPNECK = SCALE_FACTOR * 0.12   # BOT SPINE_NECK
R_tNKLINK = R_bSPNECK             # TOP NECK_LINKER
R_bNKLINK = SCALE_FACTOR * 0.15   # BOT NECK_LINKER
R_tSHLINK = R_bNKLINK             # TOP SHAFT_LINKER
R_bSHLINK = R_tSHLINK             # BOT SHAFT_LINKER
# (LAST TWO MUST BE THE SAME VALUE)



#---------------- SPINE OBJECT SHAPE B -----------------#
TOTAL_LENGTH = 1.2   # total length of spine

# Z-AXIS VALUES FOR HIGHEST POINT OF SPINE OBJECTS
Z_SYNTOP2 = SCALE_FACTOR * (TOTAL_LENGTH - 0.20)  # SYNAPSE_SURFACE
Z_PERSYN2 = SCALE_FACTOR * (TOTAL_LENGTH - 0.25)  # PERISYNAPTIC_AREA
Z_HEDTOP2 = SCALE_FACTOR * (TOTAL_LENGTH - 0.30)  # SPINE_HEAD_UPPER
Z_HEDMID2 = SCALE_FACTOR * (TOTAL_LENGTH - 0.50)  # SPINE_HEAD_MIDDLE
Z_HEDLOW2 = SCALE_FACTOR * (TOTAL_LENGTH - 0.70)  # SPINE_HEAD_LOWER
Z_SPNECK2 = SCALE_FACTOR * (TOTAL_LENGTH - 0.95)  # SPINE_NECK
Z_NKLINK2 = SCALE_FACTOR * (TOTAL_LENGTH - 1.20)  # NECK_LINKER
Z_SHLINK2 = SCALE_FACTOR * (-.1)  #Z= 0.00  SHAFT_LINKER
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

  # if (shaft_index == 1) or (shaft_index == 3) or (shaft_index == 2) or (shaft_index == 4):
  if (shaft_index == 1) or (shaft_index == 3):

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


def main(mesh_filename='meshfile', resolution=32):

  filename = '%s.%s' % (mesh_filename, 'xml')
  

  if os.path.exists(filename):
    print('Mesh file exists, loading from file')
    mesh_3d = dolfin.Mesh(filename)
  else:

    geometry_3d = get_geometry()
    dolfin.info(geometry_3d, True)
    print('Creating mesh')
    mesh_3d = dolfin.Mesh(geometry_3d, resolution)

    print('Saving to file:'); print(filename)
    data_file = dolfin.File(filename)
    data_file << mesh_3d

  # Plot in either case.
  print('Plotting in interactive mode')
  fh = dolfin.plot(mesh_3d, '3D mesh', interactive=True,
    window_width=900, window_height=700)


#########################################################################
###                   GENERATE 3D TETRAHEDRAL MESH                    ###
#########################################################################

if __name__ == '__main__':

  mesh_filename = mesh_dir+mesh_filename_base
  main(mesh_filename, resolution)
#########################################################################
