#########################################################################
###                 GENERATE PARTICLE DIFFUSION PATHS                 ###
#########################################################################

#####################      ENTER PARAMETERS       #######################

mesh_dir = 'data/'                      # mesh directory/ relative to this file
mesh_filename_base = 'dendritic_mesh'   # mesh filename without extension
mat_filename = 'vcell.mat'

k = 6.0                                 # set global diffusion rate
num_points = 50                         # set number of particles
num_steps = 100                         # set number of iterated steps

#########################################################################

##########################     PATH SETUP     ###########################
# REQUIRES:
# advance_one_step.mex
# export PATH="$PATH:~/Library/Frameworks/Python.framework/Versions/2.7/bin"
# source "~/Applications/FEniCS.app/Contents/Resources/share/fenics/fenics.conf"
#
# DISABLE:
# # source "~/Library/Enthought/Canopy_64bit/User/bin/activate"
#########################################################################

########################     PYTHON IMPORTS     #########################
from particle_diffusion_on_mesh import *
import numpy as np
import scipy.io as sio

rand = np.random.uniform
randn = np.random.randn
zeros = np.zeros
vec = np.array
#########################################################################

def savexyz(num_points, mesh_wrapper, num_steps, vcell):

  points = [Point(mesh_wrapper) for _ in xrange(num_points)]

  for nn in range(0,num_steps):
    mm=0
    for pt in (points):
      
      pt.move()
      xyz = (pt.x, pt.y, pt.z)

      vcell[nn,mm,:] = xyz
      mm = mm + 1
    if nn % 50 == 0:
      print 'step:', nn

  return vcell
# #######################################################

def _matrices_only_helper(filename):
    """Loads a previously serialized Mesh object."""
    print 'Loading mesh data from NPZ file', filename
    npzfile = np.load(filename)

    k = npzfile['k'].item()
    initial_point = npzfile['initial_point']
    initial_face_index = npzfile['initial_face_index'].item()

    all_vertices = npzfile['all_vertices']
    triangles = npzfile['triangles']
    face_local_bases = npzfile['face_local_bases']
    neighbor_faces = npzfile['neighbor_faces']

    return [k, initial_point, initial_face_index,
            all_vertices, triangles, face_local_bases, neighbor_faces]


def run_with_matrices_only(serialized_mesh_filename, num_points = 10, num_steps = 20):
    
    constructor_args = _matrices_only_helper(serialized_mesh_filename)
    constructor_args[0] = 6.0  # Override k.

    initial_point, initial_face_index = constructor_args[1:3]
    xyz_loc = np.vstack([initial_point.copy() for _ in xrange(num_points)])
    face_indices = initial_face_index * np.ones((num_points, 1), dtype=np.int64)

    print_frequency = 10
    for step_num in xrange(1, num_steps + 1):
        if print_frequency is not None and step_num % print_frequency == 0:
            print 'Step Number:', step_num

        xyz_loc, face_indices = advance_one_step(xyz_loc, face_indices, *constructor_args)
    return xyz_loc



#########################################################################
###                 GENERATE PARTICLE DIFFUSION PATHS                 ###
#########################################################################

mesh_filename = mesh_dir+mesh_filename_base+'.xml'
serialized_mesh_filename = mesh_dir+mesh_filename_base+'_serialized'+'.npz'
mesh_wrapper = Mesh.from_file(serialized_mesh_filename)
mesh_wrapper.k = k          # set global diffusion rate

vcell = zeros((num_steps, num_points, 3))
print vcell.shape

vcell = savexyz(num_points, mesh_wrapper, num_steps, vcell)
print vcell

mat_file = mesh_dir+mat_filename
sio.savemat(mat_file, {'vcell':vcell})

run_with_matrices_only(serialized_mesh_filename, num_points, num_steps)
#########################################################################
