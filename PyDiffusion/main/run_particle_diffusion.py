
from particle_diffusion_on_mesh import *
import numpy as np
import cPickle as pickle


# #######################################################
def savexyz(num_points, xyz, mesh_wrapper, num_steps):

  points = [Point(mesh_wrapper) for _ in xrange(num_points)]

  for nn in range(0,num_steps):
    for pt in (points):
      
      pt.move()
      # print pt.x, pt.y, pt.z, ' ... '

    xyz[nn] = [pt.x,pt.y,pt.z]
    # print ' ; '

  return xyz
# #######################################################


# incorporated from sample_code.py in dolfin folder
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


def run_with_matrices_only(num_points = 10, num_steps = 20, print_frequency = 10, resolution = 96):

    serialized_mesh_filename = 'serialized_mesh_res_%d.npz' % resolution
    constructor_args = _matrices_only_helper(serialized_mesh_filename)
    constructor_args[0] = 6.0  # Override k.

    initial_point, initial_face_index = constructor_args[1:3]
    xyz_loc = np.vstack([initial_point.copy() for _ in xrange(num_points)])
    face_indices = initial_face_index * np.ones((num_points, 1), dtype=np.int64)

    for step_num in xrange(1, num_steps + 1):
        if print_frequency is not None and step_num % print_frequency == 0:
            print 'Step Number:', step_num

        xyz_loc, face_indices = advance_one_step(xyz_loc, face_indices,
                                                 *constructor_args)
    return xyz_loc
# #######################################################



serialized_mesh_filename = 'serialized_mesh_res_96.npz'
mesh_wrapper = Mesh.from_file(serialized_mesh_filename)

mesh_wrapper.k = 6.0			# set global diffusion rate
num_points = 5					# set number of particles
num_steps = 50					# set number of iterated steps
xyz = np.zeros((num_steps,3))	# create n-by-3 matrix of particle locations

xyz = savexyz(num_points, xyz, mesh_wrapper, num_steps)

pickle.dump(xyz, open( "save.p", "wb" ) )

zyx = pickle.load( open( "save.p", "rb" ) )

# import pdb; pdb.set_trace()

print zyx


xyz_loc = run_with_matrices_only(num_points, num_steps)

print 'zyx_loc'
print xyz_loc

