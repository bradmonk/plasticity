
from particle_diffusion_on_mesh import *
import numpy as np
import cPickle as pickle

serialized_mesh_filename = 'serialized_mesh_res_96.npz'
mesh_wrapper = Mesh.from_file(serialized_mesh_filename)

# #######################################################
def savexyz(num_points, xyz, mesh_wrapper, num_frames):

  points = [Point(mesh_wrapper) for _ in xrange(num_points)]

  for nn in range(0,num_frames):
    for pt in (points):
      
      pt.move()
      print pt.x, pt.y, pt.z, ' ... '

    xyz[nn] = [pt.x,pt.y,pt.z]
    print ' ; '

  return xyz
# #######################################################

mesh_wrapper.k = 6.0			# set global diffusion rate
num_points = 5					# set number of particles
num_frames = 50					# set number of iterated steps
xyz = np.zeros((num_frames,3))	# create n-by-3 matrix of particle locations

xyz = savexyz(num_points, xyz, mesh_wrapper, num_frames)

pickle.dump(xyz, open( "save.p", "wb" ) )

zyx = pickle.load( open( "save.p", "rb" ) )

# import pdb; pdb.set_trace()

print zyx


