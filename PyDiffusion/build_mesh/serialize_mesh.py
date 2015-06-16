#########################################################################
###                   SERIALIZE 3D TETRAHEDRAL MESH                   ###
#########################################################################

#####################      ENTER PARAMETERS       #######################

mesh_filename_base = 'dendritic_mesh'   # mesh filename without extension
mesh_dir = 'data/'                      # mesh directory/ relative to this file

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
import numpy as np
from particle_diffusion_on_mesh import *
import dolfin
#########################################################################


def _load_serialized_mesh(filename):
    """Loads a previously serialized Mesh object."""
    print('Loading mesh data from NPZ file',filename)
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



def save_serialized_mesh(mesh_filename):
    print 'Importing dolfin, takes a bit of time...'
    import dolfin
    print 'Done importing dolfin.'

    mesh_full_filename = '%s.%s' % (mesh_filename, 'xml')
    mesh_3d = dolfin.Mesh(mesh_full_filename)

    # NOTE: This is temporary; we should package them in a different way.
    SCALE_FACTOR = 50.0
    STARTING_X = SCALE_FACTOR * 0.0
    STARTING_Y = SCALE_FACTOR * 0.0
    STARTING_Z = SCALE_FACTOR * 1.0
    STARTING_K = SCALE_FACTOR * 0.01

    initial_point = np.array((STARTING_X, STARTING_Y, STARTING_Z))
    mesh_wrapper = Mesh.from_mesh(mesh_3d, initial_point, STARTING_K)

    serialized_mesh_filename = '%s_serialized.%s' % (mesh_filename, 'npz')
    mesh_wrapper.serialize_mesh(serialized_mesh_filename)

    convert_mesh_to_matlab(serialized_mesh_filename)


def convert_mesh_to_matlab(serialized_mesh_filename):
    """Loads a previously serialized Mesh object and stores in MATLAB file."""
    (k, initial_point, initial_face_index,
     all_vertices, triangles,
     face_local_bases, neighbor_faces) = _load_serialized_mesh(serialized_mesh_filename)
    data = {
        'k': k,
        'initial_point': initial_point,
        'initial_face_index': initial_face_index,
        'all_vertices': all_vertices,
        'triangles': triangles,
        'face_local_bases': face_local_bases,
        'neighbor_faces': neighbor_faces,
    }

    root, ext = os.path.splitext(serialized_mesh_filename)
    matlab_filename = root + '.mat'
    scipy.io.savemat(matlab_filename, data)
    print 'Saved', matlab_filename



#########################################################################
###                   SERIALIZE 3D TETRAHEDRAL MESH                   ###
#########################################################################

if __name__ == '__main__':

    mesh_filename = mesh_dir+mesh_filename_base
    save_serialized_mesh(mesh_filename)
#########################################################################
