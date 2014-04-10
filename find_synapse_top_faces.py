import dolfin
import numpy as np

import full_dendrite_mesh


def main():
  resolution = 96  # 32 * 3
  mesh_filename = 'mesh_res_%d.xml' % resolution
  mesh_3d = dolfin.Mesh(mesh_filename)
  print 'Calling mesh.init(2) to compute faces'
  print '=' * 60
  mesh_3d.init(2)

  coords = mesh_3d.coordinates()
  top_coords = np.nonzero(
      coords[:, 2] >= full_dendrite_mesh.TOP_CONE_TOP_Z - 0.001)[0]
  top_coords_set = set(top_coords)

  num_faces = mesh_3d.num_faces()

  faces = dolfin.faces(mesh_3d)
  vertical_normal_faces = []
  all_top_coords_faces = []
  face_count = 0

  print 'Looping through faces to find faces on tops of synapses'
  print '=' * 60
  for face in faces:
    # The vertical normal is (0, 0, 1) and we check the angle between
    # the normal and this vector v . (0, 0, 1) = v_z (dot product).
    theta = np.arccos(face.normal().z())
    if np.allclose(theta, 0):
      vertical_normal_faces.append(face)

    # If every single index is one of the "top coordinate" indices.
    if set(face.entities(0)) <= top_coords_set:
      all_top_coords_faces.append(face)

    if face_count % 20000 == 0:
      print face_count, '/', num_faces

    face_count += 1

  if vertical_normal_faces != all_top_coords_faces:
    raise ValueError('Top face classifications disagree.')

  # Vertex indices are deterministic (since stored in file) and so is the
  # geometry (via the cells), hence the indices for the faces will be the
  # same if the data was loaded again.
  face_index_matrix = np.vstack(
      [face.entities(0) for face in vertical_normal_faces])

  faces_filename = 'faces_res_%d.npy' % resolution
  print '=' * 50
  print 'Saving to file:', faces_filename
  print '=' * 50
  np.save(faces_filename, face_index_matrix)


if __name__ == '__main__':
  main()
