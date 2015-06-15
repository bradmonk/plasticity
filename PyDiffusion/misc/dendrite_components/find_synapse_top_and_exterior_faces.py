import dolfin
import numpy as np

import full_dendrite_mesh


def main():
  resolution = 96  # 32 * 3
  mesh_full_filename = 'data/mesh_res_%d_full.xml' % resolution
  mesh_3d_full = dolfin.Mesh(mesh_full_filename)
  print 'Calling mesh.init() to compute faces / edges / etc.'
  print '=' * 60
  mesh_3d_full.init()

  coords = mesh_3d_full.coordinates()
  top_coords = np.nonzero(
      coords[:, 2] >= full_dendrite_mesh.TOP_CONE_TOP_Z - 0.001)[0]
  top_coords_set = set(top_coords)

  num_faces = mesh_3d_full.num_faces()
  num_facets = mesh_3d_full.num_facets()
  if num_faces != num_facets:
    raise ValueError('Expected tetrahedral mesh.')

  # Use facets instead of faces (though identical) because the
  # facet objects have .exterior() set.
  faces_as_facets = dolfin.facets(mesh_3d_full)
  exterior_face_indices = []
  vertical_normal_faces = []
  all_top_coords_faces = []
  face_count = 0

  print 'Looping through faces to find faces on tops of synapses'
  print '=' * 60
  for face in faces_as_facets:
    face_count += 1
    if face_count % 20000 == 0:
      print face_count, '/', num_faces

    # Only use exterior facets / faces.
    if not face.exterior():
      continue
    else:
      exterior_face_indices.append(face.index())

    # The vertical normal is (0, 0, 1) and we check the angle between
    # the normal and this vector v . (0, 0, 1) = v_z (dot product).
    theta = np.arccos(face.normal().z())
    if np.allclose(theta, 0):
      vertical_normal_faces.append(face)

    # If every single index is one of the "top coordinate" indices.
    if set(face.entities(0)) <= top_coords_set:
      all_top_coords_faces.append(face)

  # First save the exterior face indices.
  # NOTE: We are 100% sure indices are computed the same upon reloading data
  #       from XML file. This is relevant for face indices, since they are
  #       not stored in the XML file. See
  #           http://fenicsproject.org/qa/3233
  #       for details.
  exterior_face_filename = 'data/exterior_faces_res_%d_full.npz' % resolution
  print '=' * 60
  print 'Saving to file:', exterior_face_filename
  print '%d exterior faces out of %d total faces.' % (
      len(exterior_face_indices), num_faces)
  print '=' * 60
  np.savez(exterior_face_filename, np.array(exterior_face_indices))

  if vertical_normal_faces != all_top_coords_faces:
    raise ValueError('Top face classifications disagree.')

  # Indices for the faces will be the same if data loaded again. See above.
  face_index_matrix = np.vstack(
      [face.entities(0) for face in vertical_normal_faces])

  faces_full_filename = 'data/faces_top_res_%d_full.npz' % resolution
  print '=' * 60
  print 'Saving to file:', faces_full_filename
  print '=' * 60
  np.savez(faces_full_filename, face_index_matrix)

  # Save the indices of the facets on the top.
  top_indices_filename = 'data/faces_top_indices_res_%d_full.npz' % resolution
  print '=' * 60
  print 'Saving to file:', top_indices_filename
  print '=' * 60
  top_indices = [facet.index() for facet in vertical_normal_faces]
  np.savez(top_indices_filename, np.array(top_indices))


if __name__ == '__main__':
  main()
