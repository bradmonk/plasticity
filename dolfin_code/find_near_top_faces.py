import numpy as np
import scipy.io

import face_neighbors_and_top_components as utils


def get_new_edge_neighbors(boundary_vertices,
                           all_component_vertices, all_valid_indices,
                           facets, edges, vertices):
  all_edge_neighbors = set()
  for index in boundary_vertices:
    facet = facets[index]
    facet_edge_neighbors, _ = utils.get_neighbors(facet, edges, vertices,
                                                  all_valid_indices)
    all_edge_neighbors.update(facet_edge_neighbors)

  # Remove faces already known
  all_edge_neighbors = all_edge_neighbors - set(all_component_vertices)

  return sorted(all_edge_neighbors)


def get_points(face_indices, facets, vertices):
  vertices_in_faces = np.vstack([facets[index].entities(0)
                                 for index in face_indices])
  # Find the index of all possible vertices on a face.
  vertex_indices = sorted(set(vertices_in_faces.flatten()))

  # Turn the indices into actual points
  vertex_points = np.vstack([utils.get_point_in_3d(index, vertices)
                             for index in vertex_indices])

  return vertices_in_faces, vertex_indices, vertex_points


def main():
  resolution = 96  # 32 * 3
  mesh_3d_full, facets, edges, vertices = utils.get_data(resolution)
  top_facets = utils.get_top_facets(facets, resolution)

  facets_by_index, components = utils.get_connected_components(
      top_facets, edges, vertices)

  # "First" component.
  component = sorted(components[0])
  component_facets = [facets_by_index[index] for index in component]

  # All exterior indices should be allowed.
  exterior_face_filename = 'exterior_faces_res_%d_full.npy' % resolution
  exterior_face_indices = np.load(exterior_face_filename)
  exterior_face_indices = set(exterior_face_indices)

  breadth_first_extended_layers = {}
  boundary_vertices = component[:]
  all_component_vertices = component[:]
  for i in xrange(1, 4):
    new_neighbors = get_new_edge_neighbors(boundary_vertices,
                                           all_component_vertices,
                                           exterior_face_indices,
                                           facets, edges, vertices)
    # Add the new neighbors to the dictionary.
    vertices_in_faces, vertex_indices, vertex_points = get_points(
        new_neighbors, facets, vertices)
    breadth_first_extended_layers[i] = (new_neighbors, vertices_in_faces,
                                        vertex_indices, vertex_points)
    # Add the new neighbors to the component.
    all_component_vertices += new_neighbors
    # Make the boundary equal to the newly discovered layer of neighbors.
    boundary_vertices = new_neighbors[:]

  # We are forgetting about the individual values of
  # "breadth_first_extended_layers" but they can potentially by useful.

  all_component_vertices = sorted(all_component_vertices)
  vertices_in_faces, vertex_indices, vertex_points = get_points(
      all_component_vertices, facets, vertices)

  vertex_renumbering_for_matlab = {
      value: (i + 1) for i, value in enumerate(vertex_indices)}

  def renumber_method(value):
    return vertex_renumbering_for_matlab[value]
  renumber_method = np.vectorize(renumber_method)

  vertices_in_faces = renumber_method(vertices_in_faces)

  # Save these to a MATLAB file so we can use them in MATLAB.
  matlab_filename = 'trimesh_3D_data_%d.mat' % resolution
  print '=' * 60
  print 'Saving nearby first component 3D mesh data to file:',
  print matlab_filename
  data = {
      'points': vertex_points,
      'triangles': vertices_in_faces,
  }
  scipy.io.savemat(matlab_filename, data)


if __name__ == '__main__':
  main()
