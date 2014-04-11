import dolfin
import networkx
import numpy as np
import scipy.io

import full_dendrite_mesh


def get_neighbors(facet, edges, vertices, valid_face_indices):
  facet_edge_indices = facet.entities(1)
  neighbors_by_edge = set()
  for i in facet_edge_indices:
    edge = edges[i]
    neighbors_by_edge.update(
        [i for i in edge.entities(2) if i in valid_face_indices])
  # Remove the index of the current facet.
  neighbors_by_edge.remove(facet.index())

  facet_vertex_indices = facet.entities(0)
  neighbors_by_vertex = set()
  for i in facet_vertex_indices:
    vertex = vertices[i]
    neighbors_by_vertex.update(
        [i for i in vertex.entities(2) if i in valid_face_indices])
  # Remove the index of the current facet.
  neighbors_by_vertex.remove(facet.index())
  # Remove the neighbors which are already sharing an edge.
  neighbors_by_vertex = neighbors_by_vertex - neighbors_by_edge

  return sorted(neighbors_by_edge), sorted(neighbors_by_vertex)


def get_connected_components(facets, edges, vertices):
  # Get all the indices we'll use for `get_neighbors`.
  facets_by_index = {facet.index(): facet for facet in facets}
  valid_face_indices = set(facets_by_index.keys())

  # Add all the edges to a graph.
  facet_index_graph = networkx.Graph()
  for index, facet in facets_by_index.iteritems():
    # The node may already have been added as part of a previous edge.
    facet_index_graph.add_node(facet.index())

    edge_neighbors, vert_neighbors = get_neighbors(
        facet, edges, vertices, valid_face_indices)
    all_neighbors = edge_neighbors + vert_neighbors
    for neighbor_index in all_neighbors:
      facet_index_graph.add_edge(index, neighbor_index)

  return facets_by_index, networkx.connected_components(facet_index_graph)


def get_data(resolution):
  mesh_full_filename = 'mesh_res_%d_full.xml' % resolution
  mesh_3d_full = dolfin.Mesh(mesh_full_filename)
  print 'Calling mesh.init() to compute faces / edges / etc.'
  print '=' * 60
  mesh_3d_full.init()

  print 'Reading facet, edge and vertex iterators into lists'
  print '=' * 60
  facets = list(dolfin.facets(mesh_3d_full))
  edges = list(dolfin.edges(mesh_3d_full))
  vertices = list(dolfin.vertices(mesh_3d_full))

  return mesh_3d_full, facets, edges, vertices


def get_top_facets(facets, resolution):
  top_indices_filename = 'faces_top_indices_res_%d_full.npy' % resolution
  top_indices = np.load(top_indices_filename)
  top_indices = set(top_indices)

  return [facet for i, facet in enumerate(facets)
          if i in top_indices]


def get_point_in_3d(index, vertices):
  x = vertices[index].point().x()
  y = vertices[index].point().y()
  z = vertices[index].point().z()
  return np.array([x, y, z])


def cast_top_to_2d_mesh(face_indices, facets_by_index, vertices):
  # Organize the faces by the 3 vertices that define each face.
  vertices_in_faces = np.vstack([facets_by_index[index].entities(0)
                                for index in face_indices])
  # Find the index of all possible vertices on a face.
  vertex_indices = sorted(set(vertices_in_faces.flatten()))
  # Look up the (x, y, z) values for the point based on vertex index.
  vertex_points = np.vstack([get_point_in_3d(index, vertices)
                             for index in vertex_indices])

  if not np.allclose(vertex_points[:, 2], full_dendrite_mesh.TOP_CONE_TOP_Z):
    raise ValueError('On top, we expect Z-value to be same.')

  # Cast points to 2D.
  vertex_points = vertex_points[:, :2]

  # Renumber the vertices (using 1-based index for MATLAB).
  vertex_renumbering_for_matlab = {
      value: (i + 1) for i, value in enumerate(vertex_indices)}

  def renumber_method(value):
    return vertex_renumbering_for_matlab[value]
  renumber_method = np.vectorize(renumber_method)

  vertices_in_faces = renumber_method(vertices_in_faces)

  return vertex_points, vertices_in_faces


def main():
  resolution = 96  # 32 * 3
  # For some reason, this fails if mesh_3d_full is not in the
  # same scope as the values inherited from it.
  unused_mesh_3d_full, facets, edges, vertices = get_data(resolution)
  top_facets = get_top_facets(facets, resolution)

  print 'Computing connected components'
  print '=' * 60
  facets_by_index, components = get_connected_components(
      top_facets, edges, vertices)
  print 'There are %d connected components in total' % len(components)

  first_component = sorted(components[0])
  vertex_points, vertices_in_faces = cast_top_to_2d_mesh(
      first_component, facets_by_index, vertices)

  # Save these to a MATLAB file so we can use them in MATLAB.
  matlab_filename = 'trimesh_data_%d.mat' % resolution
  print '=' * 60
  print 'Saving first component 2D mesh data to file:', matlab_filename
  data = {
      'points': vertex_points,
      'triangles': vertices_in_faces,
  }
  scipy.io.savemat(matlab_filename, data)
  print '=' * 60
  print 'To use this in MATLAB perform:'
  print '>> load %s' % matlab_filename
  print '>> # triangles and points will be loaded'
  print '>> trimesh(triangles, points(:, 1), points(:, 2))'
  print '>> axis equal'
  print '=' * 60


if __name__ == '__main__':
  main()
