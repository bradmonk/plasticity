import dolfin
import networkx
import numpy as np


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

  return networkx.connected_components(facet_index_graph)


def main():
  resolution = 96  # 32 * 3
  mesh_full_filename = 'mesh_res_%d_full.xml' % resolution
  mesh_3d_full = dolfin.Mesh(mesh_full_filename)
  print 'Calling mesh.init() to compute faces / edges / etc.'
  print '=' * 60
  mesh_3d_full.init()

  top_indices_filename = 'faces_top_indices_res_%d_full.npy' % resolution
  top_indices = np.load(top_indices_filename)
  top_indices = set(top_indices)

  print 'Reading facet, edge and vertex iterators into lists'
  print '=' * 60
  facets = list(dolfin.facets(mesh_3d_full))
  top_facets = [facet for i, facet in enumerate(facets)
                if i in top_indices]
  edges = list(dolfin.edges(mesh_3d_full))
  vertices = list(dolfin.vertices(mesh_3d_full))

  # First do a data check since we aren't sure if face indexing
  # is deterministic given vertices / tetrahedral shells stored in
  # the XML file.
  for facet in top_facets:
    if not facet.exterior():
      raise ValueError('Not exterior where expected.')
    theta = np.arccos(facet.normal().z())
    if not np.allclose(theta, 0):
      raise ValueError('Not vertical normal where expected.')

  print 'Computing connected components'
  print '=' * 60
  components = get_connected_components(top_facets, edges, vertices)
  print 'There are %d connected components in total' % len(components)


if __name__ == '__main__':
  main()
