load trimesh_3D_data_96.mat
patch('Vertices', points, 'Faces', triangles, ...
      'FaceVertexCData', jet(length(triangles)), ...
      'FaceColor', 'flat');
centroids = zeros(size(triangles));
for i = 1:size(triangles, 1)
  triIndices = triangles(i, :);
  centroids(i, :) = mean(points(triIndices, :));
end

% http://www.mathworks.com/matlabcentral/newsreader/view_thread/171376
th = text(centroids(:, 1), centroids(:, 2), centroids(:, 3), ...
          num2cell(1:size(centroids, 1)), ...
          'fontsize', 10);
