load trimesh_data_96.mat
trimesh(triangles, points(:, 1), points(:, 2));
axis equal;

theta = 0:0.01:2*pi;
xCirc = 100 + 10 * cos(theta);
yCirc = 0 + 10 * sin(theta);
hold on;
plot(xCirc, yCirc);

centroids = zeros(size(triangles, 1), size(points, 2));
for i = 1:size(triangles, 1)
  triIndices = triangles(i, :);
  centroids(i, :) = mean(points(triIndices, :));
end

% http://www.mathworks.com/matlabcentral/newsreader/view_thread/171376
th = text(centroids(:, 1), centroids(:, 2), ...
          num2cell(1:size(centroids, 1)), ...
          'fontsize', 10);
