xyz_loc = [
    1, 3, 5.0;
    4.0, 1.0, 5.0
];
face_indices = int64([1; 2]);
k = 10.0;
initial_point = [11.0, 1.0, 2.0];
initial_face_index = int64(10);
all_vertices = [1.0, 2.0, 3.0];
triangles = int64([1, 1, 11]);
face_local_bases = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
neighbor_faces = int64([10, 1, 2]);
multiplier = 3;
inMatrix = [1.5, 2, 9];
outMatrix = advance_one_step(xyz_loc, face_indices, k, initial_point, ...
                             ...
                             initial_face_index, all_vertices, ...
                             triangles, face_local_bases, ...
                             neighbor_faces, multiplier, inMatrix)
