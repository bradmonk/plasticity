load('data/serialized_mesh_res_96.mat');
xyz_loc = initial_point;
face_indices = [initial_face_index];
[xyz_loc, face_indices] = advance_one_step(...
    xyz_loc, face_indices, k, initial_point, initial_face_index, ...
    all_vertices, triangles, face_local_bases, neighbor_faces);
disp('xyz_loc:');
disp(xyz_loc);
disp('face_indices:');
disp(face_indices);
