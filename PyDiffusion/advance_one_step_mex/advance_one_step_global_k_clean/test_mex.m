% IMPORTANT! DO NOT USE clear all IF U DONT WANT MATLAB TO CRASH!
clc; close all;
thisFolder = fileparts(which('test_mex.m'));
addpath(thisFolder); cd(thisFolder);


%% This must be executed to run advance_one_step from the Matlab IDE
if count(py.sys.path,'') == 0
    insert(py.sys.path,int32(0),'');
end


%% LOAD SERIALIZED MESH AND SAVED VARIABLE VALS

load('dendritic_mesh_serialized.mat');

k = 2.0;  % Override k.
xyz_loc = initial_point;                % [0 0 50]
face_indices = initial_face_index;      % 12419


%% --- MEX PYTHON LOOP advance_one_step.mexmaci64

for nn = 1:20
[xyz_loc, face_indices] = advance_one_step(...
    xyz_loc, face_indices, k, initial_point, initial_face_index, ...
    all_vertices, triangles, face_local_bases, neighbor_faces);

sp1 = sprintf('\r xyz_loc : % 5.5g % 5.5g % 5.5g \n',xyz_loc(:));
sp2 = sprintf('face_indices : % 5.5g \n',face_indices);

disp(sp1); disp(sp2);

pause(.1)
end
