% IMPORTANT! DO NOT USE clear all IF U DONT WANT MATLAB TO CRASH!
clc; close all;

% IMPORTANT! THIS FILE (run_core.m) MUST BE IN THE SAME DIRECTORY AS:
% libadvanceonestep.so & advance_one_step.mexmaci64
% AND MUST BE THE ROOT-LEVEL WORKING DIRECTORY TO RUN
this=fileparts(which('run_core.m')); addpath(this); cd(this);


%% LOAD PYTHON AND ADD TO PATH
if count(py.sys.path,'') == 0
    insert(py.sys.path,int32(0),'');
end


%% LOAD SERIALIZED MESH AND PARTICLE DIFFUSION DATA

load('dendritic_mesh_serialized.mat');


%% ---------------- USER-ENTERED PARAMETERS ----------------


Nsteps = 20;      % NUMBER OF DIFFUSION STEPS TO GENERATE
Nparticles = 4;    % NUMBER OF PARTICLES TO DIFFUSE ON MESH
k = 4.0;            % STDEV OF DIFFUSION RATE STEP-SIZE DISTRIBUTION


% ----------------------------------------------------------


%% GENERATE MULTIPLE PARTICLES USING REPMAT
xyz = initial_point;
face = initial_face_index;
xyz = repmat(xyz, Nparticles, 1);
face = repmat(face, Nparticles, 1);
k = repmat(k, Nparticles, 1);
k(3:end) = 0.1;

%% --- GENERATE STEPS USING advance_one_step.mexmaci64 AND PLOT

for nn = 1:Nsteps

    % ---------------- GENERATE PARTICLE STEPS ------------------
    [xyz, face] = advance_one_step(...
        xyz, face, k, initial_point, initial_face_index, ...
        all_vertices, triangles, face_local_bases, neighbor_faces);


    % PRINT STEP TO CONSOLE
    fprintf('\rxyz: % 4.4g % 4.4g % 4.4g\nface: % 4.4g\n',xyz(:), face);


end; % MAIN LOOP: for nn = 1:Nsteps
%%
