
% This script (run_core.m) loads a data.mat file containing a 
% serialized mesh and initial particle locations on the mesh.
% The serialized mesh was generated using python code, and 
% requires the fenics+dolfin python packages. 

% IMPORTANT! DO NOT USE 'clear all' 
% OTHERWISE MATLAB WILL CLEAR MEX FILES AND CRASH
clc; close all; clear
% this=fileparts(which('run_core.m')); addpath(this); cd(this);


%% LOAD PYTHON AND ADD TO PATH

if count(py.sys.path,'') == 0
    insert(py.sys.path,int32(0),'');
end


%% LOAD SERIALIZED MESH AND PARTICLE DIFFUSION DATA

load('serialized_mesh_res_34.mat');
% load('serialized_mesh_res_96.mat');

k = 3.0;  % Override k.
xyz_loc = initial_point;
face_indices = initial_face_index;



%% -- PROCESS VERT & TET DATA

PYdverts = all_vertices;
PYdtets = triangles + 1;
% dvertsmin = round(min(min(PYdverts)));
% dverts = PYdverts - dvertsmin;
% FBpoints = dverts;
dtets = double(PYdtets);
FBpoints = PYdverts;
FBtri = dtets;


%% -- CREATE FIGURES, AXES, PLOTS

% SET AXES LIMITS
lim.x = [-100 400]; lim.y = [-100 100];
lim.z = [-100 100]; lim.v = [-30 20];
allLims = [lim.x lim.y lim.z lim.v];

% CREATE FIGURE
f1 = figure(1);
    set(f1,'OuterPosition',[200 300 1000 700],'Color',[1 1 1]);
% CREATE AXES ONE
hax1 = axes('Position',[.05 .05 .9 .9],'Color','none');
    xlim(lim.x); ylim(lim.y); zlim(lim.z); view(lim.v);
    hold on
% CREATE AXES TWO
hax2 = axes('Position',[.05 .05 .9 .9],'Color','none');
    xlim(lim.x); ylim(lim.y); zlim(lim.z); view(lim.v);

% PLOT DENDRITIC MESH
    axes(hax1)
hts1 = trisurf(FBtri,FBpoints(:,1),FBpoints(:,2),FBpoints(:,3), ...
       'FaceColor',[.1 .9 .1],'FaceAlpha', 0.3); %axis off
        xlabel('µm (x)'); ylabel('µm (y)'); zlabel('µm (z)');
        title('trisurf of boundaryFacets from alphaShape');
        light('Position',[-193.5 10.8 -17.5]);
        set(hts1,'FaceLighting','flat','EdgeLighting','gouraud');
        shading interp; colormap('hot'); hold on;

% PLOT PARTICLE STARTING LOCATIONS
    set(f1,'CurrentAxes',hax2);
p2 = scatter3(hax2, xyz_loc(:,1), xyz_loc(:,2), xyz_loc(:,3),...
        90,'filled','MarkerEdgeColor','none','MarkerFaceColor',[.2 .5 .7]);
        xlim(lim.x); ylim(lim.y); zlim(lim.z); view(lim.v);


%% --- GENERATE STEPS USING advance_one_step.mexmaci64 AND PLOT

for nn = 1:200

    % GENERATE PARTICLE STEPS
    [xyz_loc, face_indices] = advance_one_step(...
        xyz_loc, face_indices, k, initial_point, initial_face_index, ...
        all_vertices, triangles, face_local_bases, neighbor_faces);

    % PRINT STEP TO CONSOLE
    fprintf('\r xyz_loc : % 5.5g % 5.5g % 5.5g \n face_indices % 5.5g \n',...
        xyz_loc(:), face_indices);

    % PLOT PARTICLES
    scatter3(hax2, xyz_loc(:,1), xyz_loc(:,2), xyz_loc(:,3),...
        90,'filled','MarkerEdgeColor','none','MarkerFaceColor',[.2 .5 .7]);
        xlim(lim.x); ylim(lim.y); zlim(lim.z); view(lim.v);
        drawnow; pause(.01)
end



%% NOTES
%{

sp1 = sprintf('\r %s : %s \r','k',class(k));
sp2 = sprintf('\r %s : %s \r','initial_point',class(initial_point));
sp3 = sprintf('\r %s : %s \r','xyz_loc',class(xyz_loc));
sp4 = sprintf('\r %s : %s \r','initial_face_index',class(initial_face_index));
sp5 = sprintf('\r %s : %s \r','face_indices',class(face_indices));
sp6 = sprintf('\r %s : %s \r','all_vertices',class(all_vertices));
sp7 = sprintf('\r %s : %s \r','triangles',class(triangles));
sp8 = sprintf('\r %s : %s \r','face_local_bases',class(face_local_bases));
sp9 = sprintf('\r %s : %s \r','neighbor_faces',class(neighbor_faces));
disp([sp1 sp2 sp3 sp4 sp5 sp6 sp7 sp8 sp9])
clear sp1 sp2 sp3 sp4 sp5 sp6 sp7 sp8 sp9

%}
