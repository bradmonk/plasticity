
% This script (run_core.m) loads a data.mat file containing a
% serialized mesh and initial particle locations on the mesh.
% The serialized mesh was generated using python code, and
% requires the fenics+dolfin python packages.

% IMPORTANT! DO NOT USE 'clear all'
% OTHERWISE MATLAB WILL CLEAR MEX FILES AND CRASH
clc; close all; clear


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


Nsteps = 500;      % NUMBER OF DIFFUSION STEPS TO GENERATE
Nparticles = 10;    % NUMBER OF PARTICLES TO DIFFUSE ON MESH
k = 5.0;            % STDEV OF DIFFUSION RATE STEP-SIZE DISTRIBUTION


% ----------------------------------------------------------



%% GENERATE MULTIPLE PARTICLES USING REPMAT
xyz = initial_point;
face = initial_face_index;
xyz = repmat(xyz, Nparticles, 1);
face = repmat(face, Nparticles, 1);
k = repmat(k, Nparticles, 1);
k(3:end) = 0.1;


%% -- PROCESS VERT & TET DATA

PYdverts = all_vertices;
PYdtets = triangles + 1;
% dvertsmin = round(min(min(PYdverts)));
% dverts = PYdverts - dvertsmin;
% FBpoints = dverts;
dtets = double(PYdtets);
FBpoints = PYdverts;
FBtri = dtets;



%% GET SPINE INFO AND SET STATS COLLECTORS

% size(PYdverts)
mesh_xmin = min(PYdverts(:,1));
mesh_xmax = max(PYdverts(:,1));
mesh_ymin = min(PYdverts(:,2));
mesh_ymax = max(PYdverts(:,2));
mesh_zmin = min(PYdverts(:,3));
mesh_zmax = max(PYdverts(:,3));



SPYa.Xsegmin = -50;     SPYa.Xsegmax = 50;
SPYa.Ysegmin = -50;     SPYa.Ysegmax = 49.9;
SPYa.Zsegmin = -99.9;   SPYa.Zsegmax = 50;
SPYa.Xcenter = 0;
SPYa.Ycenter = 0;
SPYa.Zcenter = 0;
SPYa.Xmin = SPYa.Xcenter - 30;
SPYa.Xmax = SPYa.Xcenter + 30;
SPYa.Ymin = SPYa.Ysegmin - 5;
SPYa.Ymax = SPYa.Ysegmax + 5;
SPYa.Zmin = 10;
SPYa.Zmax = SPYa.Zsegmax + 5;




SPYb.Xsegmin = 50;      SPYb.Xsegmax = 150;
SPYb.Ysegmin = -50;     SPYb.Ysegmax = 49.9;
SPYb.Zsegmin = -99.9;   SPYb.Zsegmax = 50;
SPYb.Xcenter = 100;
SPYb.Ycenter = 0;
SPYb.Zcenter = 0;
SPYb.Xmin = SPYb.Xcenter - 30;
SPYb.Xmax = SPYb.Xcenter + 30;
SPYb.Ymin = SPYb.Ysegmin - 5;
SPYb.Ymax = SPYb.Ysegmax + 5;
SPYb.Zmin = 10;
SPYb.Zmax = SPYb.Zsegmax + 5;



SPYc.Xsegmin = 150;     SPYc.Xsegmax = 250;
SPYc.Ysegmin = -50;     SPYc.Ysegmax = 49.9;
SPYc.Zsegmin = -99.9;   SPYc.Zsegmax = 50;
SPYc.Xcenter = 200;
SPYc.Ycenter = 0;
SPYc.Zcenter = 0;
SPYc.Xmin = SPYc.Xcenter - 30;
SPYc.Xmax = SPYc.Xcenter + 30;
SPYc.Ymin = SPYc.Ysegmin - 5;
SPYc.Ymax = SPYc.Ysegmax + 5;
SPYc.Zmin = 10;
SPYc.Zmax = SPYc.Zsegmax + 5;



SPYd.Xsegmin = 250;     SPYd.Xsegmax = 350;
SPYd.Ysegmin = -50;     SPYd.Ysegmax = 49.9;
SPYd.Zsegmin = -99.9;   SPYd.Zsegmax = 50;
SPYd.Xcenter = 300;
SPYd.Ycenter = 0;
SPYd.Zcenter = 0;
SPYd.Xmin = SPYd.Xcenter - 30;
SPYd.Xmax = SPYd.Xcenter + 30;
SPYd.Ymin = SPYd.Ysegmin - 5;
SPYd.Ymax = SPYd.Ysegmax + 5;
SPYd.Zmin = 10;
SPYd.Zmax = SPYd.Zsegmax + 5;

SPYa.N = 0;
SPYb.N = 0;
SPYc.N = 0;
SPYd.N = 0;


% SPY.N = zeros(round(Nsteps/10/2),4);
SPY.N = zeros(round(Nsteps/10),4);


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
p2 = scatter3(hax2, xyz(:,1), xyz(:,2), xyz(:,3),...
        90,'filled','MarkerEdgeColor','none','MarkerFaceColor',[.2 .5 .7]);
        xlim(lim.x); ylim(lim.y); zlim(lim.z); view(lim.v);



%% --- GENERATE STEPS USING advance_one_step.mexmaci64 AND PLOT

for nn = 1:Nsteps

    % ---------------- GENERATE PARTICLE STEPS ------------------
    [xyz, face] = advance_one_step(...
        xyz, face, k, initial_point, initial_face_index, ...
        all_vertices, triangles, face_local_bases, neighbor_faces);




    % ----------------- COLLECT PARTICLE DATA -------------------
    if mod(nn, 10) == 0 %&& nn > Nsteps/2

      for pp = 1:Nparticles
        %if xyz(pp,1)<=SPYa.Xmax && xyz(pp,3)>=10
         % SPYa.N = SPYa.N + 1;
        if xyz(pp,1)>=SPYb.Xmin && xyz(pp,1)<=SPYb.Xmax && xyz(pp,3)>=10
          SPYb.N = SPYb.N + 1;
        elseif xyz(pp,1)>=SPYc.Xmin && xyz(pp,1)<=SPYc.Xmax && xyz(pp,3)>=10
          SPYc.N = SPYc.N + 1;
        elseif xyz(pp,1)>=SPYd.Xmin && xyz(pp,3)>=10
          SPYd.N = SPYd.N + 1;
        end
      end

    SPY.N(nn/10,:) = [SPYa.N SPYb.N SPYc.N SPYd.N];
    SPYa.N=0; SPYb.N=0; SPYc.N=0; SPYd.N=0;
    end


    % ---------------- PLOT PARTICLE DIFFUSION ------------------
    % if mod(nn, 2) == 0
    scatter3(hax2, xyz(:,1), xyz(:,2), xyz(:,3),...
        90,'filled','MarkerEdgeColor','none','MarkerFaceColor',[.2 .5 .7]);
        xlim(lim.x); ylim(lim.y); zlim(lim.z); view(lim.v);
        drawnow;
    % end
    % PRINT STEP TO CONSOLE
    % fprintf('\rxyz: % 4.4g % 4.4g % 4.4g\nface: % 4.4g\n',xyz(:), face);


end; % MAIN LOOP: for nn = 1:Nsteps
%%


%% --- PLOT FINAL FIGURES

% CompressN = 20;
% subs = floor(linspace(1,CompressN+1,numel(SPY.N(:,1))));
% subs(end) = CompressN;
% for tt = 1:3
%     SPYN(:,tt) = accumarray(subs',SPY.N(:,tt),[],@mean);
% end


SMETH = {'moving','lowess','loess','sgolay','rlowess','rloess'};
for tt = 1:3
    SPYN(:,tt) = smooth( SPY.N(:,tt) , .15 , SMETH{2});
end


% CREATE FIGURE 2
f2 = figure(2); set(f2,'OuterPosition',[100 200 900 600],'Color',[1 1 1]);
hax3 = axes('Position',[.05 .05 .9 .9],'Color','none');


% PLOT SPINE COUNTS OVER TIME
    axes(hax3)
plot(SPYN)



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
