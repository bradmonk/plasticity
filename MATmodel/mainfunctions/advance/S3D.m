
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
this=fileparts(which('S3D.m')); addpath(this); cd(this);


%% LOAD PYTHON AND ADD TO PATH

if count(py.sys.path,'') == 0
    insert(py.sys.path,int32(0),'');
end


%% LOAD SERIALIZED MESH AND PARTICLE DIFFUSION DATA

load('dendritic_meshs_serialized.mat');


%% ---------------- USER-ENTERED PARAMETERS ----------------


Nsteps = 2000;      % NUMBER OF DIFFUSION STEPS TO GENERATE
Nparticles = 50;    % NUMBER OF PARTICLES TO DIFFUSE ON MESH
k = 100.0;            % STDEV OF DIFFUSION RATE STEP-SIZE DISTRIBUTION


% ----------------------------------------------------------



%% GENERATE MULTIPLE PARTICLES USING REPMAT
xyz = initial_point;
face = initial_face_index;
xyz = repmat(xyz, Nparticles, 1);
face = repmat(face, Nparticles, 1);
k = repmat(k, Nparticles, 1);
% k(10:end) = 0.5;
korig=k;


%% -- PROCESS VERT & TET DATA

PYdverts = all_vertices;
PYdtets = triangles + 1;
% dvertsmin = round(min(min(PYdverts)));
% dverts = PYdverts - dvertsmin;
% FBpoints = dverts;
dtets = double(PYdtets);
FBpoints = PYdverts;
FBtri = dtets;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%				S1 & S2 ACTIN MULTIPLEX TIP DATA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


AS1doup = 1;
AS2doup = 1;
UpdTActS1 = 10 * AS1doup;
UpdTActS2 = 10 * AS2doup;
S1TipFile = which('ATdataS1.mat');
S2TipFile = which('ATdataS2.mat');
ActS1scell = 10;
ActS2scell = 10;


load(S1TipFile);

assignin('base', 'ATs', ATs);
ATs1 = evalin('base', 'ATs');

assignin('base', 'Ax', Ax);
Ax1 = evalin('base', 'Ax');
ActinS1 = Ax1{3};

assignin('base', 'AMx', AMx);
AMx1 = evalin('base', 'AMx');

assignin('base', 'AFMx', AFMx);
AFMx1 = evalin('base', 'AFMx');

% ------------
ActinTipsP1 = ATs1;
TipCellN = numel(ActinTipsP1);

ACTnP1 = ActS1scell;
NumTipCells = TipCellN - (ACTnP1*2);
StepsPerCellP1 = ceil(Nsteps / NumTipCells);

TrimActMx = size(ActinTipsP1{1},1)+1;
ACTINp1 = ActinTipsP1{ACTnP1};
ACTINp1(:,TrimActMx:end) = [];
ACTINp1(TrimActMx:end,:) = [];


% Kernel Setup

MSK = {2.5,0,0,.18,11,.1};
hkMask = SAPkernel(MSK);

AMask=ones(AMx1{5});
Ahk = convn(ACTINp1,AMask,'same');
S = (Ahk>0).*1.0;

% ------------
S1=S;
[S1Ty,S1Tx] = find(ACTINp1);
%------------

% imagesc(S1)
PlotActin(Ax1{1},Ax1{2},Ax1{3},Ax1{4},Ax1{5},Ax1{6},Ax1{7});
pause(3);
close all



%%
% keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%				S2 ACTIN MULTIPLEX TIP DATA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% GET SPINE INFO AND SET STATS COLLECTORS
%{

Neck_Bot_X = 0;     % S1
Neck_Bot_Y = 0;     % all spines
Neck_Bot_Z = 0;     % all spines
Neck_Bot_R = 7;     % all spines

Neck_Top_X = 0;     % S1
Neck_Top_Y = 0;     % all spines
Neck_Top_Z = 30;    % all spines
Neck_Top_R = 7;     % all spines


Head_Bot_X = 0;     % S1
Head_Bot_Y = 0;     % all spines
Head_Bot_Z = 30;    % all spines
Head_Bot_R = 20;    % all spines

Head_Top_X = 0;     % S1
Head_Top_Y = 0;     % all spines
Head_Top_Z = 50;    % all spines
Head_Top_R = 20;    % all spines


Neck_Bot_X = 100;     % S2
Neck_Top_X = 100;     % S2
Head_Bot_X = 100;     % S2
Head_Top_X = 100;     % S2

Neck_Bot_X = 200;     % S3
Neck_Top_X = 200;     % S3
Head_Bot_X = 200;     % S3
Head_Top_X = 200;     % S3

Neck_Bot_X = 300;     % S4
Neck_Top_X = 300;     % S4
Head_Bot_X = 300;     % S4
Head_Top_X = 300;     % S4






%}

%{

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
SPYa.Zmin = 30;
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
SPYb.Zmin = 30;
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
SPYc.Zmin = 30;
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
SPYd.Zmin = 30;
SPYd.Zmax = SPYd.Zsegmax + 5;

SPYa.N = 0;
SPYb.N = 0;
SPYc.N = 0;
SPYd.N = 0;


% SPY.N = zeros(round(Nsteps/10/2),4);
SPY.N = zeros(round(Nsteps/10),4);


%}





mesh_xmin = min(PYdverts(:,1));
mesh_xmax = max(PYdverts(:,1));
mesh_ymin = min(PYdverts(:,2));
mesh_ymax = max(PYdverts(:,2));
mesh_zmin = min(PYdverts(:,3));
mesh_zmax = max(PYdverts(:,3));


%% -- CREATE FIGURES, AXES, PLOTS


% SET AXES LIMITS
lim.x = [mesh_xmin-100 mesh_xmax+100];
lim.y = [mesh_ymin-500 mesh_ymax+500];
lim.z = [mesh_zmin-100 mesh_zmax+100]; 
lim.v = [-15 30];
% lim.x = [-100 400]; lim.y = [-100 100];
% lim.z = [-100 100]; lim.v = [-30 20];
allLims = [lim.x lim.y lim.z lim.v];

% CREATE FIGURE
f1 = figure(1);
    set(f1,'OuterPosition',[200 200 1300 800],'Color',[1 1 1]);
% CREATE AXES ONE
hax1 = axes('Position',[.05 .05 .9 .9],'Color','none');
    xlim(lim.x); ylim(lim.y); zlim(lim.z); view(lim.v);
    hold on

% CREATE AXES TWO
hax2 = axes('Position',[.05 .05 .9 .9],'Color','none');
    xlim(lim.x); ylim(lim.y); zlim(lim.z); view(lim.v);
    hold on

% CREATE AXES FIVE SCATTER PLOT
hax5 = axes('Position',[.05 .05 .9 .9],'Color','none');
    xlim(lim.x); ylim(lim.y); zlim(lim.z); view(lim.v);



% PLOT ACTIN FILAMENTS
Actin = Ax1{3};
inPSD = Ax1{4};
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);


ph11a = scatter3(hax1, [SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',5,'ob');
ph11b = scatter3(hax1, [PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',5,'og');
ph11c = plot3(hax1, [Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]');

set(ph11a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph11b,'Marker','o','MarkerEdgeColor',[.1 .9 .1],'MarkerFaceColor',[.1 .9 .1]);
set(ph11c,'LineStyle','-','Color',[.5 .5 .5],'LineWidth',1.5);


% PLOT DENDRITIC MESH
    axes(hax1)
hts1 = trisurf(FBtri,FBpoints(:,1),FBpoints(:,2),FBpoints(:,3), ...
       'FaceColor',[.1 .9 .1],'FaceAlpha', 0.3); %axis off
        xlabel('µm (x)'); ylabel('µm (y)'); zlabel('µm (z)');
        title('trisurf of boundaryFacets from alphaShape');
        light('Position',[-193.5 10.8 -17.5]);
        set(hts1,'FaceLighting','flat','EdgeLighting','gouraud');
        shading interp; 
        colormap('hot'); % colormap('pink'); % colormap('winter'); 
        hold on;



% PLOT PARTICLE STARTING LOCATIONS
    set(f1,'CurrentAxes',hax5);
p2 = scatter3(hax5, xyz(:,1), xyz(:,2), xyz(:,3),...
        110,'filled','MarkerEdgeColor','none','MarkerFaceColor',[.95 .05 .05]);
        xlim(lim.x); ylim(lim.y); zlim(lim.z); view(lim.v);


pause(2)
set(ph11c,'LineStyle',':','Color',[.5 .5 .5],'LineWidth',1);





%% --- GENERATE STEPS USING advance_one_step.mexmaci64 AND PLOT

for nn = 1:Nsteps


    % ---------------- GENERATE PARTICLE STEPS ------------------
    [xyz, face] = advance_one_step(...
        xyz, face, k, initial_point, initial_face_index, ...
        all_vertices, triangles, face_local_bases, neighbor_faces);




    % ----------------- COLLECT PARTICLE DATA -------------------

    if mod(nn, 10) == 0 %&& nn > Nsteps/2
      k = korig;
      for pp = 1:Nparticles
        if xyz(pp,3)>=900 && xyz(pp,1)<=200
          for ti = 1:numel(PSDTips(:,1))

            difdis = PSDTips(ti,:) - xyz(pp,:); % ActinTips
            radis = sqrt(difdis.^2);

            if radis<=50
              k(pp) = 1;
            end
          end
        end
      end
    end


    % ---------------- PLOT PARTICLE DIFFUSION ------------------
    % if mod(nn, 2) == 0
    scatter3(hax5, xyz(:,1), xyz(:,2), xyz(:,3),...
        110,'filled','MarkerEdgeColor','none','MarkerFaceColor',[.95 .05 .05]);
        xlim(lim.x); ylim(lim.y); zlim(lim.z); view(lim.v);
        drawnow;
    % end
    % PRINT STEP TO CONSOLE
    % fprintf('\rxyz: % 4.4g % 4.4g % 4.4g\nface: % 4.4g\n',xyz(:), face);


end; % MAIN LOOP: for nn = 1:Nsteps
%%






return
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


%     if mod(nn, 10) == 0 %&& nn > Nsteps/2
% 
% 
%         k = korig;
%       for pp = 1:Nparticles
%         %if xyz(pp,1)<=SPYa.Xmax && xyz(pp,3)>=30
%          % SPYa.N = SPYa.N + 1;
%         if xyz(pp,1)>=SPYb.Xmin && xyz(pp,1)<=SPYb.Xmax && xyz(pp,3)>=30
%           SPYb.N = SPYb.N + 1;
%           k(pp) = .1;
%         elseif xyz(pp,1)>=SPYc.Xmin && xyz(pp,1)<=SPYc.Xmax && xyz(pp,3)>=30
%           SPYc.N = SPYc.N + 1;
%           k(pp) = .1;
%         elseif xyz(pp,1)>=SPYd.Xmin && xyz(pp,3)>=30
%           SPYd.N = SPYd.N + 1;
%           k(pp) = .1;
%         end
%       end
% 
%     SPY.N(nn/10,:) = [SPYa.N SPYb.N SPYc.N SPYd.N];
%     SPYa.N=0; SPYb.N=0; SPYc.N=0; SPYd.N=0;
%     end



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
