function [] = S3D()
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
% this=fileparts(which('S3D.m')); addpath(this); cd(this);
cd(fileparts(which('libadvanceonestep.so')));
% cd(fileparts(which(mfilename)));


%% ---------- LOAD PYTHON AND ADD TO PATH     ----------

if count(py.sys.path,'') == 0
    insert(py.sys.path,int32(0),'');
end



%% ---------- LOAD SERIALIZED MESH AND PARTICLE DIFFUSION DATA     ----------

load('dendritic_meshs_serialized.mat');





%% ---------- USER-ENTERED PARAMETERS     ----------


Nsteps = 3000;      % NUMBER OF DIFFUSION STEPS TO GENERATE
Nparticles = 50;    % NUMBER OF PARTICLES TO DIFFUSE ON MESH
k = 100.0;            % STDEV OF DIFFUSION RATE STEP-SIZE DISTRIBUTION

modAmpaPlot = 1;        % loop iterations between updates to AMPAR plot
modActinPlot = 1;      % loop iterations between updates to Actin plot
doPlotFullFilaments = 1;




% ----------     GENERATE MULTIPLE PARTICLES     ----------
xyz = initial_point;
face = initial_face_index;
xyz = repmat(xyz, Nparticles, 1);
face = repmat(face, Nparticles, 1);
k = repmat(k, Nparticles, 1);
% k(10:end) = 0.5;
korig=k;




% ----------     PROCESS VERT & TET DATA     ----------

PYdverts = all_vertices;
PYdtets = triangles + 1;
% dvertsmin = round(min(min(PYdverts)));
% dverts = PYdverts - dvertsmin;
% FBpoints = dverts;
dtets = double(PYdtets);
FBpoints = PYdverts;
FBtri = dtets;




%% ----------     S1 & S2 ACTIN MULTIPLEX TIP DATA     ----------

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
pause(.5);
close all



% ----------     GATHER AND FORMAT FILAMENT DATA     ----------

inPSD = Ax1{4};

for nn = 1:numel(AFMx1)

    Actin = AFMx1{nn};
    ACTN{nn} = [Actin(:,3) Actin(:,4) Actin(:,6) Actin(:,7) Actin(:,9) Actin(:,10)];
    ACTN2{nn} = [Actin(:,3)+1600 Actin(:,4)+1600 Actin(:,6) Actin(:,7) Actin(:,9) Actin(:,10)];
    ACTN3{nn} = [Actin(:,3)+3200 Actin(:,4)+3200 Actin(:,6) Actin(:,7) Actin(:,9) Actin(:,10)];
    ACTN4{nn} = [Actin(:,3)+4800 Actin(:,4)+4800 Actin(:,6) Actin(:,7) Actin(:,9) Actin(:,10)];



    ActinTips = [ACTN{nn}(:,2) ACTN{nn}(:,4) ACTN{nn}(:,6)];

    [Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
    SPYTips = ActinTips(Zrow2,:);
    SPYact{nn} = SPYTips;

    [Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
    PSDTips = ActinTips(Zrow1,:);
    PSDact{nn} = PSDTips;




    ActinTips2 = [ACTN2{nn}(:,2) ACTN2{nn}(:,4) ACTN2{nn}(:,6)];

    [Zrow2,Zcol2] = find(ActinTips2(:,3) < inPSD);
    SPYTips2 = ActinTips2(Zrow2,:);
    SPYact2{nn} = SPYTips2;

    [Zrow1,Zcol1] = find(ActinTips2(:,3) > inPSD);
    PSDTips2 = ActinTips2(Zrow1,:);
    PSDact2{nn} = PSDTips2;




    ActinTips3 = [ACTN3{nn}(:,2) ACTN3{nn}(:,4) ACTN3{nn}(:,6)];

    [Zrow2,Zcol2] = find(ActinTips3(:,3) < inPSD);
    SPYTips3 = ActinTips3(Zrow2,:);
    SPYact3{nn} = SPYTips3;

    [Zrow1,Zcol1] = find(ActinTips3(:,3) > inPSD);
    PSDTips3 = ActinTips3(Zrow1,:);
    PSDact3{nn} = PSDTips3;




    ActinTips4 = [ACTN4{nn}(:,2) ACTN4{nn}(:,4) ACTN4{nn}(:,6)];

    [Zrow2,Zcol2] = find(ActinTips4(:,3) < inPSD);
    SPYTips4 = ActinTips4(Zrow2,:);
    SPYact4{nn} = SPYTips4;

    [Zrow1,Zcol1] = find(ActinTips4(:,3) > inPSD);
    PSDTips4 = ActinTips4(Zrow1,:);
    PSDact4{nn} = PSDTips4;

end






%% ----------     CREATE FIGURES, AXES, PLOTS     ----------


% GET AXES LIMITS
mesh_xmin = min(PYdverts(:,1));
mesh_xmax = max(PYdverts(:,1));
mesh_ymin = min(PYdverts(:,2));
mesh_ymax = max(PYdverts(:,2));
mesh_zmin = min(PYdverts(:,3));
mesh_zmax = max(PYdverts(:,3));

% SET AXES LIMITS
lim.x = [mesh_xmin-100 mesh_xmax+100];
lim.y = [mesh_ymin-500 mesh_ymax+500];
lim.z = [mesh_zmin-100 mesh_zmax+100]; 
lim.v = [-15 30]; 

vlim = [-15 30];
axLim  = [lim.x lim.y lim.z];



% CREATE FIGURE 'f1'
f1 = figure(1);
    set(f1,'OuterPosition',[200 200 1300 800],'Color',[1 1 1]);

% CREATE AXES 'hax1'
hax1 = axes('Position',[.05 .05 .9 .9],'Color','none','NextPlot','replacechildren');
    axis(axLim); view(vlim);
    hold on

% CREATE AXES 'hax2'
hax2 = axes('Position',[.05 .05 .9 .9],'Color','none','NextPlot','replacechildren');
    axis(axLim); view(vlim);
    hold on

% CREATE AXES 'hax3'
hax3 = axes('Position',[.05 .05 .9 .9],'Color','none','NextPlot','replacechildren');
    axis(axLim); view(vlim);
    hold on

pause(1)
% ----------------




% ----------     PRE-PLOT DENDRITIC MESH     ----------
        axes(hax2); % f1.CurrentAxes = hax2; % set(f1,'CurrentAxes',hax2)
    hts1 = trisurf(FBtri,FBpoints(:,1),FBpoints(:,2),FBpoints(:,3), ...
       'FaceColor',[.1 .9 .1],'Parent',hax2); %axis off
        xlabel('µm (x)'); ylabel('µm (y)'); zlabel('µm (z)');
        title('trisurf of boundaryFacets from alphaShape');
        % light('Position',[-193 10 -17]);
        hcl = camlight(-15, 30);
        set(hcl,'Position',[-50000 -9000   8000])
        set(hts1,'FaceAlpha', 0.3,'FaceLighting','flat','EdgeLighting','gouraud');
        shading interp; colormap('hot');
        set(hax1,'HandleVisibility','off'); % since this won't change from here on out
        hold on;
        pause(.3)

pause(.3)
% ----------------




% ----------     PRE-PLOT FILAMENTS     ----------
axes(hax1)
tic
for Sn=1:10
delete(hax1.Children)
clear ph_spyTip ph_psdTip ph_fil
for Fn=1:4

    % Sn = 1; % Current Step Number
    % Fn = 2; % Filament ID

    PSDatn = {PSDact{Sn}, PSDact2{Sn}, PSDact3{Sn}, PSDact4{Sn}};
    SPYatn = {SPYact{Sn}, SPYact2{Sn}, SPYact3{Sn}, SPYact4{Sn}};

    ph_spyTip{Fn} = scatter3(hax1, SPYatn{Fn}(:,1)', SPYatn{Fn}(:,2)', SPYatn{Fn}(:,3)',5,'ob');
    ph_psdTip{Fn} = scatter3(hax1, PSDatn{Fn}(:,1)', PSDatn{Fn}(:,2)', PSDatn{Fn}(:,3)',5,'og');


        set(ph_spyTip{Fn},'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
        set(ph_psdTip{Fn},'Marker','o','MarkerEdgeColor',[.1 .9 .1],'MarkerFaceColor',[.1 .9 .1]);


    FILs = {ACTN{Sn}, ACTN2{Sn}, ACTN3{Sn}, ACTN4{Sn}};
    P3x = [FILs{Fn}(:,1) FILs{Fn}(:,2)]';
    P3y = [FILs{Fn}(:,3) FILs{Fn}(:,4)]';
    P3z = [FILs{Fn}(:,5) FILs{Fn}(:,6)]';
    P3x(3,:)=NaN; P3y(3,:)=NaN; P3z(3,:)=NaN;
    
    ph_fil{Fn} = line(P3x(:),P3y(:),P3z(:));
    %ph_fil{Fn} = plot3(hax1, P3x,P3y,P3z);
        set(ph_fil{Fn},'LineStyle','-','Color',[.5 .5 .5],'LineWidth',1.5);
    

end
pause(1e-6)
disp(Sn)
toc
end
% ----------







% ----------     PRE-PLOT PARTICLE STARTING LOCATION     ----------
    axes(hax3); % f1.CurrentAxes = hax3; % set(f1,'CurrentAxes',hax3)
    ampaPlot = scatter3(hax3, xyz(:,1), xyz(:,2), xyz(:,3),...
        110,'filled','MarkerEdgeColor','none','MarkerFaceColor',[.95 .05 .05]);
        axis(axLim); view(vlim);


pause(.3)
% ----------








% ----------     MAIN LOOP PARAMETERS     ----------

% modAmpaPlot = 1;        % loop iterations between updates to AMPAR plot
% modActinPlot = 10;      % loop iterations between updates to Actin plot
% doPlotFullFilaments = 1;

% writerObj = VideoWriter('S3Dvid');
% writerObj.FrameRate = 10;
% open(writerObj);


PSDactin = PSDact{1};

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            ==========     INITIATE MAIN LOOP     ==========               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nn=1:100



    % ---------------- GENERATE PARTICLE STEPS ------------------
    [xyz, face] = advance_one_step(...
        xyz, face, k, initial_point, initial_face_index, ...
        all_vertices, triangles, face_local_bases, neighbor_faces);




    % ---------------- PLOT FILAMENTS ------------------
    if mod(nn,modActinPlot)==0
    % delete(hax1.Children);

        Sn = nn; % Current Step Number
        PSDatn = {PSDact{Sn}, PSDact2{Sn}, PSDact3{Sn}, PSDact4{Sn}};
        SPYatn = {SPYact{Sn}, SPYact2{Sn}, SPYact3{Sn}, SPYact4{Sn}};    
        FILs = {ACTN{Sn}, ACTN2{Sn}, ACTN3{Sn}, ACTN4{Sn}};

    for Fn=1:4

        % keyboard
        set(ph_spyTip{Fn},'XData',SPYatn{Fn}(:,1)','YData',SPYatn{Fn}(:,2)','ZData',SPYatn{Fn}(:,3)');
        set(ph_psdTip{Fn},'XData',PSDatn{Fn}(:,1)','YData',PSDatn{Fn}(:,2)','ZData',PSDatn{Fn}(:,3)');

        if doPlotFullFilaments
            P3x = [FILs{Fn}(:,1) FILs{Fn}(:,2)]';
            P3y = [FILs{Fn}(:,3) FILs{Fn}(:,4)]';
            P3z = [FILs{Fn}(:,5) FILs{Fn}(:,6)]';
            P3x(3,:)=NaN; P3y(3,:)=NaN; P3z(3,:)=NaN;

        set(ph_fil{Fn},'XData',P3x(:),'YData',P3y(:),'ZData',P3z(:));
        end


    end
    % pause(1e-6)
    end



    % ----------------- COLLECT DATA -------------------
    if mod(nn, 10) == 0 %&& nn > Nsteps/2
      k = korig;
      for pp = 1:Nparticles
        if xyz(pp,3)>=900 && xyz(pp,1)<=200
          for ti = 1:numel(PSDactin(:,1))

            difdis = PSDactin(ti,:) - xyz(pp,:); % ActinTips
            radis = sqrt(difdis.^2);

            if radis<=50
              k(pp) = 1;
            end
          end
        end
      end
    end



    % ---------------- PLOT SURFACE PARTICLES ------------------ 
    if mod(nn,modAmpaPlot)==0

        set(ampaPlot,'XData',xyz(:,1),'YData',xyz(:,2),'ZData',xyz(:,3));        
            % delete(hax3.Children);
        % scatter3(hax3, xyz(:,1), xyz(:,2), xyz(:,3),...
        %    110,'filled','MarkerEdgeColor','none','MarkerFaceColor',[.95 .05 .05]);
    
    end





if mod(nn,modAmpaPlot)==0 || mod(nn,modAmpaPlot)==0
pause(1e-5); % This pause will render all graphics
% writeVideo(writerObj,getframe);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               ==========     END MAIN LOOP     ==========                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc

% close(writerObj);




end
