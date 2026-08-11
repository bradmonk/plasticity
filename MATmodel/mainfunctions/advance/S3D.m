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


ns = 5000;          % NUMBER OF DIFFUSION STEPS TO GENERATE
Nparticles = 90;    % NUMBER OF PARTICLES TO DIFFUSE ON MESH
k = 100.0;          % STDEV OF DIFFUSION RATE STEP-SIZE DISTRIBUTION

modAmpaPlot = 1;    % loop iterations between updates to AMPAR plot
modActinPlot = 1;   % loop iterations between updates to Actin plot
modBindSAP = 20;


% fprintf('% 6.31f   \r',xyz(8,:))

%% ----------     GENERATE MULTIPLE PARTICLES     ----------
xyz = [2651.8740651654261455405503511428833 ...
 585.4610340686764402562403120100498 ...
-1341.0797355336774216993944719433784];
face = int64(2110);

% xyz = initial_point;
% face = initial_face_index;



xyz = repmat(xyz, Nparticles, 1);
face = repmat(face, Nparticles, 1);
k = repmat(k, Nparticles, 1);
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

% load(which('ATdataS1.mat'));
mf = matfile('ActinMx.mat');
ActinMx = mf.ActinMx;

% AFMx = ActinMx.AFMx;
% ActXYZ = ActinMx.ActXYZ;
% PSDtips = ActinMx.PSDtips;
% SPYtips = ActinMx.SPYtips;
% Ax = ActinMx.Ax;
% inPSD = ActinMx.Ax{4}(7);



%% ----------     GATHER AND FORMAT FILAMENT DATA     ----------

as = numel(ActinMx.ActXYZ);

for nn = 1:as

    ACTN1{nn} = [ActinMx.ActXYZ{nn}(:,1:2)      ActinMx.ActXYZ{nn}(:,3:4) ActinMx.ActXYZ{nn}(:,5:6)];
    ACTN2{nn} = [ActinMx.ActXYZ{nn}(:,1:2)+1600 ActinMx.ActXYZ{nn}(:,3:4) ActinMx.ActXYZ{nn}(:,5:6)];
    ACTN3{nn} = [ActinMx.ActXYZ{nn}(:,1:2)+3200 ActinMx.ActXYZ{nn}(:,3:4) ActinMx.ActXYZ{nn}(:,5:6)];
    ACTN4{nn} = [ActinMx.ActXYZ{nn}(:,1:2)+4800 ActinMx.ActXYZ{nn}(:,3:4) ActinMx.ActXYZ{nn}(:,5:6)];

    SPYact1{nn} = [ActinMx.SPYtips{nn}(:,1)+0000 ActinMx.SPYtips{nn}(:,2:3)];
    PSDact1{nn} = [ActinMx.PSDtips{nn}(:,1)+0000 ActinMx.PSDtips{nn}(:,2:3)];
    SPYact2{nn} = [ActinMx.SPYtips{nn}(:,1)+1600 ActinMx.SPYtips{nn}(:,2:3)];
    PSDact2{nn} = [ActinMx.PSDtips{nn}(:,1)+1600 ActinMx.PSDtips{nn}(:,2:3)];
    SPYact3{nn} = [ActinMx.SPYtips{nn}(:,1)+3200 ActinMx.SPYtips{nn}(:,2:3)];
    PSDact3{nn} = [ActinMx.PSDtips{nn}(:,1)+3200 ActinMx.PSDtips{nn}(:,2:3)];
    SPYact4{nn} = [ActinMx.SPYtips{nn}(:,1)+4800 ActinMx.SPYtips{nn}(:,2:3)];
    PSDact4{nn} = [ActinMx.PSDtips{nn}(:,1)+4800 ActinMx.PSDtips{nn}(:,2:3)];
 
end

% keyboard
% %%
% %------------
% % Mask Setup
% doGMask = 1;
% MSK = {2.5, 0, 0, .18, 11, .1};
% hkMask = MaskFun(doGMask,MSK);
% 
% AFMx = ActinMx.AFMx;
% 
% AMask=ones(1);
% Ahk = convn(ACTINp1,hkMask,'same');
% S = (Ahk>0).*1.0;
% 
% S1=S;
% [S1Ty,S1Tx] = find(ACTINp1);
% %------------
% 
% %-----------------------
% %-- S1 Cluster Plots --
% %-----------------------
% axes('Position',[.31 .53 .27 .44]);
% S1Ph1 = imagesc(S1);
% colormap(clrmap);
% title('Synapse-A Scaffold Clusters');
% hold on
% %--------
% [PSDY,PSDX] = find(ACTINp1);
% S1Ph2 = scatter(PSDX,PSDY);
% set(S1Ph2,'Marker','s','SizeData',60,'LineWidth',.5,...
% 	'MarkerFaceColor',[.95 .1 .1],'MarkerEdgeColor','none')
% %set(gca,'XDir','reverse')
% set(gca,'YDir','normal')
% hold off




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



% GET SCREEN SIZE AND CALCULATE DIMS FOR FIGURE
ssz = get(groot,'ScreenSize');
pos = [(ssz(3)-1360-5) (ssz(4)-838-30) 1360 838];
% 9x6 == 1360x765+73 (toolbar is 73pixels) = 1360x838

% CREATE FIGURE 'f1'
f1 = figure(1); 
    set(f1,'OuterPosition',pos,'Color',[1 1 1]);

% CREATE AXES 'hax1'
% hax1 = axes('Position',[.05 .05 .9 .9],'Color','none','NextPlot','replacechildren');
hax1 = axes('Position',[.05 .05 .9 .9],'Color','none');
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

    PSDatn = {PSDact1{Sn}, PSDact2{Sn}, PSDact3{Sn}, PSDact4{Sn}};
    SPYatn = {SPYact1{Sn}, SPYact2{Sn}, SPYact3{Sn}, SPYact4{Sn}};

    ph_spyTip{Fn} = scatter3(hax1, SPYatn{Fn}(:,1)', SPYatn{Fn}(:,2)', SPYatn{Fn}(:,3)',5,'ob');
    ph_psdTip{Fn} = scatter3(hax1, PSDatn{Fn}(:,1)', PSDatn{Fn}(:,2)', PSDatn{Fn}(:,3)',5,'og');


        set(ph_spyTip{Fn},'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
        set(ph_psdTip{Fn},'Marker','o','MarkerEdgeColor',[.1 .9 .1],'MarkerFaceColor',[.1 .9 .1]);


    FILs = {ACTN1{Sn}, ACTN2{Sn}, ACTN3{Sn}, ACTN4{Sn}};
    P3x = [FILs{Fn}(:,1) FILs{Fn}(:,2)]';
    P3y = [FILs{Fn}(:,3) FILs{Fn}(:,4)]';
    P3z = [FILs{Fn}(:,5) FILs{Fn}(:,6)]';
    P3x(3,:)=NaN; P3y(3,:)=NaN; P3z(3,:)=NaN;
    
    ph_fil{Fn} = line('XData',P3x(:),'YData',P3y(:),'ZData',P3z(:),'Parent',hax1);
        set(ph_fil{Fn},'LineStyle','-','Color',[.5 .5 .5],'LineWidth',1.5);
    

end
pause(1e-6)
% disp(Sn)
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


PSDactin = PSDact1{1};

ActS = round(linspace(2,as-1,ns));

%%
% keyboard
%%
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            ==========     INITIATE MAIN LOOP     ==========               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nn=1:ns



    % ---------------- GENERATE PARTICLE STEPS ------------------
    [xyz, face] = advance_one_step(...
        xyz, face, k, initial_point, initial_face_index, ...
        all_vertices, triangles, face_local_bases, neighbor_faces);




    % ---------------- PLOT FILAMENTS ------------------
  if mod(nn,modActinPlot)==0

        nnAct1 = ActS(nn);
        nnAct2 = as-ActS(nn);
        nnAct3 = ActS(nn);
        nnAct4 = as-ActS(nn);

        SPYatn = {SPYact1{nnAct1}, SPYact2{nnAct2}, SPYact3{nnAct3}, SPYact4{nnAct4}}; 
        PSDatn = {PSDact1{nnAct1}, PSDact2{nnAct2}, PSDact3{nnAct3}, PSDact4{nnAct4}};
        FILs = {ACTN1{nnAct1}, ACTN2{nnAct2}, ACTN3{nnAct3}, ACTN4{nnAct4}};

    for Fn=1:4

        % update scatter plot of filament tips
        set(ph_spyTip{Fn},'XData',SPYatn{Fn}(:,1)','YData',SPYatn{Fn}(:,2)','ZData',SPYatn{Fn}(:,3)');
        set(ph_psdTip{Fn},'XData',PSDatn{Fn}(:,1)','YData',PSDatn{Fn}(:,2)','ZData',PSDatn{Fn}(:,3)');

        % update line plot of filaments
          P3x = [FILs{Fn}(:,1) FILs{Fn}(:,2)]';
          P3y = [FILs{Fn}(:,3) FILs{Fn}(:,4)]';
          P3z = [FILs{Fn}(:,5) FILs{Fn}(:,6)]';
          P3x(3,:)=NaN; P3y(3,:)=NaN; P3z(3,:)=NaN;
        set(ph_fil{Fn},'XData',P3x(:),'YData',P3y(:),'ZData',P3z(:));
    end
  end



    % 800 2500 4000
    % ----------------- COLLECT DATA -------------------
    if mod(nn, modBindSAP) == 0 %&& nn > Nsteps/2

      k = korig;
      nelA1 = numel(PSDact1{ActS(nn)}(:,1));
      nelA4 = numel(PSDact4{as-ActS(nn)}(:,1));

      for pp = 1:Nparticles

        if xyz(pp,3)>=950 && xyz(pp,1)<=800
        for ti = 1:nelA1
            difdis = PSDact1{nnAct1}(ti,:) - xyz(pp,:); % ActinTips
            radis = sqrt(difdis.^2);
            if radis<=50
              k(pp) = 1;
            end
        end
        end

        if xyz(pp,3)>=950 && xyz(pp,1)>=4000
        for ti = 1:nelA4
            difdis = PSDact4{nnAct4}(ti,:) - xyz(pp,:); % ActinTips
            radis = sqrt(difdis.^2);
            if radis<=50
              k(pp) = 1;
            end
        end
        end

      end
    end

% if nn > 100; keyboard; end




    % ---------------- PLOT SURFACE PARTICLES ------------------ 
    if mod(nn,modAmpaPlot)==0

        set(ampaPlot,'XData',xyz(:,1),'YData',xyz(:,2),'ZData',xyz(:,3));        
            % delete(hax3.Children);
        % scatter3(hax3, xyz(:,1), xyz(:,2), xyz(:,3),...
        %    110,'filled','MarkerEdgeColor','none','MarkerFaceColor',[.95 .05 .05]);
    
    end





if mod(nn,modAmpaPlot)==0 || mod(nn,modAmpaPlot)==0
pause(1e-5); % This pause will render all graphics
% drawnow
% writeVideo(writerObj,getframe);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               ==========     END MAIN LOOP     ==========                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc

% close(writerObj);


keyboard

end






%----------------------------------%
% S1_MainClusterFun
%----------------------------------%
function [S1, SMx] = S1_MainClusterFun(S1,SMx,sap,hkMask,ACTINp1,stepN)

%=========================================================%
%	TEST MASK AGAINST CLUSTER - GET CONVOLUTION MX
%=========================================================%
S=S1;
dT = sap(11);
% hkMask=[0 1 0; 1 0 1; 0 1 0];

Lon = sap(21);	% On Energy (lower = more on events)
Bon = sap(22);	% On Neighbor-Independant Rate (new growth) (lower = more on)
Ron = sap(23);	% On Neighbor-Dependant Rate (cluster fill-in) (higher = more on)

Loff = sap(24);	% Off Energy (higher = more off events)
Boff = sap(25);	% Off Neighbor-Inependant Rate (uniform off)  (lower = more off)
Roff = sap(26);	% Off Neighbor-Dependant Rate (edge off) (higher = more off)



Pmx = rand(size(S));
Soc = (S>0);
Sno = ~Soc;
hk = convn(ACTINp1,hkMask,'same');

Pon = 1 ./ (1+exp((hk-Lon).*(-Bon)));
Pkon = Sno .* ( Ron * dT * Pon );
Son = (Pkon>Pmx);

Poff = 1 ./ (1+exp(((-hk)+Loff).*(-Boff)));
Pkoff = Soc .* ( Roff * dT * Poff );
Soff = (Pkoff>Pmx);

S = (Soc-Soff) + Son;
S1=S;

% S1CP = {hk,Pkon,Sno};


%{
if stepN > 1000; keyboard; end;

figure(10)
imagesc(hk)
colormap('jet');
cbr = colorbar('location','EastOutside');
cax = caxis; caxi = floor(linspace(cax(1),cax(2),8))./10;
set(cbr(1), 'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1],...
	'YTickLabel',{caxi});
set(cbr,'YTickMode','manual');
set(get(cbr,'ylabel'),'string','P_{on}','fontsize',16)
%text(1,10,'h_k','FontSize',20,'Color',[.9 .9 .9]);

figure(11)
imagesc(Pon)
colormap('jet');
cbr = colorbar('location','EastOutside');
set(cbr(1), 'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1]);
text(1,10,'P_{on}','FontSize',20,'Color',[.9 .9 .9]);

figure(12)
imagesc(Pkon)
colormap('jet');
cbr = colorbar('location','EastOutside');
set(cbr(1), 'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1]);
text(1,10,'P_{k,on}','FontSize',20,'Color',[.9 .9 .9]);
%}


end
%----------------------------------%


%----------------------------------%
%			MaskFun
%----------------------------------%
% Creates a kernel mask using a 3D Gaussian
% that can be used for SAP clustering and
% other matrix conv() functions
%----------------------------------%
function [hkMask] = MaskFun(doGMask,MSK)

%-------------------------------%
%		Mask Setup
%-------------------------------%
hkMask=ones(10);
%----------------%
if doGMask
%----------------%

%--------------------
GNpk  = MSK{1};	% hight of peak
GNx0 = MSK{2};	% x-axis peak locations
GNy0 = MSK{3};	% y-axis peak locations
GNsd = MSK{4};	% sigma (stdev of slope)

GNnum = MSK{5};
GNres = MSK{6};
GNspr = ((GNnum-1)*GNres)/2;

a = .5/GNsd^2;
c = .5/GNsd^2;

[X, Y] = meshgrid((-GNspr):(GNres):(GNspr), (-GNspr):(GNres):(GNspr));
Z = GNpk*exp( - (a*(X-GNx0).^2 + c*(Y-GNy0).^2)) ;

%--------------------
hkMask=Z;
% hk = convn(S,hkMask,'same');
% hkor = hk(PSAsz+1,PSAsz+1);
%-----
%LBR(1) = hkor-sqrt(GNpk); LBR(2) = hkor+sqrt(GNpk);
%--------------------


%----------------%
% FIGURE: 3D Gaussian Distribution
fh5 = figure(5); %set(fh5,'OuterPosition',(scsz./[2e-3 2e-3 2 2]))
fh5op = get(fh5,'OuterPosition');
set(fh5,'OuterPosition',[200 200 500 300]);
%----------------%
figure(fh5)
subplot('Position',[.05 .55 .30 .40]); 
	ph5 = imagesc(hkMask); 
	axis equal;
	%set(gca,'XTick',[],'YTick',[])
subplot('Position',[.04 .08 .32 .42]); 
	ph7 = surf(X,Y,Z);
	axis equal; shading interp; view(90,90); 
subplot('Position',[.45 .05 .50 .90]); 
	ph7 = surf(X,Y,Z);
	axis vis3d; shading interp;
	view(-45,30); 
	xlabel('x-axis');ylabel('y-axis');zlabel('z-axis')
%----------------%
end
%-------------------------------%


end
%----------------------------------%



%----------------------------------%
%			plotMx
%----------------------------------%
% plotMx
% short description of this function
%----------------------------------%
function [] = plotMx(Fh,PhA,PhB,PhC,PhD,PhE,G1xy,G2xy,GR1c,GR1r,GR2c,GR2r,SMx)

%figure(Fh)
set(PhA,'CData',G1xy);
set(PhB,'CData',G2xy);
%---
set(PhC,'XData',GR1c,'YData',GR1r);
set(PhD,'XData',GR2c,'YData',GR2r);
%---
set(PhE,'CData',SMx);
%drawnow


% THIS PLOT MAY BE FASTER
%{
%===================================%
%				FIGURE
%-----------------------------------%
%figure(Fh)
subplot('Position',[.62 .04 .34 .93]), 
imagesc(G1xy)
hold on
imagesc(G2xy)
hold on
subplot('Position',[.62 .04 .34 .93]), 
scatter(GR1c,GR1r, 'r')
hold on
scatter(GR2c,GR2r, 'b')
hold off
%-----------------------------------%
%}
end
%----------------------------------%


%----------------------------------%
%			SPLOTS
%----------------------------------%
% plotMx
% short description of this function
%----------------------------------%
function [] = SPLOTS(SPh,S1,S2,S1Tx,S1Ty,S2Tx,S2Ty)

set(SPh.S1Ph1,'CData',S1);
set(SPh.S1Ph2,'XData',S1Tx,'YData',S1Ty);
set(SPh.S2Ph1,'CData',S2);
set(SPh.S2Ph2,'XData',S2Tx,'YData',S2Ty);
drawnow

%--------------------------------
% S1CP = {hk,Pkon,Sno};
% set(S2Ph2,'CData',hk);
% set(S2Ph3,'CData',Pkon);
% set(S2Ph4,'CData',Sno);
% % set(S2Ph5,'CData',G1oc);
%--------------------------------
	
end
%----------------------------------%

