function [varargout] = MAINBOX(dot,dr,um,sap,hr,ko,doUse,doRun,doKo,box,slt,stky,GT,GTab)
clc; format compact; format short; close all; SNSZ = get(0,'ScreenSize');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%						INPUTS FROM GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doIz = doRun(10);
doprofile = doRun(11);
if doprofile; profile on; end

dot1 = dot(1); % GluR1dots
dot2 = dot(2); % GluR2dots
dot3 = dot(3); % Steps
dot4 = dot(4); % TimeStep
dot5 = dot(5); % Scale
dot6 = dot(6); % Loops
 
dr1 = dr(1); % esDGR1
dr2 = dr(2); % spineDGR1
dr3 = dr(3); % PERIDGR1
dr4 = dr(4); % PSDDGR1
dr5 = dr(5); % esDGR2
dr6 = dr(6); % spineDGR2
dr7 = dr(7); % PERIDGR2
dr8 = dr(8); % PSDDGR2

sap1 = sap(1); % SAPdotsPSD1
sap2 = sap(2); % SAPdotsPSD2

slt1 = slt(1); % G1PSDslotN
slt2 = slt(2); % G1PERIslotN
slt3 = slt(3); % G2PSDslotN
slt4 = slt(4); % G2PERIslotN

ko1 = ko(1); % KonSpi1PSDGR2
ko2 = ko(2); % KoffSpi1PSDGR2
ko9 = ko(9); % KonSpi1PSDGR1
ko10 = ko(10); % KoffSpi1PSDGR1

doUse1 = doUse(1); % useGluR1
doUse2 = doUse(2); % useGluR2
 
doRun1 = doRun(1); % run2Dplot
doRun3 = doRun(3); % runMSDtest
doRun5 = doRun(5); % runMSDpopup      % POPUP
doRun7 = doRun(7); % runPoissonsBox

doKo1 = doKo(1); % useGluR1slots
doKo2 = doKo(2); % useGluR2slots

box1 = box(1);	% GraphTime
box2 = box(2);	% AllowedTime
box3 = box(3);	% NumLoops
box4 = box(4);	% SteadyState

esDGR1=dr1;
spineDGR1=dr2;
PERIDGR1=dr3;
PSDDGR1=dr4;
esDGR2=dr5;
spineDGR2=dr6;
PERIDGR2=dr7;
PSDDGR2=dr8;

DontDo = 0;
MSDdrop = doRun5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%						FUNCTION SWITCHES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doMainPlot = doRun1;

Nsteps = dot3;
doGluR1 = doUse1;
doGluR2 = doUse2;

G1N = dot1;		% GluR1 Particles
G2N = dot2;		% GluR2 Particles

if ~doGluR1;G1N=0;end;
if ~doGluR2;G2N=0;end;

trackMSD = doRun3;
if trackMSD; G2N=100;G1N=0;Nsteps=100;end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%				START OUTER LOOP (DATA COLLECTION LOOP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loops = dot6;
for re = 1:loops
%----------------------------------------------------------------------%
	





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%					DIFFUSION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Scale = dot5;				% scale of model
Sc = Scale;					% scale for diffusion

t = dot4/1000;				% time step
d = 2;                      % dimensions
D = esDGR1*t/Sc;			% Diffusion Rate ES (D = L² / 2d*t)
Dp = PSDDGR1*t/Sc;			% Diffusion Rate PSD
Dr = D/Dp;					% Ratio of D:Ds (1/Ls)^2;
Dn = D/Dr;					% new D after scaling L
k = sqrt(d*D);	            % stdev of D's step size distribution
MSD = 2*d*D;                % mean squared displacement
L = sqrt(2*d*D);            % average diagonal (2D) step size
Lx = L/sqrt(2);             % average linear (1D) step size
Ls = 1/sqrt(Dr);			% scales Lx values for Dn

% GLUR DIFFUSION RATES
Ds = esDGR1*t/Sc;			% ExtraSynaptic Model-Scaled Diffusion Rate 
GR1Ds = PSDDGR1*t/Sc;		% GluR1 (PSD) Model-Scaled Diffusion Rate
GR2Ds = PSDDGR2*t/Sc;		% GluR2 (PSD) Model-Scaled Diffusion Rate
Dr_GR1 = Ds/GR1Ds;			% GluR1 (PSD) Ratio Ds:GR1Ds
Dr_GR2 = Ds/GR2Ds;			% GluR2 (PSD) Ratio Ds:GR2Ds
Dn_GR1 = Ds/Dr_GR1;			% GluR1 (PSD) D value after scaling L
Dn_GR2 = Ds/Dr_GR2;			% GluR2 (PSD) D value after scaling L
LsGR1 = 1/sqrt(Dr_GR1);		% GluR1 (PSD) D Lx-Scalar (Ls scales Lx so Dn_GR1 = Ds/Dr_GR1)
LsGR2 = 1/sqrt(Dr_GR2);		% GluR2 (PSD) D Lx-Scalar (Ls scales Lx so Dn_GR2 = Ds/Dr_GR2)

% PSD DIFFUSION RATES (DEFAULT: GLUR1-BASED)
Ds = esDGR1*t/Sc;			% ExtraSynaptic Model-Scaled Diffusion Rate 
PSD1Ds = PSDDGR1*t/Sc;		% (PSD1) Model-Scaled Diffusion Rate
PSD2Ds = PSDDGR1*t/Sc;		% (PSD2) Model-Scaled Diffusion Rate
Dr_PSD1 = Ds/PSD1Ds;		% (PSD1) Ratio Ds:PSD1Ds
Dr_PSD2 = Ds/PSD2Ds;		% (PSD2) Ratio Ds:PSD2Ds
Dn_PSD1 = Ds/Dr_PSD1;		% (PSD1) D value after scaling L
Dn_PSD2 = Ds/Dr_PSD2;		% (PSD2) D value after scaling L
PSD1 = 1/sqrt(Dr_PSD1);		% (PSD1) D Lx-Scalar (Ls scales Lx so Dn_PSD1 = Ds:PSD1Ds)
PSD2 = 1/sqrt(Dr_PSD2);		% (PSD2) D Lx-Scalar (Ls scales Lx so Dn_PSD2 = Ds:PSD2Ds)


% GLUR1 DIFFUSION RATES
DGR1esm = esDGR1*t/Sc;			% GluR1 (ESM) Model-Scaled Diffusion Rate
DGR1spy = spineDGR1*t/Sc;		% GluR1 (ESM) Model-Scaled Diffusion Rate
DGR1psa = PERIDGR1*t/Sc;		% GluR1 (PSA) Model-Scaled Diffusion Rate
DGR1psd = PSDDGR1*t/Sc;			% GluR1 (PSD) Model-Scaled Diffusion Rate
kGR1 = sqrt(d*DGR1esm);			% stdev of D's step size distribution
DrGR1psa = DGR1esm/DGR1psa;		% GluR1 (PSA) Ratio DGR1esm:DGR1psa
DrGR1psd = DGR1esm/DGR1psd;		% GluR1 (PSD) Ratio DGR1esm:DGR1psd
LsGR1psa = 1/sqrt(DrGR1psa);	% GluR1 (PSA) D Lx-Scalar
LsGR1psd = 1/sqrt(DrGR1psd);	% GluR1 (PSA) D Lx-Scalar

% GLUR2 DIFFUSION RATES
DGR2esm = esDGR2*t/Sc;			% GluR2 (ESM) Model-Scaled Diffusion Rate
DGR2spy = spineDGR2*t/Sc;		% GluR2 (ESM) Model-Scaled Diffusion Rate
DGR2psa = PERIDGR2*t/Sc;		% GluR2 (PSA) Model-Scaled Diffusion Rate
DGR2psd = PSDDGR2*t/Sc;			% GluR2 (PSD) Model-Scaled Diffusion Rate
kGR2 = sqrt(d*DGR2esm);			% stdev of D's step size distribution
DrGR2psa = DGR2esm/DGR2psa;		% GluR2 (PSA) Ratio DGR2esm:DGR2psa
DrGR2psd = DGR2esm/DGR2psd;		% GluR2 (PSD) Ratio DGR2esm:DGR2psd
LsGR2psa = 1/sqrt(DrGR2psa);	% GluR2 (PSA) D Lx-Scalar
LsGR2psd = 1/sqrt(DrGR2psd);	% GluR2 (PSA) D Lx-Scalar


% DIFFUSION RATE WHEN SLOTTED
LsGR1S = LsGR1psa / 10;
LsGR2S = LsGR2psa / 1000;
% LsGR1S=1;
% LsGR2S=1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%					Synaptic Dwell-Time Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GluR2_TdwellPSD = zeros(2,G2N);
GluR1_TdwellPSD = zeros(2,G1N);
GluR2_TdwellPERI = zeros(2,G2N);
GluR1_TdwellPERI = zeros(2,G1N);
GluR2_TdwellSPYN = zeros(2,G2N);
GluR1_TdwellSPYN = zeros(2,G1N);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%				PARTICLE BOX EXIT SIMULATION (EXITBOX)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doRun7
pbox = [box1 box2 box3 box4 dot5 dot4 dot1 dot2...
    DGR1spy kGR1 DGR1psd DGR2spy kGR2 DGR2psd LsGR1psd LsGR2psd...
	doKo1 doKo2 ko9 ko10 ko1 ko2...
	sap1 sap2 slt1 slt2 slt3 slt4];
varargin = EXITBOX(pbox,um);
GluR1exT30 = [varargin{1,1}]';
GluR2exT30 = [varargin{1,2}]';
for k=1:2, varargout(k) = {eval(['GluR' int2str(k) 'exT30'])}; end
return
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%				S1 ACTIN MULTIPLEX TIP DATA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MSK = {2.5,0,0,.18,11,.1};
%savetipdir(AMX);
S1TipFile = which('ATdataP1f.mat');
load(S1TipFile);

assignin('base', 'ATs', ATs);
ATs1 = evalin('base', 'ATs');

assignin('base', 'Ax', Ax);
Ax1 = evalin('base', 'Ax');
ActinS1 = Ax1{3};

assignin('base', 'AMx', AMx);
AMx1 = evalin('base', 'AMx');

%------------
ActinTipsP1 = ATs1;
TipCellN = numel(ActinTipsP1);
ACTnP1 = TipCellN - AMx1{7};
NumTipCells = TipCellN - ACTnP1;

TrimActMx = size(ActinTipsP1{1},1)+1;
ACTINp1 = ActinTipsP1{ACTnP1};
ACTINp1(:,TrimActMx:end) = [];
ACTINp1(TrimActMx:end,:) = [];

ActUpdateP1 = AMx1{8};
%------------
% Mask Setup
doGMask = 1;
hkMask = MaskFun(doGMask,AMx1,MSK);

AMask=ones(AMx1{5});
Ahk = convn(ACTINp1,AMask,'same');
S = (Ahk>0).*1.0;

%------------
S1=S;
[S1Ty,S1Tx] = find(ACTINp1);
%------------


%ActinPlot(Ax1{1},Ax1{2},Ax1{3},Ax1{4},Ax1{5},Ax1{6},Ax1{7});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%				S2 ACTIN MULTIPLEX TIP DATA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%savetipdir(AMX);
S2TipFile = which('ATdataP1e.mat');
load(S2TipFile);

assignin('base', 'ATs', ATs);
ATs2 = evalin('base', 'ATs');

assignin('base', 'Ax', Ax);
Ax2 = evalin('base', 'Ax');
ActinS2 = Ax1{3};

assignin('base', 'AMx', AMx);
AMx2 = evalin('base', 'AMx');

%------------
ActinTipsP2 = ATs2;
TipCellN = numel(ActinTipsP2);
ACTnP2 = TipCellN - AMx2{7};
NumTipCells = TipCellN - ACTnP2;

TrimActMx = size(ActinTipsP2{1},1)+1;
ACTINp2 = ActinTipsP2{ACTnP2};
ACTINp2(:,TrimActMx:end) = [];
ACTINp2(TrimActMx:end,:) = [];

ActUpdateP2 = AMx2{8};


%------------
% Mask Setup
doGMask = 1;
hkMask = MaskFun(doGMask,AMx2,MSK);

AMask=ones(AMx2{5});
Ahk = convn(ACTINp2,AMask,'same');
S = (Ahk>0).*1.0;

%------------
S2=S;
[S2Ty,S2Tx] = find(ACTINp2);	
%------------


% ActinPlot(Ax2{1},Ax2{2},Ax2{3},Ax2{4},Ax2{5},Ax2{6},Ax2{7})



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%				FIELD SYNAPTIC MATRIX SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[DFszX,DFszY,PSD1sz,PSD2sz,PSA1sz,PSA2sz,SPY1sz,SPY2sz,...
S1sz,S2sz,SSsz,SMx,DFum,SS,DF] = FieldFun(S1,S2,um,Scale);

S1C = SS.S1C;
S2C = SS.S2C;
S1rad = SS.S1rad;
S2rad = SS.S2rad;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%				AMPAR DOT PARTICLE MATRIX SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%GLUR1
xyG1 = ones(2,G1N);
xysG1 = ones(2,G1N);

% GLUR2
xyG2 = ones(2,G2N);
xysG2 = ones(2,G2N);

% Position Dots in the middle of DF
xyG1(1,:) = xyG1(1,:) + DF.X/2;
xyG2(1,:) = xyG2(1,:) + DF.X/2;
xyG1(2,:) = xyG1(2,:) + DF.Y/2;
xyG2(2,:) = xyG2(2,:) + DF.Y/2;

% Rounded version of GluRxyl
GR1xy = round(xyG1);
GR2xy = round(xyG2);


% Radius from synapse counters
S1G1C = repmat(S1C,1,numel(xyG1(1,:)));
S2G1C = repmat(S2C,1,numel(xyG1(1,:)));
S1G2C = repmat(S1C,1,numel(xyG2(1,:)));
S2G2C = repmat(S2C,1,numel(xyG2(1,:)));


S1G1Cv = xyG1-S1G1C;
S2G1Cv = xyG1-S2G1C;
S1G2Cv = xyG2-S1G2C;
S2G2Cv = xyG2-S2G2C;

G1P1r = sqrt(sum((S1G1Cv).^2));
G1P2r = sqrt(sum((S2G1Cv).^2));
G2P1r = sqrt(sum((S1G2Cv).^2));
G2P2r = sqrt(sum((S2G2Cv).^2));

G1inP1=G1P1r<S1rad;
G1inP2=G1P2r<S2rad;
G2inP1=G2P1r<S1rad;
G2inP2=G2P2r<S2rad;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%				PARTICLE SPACE TO MATRIX SPACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	GR1x = xyG1(1,:);
GR1y = xyG1(2,:);
GR1c = round(GR1x);
GR1r = round(GR1y);
GR1c(GR1c>DF.X) = DF.X;		% any dot somehow outside, place inside
GR1r(GR1r>DF.Y) = DF.Y;		% any dot somehow outside, place inside
GR1c(GR1c<1) = 1;			% any dot somehow outside, place inside
GR1r(GR1r<1) = 1;			% any dot somehow outside, place inside
	GR1xy(1,:) = GR1c;
	GR1xy(2,:) = GR1r;

G1xy = DF.SMxE;
for xy = 1:numel(GR1c)
G1xy(GR1r(xy),GR1c(xy)) = 1;
end

%--------
GR2x = xyG2(1,:);
GR2y = xyG2(2,:);
GR2c = round(GR2x);
GR2r = round(GR2y);
GR2c(GR2c>DF.X) = DF.X;		% any dot somehow outside, place inside
GR2r(GR2r>DF.Y) = DF.Y;		% any dot somehow outside, place inside
GR2c(GR2c<1) = 1;			% any dot somehow outside, place inside
GR2r(GR2r<1) = 1;			% any dot somehow outside, place inside
	GR2xy(1,:) = GR2c;
	GR2xy(2,:) = GR2r;

G2xy = DF.SMxE;
for xy = 1:numel(GR2c)
G2xy(GR2r(xy),GR2c(xy)) = 1;
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%				TIME AND PARTICLE COUNT VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATARATE = 10;
AveOver = 10;
SaveSteps = Nsteps/DATARATE;
dataset = zeros(SaveSteps,10);
SAPdata = zeros(SaveSteps,2);
Ddata = zeros(SaveSteps,2);
AMPARdata = zeros(SaveSteps,2);
GluRdata = zeros(SaveSteps,19);
G1SLOTSDATA = zeros(SaveSteps,4);
G2SLOTSDATA = zeros(SaveSteps,4);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%				FIGURES - PLOTS - LIVE GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% fig1 = figure(1);
% fig1box(fig1);
% datafig = figure(88);
% figure('Position',[1 SCREENSIZE(4)/8 SCREENSIZE(3)/8 SCREENSIZE(4)/8])
% set(fig1,'Renderer','painters')
% set(fig1,'DoubleBuffer','off')
% set(fig1,'WindowStyle','modal')
%}
%{
% scsz = get(0,'ScreenSize'); 
% scsx=scsz(3); scsy=scsz(4);
% xlim = [0 DF.X];
% ylim = [0 DF.Y];
% 
% Fh1 = figure(1);
% set(gcf,'OuterPosition',[145,145,scsx/1.2,scsy/1.2])
% set(gcf,'Color',[1,1,1])
% 
% %================================================%
% %			Real-Space Scatter Plot
% %------------------------------------------------%
% 
% axes('Position',[.32 .25 .28 .72]);
% 
% G1Ph1 = scatter(xyG1(1,:),xyG1(2,:),10,[.9 .08 .24],'filled');
% axis([xlim, ylim]);
% %set(gca,'Color',[1,1,1])
% hold on;
% 
% G2Ph1 = scatter(xyG2(1,:),xyG2(2,:),10,[.1 .9 .1],'filled');
% %set(gca,'YDir','reverse')
% hold on;
% 
% set(gca,'xticklabel',[])
% set(gca,'yticklabel',[])
% %set(get(gca,'XLabel'),'String','GluR2:Blue \bullet \bullet GluR1:Red')
% rectangle('Position',[DF.S1L,DF.S1B,DF.S1sz,DF.S1sz],'Curvature',[1,1])
% rectangle('Position',[DF.S2L,DF.S2B,DF.S2sz,DF.S2sz],'Curvature',[1,1])
% hold on;
% 
% 
% %================================================%
% %			Actin Network Plots
% %------------------------------------------------%
% % ActinMainPlot(Fh,nT,Actin,inPSD,varargin(rot,azel,dims))
% inPSD1 = Ax1{7}(7); 
% spos1 = [.30 -.01 .18 .25];
% ActinMainPlot(Fh1,ActinS1,inPSD1,spos1)
% 
% inPSD2 = Ax2{7}(7);
% spos2 = [.48 -.01 .18 .25];
% ActinMainPlot(Fh1,ActinS2,inPSD2,spos2)
% 
% 
% %================================================%
% %			S Cluster Plots
% %------------------------------------------------%
% 
% % clrmap = 'cool';
% clrmap = [1 1 1; .1 .1 .5];
% 
% %-- S1 Cluster Plots --
% axes('Position',[.02 .53 .27 .44]);
% S1Ph1 = imagesc(S1);
% colormap(clrmap);
% title('SYNAPSE-1 CLUSTERS');
% hold on
% %--------
% [PSDY,PSDX] = find(ACTINp1);
% S1Ph2 = scatter(PSDX,PSDY,20, [.9 .08 .24],'filled');
% hold off
% %--------
% 
% 
% %-- S2 Cluster Plots --
% axes('Position',[.02 .04 .27 .44]);
% S2Ph1 = imagesc(S2);
% colormap(clrmap);
% title('SYNAPSE-2 CLUSTERS');
% hold on
% %--------
% [PSDY,PSDX] = find(ACTINp2);
% S2Ph2 = scatter(PSDX,PSDY,20, [.9 .08 .24],'filled');
% hold off
% %--------
% 
% 
% %================================================%
% %			PARTICLE MATRIX PLOT
% %------------------------------------------------%
% axes('Position',[.65 .05 .34 .92]);
% %--------
% Ph301 = imagesc(G1xy);
% set(gca,'YDir','normal')
% %colormap('bone')
% colormap(clrmap)
% hold on;
% %---
% Ph302 = imagesc(G2xy);
% set(gca,'YDir','normal')
% colormap(clrmap)
% hold on;
% %---
% Ph305 = imagesc(SMx);
% set(gca,'YDir','normal')
% colormap(clrmap)
% hold on;
% %---
% Ph303 = scatter(xyG1(1,:),xyG1(2,:),20, [.9 .08 .24],'filled');
% hold on
% %---
% Ph304 = scatter(xyG2(1,:),xyG2(2,:),20, [.1 .9 .1],'filled');
% hold off
% 
% %------------------------------------------------%
% rectangle('Position',[DF.S1L,DF.S1B,DF.S1sz,DF.S1sz],'Curvature',[1,1])
% rectangle('Position',[DF.S2L,DF.S2B,DF.S2sz,DF.S2sz],'Curvature',[1,1])
% hold on;
% 
% %----------------------%
% % SAVE PLOT HANDLES
% %----------------------%
% SPh = struct;
% SPh.S1Ph1 = S1Ph1;
% SPh.S1Ph2 = S1Ph2;
% SPh.S2Ph1 = S2Ph1;
% SPh.S2Ph2 = S2Ph2;
% %----------------------%
%{
% %-- S CLUSTER PLOTS --
% SPLXi1=.03;
% SPLXi2=.20;
% SPLYi1=.192;
% SPLYi2=.192;
% SPLHi=.14;
% SPLWi=.10;
% 
% 
% %-- S1 Additional Cluster Plots --
% % subplot(5,5,1),
% axes('Position',[(SPLXi1) (1-SPLYi1*1) SPLWi SPLHi]);
% S1Ph1 = imagesc(S1);
% colormap('hot');
% title('SYNAPSE-1 CLUSTERS');
% 
% % subplot(5,5,6),
% axes('Position',[(SPLXi1) (1-SPLYi1*2) SPLWi SPLHi]);
% S1Ph2 = imagesc(S1);
% colormap('hot');
% title('hk');
% 
% % subplot(5,5,11),
% axes('Position',[(SPLXi1) (1-SPLYi1*3) SPLWi SPLHi]);
% S1Ph3 = imagesc(S1);
% colormap('hot');
% title('Pkon');
% 
% % subplot(5,5,16),
% axes('Position',[(SPLXi1) (1-SPLYi1*4) SPLWi SPLHi]);
% S1Ph4 = imagesc(S1);
% colormap('hot');
% title('Sno');
% 
% % subplot(5,5,21),
% axes('Position',[(SPLXi1) (1-SPLYi1*5) SPLWi SPLHi]);
% S1Ph5 = imagesc(S1);
% colormap('hot');
% title('TETHERED G1');
% 
% 
% 
% 
% %-- S2 Additional Cluster Plots --
% % subplot(5,5,2); 
% axes('Position',[(SPLXi2) (1-SPLYi2*1) SPLWi SPLHi]);
% S2Ph1 = imagesc(S2);
% colormap('hot');
% title('SYNAPSE-2 CLUSTERS');
% 
% 
% % subplot(5,5,7),
% axes('Position',[(SPLXi2) (1-SPLYi2*2) SPLWi SPLHi]);
% S2Ph2 = imagesc(S2);
% colormap('hot');
% title('hk');
% 
% 
% % subplot(5,5,12),
% axes('Position',[(SPLXi2) (1-SPLYi2*3) SPLWi SPLHi]);
% S2Ph3 = imagesc(S2);
% colormap('hot');
% title('Pkon');
% 
% 
% % subplot(5,5,17),
% axes('Position',[(SPLXi2) (1-SPLYi2*4) SPLWi SPLHi]);
% S2Ph4 = imagesc(S2);
% colormap('hot');
% title('Sno');
% 
% 
% % subplot(5,5,22),
% axes('Position',[(SPLXi2) (1-SPLYi2*5) SPLWi SPLHi]);
% S2Ph5 = imagesc(S2);
% colormap('hot');
% title('TETHERED G1');
%}
%------------------------------------------------%
%%
%}
%{
%================================================%
%			Real-Space Scatter Plot
%------------------------------------------------%
% axes('Position',[.32 .25 .28 .72]);
% 
% G1Ph1 = scatter(xyG1(1,:),xyG1(2,:),10,[.9 .08 .24],'filled');
% axis([xlim, ylim]);
% %set(gca,'Color',[1,1,1])
% hold on;
% 
% G2Ph1 = scatter(xyG2(1,:),xyG2(2,:),10,[.1 .9 .1],'filled');
% %set(gca,'YDir','reverse')
% hold on;
% 
% set(gca,'xticklabel',[])
% set(gca,'yticklabel',[])
% %set(get(gca,'XLabel'),'String','GluR2:Blue \bullet \bullet GluR1:Red')
% rectangle('Position',[DF.S1L,DF.S1B,DF.S1sz,DF.S1sz],'Curvature',[1,1])
% rectangle('Position',[DF.S2L,DF.S2B,DF.S2sz,DF.S2sz],'Curvature',[1,1])
% hold on;
%}

%-----------------------------------%
%	PLOTS: SETUP
%-----------------------------------%
scsz = get(0,'ScreenSize'); 
scsx=scsz(3); scsy=scsz(4);
xlim = [0 DF.X];
ylim = [0 DF.Y];

Fh1 = figure(1);
set(gcf,'OuterPosition',[145,145,scsx/1.2,scsy/1.2])
set(gcf,'Color',[1,1,1])



%-----------------------------------%
%	PLOTS: ACTIN NETWORK
%-----------------------------------%
% ActinMainPlot(Fh,nT,Actin,inPSD,varargin(rot,azel,dims))
inPSD1 = Ax1{7}(7); 
spos1 = [.02 .50 .26 .47];
ActinMainPlot(Fh1,ActinS1,inPSD1,spos1)

inPSD2 = Ax2{7}(7);
spos2 = [.02 .02 .26 .47];
ActinMainPlot(Fh1,ActinS2,inPSD2,spos2)


%-----------------------------------%
%	PLOTS: CLUSTERS (S1 & S2)
%-----------------------------------%
clrmap = [1 1 1; .55 .55 .55];
%-- S1 Cluster Plots --
axes('Position',[.31 .53 .27 .44]);
S1Ph1 = imagesc(S1);
colormap(clrmap);
title('Synapse-A Scaffold Clusters');
hold on
%--------
[PSDY,PSDX] = find(ACTINp1);
S1Ph2 = scatter(PSDX,PSDY);
set(S1Ph2,'Marker','s','SizeData',60,'LineWidth',.5,...
	'MarkerFaceColor',[.95 .1 .1],'MarkerEdgeColor','none')
hold off


%-- S2 Cluster Plots --
axes('Position',[.31 .04 .27 .44]);
S2Ph1 = imagesc(S2);
colormap(clrmap);
title('Synapse-B Scaffold Clusters');
hold on
%--------
[PSDY,PSDX] = find(ACTINp2);
S2Ph2 = scatter(PSDX,PSDY);
set(S2Ph2,'Marker','s','SizeData',60,'LineWidth',.5,...
	'MarkerFaceColor',[.95 .1 .1],'MarkerEdgeColor','none')
hold off
%--------


%-----------------------------------%
%	PLOTS: PARTICLE MATRIX
%-----------------------------------%
axes('Position',[.62 .04 .34 .93]);
%--------
Ph301 = imagesc(G1xy);
set(gca,'YDir','normal')
%colormap('bone')
colormap(clrmap)
hold on;
%---
Ph302 = imagesc(G2xy);
set(gca,'YDir','normal')
colormap(clrmap)
hold on;
%---
Ph305 = imagesc(SMx);
set(gca,'YDir','normal')
colormap(clrmap)
hold on;
%---
Ph303 = scatter(xyG1(1,:),xyG1(2,:),50, [.0 .5 .8],'filled');
hold on
%---
%xyG2 = xyG1.*rand(2,150)
Ph304 = scatter(xyG2(1,:),xyG2(2,:),50, [.2 .7 .0],'filled');
set(Ph304,'Marker','s')
hold off

%--------
rectangle('Position',[DF.S1L,DF.S1B,DF.S1sz,DF.S1sz],'Curvature',[1,1])
rectangle('Position',[DF.S2L,DF.S2B,DF.S2sz,DF.S2sz],'Curvature',[1,1])
rectangle('Position',[DF.PS1L,DF.PS1B,DF.PS1sz,DF.PS1sz],'Curvature',[1,1])
rectangle('Position',[DF.PS2L,DF.PS2B,DF.PS2sz,DF.PS2sz],'Curvature',[1,1])
hold on;
%--------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%				PLOT HANDLES AND UPDATE RATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SPh = struct;
SPh.S1Ph1 = S1Ph1;
SPh.S1Ph2 = S1Ph2;
SPh.S2Ph1 = S2Ph1;
SPh.S2Ph2 = S2Ph2;
%----
SPLOTr=100;
MPLOTr=200;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%					START Main Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for stepN = 1:Nsteps 
%------------------------------------------------------------------------------%



	%------------------------------------------%
    %		STEP DIRECTION & SIZE
    %------------------------------------------%
	if doGluR1; xysG1 = STEPxyds(G1N, kGR1); end
	if doGluR2; xysG2 = STEPxyds(G2N, kGR2); end
    
  
	
	%------------------------------------------%
    %		SAP MATRIX UPDATE
    %------------------------------------------%
	% THIS SHIT SEEMS BACKWARDS BUT IT WORKS
	SMx(SS.S1B:SS.S1T,SS.S1L:SS.S1R) = S2;
	SMx(SS.S2B:SS.S2T,SS.S2L:SS.S2R) = S1;
	SS.SMx = SMx;
	
	
	
	%------------------------------------------%
    %		MAKE PARTICLE MATRIX
	%------------------------------------------%
	GR1x = xyG1(1,:);
	GR1y = xyG1(2,:);
	GR1c = round(GR1x);
	GR1r = round(GR1y);
	GR1c(GR1c>DF.X) = DF.X;		% any dot somehow outside, place inside
	GR1r(GR1r>DF.Y) = DF.Y;		% any dot somehow outside, place inside
	GR1c(GR1c<1) = 1;			% any dot somehow outside, place inside
	GR1r(GR1r<1) = 1;			% any dot somehow outside, place inside
	GR1xy(1,:) = GR1c;
	GR1xy(2,:) = GR1r;
	%---
	G1xy = DF.SMxE;
	for xy = 1:numel(GR1c)
	G1xy(GR1r(xy),GR1c(xy)) = 1;
	end
	%--------
	GR2x = xyG2(1,:);
	GR2y = xyG2(2,:);
	GR2c = round(GR2x);
	GR2r = round(GR2y);
	GR2c(GR2c>DF.X) = DF.X;		% any dot somehow outside, place inside
	GR2r(GR2r>DF.Y) = DF.Y;		% any dot somehow outside, place inside
	GR2c(GR2c<1) = 1;			% any dot somehow outside, place inside
	GR2r(GR2r<1) = 1;			% any dot somehow outside, place inside
	GR2xy(1,:) = GR2c;
	GR2xy(2,:) = GR2r;
	%---
	G2xy = DF.SMxE;
	for xy = 1:numel(GR2c)
	G2xy(GR2r(xy),GR2c(xy)) = 1;
	end
	
	

	
	
	%------------------------------------------%
    %		SUPERSLOT FUNCTION
    %------------------------------------------%
	[xyG1 xyG2 xysG1 xysG2...
	G1FSLOTS G2FSLOTS...
	GluR1_TdwellSPYN GluR1_TdwellPSD GluR1_TdwellPERI...
	GluR2_TdwellSPYN GluR2_TdwellPSD GluR2_TdwellPERI...
	G1INSPA G2INSPA SMx G1v G2v...
	G1P1r G1P2r G2P1r G2P2r G1inP1 G1inP2 G2inP1 G2inP2...
	S1G1Cv S2G1Cv S1G2Cv S2G2Cv]...
	= SUPERSLOT(stepN,t,xyG1,xyG2,xysG1,xysG2,...
	GluR1_TdwellPSD,GluR1_TdwellPERI,GluR1_TdwellSPYN,LsGR1psa,LsGR1psd,LsGR1S,...
	GluR2_TdwellPSD,GluR2_TdwellPERI,GluR2_TdwellSPYN,LsGR2psa,LsGR2psd,LsGR2S,...	
	SSsz,G1xy,G2xy,SMx,SS,DF,GR1xy,GR2xy,...
	S1G1C,S2G1C,S1G2C,S2G2C,S1rad,S2rad);

	

    %------------------------------------------%
    %		SAP CLUSTERS
    %------------------------------------------%
	if mod(stepN,ActUpdateP1) == 0
	ACTINp1 = ActinTipsP1{ACTnP1};
	ACTINp1 = convn(ACTINp1,AMask,'same');
	[S1Ty,S1Tx] = find(ACTINp1);
	
	ACTnP1 = ACTnP1+1;
	end
	
	if mod(stepN,ActUpdateP2) == 0
	ACTINp2 = ActinTipsP2{ACTnP2};
	ACTINp2 = convn(ACTINp2,AMask,'same');
	[S2Ty,S2Tx] = find(ACTINp2);
	
	ACTnP2 = ACTnP2+1;
	end

	
	
	
	%------------------------------------------%
	%		S1 & S2  MainClusterFun
	%------------------------------------------%
	[S1 SMx] = S1_MainClusterFun(S1,SMx,sap,hkMask,ACTINp1,stepN);
	[S2 SMx] = S2_MainClusterFun(S2,SMx,sap,hkMask,ACTINp2,stepN);
	
	S1sum = sum(S1(:))/2;
	S2sum = sum(S2(:))/2;
		
	
	
	
	%------------------------------------------%
    %		MOVE PARTICLES MAIN FUNCTION
	%------------------------------------------%
	[xysG1 xyG1 xysG2 xyG2] = MOVEGLUR(stepN,DF,G1N,xysG1,xyG1, G2N,xysG2,xyG2);
	

	
	
	
	%------------------------------------------%
    %		CLUSTER PLOT
    %------------------------------------------%
	if doMainPlot;if mod(stepN,SPLOTr)==0;
		SPLOTS(SPh,S1,S2,S1Tx,S1Ty,S2Tx,S2Ty);
	end;end;
	
    %------------------------------------------%
    %		MAIN PLOT XYL SIMULATION
    %------------------------------------------%
	if doMainPlot;if mod(stepN,MPLOTr)==0;
	plotMx(Fh1,Ph301,Ph302,Ph303,Ph304,Ph305,G1xy,G2xy,GR1c,GR1r,GR2c,GR2r,SMx);
	end;end;


 	
 


	%==========================================================%
	%				START PARTICLE COUNTERS
	%----------------------------------------------------------%
	if mod(stepN, DATARATE) == 0
	  
    %-------------------------------%
    %       PARTICLE COUNTERS
    %-------------------------------%
	[GR1SPY1N GR1SPY2N GR1PER1N GR1PER2N GR1PSD1N GR1PSD2N...
		GR2SPY1N GR2SPY2N GR2PER1N GR2PER2N GR2PSD1N GR2PSD2N...
		GR1PERN GR1PSDN GR2PERN GR2PSDN GR1SPYN GR2SPYN...
		SPY1N SPY2N SPYNT ESNT EST G1INSP G2INSP GINSP...
		G1xyP1 G1xyP2 G2xyP1 G2xyP2]...
		= dotCount(stepN,G2N, G1N,...
		G1P1r, G1P2r, G2P1r, G2P2r, G1inP1, G1inP2, G2inP1, G2inP2,...
		S1G1Cv,S2G1Cv,S1G2Cv,S2G2Cv,G1INSPA,G2INSPA,xyG1,xyG2);
	
	%-------------------------------%
    %	SAVE DATA
    %-------------------------------%
       
       data = [stepN/Nsteps ESNT SPY1N SPY2N SPYNT;... 
			   (Nsteps-stepN) D/10 Dn_PSD1/10 Dn_PSD2/10 k];
       dataset(stepN/DATARATE,:) = [data(1,:) data(2,:)];
	   SAPdata(stepN/DATARATE,:) = [S1sum S2sum];
	   Ddata(stepN/DATARATE,:) = [Dn_PSD1/10 Dn_PSD2/10];
	   AMPARdata(stepN/DATARATE,:) = [GR1SPYN GR2SPYN];
	   
	   
	   GluRdata(stepN/DATARATE,:) = [GR1PSDN GR1PERN GR2PSDN GR2PERN...
								GR1SPY1N GR1SPY2N GR2SPY1N GR2SPY2N...
								G1INSP G2INSP GINSP...
								GR1PER1N GR1PER2N GR1PSD1N GR1PSD2N...
								GR2PER1N GR2PER2N GR2PSD1N GR2PSD2N];
							
							
	   G1SLOTSDATA(stepN/DATARATE,:) = G1FSLOTS;
	   G2SLOTSDATA(stepN/DATARATE,:) = G2FSLOTS; 
	   
	   SGr{stepN/DATARATE}={G1P1r,G1P2r,G2P1r,G2P2r,G1inP1,G1inP2,G2inP1,G2inP2,...
							G1xyP1 G1xyP2 G2xyP1 G2xyP2};
	   
	   if mod(stepN, DATARATE*10) == 0
	   head1 = ('Time%2go    CntDwn');
	   disp(head1); disp([data(1,1) data(2,1)]); %snapnow;
	   end
       
	end
	%----------------------------------------------------------%
	%				END PARTICLE COUNTERS
	%==========================================================%
	

 
 
 
 
 
 
 
 

%------------------------------------------------------------------------------%
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%							END SIMULATION LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doprofile; profile viewer; end


DATAdataset(:,1) = mean(reshape(dataset(:,1), [AveOver, length(dataset(:,1))/AveOver])); 
DATAdataset(:,2) = mean(reshape(dataset(:,2), [AveOver, length(dataset(:,2))/AveOver])); 
DATAdataset(:,3) = mean(reshape(dataset(:,3), [AveOver, length(dataset(:,3))/AveOver])); 
DATAdataset(:,4) = mean(reshape(dataset(:,4), [AveOver, length(dataset(:,4))/AveOver])); 
DATAdataset(:,5) = mean(reshape(dataset(:,5), [AveOver, length(dataset(:,5))/AveOver])); 
DATAdataset(:,6) = mean(reshape(dataset(:,6), [AveOver, length(dataset(:,6))/AveOver])); 
DATAdataset(:,7) = mean(reshape(dataset(:,7), [AveOver, length(dataset(:,7))/AveOver])); 
DATAdataset(:,8) = mean(reshape(dataset(:,8), [AveOver, length(dataset(:,8))/AveOver])); 
DATAdataset(:,9) = mean(reshape(dataset(:,9), [AveOver, length(dataset(:,9))/AveOver])); 
DATAdataset(:,10) = mean(reshape(dataset(:,10), [AveOver, length(dataset(:,10))/AveOver])); 


DATASAPdata(:,1) = mean(reshape(SAPdata(:,1), [AveOver, length(SAPdata(:,1))/AveOver])); 
DATASAPdata(:,2) = mean(reshape(SAPdata(:,2), [AveOver, length(SAPdata(:,2))/AveOver])); 

DATADdata(:,1) = mean(reshape(Ddata(:,1), [AveOver, length(Ddata(:,1))/AveOver])); 
DATADdata(:,2) = mean(reshape(Ddata(:,2), [AveOver, length(Ddata(:,2))/AveOver])); 

DATAAMPARdata(:,1) = mean(reshape(AMPARdata(:,1), [AveOver, length(AMPARdata(:,1))/AveOver])); 
DATAAMPARdata(:,2) = mean(reshape(AMPARdata(:,2), [AveOver, length(AMPARdata(:,2))/AveOver])); 

DATAGluRdata(:,1) = mean(reshape(GluRdata(:,1), [AveOver, length(GluRdata(:,1))/AveOver])); 
DATAGluRdata(:,2) = mean(reshape(GluRdata(:,2), [AveOver, length(GluRdata(:,2))/AveOver])); 
DATAGluRdata(:,3) = mean(reshape(GluRdata(:,3), [AveOver, length(GluRdata(:,3))/AveOver])); 
DATAGluRdata(:,4) = mean(reshape(GluRdata(:,4), [AveOver, length(GluRdata(:,4))/AveOver])); 
DATAGluRdata(:,5) = mean(reshape(GluRdata(:,5), [AveOver, length(GluRdata(:,5))/AveOver])); 
DATAGluRdata(:,6) = mean(reshape(GluRdata(:,6), [AveOver, length(GluRdata(:,6))/AveOver])); 
DATAGluRdata(:,7) = mean(reshape(GluRdata(:,7), [AveOver, length(GluRdata(:,7))/AveOver])); 
DATAGluRdata(:,8) = mean(reshape(GluRdata(:,8), [AveOver, length(GluRdata(:,8))/AveOver])); 
DATAGluRdata(:,9) = mean(reshape(GluRdata(:,9), [AveOver, length(GluRdata(:,9))/AveOver])); 
DATAGluRdata(:,10) = mean(reshape(GluRdata(:,10), [AveOver, length(GluRdata(:,10))/AveOver])); 
DATAGluRdata(:,11) = mean(reshape(GluRdata(:,11), [AveOver, length(GluRdata(:,11))/AveOver])); 
DATAGluRdata(:,12) = mean(reshape(GluRdata(:,12), [AveOver, length(GluRdata(:,12))/AveOver])); 
DATAGluRdata(:,13) = mean(reshape(GluRdata(:,13), [AveOver, length(GluRdata(:,13))/AveOver])); 
DATAGluRdata(:,14) = mean(reshape(GluRdata(:,14), [AveOver, length(GluRdata(:,14))/AveOver])); 
DATAGluRdata(:,15) = mean(reshape(GluRdata(:,15), [AveOver, length(GluRdata(:,15))/AveOver])); 
DATAGluRdata(:,16) = mean(reshape(GluRdata(:,16), [AveOver, length(GluRdata(:,16))/AveOver])); 
DATAGluRdata(:,17) = mean(reshape(GluRdata(:,17), [AveOver, length(GluRdata(:,17))/AveOver])); 
DATAGluRdata(:,18) = mean(reshape(GluRdata(:,18), [AveOver, length(GluRdata(:,18))/AveOver])); 
DATAGluRdata(:,19) = mean(reshape(GluRdata(:,19), [AveOver, length(GluRdata(:,19))/AveOver])); 


DATAG1SLOTSDATA(:,1) = mean(reshape(G1SLOTSDATA(:,1), [AveOver, length(G1SLOTSDATA(:,1))/AveOver])); 
DATAG1SLOTSDATA(:,2) = mean(reshape(G1SLOTSDATA(:,2), [AveOver, length(G1SLOTSDATA(:,2))/AveOver])); 
DATAG1SLOTSDATA(:,3) = mean(reshape(G1SLOTSDATA(:,3), [AveOver, length(G1SLOTSDATA(:,3))/AveOver])); 
DATAG1SLOTSDATA(:,4) = mean(reshape(G1SLOTSDATA(:,4), [AveOver, length(G1SLOTSDATA(:,4))/AveOver])); 

DATAG2SLOTSDATA(:,1) = mean(reshape(G2SLOTSDATA(:,1), [AveOver, length(G2SLOTSDATA(:,1))/AveOver])); 
DATAG2SLOTSDATA(:,2) = mean(reshape(G2SLOTSDATA(:,2), [AveOver, length(G2SLOTSDATA(:,2))/AveOver])); 
DATAG2SLOTSDATA(:,3) = mean(reshape(G2SLOTSDATA(:,3), [AveOver, length(G2SLOTSDATA(:,3))/AveOver])); 
DATAG2SLOTSDATA(:,4) = mean(reshape(G2SLOTSDATA(:,4), [AveOver, length(G2SLOTSDATA(:,4))/AveOver])); 


reDATAdataset{re} = DATAdataset;
reDATASAPdata{re} = DATASAPdata;
reDATADdata{re} = DATADdata;
reDATAAMPARdata{re} = DATAAMPARdata;
reDATAGluRdata{re} = DATAGluRdata;
reDATAG1SLOTSDATA{re} = DATAG1SLOTSDATA;
reDATAG2SLOTSDATA{re} = DATAG2SLOTSDATA;



reSGr{re} = SGr;





%==============================================================================%
end % for re = 1:loops
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%								END DATA LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%=============================================%
if loops>1
%=============================================%
redDATAdataset = 0;
redDATASAPdata  = 0;
redDATADdata  = 0;
redDATAAMPARdata  = 0;
redDATAGluRdata  = 0;
redDATAG1SLOTSDATA  = 0;
redDATAG2SLOTSDATA  = 0;

for red = 1:loops
redDATAdataset = redDATAdataset + reDATAdataset{red};
redDATASAPdata  = redDATASAPdata + reDATASAPdata{red};
redDATADdata  = redDATADdata + reDATADdata{red};
redDATAAMPARdata  = redDATAAMPARdata + reDATAAMPARdata{red};
redDATAGluRdata  = redDATAGluRdata + reDATAGluRdata{red};
redDATAG1SLOTSDATA  = redDATAG1SLOTSDATA + reDATAG1SLOTSDATA{red};
redDATAG2SLOTSDATA  = redDATAG2SLOTSDATA + reDATAG2SLOTSDATA{red};
end

redoDATAdataset = redDATAdataset ./loops;
redoDATASAPdata  = redDATASAPdata ./loops;
redoDATADdata  = redDATADdata ./loops;
redoDATAAMPARdata  = redDATAAMPARdata ./loops;
redoDATAGluRdata  = redDATAGluRdata ./loops;
redoDATAG1SLOTSDATA  = redDATAG1SLOTSDATA ./loops;
redoDATAG2SLOTSDATA  = redDATAG2SLOTSDATA ./loops;

DATAdataset = redoDATAdataset;
DATASAPdata = redoDATASAPdata;
DATADdata  = redoDATADdata;
DATAAMPARdata  = redoDATAAMPARdata;
DATAGluRdata  = redoDATAGluRdata;
DATAG1SLOTSDATA  = redoDATAG1SLOTSDATA;
DATAG2SLOTSDATA  = redoDATAG2SLOTSDATA;
%=============================================%
end %if loops>1
%=============================================%



%%
reDATAdatasetFull = reDATAdataset;
reDATASAPdataFull = reDATASAPdata;
reDATADdataFull = reDATADdata;
reDATAAMPARdataFull = reDATAAMPARdata;
reDATAGluRdataFull = reDATAGluRdata;
reDATAG1SLOTSDATAFull = reDATAG1SLOTSDATA;
reDATAG2SLOTSDATAFull = reDATAG2SLOTSDATA;

assignin('base', 'BMData1', reDATAdatasetFull)
assignin('base', 'BMData2', reDATASAPdataFull)
assignin('base', 'BMData3', reDATADdataFull)
assignin('base', 'BMData4', reDATAAMPARdataFull)
assignin('base', 'BMData5', reDATAGluRdataFull)
assignin('base', 'BMData6', reDATAG1SLOTSDATAFull)
assignin('base', 'BMData7', reDATAG2SLOTSDATAFull)
assignin('base', 'BMData8', reSGr)
BMData1= reDATAdatasetFull;
BMData2= reDATASAPdataFull;
BMData3= reDATADdataFull;
BMData4= reDATAAMPARdataFull;
BMData5= reDATAGluRdataFull;
BMData6= reDATAG1SLOTSDATAFull;
BMData7= reDATAG2SLOTSDATAFull;
BMData8= reSGr;





%=============================================%
% UNCOMMENT TO REINITIALIZE DATA
%{
AveOver = 10;
DATARATE = 100;
t = .1;
%}
reDATAdataset = BMData1;
reDATASAPdata = BMData2;
reDATADdata = BMData3;
reDATAAMPARdata = BMData4;
reDATAGluRdata = BMData5;
reDATAG1SLOTSDATA = BMData6;
reDATAG2SLOTSDATA = BMData7;
reSGr = BMData8;
%=============================================%




%%
%=============================================%
if ~doRun(3);
ChopStart = 1; 
if ChopStart
%------------------
chop = 4;
	
for p = 1:numel(reDATAdataset)
reDATAdataset{p}(1:chop,:) = [];
end

for p = 1:numel(reDATASAPdata)
reDATASAPdata{p}(1:chop,:) = [];
end

for p = 1:numel(reDATADdata)
reDATADdata{p}(1:chop,:) = [];
end

for p = 1:numel(reDATAAMPARdata)
reDATAAMPARdata{p}(1:chop,:) = [];
end

for p = 1:numel(reDATAGluRdata)
reDATAGluRdata{p}(1:chop,:) = [];
end

for p = 1:numel(reDATAG1SLOTSDATA)
reDATAG1SLOTSDATA{p}(1:chop,:) = [];
end

for p = 1:numel(reDATAG2SLOTSDATA)
reDATAG2SLOTSDATA{p}(1:chop,:) = [];
end


%------------------
end; % if ChopStart
%------------------
% for p = 1:numel(reDATAdataset)
% reDATAdataset = smooth(reDATAdataset{p},5);
% end





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					SPINE RADIAL DISTANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SGr={  G1P1r,	G1P2r,	G2P1r,	G2P2r,
%       G1inP1,	G1inP2,	G2inP1,	G2inP2,
%       G1xyP1,	G1xyP2,	G2xyP1,	G2xyP2  };
doRadDist = 1;
doRadDistDPlots=1;
if doRadDist



%=================================================
for reN = 1:loops
%=================================================

SGr = reSGr{reN};

%--------------------------------------
% Make concentric circle annulus
%--------------------------------------
S1r = round(S1rad);
Qhr = S1r/2;
Qx = linspace(0,S1r,6);
QxL = Qx - Qhr;
yL = sqrt(Qhr.^2 - QxL.^2);
Qr = sqrt(Qx.^2 + yL.^2);	% Concentric circles

%--------------------------------------
% Tally GluRs in each spine annulus
%--------------------------------------
for qq = 1:4
Q0=0;Q1=0;Q2=0;Q3=0;Q4=0;Q5=0;Qnan=0;
for mm = 1:numel(SGr)
	QD1 = (SGr{mm}{qq}(SGr{mm}{qq+4}(:)));
	for Nin=1:numel(QD1)
		if	   QD1(Nin) <= Qr(1)
			Q0 = Q0+1;
		elseif QD1(Nin) <= Qr(2)
			Q1 = Q1+1;
		elseif QD1(Nin) <= Qr(3)
			Q2 = Q2+1;
		elseif QD1(Nin) <= Qr(4)
			Q3 = Q3+1;
		elseif QD1(Nin) <= Qr(5)
			Q4 = Q4+1;
		elseif QD1(Nin) <= Qr(6)
			Q5 = Q5+1;
		elseif QD1(Nin) > Qr(7)
			Qnan = Qnan+1;
		end
	end
end
Qqs(qq,:) = [Q1 Q2 Q3 Q4 Q5] ./(stepN/DATARATE);
end

%--------------------------------------
% Parse out G1 G2 in SPY1 SPY2
%--------------------------------------
% Qqs =   
%  anullus 1  2  3  4  5
% G1-SPY1: N  N  N  N  N
% G1-SPY2: N  N  N  N  N
% G2-SPY1: N  N  N  N  N 
% G2-SPY2: N  N  N  N  N

% Raw Single
G1SP1 = Qqs(1,:);
G1SP2 = Qqs(2,:);
G2SP1 = Qqs(3,:);
G2SP2 = Qqs(4,:);

% Raw Multi					% Means						% Sums
G12SPY12 = Qqs;				G12SPY12m = mean(G12SPY12);	G12SPY12s = sum(G12SPY12);
G12SPY1 = [G1SP1 ; G2SP1];	G12SPY1m = mean(G12SPY1);	G12SPY1s = sum(G12SPY1);
G12SPY2 = [G1SP2 ; G2SP2];	G12SPY2m = mean(G12SPY2);	G12SPY2s = sum(G12SPY2);
G1SPY12 = [G1SP1 ; G1SP2];	G1SPY12m = mean(G1SPY12);	G1SPY12s = sum(G1SPY12);
G2SPY12 = [G2SP1 ; G2SP2];	G2SPY12m = mean(G2SPY12);	G2SPY12s = sum(G2SPY12);

%-------------------------------%

reG1SP1(reN,:) = G1SP1;
reG1SP2(reN,:) = G1SP2;
reG2SP1(reN,:) = G2SP1;
reG2SP2(reN,:) = G2SP2;

%=================================================
end; %for reN = 1:loops
%=================================================



mG1SP1 = mean(reG1SP1);
mG1SP2 = mean(reG1SP2);
mG2SP1 = mean(reG2SP1);
mG2SP2 = mean(reG2SP2);

sG1SP1 = std(reG1SP1);
sG1SP2 = std(reG1SP2);
sG2SP1 = std(reG2SP1);
sG2SP2 = std(reG2SP2);

semG1SP1 = sG1SP1 ./ sqrt(loops);
semG1SP2 = sG1SP2 ./ sqrt(loops);
semG2SP1 = sG2SP1 ./ sqrt(loops);
semG2SP2 = sG2SP2 ./ sqrt(loops);









%=====================================================
if doRadDistDPlots
for GSPswitch = 1:4;	
%=====================================================

switch GSPswitch
    case 1
		Qs = mean(reG1SP1);
	case 2
		Qs = mean(reG1SP2);
	case 3
		Qs = mean(reG2SP1);
	case 4
		Qs = mean(reG2SP2);
end

%--------------------------------------
% Interpolate Annulus Patch Densities
%--------------------------------------
[clrmap1 CRank] = ClrMap2;
[QsortVal, QsortInd] = sort(Qs, 2);
QcRankA = [fliplr(QsortInd) 6];
QcRankB =  circshift(QcRankA,[0 -1]);

for Nx = 1:5
clear xverts yverts verts faces Adatc Bdatc fvals xdata ydata cdata t xdat ydat cdat
	QcA = QcRankA(Nx);
	QcB = QcRankB(Nx);
	
	nverts = 50;
	t = linspace(0,2*pi,nverts);	

	xverts = [(Qr(Nx+1)*sin(t)) (Qr(Nx)*sin(t))]';
	yverts = [(Qr(Nx+1)*cos(t)) (Qr(Nx)*cos(t))]';
	verts = [xverts yverts];
	
	fvals = [1 2 (nverts+2) (nverts+1)];
	fvals = repmat(fvals,nverts-1,1);
	for mm = 1:(nverts-2); fvals(mm+1,:) = fvals(mm+1,:)+mm; end
	faces = fvals;
	
	Adatc = repmat(CRank{QcB}, nverts, 1);
	Bdatc = repmat(CRank{QcA}, nverts, 1);
	cdata = [Adatc; Bdatc];

	Qfaces(Nx) = {faces};
	Qverts(Nx) = {verts};
	Qcdata(Nx) = {cdata};
end


%--------------------------------------
% Plot Annulus Patch Densities
%--------------------------------------
Ax1LBWH = [-.02 .53 .40 .40	;...
			.47	.53 .40 .40	;...
		   -.02 .05 .40 .40	;...
			.47 .05 .40 .40 ;...
			.38 .60 .09 .35	;...
			.87	.60 .09 .35	;...
		    .38 .12 .09 .35	;...
			.87 .12 .09 .35	];
AxCB = [ 0.322221476510067 0.573033707865168 0.0168 0.320848938826467 ;...
		 0.813060402684563 0.576779026217228 0.0168 0.31585518102372 ;...
		 0.32289932885906 0.0898876404494382 0.0168 0.32334581772784 ;...
		 0.81389932885906 0.0898876404494382 0.0168 0.318352059925094 ];
figure1 = figure(44);
set(figure1,'Colormap',clrmap1,'OuterPosition',[400 150 1100 800],'Color',[1 1 1])
pause(.1)
%figure(figure1)
axes1 = subplot('Position',Ax1LBWH(GSPswitch,:),'Visible','off');
hold(axes1,'all');
pause(.1)
%--------
for Ntx = 1:5
p = patch('Faces',Qfaces{Ntx},'Vertices',Qverts{Ntx},'FaceColor','b');
set(p,'FaceColor','interp','FaceVertexCData',Qcdata{Ntx},...
'EdgeColor','interp','LineWidth',.1,...
'EdgeLighting','flat','FaceLighting','flat',...
'BackFaceLighting','lit')
hold on; material shiny; axis vis3d off;
end
pause(.1)
%--------
colorbar('peer',axes1,AxCB(GSPswitch,:),'YTickLabel',{num2str(QsortVal(1), '% 10.2f'),...
'-','-','-','-',num2str(QsortVal(3), '% 10.2f'),...
'-','-','-','-',num2str(QsortVal(end), '% 10.2f')},...
'LineWidth',1);
hold on;			  
%--------

switch GSPswitch

    case 1 % G1xyP1
		CQRall = [];CQG1P1=[];
		figure(figure1)
		for mm = 340:1:numel(SGr)
		QPh2 = scatter(SGr{mm}{9}(1,:),SGr{mm}{9}(2,:));
		hold on
		set(QPh2,'Marker','o','SizeData',20,'LineWidth',.5,...
			'MarkerFaceColor',[.2 1 1],'MarkerEdgeColor',[0 .4 .4])
		CQRall = cat(2,CQRall,[SGr{mm}{9}(1,:);SGr{mm}{9}(2,:)]);
		CQG1P1 = cat(2,CQG1P1,[SGr{mm}{9}(1,:);SGr{mm}{9}(2,:)]);
		end
		hTit1  = title('G1-S1');
		set([hTit1],'FontName','Euclid Extra','FontSize',16);
		hold on
		%--------
		[PSDY,PSDX] = find(ACTINp1);
		S1Ph2 = scatter(PSDX-50,PSDY-50);
		set(S1Ph2,'Marker','o','SizeData',60,'LineWidth',.5,...
			'MarkerFaceColor',[.9 .2 .2],'MarkerEdgeColor',[.4 0 0])
		hold on
		%--------
		axes1 = subplot('Position',Ax1LBWH(GSPswitch+4,:));
		axis off
		hold(axes1,'all');
		bpcols = [1 0 0;1 0 1;0 0 1;0 1 1];
		bplabs = {num2str(mG1SP1(1),'% 10.1f');...
				num2str(mG1SP1(2),'% 10.1f');...
				num2str(mG1SP1(3),'% 10.1f');...
				num2str(mG1SP1(4),'% 10.1f');...
				num2str(mG1SP1(5),'% 10.1f')};
		hbp = boxplot(reG1SP1 ...
			,'whisker',1 ...
			,'boxstyle','filled' ...
			,'colors',bpcols ...
			,'fullfactors','off'...
			,'labels',bplabs ...
			,'widths',.8 ...
			,'factorgap',[0] ...
			,'medianstyle','target');
		hold on
		annotation(figure1,'textarrow',[0.355 0.455],...
	[0.941 0.942],'TextEdgeColor','none','FontSize',14,'String','Center');
	case 2 % G1xyP2
		CQRall = [];CQG1P2=[];
		figure(figure1)
		for mm = 340:1:numel(SGr)
		QPh2 = scatter(SGr{mm}{10}(1,:),SGr{mm}{10}(2,:));
		set(QPh2,'Marker','o','SizeData',20,'LineWidth',.5,...
			'MarkerFaceColor',[.2 1 1],'MarkerEdgeColor',[0 .4 .4])
		CQRall = cat(2,CQRall,[SGr{mm}{10}(1,:);SGr{mm}{10}(2,:)]);
		CQG1P2 = cat(2,CQG1P2,[SGr{mm}{9}(1,:);SGr{mm}{9}(2,:)]);
		end
		hTit1  = title('G1-S2');
		set([hTit1],'FontName','Euclid Extra','FontSize',16);
		hold on
		%--------
		[PSDY,PSDX] = find(ACTINp2);
		S1Ph2 = scatter(PSDX-50,PSDY-50);
		set(S1Ph2,'Marker','o','SizeData',60,'LineWidth',.5,...
			'MarkerFaceColor',[.9 .2 .2],'MarkerEdgeColor',[.4 0 0])
		hold on
		%--------
		axes1 = subplot('Position',Ax1LBWH(GSPswitch+4,:));
		axis off
		hold(axes1,'all');
		bpcols = [1 0 0;1 0 1;0 0 1;0 1 1];
		bplabs = {num2str(mG1SP2(1),'% 10.1f');...
				  num2str(mG1SP2(2),'% 10.1f');...
				  num2str(mG1SP2(3),'% 10.1f');...
				  num2str(mG1SP2(4),'% 10.1f');...
				  num2str(mG1SP2(5),'% 10.1f')};
		hbp = boxplot(reG1SP2 ...
			,'whisker',1 ...
			,'boxstyle','filled' ...
			,'colors',bpcols ...
			,'fullfactors','off'...
			,'labels',bplabs ...
			,'widths',.8 ...
			,'factorgap',[0] ...
			,'medianstyle','target');
		hold on
		annotation(figure1,'textarrow',[0.845637583892617 0.942114093959731],...
	[0.941323345817728 0.941323345817728],'TextEdgeColor','none','FontSize',14,'String','Center');
	case 3 % G2xyP1
		CQRall = [];CQG2P1=[];
		figure(figure1)
		for mm = 280:1:numel(SGr)
		QPh2 = scatter(SGr{mm}{11}(1,:),SGr{mm}{11}(2,:));
		set(QPh2,'Marker','o','SizeData',20,'LineWidth',.5,...
			'MarkerFaceColor',[.2 1 1],'MarkerEdgeColor',[0 .4 .4])
		CQRall = cat(2,CQRall,[SGr{mm}{11}(1,:);SGr{mm}{11}(2,:)]);
		CQG2P1 = cat(2,CQG2P1,[SGr{mm}{9}(1,:);SGr{mm}{9}(2,:)]);
		end
		hTit1  = title('G2-S1');
		set([hTit1],'FontName','Euclid Extra','FontSize',16);
		hold on
		%--------
		[PSDY,PSDX] = find(ACTINp1);
		S1Ph2 = scatter(PSDX-50,PSDY-50);
		set(S1Ph2,'Marker','o','SizeData',60,'LineWidth',.5,...
			'MarkerFaceColor',[.9 .2 .2],'MarkerEdgeColor',[.4 0 0])
		hold on
		%--------
		axes1 = subplot('Position',Ax1LBWH(GSPswitch+4,:));
		axis off
		hold(axes1,'all');
		bpcols = [1 0 0;1 0 1;0 0 1;0 1 1];
		bplabs = {num2str(mG2SP1(1),'% 10.1f');...
				  num2str(mG2SP1(2),'% 10.1f');...
				  num2str(mG2SP1(3),'% 10.1f');...
				  num2str(mG2SP1(4),'% 10.1f');...
				  num2str(mG2SP1(5),'% 10.1f')};
		hbp = boxplot(reG2SP1 ...
			,'whisker',1 ...
			,'boxstyle','filled' ...
			,'colors',bpcols ...
			,'fullfactors','off'...
			,'labels',bplabs ...
			,'widths',.8 ...
			,'factorgap',[0] ...
			,'medianstyle','target');
		hold on
		annotation(figure1,'textarrow',[0.359060402684564 0.462248322147651],...
	[0.461922596754057 0.461922596754057],'TextEdgeColor','none','FontSize',14,'String','Center');
	case 4 % G2xyP2
		CQRall = [];CQG2P2=[];
		figure(figure1)
		for mm = 280:1:numel(SGr)
		QPh2 = scatter(SGr{mm}{12}(1,:),SGr{mm}{12}(2,:));
		set(QPh2,'Marker','o','SizeData',20,'LineWidth',.5,...
			'MarkerFaceColor',[.2 1 1],'MarkerEdgeColor',[0 .4 .4])
		CQRall = cat(2,CQRall,[SGr{mm}{12}(1,:);SGr{mm}{12}(2,:)]);
		CQG2P2 = cat(2,CQG2P2,[SGr{mm}{9}(1,:);SGr{mm}{9}(2,:)]);
		end
		hTit1  = title('G2-S2');
		set([hTit1],'FontName','Euclid Extra','FontSize',16);
		hold on
		%--------
		[PSDY,PSDX] = find(ACTINp2);
		S1Ph2 = scatter(PSDX-50,PSDY-50);
		set(S1Ph2,'Marker','o','SizeData',60,'LineWidth',.5,...
			'MarkerFaceColor',[.9 .2 .2],'MarkerEdgeColor',[.4 0 0])
		hold on
		%--------
		axes1 = subplot('Position',Ax1LBWH(GSPswitch+4,:));
		axis off
		hold(axes1,'all');
		bpcols = [1 0 0;1 0 1;0 0 1;0 1 1];
		bplabs = {num2str(mG2SP2(1),'% 10.1f');...
				  num2str(mG2SP2(2),'% 10.1f');...
				  num2str(mG2SP2(3),'% 10.1f');...
				  num2str(mG2SP2(4),'% 10.1f');...
				  num2str(mG2SP2(5),'% 10.1f')};
		hbp = boxplot(reG2SP2 ...
			,'whisker',1 ...
			,'boxstyle','filled' ...
			,'colors',bpcols ...
			,'fullfactors','off'...
			,'labels',bplabs ...
			,'widths',.8 ...
			,'factorgap',[0] ...
			,'medianstyle','target');
		hold on
		annotation(figure1,'textarrow',[0.849 0.944],...
	[0.465 0.465],'TextEdgeColor','none','FontSize',14,'String','Center');
	otherwise
        warning('Unexpected plot type.');
end

%--------
for mm = 2:5
rectangle('Position',[-Qr(mm),-Qr(mm),Qr(mm)*2,Qr(mm)*2],'Curvature',[1,1])
hold on;
end
%--------












%{
X = CQG1P1';
opts = statset('Display','final');
[idx,ctrs] = kmeans(X,15,'Distance','city','Replicates',10,'Options',opts);

for pidx = 1:15
plot(X(idx==pidx,1),X(idx==pidx,2),'o','MarkerSize',5,'MarkerEdgeColor',[rand rand rand])
hold on
end
plot(ctrs(:,1),ctrs(:,2),'kx','MarkerSize',10,'LineWidth',1)
plot(ctrs(:,1),ctrs(:,2),'ko','MarkerSize',10,'LineWidth',1)
hold on
[PSDY,PSDX] = find(ACTINp1);
S1Ph2 = scatter(PSDX-50,PSDY-50);
set(S1Ph2,'Marker','s','SizeData',60,'LineWidth',1,...
			'MarkerFaceColor',[.1 .9 .9],'MarkerEdgeColor',[.0 .7 .7])
hold on
pause(.5)
options = statset('Display','final','MaxIter',300);
obj = gmdistribution.fit(X,12,'Options',options);
h = ezcontour(@(x,y)pdf(obj,[x y]),[-50 50],[-50 50]);
%}


%=====================================================
end %for GSPswitch = 1:4;	
end; %if doRadDistDPlots
%=====================================================





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end; %if doRadDist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%%
keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					FINAL OUTPUT (SPLINE) FIGURE 1 OF 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveoutputfigs('STARShiP0')
%----------------------------------------------------------------------%
fig21 = figure(21);
set(21,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig21,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%

c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6]; c5= [.01 .9 .01];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7]; c55=[.01 .9 .01];
cOLOR = [c1; c2; c3; c4];
pause(1);

%===========================================================%
% FIG1 TOP LEFT: GluR Subtypes in Synapses
%===========================================================%
sbpos = [.05 .57 .4 .38]; ptype = 4;
cOLOR = [c1; c2; c3; c4; c1; c2; c3; c4];
itemN = 1; 
[ph1 hax1] = CIenvFun(reDATAAMPARdata,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'GluR1/2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 2;
[ph2 hax2] = CIenvFun(reDATAAMPARdata,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'GluR2/3');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%------------------------------------------%
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','Occupied Slots')
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 2;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
hTitle  = title('GluR Subtypes in Synapses');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles (+/- SEM)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
% ylim([0 haxes(4)*1.2 ]);
% xlim([0 (haxes(2)*.9)]);
%===========================================================%


%===========================================================%
% FIG1 TOP RIGHT: GluR Subtypes in PSD vs PSA
%===========================================================%
sbpos = [.55 .77 .4 .18]; ptype = 4;
cOLOR = [c1; c2; c3; c4];
itemN = 1; 
[ph1 hax1] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'PSD');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 2;
[ph2 hax2] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'PSA');
[LEGH,OBJH,OUTH,OUTM] = legend;
%------------------------------------------%
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','Occupied Slots')
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 5;
set(ph1,'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
hTitle  = title ('GluR Subtypes in PSD vs PSA');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('GluR1/2');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
% ylim([0 haxes(4)*1.2 ]);
% xlim([0 (haxes(2)*.9)]);
%===========================================================%
sbpos = [.55 .57 .4 .18]; ptype = 4;
cOLOR = [c1; c2; c3; c4];
itemN = 3; 
[ph1 hax1] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'PSD');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 4;
[ph2 hax2] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'PSA');
[LEGH,OBJH,OUTH,OUTM] = legend;
%------------------------------------------%
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','Occupied Slots')
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 5;
set(ph1,'LineStyle','-','Color', c3,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
set(ph2,'LineStyle','-','Color',c4,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('-');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('GluR2/3');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
% ylim([0 haxes(4)*1.2 ]);
% xlim([0 (haxes(2)*.9)]);
%===========================================================%



%===========================================================%
% FIG1 BOTTOM LEFT: GluR Subtypes in Synapse 1v2
%===========================================================%
sbpos = [.05 .28 .4 .18]; ptype = 4;
cOLOR = [c1; c2; c3; c4; c1; c2; c3; c4];
itemN = 5; 
[ph1 hax1] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'Synapse 1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 6;
[ph2 hax2] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'Synapse 2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 2;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
hTitle  = title ('GluR Subtypes in Synapse 1v2');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('GluR1/2');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
% ylim([0 haxes(4)]);
% xlim([0 (haxes(2)*.9)]);
%===========================================================%
sbpos = [.05 .09 .4 .18]; ptype = 4;
cOLOR = [c1; c2; c3; c4; c1; c2; c1; c2];
itemN = 7; 
[ph1 hax1] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'Synapse 1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 8;
[ph2 hax2] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'Synapse 2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 2;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
hTitle  = title ('-');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('GluR2/3');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
% ylim([0 haxes(4)]);
% xlim([0 (haxes(2)*.9)]);
%===========================================================%





%===========================================================%
% FIG1 BOTTOM RIGHT: Synapse Particle Counts
%===========================================================%
sbpos = [.55 .09 .4 .38]; ptype = 4;
cOLOR = [c1; c2; c1; c2; c3; c2; c3; c4; c5;c5;c5;c5;c5;c5;c5;c5;];
itemN = 3; 
[ph1 hax1] = CIenvFun(reDATAdataset,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'Synapse 1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 4;
[ph2 hax2] = CIenvFun(reDATAdataset,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'Synapse 2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 5; 
[ph3 hax3] = CIenvFun(reDATAdataset,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph3],OUTM{:},'Synapse Total');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on

itemN = 11;
[ph4 hax4] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph4],OUTM{:},'DFE');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%------------------------------------------%
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','Particles')
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 5;
set(ph1,'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(ph3,'LineStyle','-','Color',c3,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
set(ph4,'LineStyle','-','Color',c5,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c5,'MarkerFaceColor',c55);
hTitle  = title ('Synapse Particle Counts');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles (+/- SEM)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
% ylim([0 haxes(4)*1.2 ]);
% xlim([0 (haxes(2)*.9)]);
%===========================================================%







%===========================================================%
saveoutputfigs('STARShiP1')
%===========================================================%









%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					FINAL OUTPUT (SPLINE) FIGURE 2 OF 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig12 = figure(12);
set(12,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig12,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
cOLOR = [c1; c2; c3; c4; c1; c2; c3; c4];



%===========================================================%
% FIG2 TOP LEFT: Synaptic AMPARs
%===========================================================%
sbpos = [.055 .57 .4 .38]; ptype = 4;
cOLOR = [c1; c4; c1; c2; c3; c3];
itemN = 3; 
[ph1 hax1] = CIenvFun(reDATAdataset,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'Synapse 1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 4;
[ph2 hax2] = CIenvFun(reDATAdataset,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'Synapse 2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 5; 
%------------------------------------------%
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','Occupied Slots')
% xt = (get(gca,'XTick'));
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 2;
set(ph1,'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
% set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(ph2,'LineStyle','-.','Color',c2,'LineWidth',3,...
'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c2);
hTitle  = title ('Synaptic AMPARs');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles (+/- SEM)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
% ylim([haxes(3)/1.1 haxes(4)*1.1]);
% xlim([0 (haxes(2)*.9)]);
%======================================================================%





%===========================================================%
% FIG2 TOP RIGHT: GluR Subtypes in Synapse 1v2
%===========================================================%
% sbpos = [.05 .09 .4 .38]; ptype = 4;
sbpos = [.55 .57 .4 .38]; ptype = 4;
cOLOR = [c1; c2; c3; c4; c1; c2; c3; c4];
itemN = 5; 
[ph1 hax1] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'G1/2 Synapse 1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 6;
[ph2 hax2] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'G1/2 Synapse 2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 7; 
[ph3 hax3] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph3],OUTM{:},'G2/3 Synapse 1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 8;
[ph4 hax4] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph4],OUTM{:},'G2/3 Synapse 2');
hold on
%------------------------------------------%
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','Occupied Slots')
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 5;
set(ph1,'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(ph3,'LineStyle','-','Color',c3,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
set(ph4,'LineStyle','-','Color',c4,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('GluR Subtypes in Synapse 1v2');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles (+/- SEM)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
% ylim([0 haxes(4)]);
% xlim([0 (haxes(2)*.9)]);
%======================================================================%



%===========================================================%
% FIG2 BOTTOM RIGHT: Occupied Slots
%===========================================================%
sbpos = [.55 .09 .4 .38]; ptype = 4;
cOLOR = [c1; c2; c3; c4; c1; c2; c3; c4];
itemN = 1; 
[ph1 hax1] = CIenvFun(reDATAG1SLOTSDATA,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'G1/2 PSD 1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 2;
[ph2 hax2] = CIenvFun(reDATAG1SLOTSDATA,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'G1/2 PSD 2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
cOLOR = [c3; c4; c3; c4; c1; c2; c3; c4];
itemN = 1; 
[ph3 hax3] = CIenvFun(reDATAG2SLOTSDATA,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph3],OUTM{:},'G2/3 PSD 1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 2;
[ph4 hax4] = CIenvFun(reDATAG2SLOTSDATA,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph4],OUTM{:},'G2/3 PSD 2');
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 2;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(ph3,'LineStyle','-','Color',c3,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
set(ph4,'LineStyle','-','Color',c4,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('Occupied Slots');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles (+/- SEM)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
% ylim([0 haxes(4)*1.1 ]);
% xlim([0 (haxes(2)*.9)]);
%======================================================================%



%===========================================================%
% FIG2 BOTTOM LEFT: Synaptic SAP Cluster Size
%===========================================================%
sbpos = [.055 .09 .4 .38]; ptype = 4;
cOLOR = [c1; c2; c3; c4; c1; c2; c3; c4];
itemN = 1; 
[ph1 hax1] = CIenvFun(reDATASAPdata,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'Synapse 1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 2;
[ph2 hax2] = CIenvFun(reDATASAPdata,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'Synapse 2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 2;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
% set(ph2,'LineStyle',':','Color',c2,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(ph2,'LineStyle','-.','Color',c2,'LineWidth',3,...
'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c2);
hTitle  = title ('Synaptic SAP Cluster Size');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles (+/- SEM)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
% ylim([haxes(3)/1.1 haxes(4)*1.1 ]);
% xlim([0 (haxes(2)-3)]);
% ylim([min(min(reDATASAPdata{2}))/2 max(max(reDATASAPdata{2}))*1.2]);
% xlim([0 (haxes(2)*.9)]);
%======================================================================%


%======================================================================%
saveoutputfigs('STARShiP2')
%======================================================================%
%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					FINAL OUTPUT (SPLINE) FIGURE 3 OF 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig23 = figure(23);
set(23,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig23,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6]; c5= [.01 .9 .01];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7]; c55=[.01 .9 .01];


applered= [.9 .2 .2];
oceanblue= [.2 .4 .6];
neongreen = [.1 .9 .1];
liteblue = [.2 .9 .9];
hotpink=[.9 .1 .9];

cOLOR = [c1; c2; c3; c4];
%----------------------------------------------------------------------%



%===========================================================%
% FIG3 TOP LEFT: GluR Subtypes in Synapses
%===========================================================%
sbpos = [.05 .57 .4 .38]; ptype = 4;
cOLOR = [c1; c2; c3; c4];
cOLOR = repmat(cOLOR,15,1);
%------------------------------------------%
itemN = 12; 
[ph1 hax1 mu1] = CIplot(reDATAGluRdata,applered,5,sbpos,itemN,1);
legend(ph1,'G1-PSA1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%ph1 = plot(mu1);
%hold on
%---
itemN = 13; 
[ph2 hax2 mu2] = CIplot(reDATAGluRdata,oceanblue,5,sbpos,itemN,1);
legend([OUTH;ph2],OUTM{:},'G1-PSA2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%ph2 = plot(mu2);
%hold on
%---
itemN = 14; 
[ph3 hax3 mu3] = CIplot(reDATAGluRdata,neongreen,5,sbpos,itemN,1);
legend([OUTH;ph3],OUTM{:},'G1-PSD1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%ph3 = plot(mu3);
%hold on
%---
itemN = 15; 
[ph4 hax4 mu4] = CIplot(reDATAGluRdata,hotpink,5,sbpos,itemN,1);
legend([OUTH;ph4],OUTM{:},'G1-PSD2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%ph4 = plot(mu4);
%hold on
%------------------------------------------%
xt = (get(gca,'XTick'))*100;
set(gca,'XTickLabel', sprintf('%.0f|',xt))
%------------------------------------------%
MS1 = 5; LnW=1;
set(ph1,'LineStyle','-','Color',applered,'LineWidth',LnW,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',applered,'MarkerFaceColor',applered);
set(ph2,'LineStyle','-','Color',oceanblue,'LineWidth',LnW,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',oceanblue,'MarkerFaceColor',oceanblue);
set(ph3,'LineStyle',':','Color',neongreen,'LineWidth',LnW,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',neongreen,'MarkerFaceColor',neongreen);
set(ph4,'LineStyle',':','Color',hotpink,'LineWidth',LnW,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',hotpink,'MarkerFaceColor',hotpink);
hTitle  = title('Distribution of Particles with Brownian Motion');
hXLabel = xlabel('Time');
hYLabel = ylabel('Particles (+/- SEM)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12);
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
%ylim([0 haxes(4)*1.2 ]);
% xlim([0 (haxes(2)*.9)]);
%======================================================================%






%===========================================================%
% FIG3 TOP RIGHT: Synapse Particle Counts
%===========================================================%
sbpos = [.55 .57 .4 .38]; ptype = 4;
cOLOR = [c1; c2; c3; c4];
cOLOR = repmat(cOLOR,15,1);
%------------------------------------------%
itemN = 16; 
[ph1 hax1 mu1] = CIplot(reDATAGluRdata,applered,5,sbpos,itemN,1);
legend(ph1,'G2-PSA1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%ph1 = plot(mu1);
%hold on
%---
itemN = 17; 
[ph2 hax2 mu2] = CIplot(reDATAGluRdata,oceanblue,5,sbpos,itemN,1);
legend([OUTH;ph2],OUTM{:},'G2-PSA2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%ph2 = plot(mu2);
%hold on
%---
itemN = 18; 
[ph3 hax3 mu3] = CIplot(reDATAGluRdata,neongreen,5,sbpos,itemN,1);
legend([OUTH;ph3],OUTM{:},'G2-PSD1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%ph3 = plot(mu3);
%hold on
%---
itemN = 19; 
[ph4 hax4 mu4] = CIplot(reDATAGluRdata,hotpink,5,sbpos,itemN,1);
legend([OUTH;ph4],OUTM{:},'G2-PSD2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%ph4 = plot(mu4);
%hold on
%------------------------------------------%
xt = (get(gca,'XTick'))*100;
set(gca,'XTickLabel', sprintf('%.0f|',xt))
%------------------------------------------%
MS1 = 5; LnW=1;
set(ph1,'LineStyle','-','Color',applered,'LineWidth',LnW,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',applered,'MarkerFaceColor',applered);
set(ph2,'LineStyle','-','Color',oceanblue,'LineWidth',LnW,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',oceanblue,'MarkerFaceColor',oceanblue);
set(ph3,'LineStyle',':','Color',neongreen,'LineWidth',LnW,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',neongreen,'MarkerFaceColor',neongreen);
set(ph4,'LineStyle',':','Color',hotpink,'LineWidth',LnW,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',hotpink,'MarkerFaceColor',hotpink);

hTitle  = title('Distribution of Particles with Brownian Motion');
hXLabel = xlabel('Time');
hYLabel = ylabel('Particles (+/- SEM)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12);
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
%ylim([0 haxes(4)*1.2 ]);
% xlim([0 (haxes(2)*.9)]);
%======================================================================%





%===========================================================%
% FIG3 BOTTOM LEFT: GluR Subtypes in Synapse 1v2
%===========================================================%
sbpos = [.05 .09 .4 .38]; ptype = 4;
cOLOR = [c1; c2; c3; c4];
cOLOR = repmat(cOLOR,15,1);
%------------------------------------------%
itemN = 12; 
[ph1 hax1 mu1] = CIplot(reDATAGluRdata,applered,5,sbpos,itemN,1);
legend(ph1,'G1-PSA1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%ph1 = plot(mu1);
%hold on
%---
itemN = 13; 
[ph2 hax2 mu2] = CIplot(reDATAGluRdata,oceanblue,5,sbpos,itemN,1);
legend([OUTH;ph2],OUTM{:},'G1-PSA2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%ph2 = plot(mu2);
%hold on
%---
itemN = 16; 
[ph3 hax3 mu3] = CIplot(reDATAGluRdata,neongreen,5,sbpos,itemN,1);
legend([OUTH;ph3],OUTM{:},'G2-PSA1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%ph3 = plot(mu3);
%hold on
%---
itemN = 17; 
[ph4 hax4 mu4] = CIplot(reDATAGluRdata,hotpink,5,sbpos,itemN,1);
legend([OUTH;ph4],OUTM{:},'G2-PSA2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%ph4 = plot(mu4);
%hold on
%------------------------------------------%
xt = (get(gca,'XTick'))*100;
set(gca,'XTickLabel', sprintf('%.0f|',xt))
%------------------------------------------%
MS1 = 5; LnW=1;
set(ph1,'LineStyle','-','Color',applered,'LineWidth',LnW,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',applered,'MarkerFaceColor',applered);
set(ph2,'LineStyle','-','Color',oceanblue,'LineWidth',LnW,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',oceanblue,'MarkerFaceColor',oceanblue);
set(ph3,'LineStyle',':','Color',neongreen,'LineWidth',LnW,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',neongreen,'MarkerFaceColor',neongreen);
set(ph4,'LineStyle',':','Color',hotpink,'LineWidth',LnW,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',hotpink,'MarkerFaceColor',hotpink);
hTitle  = title('Distribution of Particles with Brownian Motion');
hXLabel = xlabel('Time');
hYLabel = ylabel('Particles (+/- SEM)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12);
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
%ylim([0 haxes(4)*1.2 ]);
% xlim([0 (haxes(2)*.9)]);
%======================================================================%





%===========================================================%
% FIG3 BOTTOM RIGHT: Synapse Particle Counts
%===========================================================%
sbpos = [.55 .09 .4 .38]; ptype = 4;
cOLOR = [c1; c2; c3; c4];
cOLOR = repmat(cOLOR,15,1);
%------------------------------------------%
itemN = 14; 
[ph1 hax1 mu1] = CIplot(reDATAGluRdata,applered,5,sbpos,itemN,1);
legend(ph1,'G1-PSD1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%ph1 = plot(mu1);
%hold on
%---
itemN = 15; 
[ph2 hax2 mu2] = CIplot(reDATAGluRdata,oceanblue,5,sbpos,itemN,1);
legend([OUTH;ph2],OUTM{:},'G1-PSD2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%ph2 = plot(mu2);
%hold on
%---
itemN = 18; 
[ph3 hax3 mu3] = CIplot(reDATAGluRdata,neongreen,5,sbpos,itemN,1);
legend([OUTH;ph3],OUTM{:},'G2-PSD1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%ph3 = plot(mu3);
%hold on
%---
itemN = 19; 
[ph4 hax4 mu4] = CIplot(reDATAGluRdata,hotpink,5,sbpos,itemN,1);
legend([OUTH;ph4],OUTM{:},'G2-PSD2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%ph4 = plot(mu4);
%hold on
%------------------------------------------%
xt = (get(gca,'XTick'))*100;
set(gca,'XTickLabel', sprintf('%.0f|',xt))
%------------------------------------------%
MS1 = 5; LnW=1;
set(ph1,'LineStyle','-','Color',applered,'LineWidth',LnW,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',applered,'MarkerFaceColor',applered);
set(ph2,'LineStyle','-','Color',oceanblue,'LineWidth',LnW,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',oceanblue,'MarkerFaceColor',oceanblue);
set(ph3,'LineStyle',':','Color',neongreen,'LineWidth',LnW,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',neongreen,'MarkerFaceColor',neongreen);
set(ph4,'LineStyle',':','Color',hotpink,'LineWidth',LnW,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',hotpink,'MarkerFaceColor',hotpink);
hTitle  = title('Distribution of Particles with Brownian Motion');
hXLabel = xlabel('Time');
hYLabel = ylabel('Particles (+/- SEM)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12);
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
%ylim([0 haxes(4)*1.2 ]);
% xlim([0 (haxes(2)*.9)]);
%======================================================================%


%======================================================================%
saveoutputfigs('STARShiP3')
%======================================================================%


















%%
%---------------
end % if ~doRun(3);
%---------------


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%						MSD REVERSE TRACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRACK MSD AND STEP SIZE MEANS (0:no 1:yes)
trackMSD = doRun(3);
if trackMSD

dT = t;
G2N = 100;
Nsteps = G2N;
tracks = cell(G2N, 1);
lims = ((D+1)^2)*10;

trackStepMeans = DontDo;
MSDdrop = MSDdrop;
strdropES = 'D-ES';
strdropPSD1 = 'D-PSD1';
strdropPSD2 = 'D-PSD2';
testESMSD = strcmp(MSDdrop,strdropES);
testPSD1MSD = strcmp(MSDdrop,strdropPSD1);
testPSD2MSD = strcmp(MSDdrop,strdropPSD2);
MSDtest = [testESMSD testPSD1MSD testPSD2MSD];
MANstepsize = 0;

stepN = 1;
for Nt = 1:Nsteps 
%-------------------------------%
    %-------------------------------%
    %     STEP DIRECTION & SIZE
    %-------------------------------%
    xysG2 = STEPxyds(G2N, k);
	%-------------------------------%
    %     DO MANUAL STEP SIZES
    %-------------------------------%
	if MANstepsize
		xysG2 = CIRCSTEPxyds(G2N, k);
	end % xyd = DIRxyd(xyd);
	
		
	[xyG2 xysG2] = MSDAMPARSTEP(G2N, xysG2, xyG2,... 
	testPSD1MSD, testPSD2MSD, testESMSD, PSD1, PSD2, D, Dn_PSD1, Dn_PSD2);


    [tracks] = MSDfun(stepN, Nsteps, tracks, xysG2);

    
    if trackStepMeans
     if mod(stepN, 10) == 0
        StepMeans = meanStep(xysG2, k, D);
        head = ['MeanAbsStep CalcHalfNorm MeanPythStep CalcPythStep CalcHalfPythStep'];
       disp(head)
       disp(StepMeans)
     end
    end


 stepN = stepN+1;
end
output1 = tracks;
[rawMSDout] = MSDfunction(tracks, G2N, Nsteps, D, Dn_PSD1, Dn_PSD2, dT, k, MSDtest);


%===========================================================%
%               LIVE PARTICLE DIFFUSION
%-----------------------------------------------------------%
for Nt = 1:Nsteps

    xysG2 = STEPxyds2(G2N, k);
	[xyG2] = AMPARSTEP2(G2N, xysG2, xyG2);
	MAINPLOT2(xyG2, lims);

end
%===========================================================%


%===========================================================%
%               MSD RANDOM STEPS ANALYSIS
%-----------------------------------------------------------%
tracks = cell(G2N, 1);

stepN = 1;
for Nt = 1:Nsteps 

	xysG2 = STEPxyds2(G2N, k);
	[xyG2] = AMPARSTEP2(G2N, xysG2, xyG2);
    [tracks] = MSDfun2(stepN, Nsteps, tracks, xysG2);

stepN = stepN+1;
end
[ESMSDout] = MSDfunction2(tracks,G2N,Nsteps,D,Dn,L,dT,k,Scale,MSDtest);

%===========================================================%
%               MSD UNIFORM STEPS ANALYSIS
%-----------------------------------------------------------%
stepN = 1;
for t = 1:Nsteps 

	xysG2 = UNIFORMSTEPS2(G2N, Lx);
	[xyG2 xysG2] = SCALEUNIFORMSTEPS2(G2N, xysG2, xyG2, Ls, MSDtest);
    [tracks] = MSDfun2(stepN, Nsteps, tracks, xysG2);

stepN = stepN+1;
end
[PSDMSDout] = MSDfunction2(tracks,G2N,Nsteps,D,Dn,L,dT,k,Scale,MSDtest);

figure(80)
set(gcf,'OuterPosition',[600,400,300,300])
MSDfig = gcf;
figure(MSDfig)
subplot(2,1,1),text(0.5,0.5,strcat('PSD:',' \bullet  ', num2str(ESMSDout), '  µm²/s'),...
	'FontSize',24,'HorizontalAlignment','center','BackgroundColor',[.7 .9 .7]);
subplot(2,1,2),text(0.5,0.5,strcat('PSD:',' \bullet  ', num2str(PSDMSDout), '  µm²/s'),...
	'FontSize',24,'HorizontalAlignment','center','BackgroundColor',[.7 .9 .7]);

end % if trackMSD
%===========================================%









%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Do Izhikevich Spiking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doIz
spikes=DATAdataset(:,3:4);
IzhFun(spikes);
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%				OUTPUT RETURNS FOR TOP LEVEL FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% varargout = {GluRdata;SAPdata;Ddata;AMPARdata};
varargout = {DATAdataset;DATASAPdata;DATADdata;DATAAMPARdata;DATAGluRdata;...
	DATAG1SLOTSDATA;DATAG2SLOTSDATA};
assignin('base', 'varargout', varargout)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %			##    END MAIN FUNCTION   ##
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 













%%						 SUBROUTINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					 ##		INDEX		##
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%				1. BROWNIAN MOTION FUNCTIONS
%				2. SAPS & SLOTS
%				3. DENDRITIC FIELD MAPPING TOOLS
%				4. PLOTTING AND LIVE SIMULATION
%				5. PARTICLE COUNTERS
%				6. MSD BROWNIAN MOTION ANALYSIS TOOLS
%=================############################=================%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%					BROWNIAN MOTION FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%===================================%
% STEP SIZE GENERATOR
%===================================%
function xyds = STEPxyds(Ndots, k)
    xyds = (k * randn(2,Ndots));
end
%===================================%



%===================================%
% MOVE PARTICLES MAIN FUNCTION
%===================================%
function [xysG1 xyG1 xysG2 xyG2]...
    = MOVEGLUR(stepN,DF,...
	G1N,xysG1,xyG1,...
    G2N,xysG2,xyG2)

%------------------

% This function will add step-sizes to the particles current location.
% This function has been updated so particles located in a position
% less than 1 (ie. xyG1<1) are out-of-bounds (instead of particles
% located in a position less than zero). This will make it so all particles
% will have a rounded particle location of 1 or greater, which can be used
% as a matrix index ID value (zero cannot be used as a Mx index number). 


%------------------
% add step to current location
xyG1 = xyG1+xysG1;
xyG2 = xyG2+xysG2;
%------------------



%------------------
% G1 & G2 Containment

%LOWVAL=0;
LOWVAL=1;

G1XLO = xyG1(1,:) < LOWVAL;
G1XHI = xyG1(1,:) > DF.X;
G1YLO = xyG1(2,:) < LOWVAL;
G1YHI = xyG1(2,:) > DF.Y;

G2XLO = xyG2(1,:) < LOWVAL;
G2XHI = xyG2(1,:) > DF.X;
G2YLO = xyG2(2,:) < LOWVAL;
G2YHI = xyG2(2,:) > DF.Y;


xyG1(1,G1XLO) = LOWVAL;
xyG1(1,G1XHI) = DF.X;
xyG1(2,G1YLO) = LOWVAL;
xyG1(2,G1YHI) = DF.Y;


xyG2(1,G2XLO) = LOWVAL;
xyG2(1,G2XHI) = DF.X;
xyG2(2,G2YLO) = LOWVAL;
xyG2(2,G2YHI) = DF.Y;



%------------------
% OLD
%------------------
%{
% GluR1

% add step to current location
xyG1 = xyG1+xysG1;

% containment
for j = 1:G1N 
    if xyG1(1,j)>(DF.X) || xyG1(1,j)<LOWVAL
            xyG1(1,j) = uint16(sign(xyG1(1,j)))*(DF.X);
    elseif xyG1(2,j)>(DF.Y) || xyG1(2,j)<LOWVAL
            xyG1(2,j) = uint16(sign(xyG1(2,j)))*(DF.Y);
    end    
        
end


%------------------
% GluR2

% add step to current location
xyG2 = xyG2+xysG2;

% containment
for j = 1:G2N 
	if xyG2(1,j)>(DF.X) || xyG2(1,j)<LOWVAL
            xyG2(1,j) = uint16(sign(xyG2(1,j)))*(DF.X);
	elseif xyG2(2,j)>(DF.Y) || xyG2(2,j)<LOWVAL
            xyG2(2,j) = uint16(sign(xyG2(2,j)))*(DF.Y);
	end	   
        
end
%}

%===================================%
% if stepN == 500; keyboard; end
end
%===================================%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%						SAPS & SLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%----------------------------------%
%		SUPERSLOT
%----------------------------------%
function [xyG1 xyG2 xysG1 xysG2...
	G1FSLOTS G2FSLOTS...
	GluR1_TdwellSPYN GluR1_TdwellPSD GluR1_TdwellPERI...
	GluR2_TdwellSPYN GluR2_TdwellPSD GluR2_TdwellPERI...
	G1INSPA G2INSPA SMx G1v G2v...
	G1P1r G1P2r G2P1r G2P2r G1inP1 G1inP2 G2inP1 G2inP2...
	S1G1Cv S2G1Cv S1G2Cv S2G2Cv]...
	= SUPERSLOT(stepN,t,xyG1,xyG2,xysG1,xysG2,...
	GluR1_TdwellPSD,GluR1_TdwellPERI,GluR1_TdwellSPYN,LsGR1psa,LsGR1psd,LsGR1S,...
	GluR2_TdwellPSD,GluR2_TdwellPERI,GluR2_TdwellSPYN,LsGR2psa,LsGR2psd,LsGR2S,...	
	SSsz,G1xy,G2xy,SMx,SS,DF,GR1xy,GR2xy,...
	S1G1C,S2G1C,S1G2C,S2G2C,S1rad,S2rad)



%---------------------------
% EXPAND SMx to MASK
%---------------------------
%SMxMask = [0 0 0; 0 1 0; 0 0 0];
%SMxMask = ones(10);
%cSMx = convn(SMx,SMxMask,'same');

cSMx = SMx;


%---------------------------
% GLUR1 PSD STUFF
%---------------------------
G1SMx = G1xy.*cSMx;

[yMx xMx] = find(G1SMx);
G1SMxy = [xMx'; yMx'];
% ALL PARTICLES AT G1SMxy SHOULD BE SLOTTED

G1v=[];
doSLT = numel(G1SMxy);
if doSLT
mnm=1;
for mm = 1:numel(GR1xy(1,:))
	for nn = 1:numel(G1SMxy(1,:))
	matchd = (GR1xy(:,mm) == G1SMxy(:,nn));
	if matchd
		G1v(mnm) = mm;
		mnm=mnm+1;
	end
	end
end
end



%---------------------------
% GLUR2 PSD STUFF
%---------------------------
G2SMx = G2xy.*cSMx;

[yMx xMx] = find(G2SMx);
G2SMxy = [xMx'; yMx'];
% ALL PARTICLES IN G2SMxy SHOULD BE SLOTTED

G2v=[];
doSLT = numel(G2SMxy);
if doSLT
mnm=1;
for mm = 1:numel(GR2xy(1,:))
	for nn = 1:numel(G2SMxy(1,:))
	matchd = (GR2xy(:,mm) == G2SMxy(:,nn));
	if matchd
		G2v(mnm) = mm;
		mnm=mnm+1;
	end
	end
end
end


%if stepN == 1000; keyboard; end




% ADJUST STEP SIZE OF SLOTTED GLUR
%-------------
xysG1(:,G1v) = xysG1(:,G1v)*(LsGR1S);
xysG2(:,G2v) = xysG2(:,G2v)*(LsGR2S);




% SAVE SLOT INFO
%-------------
% GluR1
G1Slot = xyG1(:,G1v);
G1inS1 = G1Slot(2,:) > (DF.Y/2);
G1inS2 = ~G1inS1;
G1PSD1FSLOTS = sum(G1inS1);
G1PSD2FSLOTS = sum(G1inS2);
G1FSLOTS = [G1PSD1FSLOTS G1PSD2FSLOTS 0 0];

% GluR2
G2Slot = xyG2(:,G2v);
G2inS1 = G2Slot(2,:) > (DF.Y/2);
G2inS2 = ~G2inS1;
G2PSD1FSLOTS = sum(G2inS1);
G2PSD2FSLOTS = sum(G2inS2);
G2FSLOTS = [G2PSD1FSLOTS G2PSD2FSLOTS 0 0];




%############################################################


S1G1Cv = xyG1-S1G1C;
S2G1Cv = xyG1-S2G1C;
S1G2Cv = xyG2-S1G2C;
S2G2Cv = xyG2-S2G2C;

G1P1r = sqrt(sum((S1G1Cv).^2));
G1P2r = sqrt(sum((S2G1Cv).^2));
G2P1r = sqrt(sum((S1G2Cv).^2));
G2P2r = sqrt(sum((S2G2Cv).^2));

G1inP1=G1P1r<S1rad;
G1inP2=G1P2r<S2rad;
G2inP1=G2P1r<S1rad;
G2inP2=G2P2r<S2rad;



%----------------------------------------
% Step Size adjustments (location based)
%----------------------------------------

xysG1(:,G1inP1) = xysG1(:,G1inP1) .* LsGR1psa;
xysG1(:,G1inP2) = xysG1(:,G1inP2) .* LsGR1psa;
xysG2(:,G2inP1) = xysG2(:,G2inP1) .* LsGR2psa;
xysG2(:,G2inP2) = xysG2(:,G2inP2) .* LsGR2psa;



% IN EXTRASYNAPSE THE SIZE OF SPINES
G1INSPA = inboxfun([101 250],[201 350],xyG1);
G2INSPA = inboxfun([101 250],[201 350],xyG2);



% % GluR1 
% %-------------
% % SPINE %
% %-------%
% G1INSPYN1 = inboxfun(DF.S1LB,DF.S1RT,xyG1);
% G1INSPYN2 = inboxfun(DF.S2LB,DF.S2RT,xyG1);
% 
% xysG1(:,G1INSPYN1) = xysG1(:,G1INSPYN1)*(LsGR1psa);
% xysG1(:,G1INSPYN2) = xysG1(:,G1INSPYN2)*(LsGR1psa);
% 
% GluR1_TdwellSPYN(1,G1INSPYN1) = GluR1_TdwellSPYN(1,G1INSPYN1)+t;
% GluR1_TdwellSPYN(1,~G1INSPYN1) = 0;
% GluR1_TdwellSPYN(2,G1INSPYN2) = GluR1_TdwellSPYN(2,G1INSPYN2)+t;
% GluR1_TdwellSPYN(2,~G1INSPYN2) = 0;
% 
% %-------------
% % PSD   %
% %-------%
% G1INPSD1 = inboxfun(DF.PS1LB,DF.PS1RT,xyG1);
% G1INPSD2 = inboxfun(DF.PS2LB,DF.PS2RT,xyG1);
%   
% GluR1_TdwellPSD(1,G1INPSD1) = GluR1_TdwellPSD(1,G1INPSD1)+t;
% GluR1_TdwellPSD(1,~G1INPSD1) = 0;
% GluR1_TdwellPSD(2,G1INPSD2) = GluR1_TdwellPSD(2,G1INPSD2)+t;
% GluR1_TdwellPSD(2,~G1INPSD2) = 0;
% 
% %-------------
% % PERI  %
% %-------%
% G1INPERI1 = G1INSPYN1 - G1INPSD1; G1INPERI1 = G1INPERI1>0;
% G1INPERI2 = G1INSPYN2 - G1INPSD2; G1INPERI2 = G1INPERI2>0;
% 
% GluR1_TdwellPERI(1,G1INPERI1) = GluR1_TdwellPERI(1,G1INPERI1)+t;
% GluR1_TdwellPERI(1,~G1INPERI1) = 0;
% GluR1_TdwellPERI(2,G1INPERI2) = GluR1_TdwellPERI(2,G1INPERI2)+t;
% GluR1_TdwellPERI(2,~G1INPERI2) = 0;
% %-------------





% % GluR2
% %-------------
% % SPINE %
% %-------%
% G2INSPYN1 = inboxfun(DF.S1LB,DF.S1RT,xyG2);
% G2INSPYN2 = inboxfun(DF.S2LB,DF.S2RT,xyG2);
% 
% xysG2(:,G2INSPYN1) = xysG2(:,G2INSPYN1)*(LsGR2psa);
% xysG2(:,G2INSPYN2) = xysG2(:,G2INSPYN2)*(LsGR2psa);
% 
% GluR2_TdwellSPYN(1,G2INSPYN1) = GluR2_TdwellSPYN(1,G2INSPYN1)+t;
% GluR2_TdwellSPYN(1,~G2INSPYN1) = 0;
% GluR2_TdwellSPYN(2,G2INSPYN2) = GluR2_TdwellSPYN(2,G2INSPYN2)+t;
% GluR2_TdwellSPYN(2,~G2INSPYN2) = 0;
% 
% %-------------
% % PSD   %
% %-------%
% G2INPSD1 = inboxfun(DF.PS1LB,DF.PS1RT,xyG2);
% G2INPSD2 = inboxfun(DF.PS2LB,DF.PS2RT,xyG2);
% 
% GluR2_TdwellPSD(1,G2INPSD1) = GluR2_TdwellPSD(1,G2INPSD1)+t;
% GluR2_TdwellPSD(1,~G2INPSD1) = 0;
% GluR2_TdwellPSD(2,G2INPSD2) = GluR2_TdwellPSD(2,G2INPSD2)+t;
% GluR2_TdwellPSD(2,~G2INPSD2) = 0;
% 
% %-------------
% % PERI  %
% %-------%
% G2INPERI1 = G2INSPYN1 - G2INPSD1; G2INPERI1 = G2INPERI1>0;
% G2INPERI2 = G2INSPYN2 - G2INPSD2; G2INPERI2 = G2INPERI2>0;
% 
% GluR2_TdwellPERI(1,G2INPERI1) = GluR2_TdwellPERI(1,G2INPERI1)+t;
% GluR2_TdwellPERI(1,~G2INPERI1) = 0;
% GluR2_TdwellPERI(2,G2INPERI2) = GluR2_TdwellPERI(2,G2INPERI2)+t;
% GluR2_TdwellPERI(2,~G2INPERI2) = 0;
% %-------------









%if stepN == 2000; keyboard; end


end
%----------------------------------%




%----------------------------------%
% S1_MainClusterFun
%----------------------------------%
function [S1 SMx] = S1_MainClusterFun(S1,SMx,sap,hkMask,ACTINp1,stepN)

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


%------------
Pmx = rand(size(S));
Soc = (S>0);
Sno = ~Soc;
%---
hk = convn(ACTINp1,hkMask,'same');
%------------
Pon = 1 ./ (1+exp((hk-Lon).*(-Bon)));
Pkon = Sno .* ( Ron * dT * Pon );
%---
Son = (Pkon>Pmx);
%------------
Poff = 1 ./ (1+exp(((-hk)+Loff).*(-Boff)));
Pkoff = Soc .* ( Roff * dT * Poff );
%---
Soff = (Pkoff>Pmx);
%------------

S = (Soc-Soff) + Son;
S1=S;

% S1CP = {hk,Pkon,Sno};

end
%----------------------------------%



%----------------------------------%
% S2_MainClusterFun
%----------------------------------%
function [S2 SMx] = S2_MainClusterFun(S2,SMx,sap,hkMask,ACTINp2,stepN)
 

%=========================================================%
%   TEST MASK AGAINST CLUSTER - GET CONVOLUTION MX
%=========================================================%
S=S2;
dT = sap(12);

Lon = sap(27);	% On Energy (lower = more on events)
Bon = sap(28);	% On Neighbor-Independant Rate (new growth) (lower = more on)
Ron = sap(29);	% On Neighbor-Dependant Rate (cluster fill-in) (higher = more on)

Loff = sap(30);	% Off Energy (higher = more off events)
Boff = sap(31);	% Off Neighbor-Inependant Rate (uniform off)  (lower = more off)
Roff = sap(32);	% Off Neighbor-Dependant Rate (edge off) (higher = more off)

%------------
Pmx = rand(size(S));
Soc = (S>0);
Sno = ~Soc;
%---
hk = convn(ACTINp2,hkMask,'same');
%------------
Pon = 1 ./ (1+exp((hk-Lon).*(-Bon)));
Pkon = Sno .* ( Ron * dT * Pon );
%---
Son = (Pkon>Pmx);
%------------
Poff = 1 ./ (1+exp(((-hk)+Loff).*(-Boff)));
Pkoff = Soc .* ( Roff * dT * Poff );
%---
Soff = (Pkoff>Pmx);
%------------

S = (Soc-Soff) + Son;
S2=S;

% S2CP = {hk,Pkon,Sno};

end
%----------------------------------%




%----------------------------------%
%			MaskFun
%----------------------------------%
% Creates a kernel mask using a 3D Gaussian
% that can be used for SAP clustering and
% other matrix conv() functions
%----------------------------------%
function [hkMask] = MaskFun(doGMask,AMX,MSK)

%-------------------------------%
%		Mask Setup
%-------------------------------%
hkMask=ones(AMX{6});
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
%		inboxfun
%----------------------------------%
% Tests whether particles are in a box polygon
% and returns a logical vector
%----------------------------------%
function [inbox] = inboxfun(LB,RT,xyl)

if LB(1)>RT(1)
	LBt=LB;
	RTt=RT;
	LB(1)=RTt(1);
	RT(1)=LBt(1);
end
if LB(2)>RT(2)
	LBt=LB;
	RTt=RT;
	LB(2)=RTt(2);
	RT(2)=LBt(2);
end

xylLB1 = xyl(1,:) > LB(1);
xylRT1 = xyl(1,:) < RT(1);
xylLB2 = xyl(2,:) > LB(2);
xylRT2 = xyl(2,:) < RT(2);

inbox = xylLB1 & xylRT1 & xylLB2 & xylRT2;

end
%----------------------------------%



%----------------------------------%
%			CircMask
%----------------------------------%
% Creates a kernel mask using a 3D Gaussian
% that can be used for SAP clustering and
% other matrix conv() functions
%----------------------------------%
function [hkMask] = CircMask(GN)


%--------------------
GNpk = GN(1);	% hight of peak
GNx0 = 0;	% x-axis peak locations
GNy0 = 0;	% y-axis peak locations
GNsd = GN(2);	% sigma (stdev of slope)

GNnum = GN(3);
GNres = GN(4);
GNspr = ((GNnum-1)*GNres)/2;

a = .5/GNsd^2;
c = .5/GNsd^2;

[X, Y] = meshgrid((-GNspr):(GNres):(GNspr), (-GNspr):(GNres):(GNspr));
Z = GNpk*exp( - (a*(X-GNx0).^2 + c*(Y-GNy0).^2)) ;


hkMask=Z;

%----------------%
% FIGURE: 3D Gaussian Distribution
fh65 = figure(65); %set(fh5,'OuterPosition',(scsz./[2e-3 2e-3 2 2]))
set(fh65,'OuterPosition',[200 200 500 300]);
%----------------%
figure(fh65)
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
%----------------------------------%








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%				DENDRITIC FIELD MAPPING TOOLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%----------------------------------%
%			FieldFun
%----------------------------------%
% FieldFun
% short description of this function
%----------------------------------%
function [DFszX,DFszY,PSD1sz,PSD2sz,PSA1sz,PSA2sz,SPY1sz,SPY2sz,...
		 S1sz,S2sz,SSsz,SMx,DFum,SS,DF] = FieldFun(S1,S2,um,Scale)


% DENDRITIC FIELD
DFszX = round(um(1)/Scale);
DFszY = round(um(2)/Scale);
PSD1sz = round(um(3)/Scale);
PSD2sz = round(um(4)/Scale);
PSA1sz = round(um(5)/Scale);
PSA2sz = round(um(6)/Scale);
SPY1sz = PSD1sz+(PSA1sz*2);
SPY2sz = PSD2sz+(PSA2sz*2);


% SAP/ACTIN MATRIX DIMS
S1sz = size(S1,1);		% Full Length of S1 Mx along 1 dimension
S2sz = size(S2,1);		% Full Length of S2 Mx along 1 dimension

S1num = numel(S1);		% Full Number of elements in the S1 Mx
S2num = numel(S2);		% Full Number of elements in the S2 Mx

S1um = um(3)/Scale + (um(5)/Scale*2); % S1 size based on synapse size inputs
S2um = um(4)/Scale + (um(6)/Scale*2); % S2 size based on synapse size inputs






%---
% Here's where it gets a little tricky. 
% When converting between real-space and matrix-space you must consider that,
% if real space extends between 1 to 100, the total distance is 99 units (100-1=99).
% Meanwhile, a  Vx=[1:100]  has a  length(Vx)==100  instead of 99. 
% Here I'm choosing to add +1 to the synaptic (S) starting position to mitigate 
% any potential positioning errors. Therefore the universal convention will be
% such that, if a real-space dendritic width of 300 is chosen, S will be positioned
% at SMx(101:200,101:200), and the real-space a particle can diffuse freely within
% will be from 1-300 (1 unit smaller than GUI-entered dendritic width).

S1xcL = round(DFszX/2 - S1sz/2)   +1;
S1yrB = round(DFszY/4*1 - S1sz/2) +1;
S1xcR = S1xcL-1 + S1sz;
S1yrT = S1yrB-1 + S1sz;

S2xcL = round(DFszX/2 - S1sz/2)   +1;
S2yrB = round(DFszY/4*3 - S1sz/2) +1;
S2xcR = S2xcL-1 + S1sz;
S2yrT = S2yrB-1 + S1sz;

S1L = S1xcL;
S1B = S1yrB;
S1R = S1xcR;
S1T = S1yrT;
S2L = S2xcL;
S2B = S2yrB;
S2R = S2xcR;
S2T = S2yrT;

S1C = [((S1L + S1R) / 2); ((S1T + S1B) / 2)];
S2C = [((S2L + S2R) / 2); ((S2T + S2B) / 2)];
S1rad = ((S1R - S1L) / 2);
S2rad = ((S2R - S2L) / 2);


PS1sz = round(S1sz/2);
PS2sz = round(S2sz/2);
PS1xy = round(S1sz/4);
PS2xy = round(S2sz/4);
PS1L = S1L + PS1xy;
PS1B = S1B + PS1xy;
PS1R = S1R - PS1xy;
PS1T = S1T - PS1xy;
PS2L = S2L + PS2xy;
PS2B = S2B + PS2xy;
PS2R = S2R - PS2xy;
PS2T = S2T - PS2xy;



SMxE = zeros(DFszY,DFszX);
SMx = SMxE;
SMx(S1B:S1T,S1L:S1R) = S1;
SMx(S2B:S2T,S2L:S2R) = S2;


%-----------------------
% Function Return Cells & Structure Arrays
%-----------------------
% DFum = {DFszX,DFszY,PSD1sz,PSD2sz,PSA1sz,PSA2sz,SPY1sz,SPY2sz};
SSsz = {S1sz,S2sz,S1num,S2num,S1um,S2um};



DFum = struct;
%---
DFum.X = DFszX;
DFum.Y = DFszY;
DFum.PSD1 = PSD1sz;
DFum.PSD2 = PSD2sz;
DFum.PSA1 = PSA1sz;
DFum.PSA2 = PSA2sz;
DFum.SPY1 = SPY1sz;
DFum.SPY2 = SPY2sz;

% Matrix Space
SS = struct;
%---
SS.SMx = SMx;
SS.SMxE = SMxE;
SS.S1 = S2;
SS.S1 = S2;
SS.S1L = S1L;
SS.S1B = S1B;
SS.S1R = S1R;
SS.S1T = S1T;
SS.S2L = S2L;
SS.S2B = S2B;
SS.S2R = S2R;
SS.S2T = S2T;
SS.S1LB = [SS.S1L SS.S1B];
SS.S1RT = [SS.S1R SS.S1T];
SS.S2LB = [SS.S2L SS.S2B];
SS.S2RT = [SS.S2R SS.S2T];
SS.S1sz = S1sz;
SS.S2sz = S2sz;
SS.S1num = S1num;
SS.S2num = S2num;
SS.S1um = S1um;
SS.S2um = S2um;
SS.PS1L = PS1L;
SS.PS1B = PS1B;
SS.PS1R = PS1R;
SS.PS1T = PS1T;
SS.PS2L = PS2L;
SS.PS2B = PS2B;
SS.PS2R = PS2R;
SS.PS2T = PS2T;
SS.PS1sz = PS1sz;
SS.PS2sz = PS2sz;
SS.S1C = S1C;
SS.S2C = S2C;
SS.S1rad = S1rad;
SS.S2rad = S2rad;

% Real Space
DF = struct;
%---
DF.X = size(SMx,2);
DF.Y = size(SMx,1);
DF.S1sz = S1sz;
DF.S2sz = S2sz;
DF.S1L = S2L;
DF.S1B = S2B;
DF.S1R = S2R;
DF.S1T = S2T;
DF.S2L = S1L;
DF.S2B = S1B;
DF.S2R = S1R;
DF.S2T = S1T;
DF.S1LB = [DF.S1L DF.S1B];
DF.S1RT = [DF.S1R DF.S1T];
DF.S2LB = [DF.S2L DF.S2B];
DF.S2RT = [DF.S2R DF.S2T];
DF.SMxE = SMxE;
DF.SMx = SMx;
DF.S1 = S1;
DF.S2 = S2;

DF.PS1L = PS2L;
DF.PS1B = PS2B;
DF.PS1R = PS2R;
DF.PS1T = PS2T;
DF.PS2L = PS1L;
DF.PS2B = PS1B;
DF.PS2R = PS1R;
DF.PS2T = PS1T;
DF.PS1LB = [DF.PS1L DF.PS1B];
DF.PS1RT = [DF.PS1R DF.PS1T];
DF.PS2LB = [DF.PS2L DF.PS2B];
DF.PS2RT = [DF.PS2R DF.PS2T];
DF.PS1sz = PS1sz;
DF.PS2sz = PS2sz;
DF.PS1C = S2C;
DF.PS2C = S1C;
DF.PS1rad = S2rad;
DF.PS2rad = S1rad;



%------------
%{


S1Hsz = size(S1,1)/2;	% Half Length of S1 Mx along 1 dimension
S2Hsz = size(S2,1)/2;	% Half Length of S1 Mx along 1 dimension

S1Hnum = numel(S1)/2;	% Half Number of elements in the S1 Mx
S2Hnum = numel(S2)/2;	% Half Number of elements in the S1 Mx

% SS.S1sz = S1sz;
% SS.S2sz = S2sz;
% SS.S1numel = S1num;
% SS.S2numel = S2num;
% SS.S1Hsz = S1Hsz;
% SS.S2Hsz = S2Hsz;
% SS.S1Hnumel = S1Hnum;
% SS.S2Hnumel = S2Hnum;
% SS.S1um = S1um;
% SS.S2um = S2um;

% % ZMx1 = S1um;
% % ZMx2 = S2um;
% ZMx1=S1Hsz;
% ZMx2=S2Hsz;
% 
% oVx1 = round(.5:.5:ZMx1);
% oVx1 = [oVx1;oVx1];
% zx = 0;
% for ox = 1:2:ZMx1*2
% 	rVx1(ox:ox+1,:) = oVx1 + (ZMx1*zx);
% 	zx = zx+1;
% end
% 
% oVx2 = round(.5:.5:ZMx2);
% oVx2 = [oVx2;oVx2];
% zx = 0;
% for ox = 1:2:ZMx2*2
% 	rVx2(ox:ox+1,:) = oVx2 + (ZMx2*zx);
% 	zx = zx+1;
% end



%---
% Dot Matrix Preallocate
% PG1xy = zeros(DFszY);
% PG1xy((DFszX+1):end,:) = [];
% PG1xy=PG1xy';
% PG2xy = PG1xy;
% SMx=PG1xy;
% SMx(401:500,101:200) = S1;
% SMx(101:200,101:200) = S2;
%------


%}
%------------
% OLD
%{

PSD1size=PSD1sz;
PSD2size=PSA1sz;
periPSD1size=PSD2sz;
periPSD2size=PSA2sz;


Syn1Size = PSD1size+(periPSD1size*2);
Syn2Size = PSD2size+(periPSD2size*2);

PSD1padX = round((fsizeX - Syn1Size)/2);
PSD1padY = round(((fsizeY - Syn1Size)/2)/2);
PSD2padX = round((fsizeX - Syn2Size)/2);
PSD2padY = round(((fsizeY - Syn2Size)/2)/2);


fPSD1 = ones(PSD1size+1);						% PSD1 SIZE
fPSD2 = ones(PSD2size+1);						% PSD2 SIZE
fPSD1 = padarray(fPSD1,[periPSD1size periPSD1size], 2);	% PSD1 SIZE
fPSD2 = padarray(fPSD2,[periPSD2size periPSD2size], 2);	% PSD2 SIZE
pfPSD1 = padarray(fPSD1,[PSD1padY PSD1padX], 0);	% PAD PSD1 [Y-rows X-cols]
pfPSD2 = padarray(fPSD2,[PSD2padY PSD2padX], 0);	% PAD PSD2 [Y-rows X-cols]

[pfPSD1R pfPSD1C] = size(pfPSD1);
[pfPSD2R pfPSD2C] = size(pfPSD2);


if pfPSD1C > pfPSD2C
		pfPSD2 = padarray(pfPSD2,[0 (pfPSD1C - pfPSD2C)], 0,'post');
elseif pfPSD1C < pfPSD2C
		pfPSD1 = padarray(pfPSD1,[0 (pfPSD2C - pfPSD1C)], 0,'post');
end


%=========================================================%
PSDfield = cat(1, pfPSD1, pfPSD2);				% CONCAT PSD FIELDS
[Yrows Xcols] = size(PSDfield);					% Even value, total field
[YrPr1top XcPr1lft] = find(pfPSD1,1,'first');	% PSD1 1st row & col
[YrPr1bot XcPr1rit] = find(pfPSD1,1,'last');	% PSD1 Lst row & col
[YrPr2top XcPr2lft] = find(pfPSD2,1,'first');	% PSD2 1st row & col
[YrPr2bot XcPr2rit] = find(pfPSD2,1,'last');	% PSD2 Lst row & col
YrowsHIGH = Yrows/2;							% Rows high from Y=0
XcolsWIDE = Xcols;								% Cols wide from X=0
YHIGH = YrowsHIGH;								% Rows high from Y=0
XWIDE = XcolsWIDE;								% Cols wide from X=0
%=========================================================%
Y1RowsUpTop = (YrowsHIGH-YrPr1top);
Y1RowsUpBot = (YrowsHIGH-YrPr1bot);
X1ColsOvrL = XcPr1lft;
X1ColsOvrR = XcPr1rit;

Y2RowsDnTop = (YrowsHIGH-YrPr2bot)*-1;
Y2RowsDnBot = (YrowsHIGH-YrPr2top)*-1;
X2ColsOvrL = XcPr2lft;
X2ColsOvrR = XcPr2rit;
%=========================================================%
PSD1WH = PSD1size;
PSD2WH = PSD2size;
PERI1WH = periPSD1size;
PERI2WH = periPSD2size;
SPYN1WH = PSD1WH+(PERI1WH*2);
SPYN2WH = PSD2WH+(PERI2WH*2);
%=========================================================%
XYLTpr1 = [X1ColsOvrL Y1RowsUpTop];
XYRTpr1 = [X1ColsOvrR Y1RowsUpTop];
XYLBpr1 = [X1ColsOvrL Y1RowsUpBot];
XYRBpr1 = [X1ColsOvrR Y1RowsUpBot];
XYLTpr2 = [X2ColsOvrL Y2RowsDnTop];
XYRTpr2 = [X2ColsOvrR Y2RowsDnTop];
XYLBpr2 = [X2ColsOvrL Y2RowsDnBot];
XYRBpr2 = [X2ColsOvrR Y2RowsDnBot];
XYLTp1 = [(XYLTpr1(1,1)+PERI1WH) (XYLTpr1(1,2)-PERI1WH)];
XYRTp1 = [(XYRTpr1(1,1)-PERI1WH) (XYRTpr1(1,2)-PERI1WH)];
XYLBp1 = [(XYLBpr1(1,1)+PERI1WH) (XYLBpr1(1,2)+PERI1WH)];
XYRBp1 = [(XYRBpr1(1,1)-PERI1WH) (XYRBpr1(1,2)+PERI1WH)];
XYLTp2 = [(XYLTpr2(1,1)+PERI2WH) (XYLTpr2(1,2)-PERI2WH)];
XYRTp2 = [(XYRTpr2(1,1)-PERI2WH) (XYRTpr2(1,2)-PERI2WH)];
XYLBp2 = [(XYLBpr2(1,1)+PERI2WH) (XYLBpr2(1,2)+PERI2WH)];
XYRBp2 = [(XYRBpr2(1,1)-PERI2WH) (XYRBpr2(1,2)+PERI2WH)];
%=========================================================%
% rectangle('Position',[x,y,w,h])
XYBOXpr1 = [XYLBpr1 SPYN1WH SPYN1WH];
XYBOXpr2 = [XYLBpr2 SPYN2WH SPYN2WH];
XYBOXp1 = [XYLBp1 PSD1WH PSD1WH];
XYBOXp2 = [XYLBp2 PSD2WH PSD2WH];
%=========================================================%
Zfield = changem(PSDfield,[1 2 3],[0 2 1]);
%=========================================================%



%%
%=========================================================%
if doRun(8)
%=========================================================%
xlim = [0 XcolsWIDE];
ylim = [-YrowsHIGH YrowsHIGH];
%---
scsz = get(0,'ScreenSize'); 
scsx=scsz(3); scsy=scsz(4);

Fh96 = figure(96);
set(gcf,'OuterPosition',[145,145,scsx/2,scsy/1.2])

axes('Position',[.25 .05 .7 .9]);
scatter(1:4,1:4,5,[1 0 0]);
axis([xlim, ylim]);
grid on
set(get(gca,'XLabel'),'String','GluR2:Blue \bullet \bullet GluR1:Red')
rectangle('Position',[XYBOXpr1(1),XYBOXpr1(2),XYBOXpr1(3),XYBOXpr1(4)])
rectangle('Position',[XYBOXpr2(1),XYBOXpr2(2),XYBOXpr2(3),XYBOXpr2(4)])
rectangle('Position',[XYBOXp1(1),XYBOXp1(2),XYBOXp1(3),XYBOXp1(4)])
rectangle('Position',[XYBOXp2(1),XYBOXp2(2),XYBOXp2(3),XYBOXp2(4)])
%---
text(XYLTpr1(1),XYLTpr1(2),...
strcat(num2str(XYLTpr1(1)), '\bullet',num2str(XYLTpr1(2)),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(XYRTpr1(1),XYRTpr1(2),...
strcat('\leftarrow ', num2str(XYRTpr1(1)), '\bullet',num2str(XYRTpr1(2))),...
'FontSize',12,'HorizontalAlignment','left');
text(XYLBpr1(1),XYLBpr1(2),...
strcat(num2str(XYLBpr1(1)), '\bullet',num2str(XYLBpr1(2)),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(XYRBpr1(1),XYRBpr1(2),...
strcat('\leftarrow ', num2str(XYRBpr1(1)), '\bullet',num2str(XYRBpr1(2))),...
'FontSize',12,'HorizontalAlignment','left');
%---
text(XYLTpr2(1),XYLTpr2(2),...
strcat(num2str(XYLTpr2(1)), '\bullet',num2str(XYLTpr2(2)),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(XYRTpr2(1),XYRTpr2(2),...
strcat('\leftarrow ', num2str(XYRTpr2(1)), '\bullet',num2str(XYRTpr2(2))),...
'FontSize',12,'HorizontalAlignment','left');
text(XYLBpr2(1),XYLBpr2(2),...
strcat(num2str(XYLBpr2(1)), '\bullet',num2str(XYLBpr2(2)),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(XYRBpr2(1),XYRBpr2(2),...
strcat('\leftarrow ', num2str(XYRBpr2(1)), '\bullet',num2str(XYRBpr2(2))),...
'FontSize',12,'HorizontalAlignment','left');
%---
text(XYLTp1(1),XYLTp1(2),...
strcat(num2str(XYLTp1(1)), '\bullet',num2str(XYLTp1(2)),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(XYRTp1(1),XYRTp1(2),...
strcat('\leftarrow ', num2str(XYRTp1(1)), '\bullet',num2str(XYRTp1(2))),...
'FontSize',12,'HorizontalAlignment','left');
text(XYLBp1(1),XYLBp1(2),...
strcat(num2str(XYLBp1(1)), '\bullet',num2str(XYLBp1(2)),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(XYRBp1(1),XYRBp1(2),...
strcat('\leftarrow ', num2str(XYRBp1(1)), '\bullet',num2str(XYRBp1(2))),...
'FontSize',12,'HorizontalAlignment','left');
%---
text(XYLTp2(1),XYLTp2(2),...
strcat(num2str(XYLTp2(1)), '\bullet',num2str(XYLTp2(2)),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(XYRTp2(1),XYRTp2(2),...
strcat('\leftarrow ', num2str(XYRTp2(1)), '\bullet',num2str(XYRTp2(2))),...
'FontSize',12,'HorizontalAlignment','left');
text(XYLBp2(1),XYLBp2(2),...
strcat(num2str(XYLBp2(1)), '\bullet',num2str(XYLBp2(2)),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(XYRBp2(1),XYRBp2(2),...
strcat('\leftarrow ', num2str(XYRBp2(1)), '\bullet',num2str(XYRBp2(2))),...
'FontSize',12,'HorizontalAlignment','left');
%---
text(XcolsWIDE,0,...
strcat(num2str(XcolsWIDE), '\bullet',num2str(0),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(0,0,...
strcat('\leftarrow ', num2str(0), '\bullet',num2str(0)),...
'FontSize',12,'HorizontalAlignment','left');
text(XcolsWIDE,YrowsHIGH,...
strcat(num2str(XcolsWIDE), '\bullet',num2str(YrowsHIGH),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(0,-YrowsHIGH,...
strcat('\leftarrow ', num2str(0), '\bullet',num2str(-YrowsHIGH)),...
'FontSize',12,'HorizontalAlignment','left');
%---

axes('Position',[.05 .5 .15 .4]);
imagesc(Zfield);     % <--FIG--##
colormap('bone')
title('PSD1');
axes('Position',[.05 .05 .15 .4]); 
imagesc(PSDfield); % <--FIG--##
title('PSD2');


%%
%=========================================================%
end %if doRun(8)
%=========================================================%

SPYN1xv = [XYLTpr1(1) XYRTpr1(1) XYRBpr1(1) XYLBpr1(1) XYLTpr1(1)]';
SPYN1yv = [XYLTpr1(2) XYRTpr1(2) XYRBpr1(2) XYLBpr1(2) XYLTpr1(2)]';
SPYN2xv = [XYLTpr2(1) XYRTpr2(1) XYRBpr2(1) XYLBpr2(1) XYLTpr2(1)]';
SPYN2yv = [XYLTpr2(2) XYRTpr2(2) XYRBpr2(2) XYLBpr2(2) XYLTpr2(2)]';

PSD1xv = [XYLTp1(1) XYRTp1(1) XYRBp1(1) XYLBp1(1) XYLTp1(1)]';
PSD1yv = [XYLTp1(2) XYRTp1(2) XYRBp1(2) XYLBp1(2) XYLTp1(2)]';
PSD2xv = [XYLTp2(1) XYRTp2(1) XYRBp2(1) XYLBp2(1) XYLTp2(1)]';
PSD2yv = [XYLTp2(2) XYRTp2(2) XYRBp2(2) XYLBp2(2) XYLTp2(2)]';
 
[PERI1xv, PERI1yv] = polybool('xor', PSD1xv, PSD1yv, SPYN1xv, SPYN1yv);
[PERI2xv, PERI2yv] = polybool('xor', PSD2xv, PSD2yv, SPYN2xv, SPYN2yv);
%}
%------------


end
%----------------------------------%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%					PLOTTING AND LIVE SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





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




%----------------------------------%
%			MAINPLOT
%----------------------------------%
% FieldFun
% short description of this function
%----------------------------------%
function [] = MAINPLOT(G1Ph1,G2Ph1,xyG2,xyG1)
%-------------------------------------------%

set(G2Ph1,'XData',xyG2(1,:),'YData',xyG2(2,:));
%drawnow

set(G1Ph1,'XData',xyG1(1,:),'YData',xyG1(2,:));
drawnow


end
%----------------------------------%



%----------------------------------%
%			ActinPlot
%----------------------------------%
function [varargout] = ActinMainPlot(Fh,Actin,inPSD,spos,varargin)


scsz = get(0,'ScreenSize');
basepos = scsz./[2.5e-3 2.5e-3 1.5 2];
baseazel = [-32 12];
baserot=[0 0];
% spos = [.03 .53 .27 .44];

if nargin == 7
	rot=varargin{1};
	azel=varargin{2};
	dims=varargin{3};
elseif nargin == 6 
	rot=varargin{1};
	azel=varargin{2};
	dims = [25   300   200    50    50    25   275];
elseif nargin == 5 
	rot=varargin{1};
	azel = baseazel;
	dims = [25   300   200    50    50    25   275];
else
	rot=baserot;
	azel = baseazel;
	dims = [25   300   200    50    50    25   275];
end

AxLims = [-dims(4) dims(4) -dims(5) dims(5) 0 dims(2)].*1.2;

%--------------------
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zin,Zic] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zin,:);
[Zout,Zoc] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zout,:);

% [Zin,Zic] = find(ActinTips(:,3) > SPYz);
% XYsq = sqrt(ActinTips(:,1).^2 + ActinTips(:,2).^2);
% [XYin,XYic] = find(XYsq < (HX-SPYxy));
% XYZin = intersect(XYin, Zin);
% % XYt = intersect(Xrow1, Yrow1);
% %[tf, loc] = ismember(Xrow1, Yrow1)
%--------------------
figure(Fh)
hold on;
subplot('Position',spos), 
ph11c = plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]');
	axis(AxLims); 
	%axis vis3d;
	%xlabel('X');ylabel('Y');zlabel('Z');
	set(gca,'xticklabel',[])
	set(gca,'yticklabel',[])
	set(gca,'zticklabel',[])
	view(azel); grid off;
	hold on;
ph11a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
	hold on;
ph11b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
	axis(AxLims)
	%xlabel('X');ylabel('Y');zlabel('Z');
	view(azel+rot)
	grid off
	set(ph11a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
	set(ph11b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
	set(ph11c,'LineStyle','-','Color',[.7 .7 .7],'LineWidth',.1);
	%hold off;
	set(gca,'xticklabel',[])
	set(gca,'yticklabel',[])
	set(gca,'zticklabel',[])
%--------------------
% subplot('Position',[.13 .04 .27 .44]), 
% ph12a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
% 	hold on;
% ph12b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
% 	axis(AxLims)
% 	view([0 90])
% 	grid off
% 	set(gca,'Color',[1,1,1])
% 	set(ph12a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
% 	set(ph12b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
% 	%--------------------
% 	%set(gca,'XTickLabel', sprintf('%.1f|',nT),'FontSize',10)
% 	%hold off;
% %--------------------




varargout = {Fh};
end
%----------------------------------%





%----------------------------------%
%			ActinPlot
%----------------------------------%
function [varargout] = ActinPlot(Fh,nT,Actin,inPSD,varargin)


scsz = get(0,'ScreenSize');
basepos = scsz./[2.5e-3 2.5e-3 1.5 2];
baseazel = [-32 12];
baserot=[0 0];

if nargin >= 3
	rot=varargin{1};
	azel=varargin{2};
	dims=varargin{3};
elseif nargin == 2 
	rot=varargin{1};
	azel=varargin{2};
	dims = [60   400   280   160   160    40   360];
elseif nargin == 1 
	rot=varargin{1};
	azel = baseazel;
	dims = [60   400   280   160   160    40   360];
else
	rot=baserot;
	azel = baseazel;
	dims = [60   400   280   160   160    40   360];
end

AxLims = [-dims(4) dims(4) -dims(5) dims(5) 0 dims(2)].*1.2;

%--------------------
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%--------------------
figure(Fh)
subplot('Position',[.03 .03 .65 .95]), 
ph11c = plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]');
axis(AxLims); axis vis3d;
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
view(azel)
grid off
set(gca,'Color',[1,1,1])
hold on;
ph11a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph11b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis(AxLims)
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
view(azel+rot)
grid off
set(gca,'Color',[1,1,1])
set(ph11a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph11b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
set(ph11c,'LineStyle','-','Color',[.7 .7 .7],'LineWidth',.1);
hold off;
%--------------------
figure(Fh)
subplot('Position',[.7 .55 .28 .38]), 
ph12a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph12b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis(AxLims)
view([0 90])
grid off
set(gca,'Color',[1,1,1])
set(ph12a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph12b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
%--------------------
set(gca,'XTickLabel', sprintf('%.1f|',nT),'FontSize',10)
hold off;
%--------------------


%{
%==================================================%
if nargin >= 3
dims=varargin{3};
%-----------------------------------%
% Spine Dimensions
SPYneckXY = dims(1);
SPYheadZN = dims(2);
SPYheadZS = dims(3);
SPYheadX = dims(4);
SPYheadY = dims(5);
PSDproxy = dims(6);
inPSD = dims(7);
%-----------------------------------%
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%-----------------------------------%
PSDXYZ = [PSDTips(:,1) PSDTips(:,2) PSDTips(:,3)];
PSDXY = round([PSDTips(:,1) PSDTips(:,2)]);
PSDactMx = zeros(SPYheadY+100,SPYheadX+100);
for mxp = 1:numel(PSDXY(:,1))
PSDactMx(PSDXY(mxp,2)+SPYheadY+10, PSDXY(mxp,1)+SPYheadX+10) = 1;
end
ActMask=[1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1];
ActMx = convn(PSDactMx,ActMask,'same');
ActMx = (ActMx>0).*1.0;
%===================================%
%				FIGURE
%-----------------------------------%
figure(Fh)
subplot('Position',[.6 .1 .38 .38]), 
imagesc(ActMx)
colormap(bone)
hold on
% subplot('Position',[.55 .05 .40 .40]), 
subplot('Position',[.6 .1 .38 .38]), 
scatter(PSDXY(:,1),PSDXY(:,2), 'r')
hold off
%-----------------------------------%

end
%==================================================%
if nT >40000; keyboard; end
%}

varargout = {Fh};
end
%----------------------------------%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%					PARTICLE COUNTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-------------------------------%
%       MAIN PARTICLE COUNTER
%-------------------------------%
function [GR1SPY1N GR1SPY2N GR1PER1N GR1PER2N GR1PSD1N GR1PSD2N...
		GR2SPY1N GR2SPY2N GR2PER1N GR2PER2N GR2PSD1N GR2PSD2N...
		GR1PERN GR1PSDN GR2PERN GR2PSDN GR1SPYN GR2SPYN...
		SPY1N SPY2N SPYNT ESNT EST G1INSP G2INSP GINSP...
		G1xyP1 G1xyP2 G2xyP1 G2xyP2]...
		= dotCount(stepN,G2N, G1N,...
		G1P1r, G1P2r, G2P1r, G2P2r, G1inP1, G1inP2, G2inP1, G2inP2,...
		S1G1Cv,S2G1Cv,S1G2Cv,S2G2Cv,G1INSPA,G2INSPA,xyG1,xyG2)


%================%
GR1SPY1N = sum(G1inP1);	% GR1PSD1
GR1SPY2N = sum(G1inP2);	% GR1PSD2
GR1PER1N = sum(G1inP1);	% GR1PS1T
GR1PER2N = sum(G1inP2);	% GR1PS2T
GR1PSD1N = sum(G1inP1);	% GR1PSD1T
GR1PSD2N = sum(G1inP2);	% GR1PSD2T
%==========%
GR2SPY1N = sum(G2inP1);	% GR2SPY1N
GR2SPY2N = sum(G2inP2);	% GR2SPY2N
GR2PER1N = sum(G2inP1);	% GR2PER1N
GR2PER2N = sum(G2inP2);	% GR2PER2N
GR2PSD1N = sum(G2inP1);	% GR2PSD1N
GR2PSD2N = sum(G2inP2);	% GR2PSD2N
%================%

GR1PERN = GR1PER1N + GR1PER2N;	% GR1PERN
GR1PSDN = GR1PSD1N + GR1PSD2N;	% GR1PSDN
GR2PERN = GR2PER1N + GR2PER2N;	% GR2PERN
GR2PSDN = GR2PSD1N + GR2PSD2N;	% GR2PSDN

GR1SPYN = GR1SPY1N + GR1SPY2N;	% GR1SPYN
GR2SPYN = GR2SPY1N + GR2SPY2N;	% GR2SPYN

SPY1N = GR1SPY1N + GR2SPY1N;	% SPY1N SPY1N
SPY2N = GR1SPY2N + GR2SPY2N;	% SPY2N SPY2N

SPYNT = SPY1N + SPY2N;			% SPYNT SPYNT
ESNT = (G2N+G1N)-(SPY1N+SPY2N);	% ESNT
EST = [ESNT SPY1N SPY2N SPYNT];	% EST

G1INSP = sum(G1INSPA);
G2INSP = sum(G2INSPA);
GINSP = G1INSP + G2INSP;




G1xyP1 = S1G1Cv(:,G1inP1);
G1xyP2 = S2G1Cv(:,G1inP2);
G2xyP1 = S1G2Cv(:,G2inP1);
G2xyP2 = S2G2Cv(:,G2inP2);



end
%================================%



%-------------------------------%
%   GRAPH PARTICLE COUNTERS
%-------------------------------%
function [] = GRAPHCOUNTS(stepN,G2N,doGluR1,doGluR2,doRun1,doRun2,doRun6,...
	HSlow, HShi, G1FSLOTS,G2FSLOTS,SAP5,S1sum,S2sum,doKo,...
	EST,SPY1N,SPY2N)



doPLOT2D = doRun1;
doPLOT3D = doRun2;
doHOMEO = doRun6;

if G2N == 0
	G2N = 50;
end


%===============================%
if doPLOT2D %MAIN 2D PLOT IS ON
%===============================%
	
%===============================%
% Homeostatic Bar Graph
%-------------------------------%
if doPLOT3D == 0 && doHOMEO ==1
if mod(stepN, 100) == 0
PSDHSn = [SPY1N (SPY1N+SPY2N) SPY2N];
% if stepN == 10
HSline = [HSlow HSlow HSlow HShi HShi HShi];
HSshape = [0 4 0 0 4 4]; 
figure(1);
subplot(5,5,[1 7]);
plot(HSshape,HSline,'LineWidth',3)
hold on
subplot(5,5,[1 7]);
bar(1:3, [PSDHSn']);
set(gca,'Ylim',[0 (HShi+20)])
hold off
end
end
%===============================%
	

%===============================%
% Live Particle Totals Bar Graph
%-------------------------------%
if doPLOT3D == 0 && doHOMEO == 0
	if mod(stepN, 100) == 0
	figure(1)
	subplot(5,5,[1 7]);
	bar(1:4, [EST'], .8, 'EdgeColor',[1 0.5 0.5]);
	set(gca,'Ylim',[0 G2N],'XLim',[.5 4.5])
	title('AMPAR#:   EC    -   SPINE1   -     SPINE2  -    SPINES');
	end
end

if doPLOT3D == 1
	if mod(stepN, 100) == 0
	figure(1)
	subplot(5,5,[21 22]);
	bar(1:4, [EST'], .8, 'EdgeColor',[1 0.5 0.5]);
	set(gca,'Ylim',[0 G2N],'XLim',[.5 4.5])
	title('AMPAR#:   EC    -   SPINE1   -     SPINE2  -    SPINES');
	end
end
%===============================%

	

%===============================%
if mod(stepN, 50) == 0
G1PSD1slotN = fix(S1sum/4*SAP5(1))*doKo(3);
G1PSD2slotN = fix(S2sum/4*SAP5(1))*doKo(3);
G1PERI1slotN = fix(S1sum/4*SAP5(2))*doKo(3);
G1PERI2slotN = fix(S2sum/4*SAP5(2))*doKo(3);
G2PSD1slotN = fix(S1sum/4*SAP5(1))*doKo(4);
G2PSD2slotN = fix(S2sum/4*SAP5(1))*doKo(4);
G2PERI1slotN = fix(S1sum/4*SAP5(2))*doKo(4);
G2PERI2slotN = fix(S2sum/4*SAP5(2))*doKo(4);
SLOTMX = [G1PSD1slotN G1PSD2slotN G1PERI1slotN G1PERI2slotN;
	G2PSD1slotN G2PSD2slotN G2PERI1slotN G2PERI2slotN];
SLOTMAX = max(max(SLOTMX));
if SLOTMAX < 11
	ymax = 10;
elseif SLOTMAX >= 11 && SLOTMAX < 21
	ymax = 20;
else
	ymax = SLOTMAX;
end
if doGluR1	
	subplot(5,5,[16 17]), bar(1:4, [G1FSLOTS'], .8, 'EdgeColor',[1 0.5 0.5]);
		set(gca,'Ylim',[0 ymax],'XLim',[.5 4.5])
        title('GLUR1 IN SLOTS');
		xt = [G1PSD1slotN G1PSD2slotN G1PERI1slotN G1PERI2slotN];
		set(gca,'XTickLabel', sprintf('%i|',xt))
		%set(get(gca,'XLabel'),'String','\bullet \bullet  PSD  \bullet  \bullet PERI')
end
if doGluR2
	subplot(5,5,[21 22]), bar(1:4, [G2FSLOTS'], .8, 'EdgeColor',[1 0.5 0.5]);
		set(gca,'Ylim',[0 ymax],'XLim',[.5 4.5])
        title('GLUR2 IN SLOTS');
		xt = [G2PSD1slotN G2PSD2slotN G2PERI1slotN G2PERI2slotN];
		set(gca,'XTickLabel', sprintf('%i|',xt))
		set(get(gca,'XLabel'),'String',...
			strcat('PSD \bullet \bullet \bullet \bullet \bullet \bullet',...
			'\bullet \bullet \bullet \bullet \bullet \bullet PERI'))
end
end
%===============================%
%{
text(0,-5, strcat(['slots: '...
		num2str(G1PSD1slotN),' \bullet ' num2str(G1PSD2slotN) ' \bullet \bullet '...
		' \bullet ' num2str(G1PERI1slotN),' \bullet ' num2str(G1PERI2slotN)]),...
		'FontSize',12, 'HorizontalAlignment','left',... 
		'BackgroundColor',[.8 .8 .8]);

text(0,-5, strcat(['slots: '...
num2str(G2PSD1slotN),' \bullet ' num2str(G2PSD2slotN) ' \bullet \bullet '...
' \bullet ' num2str(G2PERI1slotN),' \bullet ' num2str(G2PERI2slotN)]),...
'FontSize',12, 'HorizontalAlignment','left',... 
'BackgroundColor',[.8 .8 .8]);
%}

%===============================%
end %if doRun1 %MAIN 2D PLOT IS ON
%===============================%

% if stepN == 200; keyboard; end
	
	
end
%-------------------------------%










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%				   IzhFun NEURAL SPKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%-------------------------------------------%
% IzhFun SPIKE GENERATOR
%-------------------------------------------%
function [] = IzhFun(spikes)
global v
izhikevich(spikes)
end
%-------------------------------------------%











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%				MSD BROWNIAN MOTION ANALYSIS TOOLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------%
% MSD PARTICLE MOVEMENT GENERATION
%-------------------------------------------%
function [xyG2 xysG2] = MSDAMPARSTEP(G2N, xysG2, xyG2,... 
	testPSD1MSD, testPSD2MSD, testESMSD, PSD1, PSD2, D, Dn_PSD1, Dn_PSD2)
	
	for j = 1:G2N
        
		if testPSD1MSD
           xysG2(:,j) = xysG2(:,j)*PSD1;
        elseif testPSD2MSD
            xysG2(:,j) = xysG2(:,j)*PSD2;
		elseif testESMSD
            xysG2(:,j) = xysG2(:,j);
		end
		
	   
        xyG2(:,j) = xyG2(:,j)+xysG2(:,j);
	end

end

%-------------------------------------------%
% MEAN SQUARED DISPLACEMENT FUNCTION
%-------------------------------------------%
function [tracks] = MSDfun(stepN, Nsteps, tracks, xysG2)
    time = (0:Nsteps-1)';
    xymsd = xysG2';
    xymsd = cumsum(xymsd,1);
    tracks{stepN} = [time xymsd];

end


%-------------------------------------------%
% MSDfunction - FINAL TRACKS ANALYSIS
%-------------------------------------------%
function [MSDout] = MSDfunction(tracks, G2N, Nsteps, D, Dn_PSD1, Dn_PSD2, dT, k, MSDtest)

SPACE_UNITS = 'µm';
TIME_UNITS = 's';
N_PARTICLES = G2N;
N_TIME_STEPS = Nsteps;
N_DIM = 2;

oDes = D;				% raw		µm^2/s
D  = D*.1;              % to-scale	µm^2/s

oDpsd1 = Dn_PSD1;		% raw		µm^2/s
Dpsd1 = Dn_PSD1*.1;		% to-scale	µm^2/s

oDpsd2 = Dn_PSD2;		% raw		µm^2/s
Dpsd2 = Dn_PSD2*.1;		% to-scale	µm^2/s

dTbase = dT;			% raw time-step 
dT = dT*1;				% to-scale time-step
k = k;					% stdv of step distribution

ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
ma = ma.addAll(tracks);
disp(ma)

figure
ma.plotTracks;
ma.labelPlotTracks;

ma = ma.computeMSD;
ma.msd;

t = (0 : N_TIME_STEPS)' * dT;
[T1, T2] = meshgrid(t, t);
all_delays = unique( abs(T1 - T2) );

figure
ma.plotMSD;


cla
ma.plotMeanMSD(gca, true)

mmsd = ma.getMeanMSD;
t = mmsd(:,1);
x = mmsd(:,2);
dx = mmsd(:,3) ./ sqrt(mmsd(:,4));
errorbar(t, x, dx, 'k')

[fo, gof] = ma.fitMeanMSD;
plot(fo)
ma.labelPlotMSD;
legend off


ma = ma.fitMSD;

good_enough_fit = ma.lfit.r2fit > 0.8;
Dmean = mean( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;
Dstd  =  std( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;

Dheader1 = ['Raw Unscaled Values'];
Dhead1 = ['     D.es	D.psd1	D.psd2'];
Ddat1 = [oDes oDpsd1 oDpsd2];
disp(' ')
disp(Dheader1)
disp(Dhead1)
disp(Ddat1)


yourtesthead = ['YOU ARE TESTING DIFFUSION FOR:'];
if MSDtest(1)
	yourtest = ['   Des:   extrasynaptic diffusion rate'];
elseif MSDtest(2)
	yourtest = ['   Dpsd1:  PSD-1 diffusion rate'];
elseif MSDtest(3)
	yourtest = ['   Dpsd2:  PSD-2 diffusion rate'];
else
	yourtest = ['   generic diffusion rate'];
end
disp(yourtesthead)
disp(yourtest)

disp(' ')
fprintf('Estimation of raw D coefficient from MSD:\n')
fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', ...
    Dmean, Dstd, sum(good_enough_fit));




% Retrieve instantaneous velocities, per track
 trackV = ma.getVelocities;

 % Pool track data together
 TV = vertcat( trackV{:} );

 % Velocities are returned in a N x (nDim+1) array: [ T Vx Vy ...]. So the
 % velocity vector in 2D is:
 V = TV(:, 2:3);

 % Compute diffusion coefficient
varV = var(V);
mVarV = mean(varV); % Take the mean of the two estimates
Dest = mVarV / 2 * dT;



Dheader2 = ['Scaling to model (10units = 1µm)...'];
Dhead2 = ['     D.es	D.psd1	D.psd2'];
Ddat2 = [D Dpsd1 Dpsd2];

disp(' ')
disp(Dheader2)
disp(Dhead2)
disp(Ddat2)
fprintf('Estimation from velocities histogram:\n')
fprintf('Tested D = %.3g %s, compare to scaled Des value of %.3g %s\n', ...
    Dest, [SPACE_UNITS '²/' TIME_UNITS], D, [SPACE_UNITS '²/' TIME_UNITS]);

% printf('D.psd target value was %.3g %s\n', ...
%     Dest, msdDpsd, [SPACE_UNITS '²/' TIME_UNITS]);

MSDout = D;

end



%-----------------------------------------------------------%
%			  ##   MSD2 SUBFUNCTIONS    ##
%-----------------------------------------------------------% 

%-------------------------------------------%
% STEP SIZE GENERATOR
%-------------------------------------------%
function xysG2 = STEPxyds2(G2N, k)

    xysG2 = (k * randn(2,G2N));

end


%-------------------------------------------%
% MOVE PARTICLES MAIN FUNCTION
%-------------------------------------------%
function [xyG2] = AMPARSTEP2(G2N, xysG2, xyG2)
	
	for j = 1:G2N
        xyG2(:,j) = xyG2(:,j)+xysG2(:,j);
	end
	
end


%-------------------------------------------%
% LIVE DIFFUSION PLOT
%-------------------------------------------%
function [] = MAINPLOT2(xyG2, lims)
%-------------------------------------------%
xlim = [-lims lims];
ylim = [-lims lims];
zlim = [-5 5];

%=================================%
%       MAIN 2D PLOT
%---------------------------------%
figure(1)
subplot(2,1,1), 
AMPARPlot = gscatter(xyG2(1,:),xyG2(2,:));
axis([xlim, ylim]);
set(AMPARPlot,'marker','.','markersize',[6],'color',[1 0 0])


%=================================%
%           3D PLOT
%---------------------------------%
figure(1);
subplot(2,1,2), 
gscatter(xyG2(1,:),xyG2(2,:)); view(20, 30);
axis normal;
grid off
axis([xlim, ylim, zlim]);
set(gca, 'Box', 'on');

end


%-------------------------------------------%
% MANUAL STEP SIZE FUNCTION
%-------------------------------------------%
function xysG2 = UNIFORMSTEPS2(G2N, Lx)
%-------------------------------------------%

   Lx(1:2,1:G2N) = Lx;
   xyd = randi([0 1],G2N,2)';
   xyd(xyd == 0) = -1;
   xysG2 = (Lx.*xyd);
   
end


%-------------------------------------------%
% MSD SCALED STEPS FUNCTION
%-------------------------------------------%
function [xyG2 xysG2] = SCALEUNIFORMSTEPS2(G2N, xysG2, xyG2, Ls, MSDtest)


if MSDtest(1)
	Ls = 1;	
end
	
	for j = 1:G2N
        xysG2(:,j) = xysG2(:,j)*Ls;
        xyG2(:,j) = xyG2(:,j)+xysG2(:,j);
	end	
	
end


%-------------------------------------------%
% MSD TRACKS GENERATOR
%-------------------------------------------%
function [tracks] = MSDfun2(stepN, Nsteps, tracks, xysG2)
    time = (0:Nsteps-1)';
    xymsd = xysG2';
    xymsd = cumsum(xymsd,1);
    tracks{stepN} = [time xymsd];

end


%-------------------------------------------%
% MSD TRACKS ANALYSIS
%-------------------------------------------%
function [Dest] = MSDfunction2(tracks,G2N,Nsteps,D,Dn,L,dT,k,Scale,MSDtest)


% printf('D.psd target value was %.3g %s\n', ...
%     Dest, msdDpsd, [SPACE_UNITS '²/' TIME_UNITS]);

yourtesthead = ['-------TESTING DIFFUSION---------'];
if MSDtest(1)
	yourtest = ['   D:   original diffusion rate'];
elseif MSDtest(2)
	yourtest = ['   Dn:  PSD diffusion rate'];
elseif MSDtest(3)
	yourtest = ['   L:  step length'];
else
	yourtest = ['   generic diffusion rate'];
end
disp(yourtesthead)
disp(yourtest)



SPACE_UNITS = 'µm';
TIME_UNITS = 's';
N_PARTICLES = G2N;
N_TIME_STEPS = Nsteps;
N_DIM = 2;

oD = D;				% raw		µm^2/s
D  = D*Scale;       % to-scale	µm^2/s

oDn = Dn;			% raw		µm^2/s
Dn = Dn*Scale;		% to-scale	µm^2/s

oL = L;				% raw		µm
L = L*Scale;		% to-scale	µm

dTbase = dT;		% raw time-step 
dT = dT*Scale;		% to-scale time-step
k = k;				% stdv of step distribution

ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
ma = ma.addAll(tracks);
disp(ma)

figure
ma.plotTracks;
ma.labelPlotTracks;

ma = ma.computeMSD;
ma.msd;

t = (0 : N_TIME_STEPS)' * dT;
[T1, T2] = meshgrid(t, t);
all_delays = unique( abs(T1 - T2) );

figure
ma.plotMSD;


cla
ma.plotMeanMSD(gca, true)

mmsd = ma.getMeanMSD;
t = mmsd(:,1);
x = mmsd(:,2);
dx = mmsd(:,3) ./ sqrt(mmsd(:,4));
errorbar(t, x, dx, 'k')

[fo, gof] = ma.fitMeanMSD;
plot(fo)
ma.labelPlotMSD;
legend off


ma = ma.fitMSD;

good_enough_fit = ma.lfit.r2fit > 0.8;
Dmean = mean( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;
Dstd  =  std( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;

Dheader1 = ['Raw Unscaled Values'];
Dhead1 = ['    D        Dn        L'];
Ddat1 = [oD oDn oL];
disp(' ')
disp(Dheader1)
disp(Dhead1)
disp(Ddat1)


yourtesthead = ['YOU ARE TESTING DIFFUSION FOR:'];
if MSDtest(1)
	yourtest = ['   D:   original diffusion rate'];
elseif MSDtest(2)
	yourtest = ['   Dn:  PSD diffusion rate'];
elseif MSDtest(3)
	yourtest = ['   L:  step length'];
else
	yourtest = ['   generic diffusion rate'];
end
disp(yourtesthead)
disp(yourtest)

disp(' ')
fprintf('Estimation of raw D coefficient from MSD:\n')
fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', ...
    Dmean, Dstd, sum(good_enough_fit));




% Retrieve instantaneous velocities, per track
 trackV = ma.getVelocities;

 % Pool track data together
 TV = vertcat( trackV{:} );

 % Velocities are returned in a N x (nDim+1) array: [ T Vx Vy ...]. So the
 % velocity vector in 2D is:
 V = TV(:, 2:3);

 % Compute diffusion coefficient
varV = var(V);
mVarV = mean(varV); % Take the mean of the two estimates
Dest = mVarV / 2 * dT;



Dheader2 = ['Scaling to model...'];
Dhead2 = ['    D        Dn        L'];
Ddat2 = [D Dn L];

disp(' ')
disp(Dheader2)
disp(Dhead2)
disp(Ddat2)
fprintf('Estimation from velocities histogram:\n')
fprintf('Computed scaled velocity D = %.3g %s, generated from set D = %.3g %s\n', ...
    Dest, [SPACE_UNITS '²/' TIME_UNITS], D, [SPACE_UNITS '²/' TIME_UNITS]);


end










					
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%						CODE GRAVE YARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------
% Alternative Anullus Plots:
%{
TQs = [PSDY PSDX]'-50;
TQRz = sqrt(sum(TQs)).*0+1;
CQRz = sqrt(sum(CQRall.^2)).*0+2;
figure2 = figure('XVisual','');
axes2 = axes('Parent',figure2,'PlotBoxAspectRatio',[0.555555555555556 0.555555555555556 1],...
	'OuterPosition',[0 0 1 1],...
	'DataAspectRatio',[10 10 1],...
	'CameraViewAngle',4.36580191767115);
view(axes2,[50.5 40]);
grid(axes2,'on');
hold(axes2,'all');
PhS3a = scatter3(TQs(1,:),TQs(2,:),TQRz,50, [.8 .2 .2],'filled');
set(PhS3a,'Marker','s')
hold on
PhS3b = scatter3(CQRall(1,:),CQRall(2,:),CQRz,10, [.2 .7 .0],'filled');
set(PhS3b,'Marker','o')
hold off
%}
%{

QallTG = cat(2,CQRall,TQs);
QallTGz = cat(2,CQRz,TQRz);

Qallz = cat(1,QallTG,QallTGz);
spinesurf2(Qallz(1,:)', Qallz(2,:)', Qallz(3,:)')
sftool(Qallz(1,:)', Qallz(2,:)', Qallz(3,:)')

%--------------------------------------
CQRz = sqrt(sum(CQRall.^2));
sja=0;
sjb=max(CQRz);
sjc=2;
sjd=1;
linsc = @(jx) (sjc*(1-(jx-sja)/(sjb-sja)) + sjd*((jx-sja)/(sjb-sja)));
for mm = 1:numel(CQRz)
	CQRz(mm) = linsc(CQRz(mm));
end
Qallz = cat(1,CQRall,CQRz);
spinesurf2(Qallz(1,:)', Qallz(2,:)', Qallz(3,:)')
sftool(Qallz(1,:)', Qallz(2,:)', Qallz(3,:)')
%}
%{

%-------------------------------------
%Generate values from the uniform distribution on the interval [a, b]:
% r = a + (b-a).*rand(100,1);
clear pxr pyr Qpr
for mm = 1:(numel(Qr)-1)
clear pxr pyr
	Qr1 = Qr(mm+1); Qr2 = Qr(mm); 
	for nn = 1:Qs(mm)*8e3 % <------- multiplier
		theta = 2*pi*rand;
		r = sqrt(rand*(Qr1^2-Qr2^2)+Qr2^2);
		pxr(nn) = r*cos(theta); 
		pyr(nn) = r*sin(theta);
	end
	Qpr{mm} = [pxr; pyr];
end

% figure
% scatter(Qpr{1}(1,:), Qpr{1}(2,:),3,[1 0 0],'filled')
% hold on
% scatter(Qpr{2}(1,:), Qpr{2}(2,:),3,[.8 .9 .1],'filled')
% hold on
% scatter(Qpr{3}(1,:), Qpr{3}(2,:),3,[.6 .5 0],'filled')
% hold on
% scatter(Qpr{4}(1,:), Qpr{4}(2,:),3,[.4 .9 .9],'filled')
% hold on
% scatter(Qpr{5}(1,:), Qpr{5}(2,:),3,[.2 0 1],'filled')
% hold on

%-------------------------------------
QA1 = round([Qpr{1}(1,:); Qpr{1}(2,:)]+51);
QA2 = round([Qpr{2}(1,:); Qpr{2}(2,:)]+51);
QA3 = round([Qpr{3}(1,:); Qpr{3}(2,:)]+51);
QA4 = round([Qpr{4}(1,:); Qpr{4}(2,:)]+51);
QA5 = round([Qpr{5}(1,:); Qpr{5}(2,:)]+51);
G1xy = zeros(100,100);
for xy = 1:numel(QA1(1,:));G1xy(QA1(1,xy),QA1(2,xy)) = Qs(1);end
for xy = 1:numel(QA2(1,:));G1xy(QA2(1,xy),QA2(2,xy)) = Qs(2);end
for xy = 1:numel(QA3(1,:));G1xy(QA3(1,xy),QA3(2,xy)) = Qs(3);end
for xy = 1:numel(QA4(1,:));G1xy(QA4(1,xy),QA4(2,xy)) = Qs(4);end
for xy = 1:numel(QA5(1,:));G1xy(QA5(1,xy),QA5(2,xy)) = Qs(5);end

G1xy(101:end,:) = [];
G1xy(:,101:end) = [];

GN=[4, .18, 9, .1];
ccMask = CircMask(GN);
ccMask = ccMask -.2;
G1xy(G1xy<min(Qs))=min(Qs)*.3;
G1Cv = convn(G1xy,ccMask,'same');

figure
imagesc(G1xy)
colormap(jet)
colorbar

figure
imagesc(G1Cv)
colormap(jet)
colorbar

%-------------------------------------
pxr=0; pyr=0; Qpr={0};
for mm = 1:(numel(Qr)-1)
clear pxr pyr
	Qr1 = Qr(mm+1); Qr2 = Qr(mm); 
	for nn = 1:round(20^(sqrt(Qs(mm))))
		theta = 2*pi*rand;
		r = sqrt(rand*(Qr1^2-Qr2^2)+Qr2^2);
		pxr(nn) = r*cos(theta); 
		pyr(nn) = r*sin(theta);
	end
	Qpr{mm} = [pxr; pyr];
end

Qall = [Qpr{1}(1,:) Qpr{2}(1,:) Qpr{3}(1,:) Qpr{4}(1,:) Qpr{5}(1,:);
		Qpr{1}(2,:) Qpr{2}(2,:) Qpr{3}(2,:) Qpr{4}(2,:) Qpr{5}(2,:)];

Qzz =	[(.8 + (.3).*rand(1,numel(Qpr{1}(1,:)))) ...
		(.6 + (.3).*rand(1,numel(Qpr{2}(1,:)))) ...
		(.4 + (.3).*rand(1,numel(Qpr{3}(1,:)))) ...
		(.2 + (.3).*rand(1,numel(Qpr{4}(1,:)))) ...
		(0 + (.3).*rand(1,numel(Qpr{5}(1,:))))];
Qallz = cat(1,Qall,Qzz);
Qallz(isinf(Qallz)) = 0;
Qallz(isnan(Qallz)) = 0;

spinesurf2(Qallz(1,:)', Qallz(2,:)', Qallz(3,:)')
%sftool(Qallz(1,:)', Qallz(2,:)', Qallz(3,:)')
%-------------------------------------
%}
%--------------------------------------

%{


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%					SUPPLEMENTAL GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------%
% PLOT3DS
%-------------------------------------------%
function [] = PLOT3DS(stepN,xyG2,xyG1,XWIDE,YHIGH,G1Z,G2Z,Zfield)
%-------------------------------------------%
% xlim = [-10 (XWIDE+10)];
% ylim = [(-XWIDE-10) (XWIDE+10)];
% zlim = [0 4];
%---

xlim = [0 (XWIDE)];
ylim = [(-YHIGH) (YHIGH)];
zlim = [0 4];

%{
if mod(stepN, 100) == 0 && stepN>=100
figure(1);
subplot('Position',[.02 .61 .4 .35]),...
axes('position',[.02 .61 .4 .35])
surf(Zfield),view(70, 10),
axis([xlim, ylim, zlim]),
axis tight
grid off,
set(gca, 'Box', 'on'),
set(gca,'xticklabel',[]),
set(gca,'yticklabel',[]),
set(gca,'zticklabel',[]),
hold on
end
%}

if mod(stepN, 20) == 0 && stepN>=100
ti = -XWIDE:2:XWIDE; 
[XI,YI] = meshgrid(ti,ti);
ZI = griddata(xyG2(1,:),xyG2(2,:),G2Z(:),XI,YI,'v4');
end


%{
if mod(stepN, 20) == 0 && stepN>=100
figure(1);
subplot('Position',[.02 .61 .4 .35]),...
hold on,
grid off,
axis([xlim, ylim, zlim]),
set(gca, 'Box', 'on'),
set(gca,'xticklabel',[]),
set(gca,'yticklabel',[]),
set(gca,'zticklabel',[]),
mesh(XI,YI,ZI),view(70, 10)

end
%}

%{
if mod(stepN, 20) == 0 && stepN>=100
figure(36);
ti = -XWIDE:2:XWIDE; 
[XI,YI] = meshgrid(ti,ti);
ZI = griddata(xyG2(1,:),xyG2(2,:),G2Z(:),XI,YI,'v4');
hold on,
grid off,
axis([xlim, ylim, zlim]),
set(gca, 'Box', 'on'),
set(gca,'xticklabel',[]),
set(gca,'yticklabel',[]),
set(gca,'zticklabel',[]),
mesh(XI,YI,ZI),view(70, 10)
scatter3(xyG2(1,:),xyG2(2,:),G2Z(:),'.','MarkerEdgeColor',[0 0 1]),...
view(70, 10)
hold on;
scatter3(xyG1(1,:),xyG1(2,:),G1Z(:),'.','MarkerEdgeColor',[1 0 0]),...
	view(70, 10)
grid off
axis([xlim, ylim, zlim]);
set(gca, 'Box', 'on');
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])



set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 -tiff CIplots1
end
%}


%======================%
% GLUR2 DOTS
%----------------------%
figure(1);
subplot('Position',[.02 .61 .4 .35]),...
scatter3(xyG2(1,:),xyG2(2,:),G2Z(:),'.','MarkerEdgeColor',[0 0 1]),...
view(70, 10)
hold on;
subplot('Position',[.02 .61 .4 .35]),...
scatter3(xyG1(1,:),xyG1(2,:),G1Z(:),'.','MarkerEdgeColor',[1 0 0]),...
	view(70, 10)
grid off
axis([xlim, ylim, zlim]);
set(gca, 'Box', 'on');
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
%-------%
hold off;



%{
%======================%
% GLUR2 DOTS
%----------------------%
figure(1);
subplot(5,5,[1 7]), subplot('Position',[.05 .65 .36 .32]),...
scatter3(xyG2(1,:),xyG2(2,:),G2Z(:),'.'), view(60, 10)
axis normal;
grid off
axis([xlim, ylim, zlim]);
set(gca, 'Box', 'on');
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
%-------%
hold on;


%======================%
% GLUR1 DOTS
%----------------------%
figure(1);
subplot(5,5,[1 7]), subplot('Position',[.05 .65 .36 .32]),...
scatter3(xyG1(1,:),xyG1(2,:),G1Z(:),'.'), view(60, 10)
axis normal;
grid off
axis([xlim, ylim, zlim]);
set(gca, 'Box', 'on');
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
%-------%
hold off;
%}

%================%
% if stepN == 50; keyboard; end
%================%
end



%-------------------------------------------%
% ONEDOTPLOT
%-------------------------------------------%
function [] = ONEDOTPLOT(stepN,Nsteps,D, xyG2, xyl2, XWIDE,YHIGH,...
				XYBOXpr1,XYBOXpr2,XYBOXp1,XYBOXp2)
%-------------------------------------------%
xlim = [0 XWIDE];
ylim = [-YHIGH YHIGH];
%---

xylcat = cat(2, xyl2,xyG2)';
xl = xylcat(:,1)';
yl = xylcat(:,2)';

figure(1)
subplot(5,5,[3 25]), plot(xl, yl);
axis([xlim, ylim]);
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(get(gca,'XLabel'),'String','GluR2:Blue \bullet \bullet GluR1:Red')
rectangle('Position',[XYBOXpr1(1),XYBOXpr1(2),XYBOXpr1(3),XYBOXpr1(4)])
rectangle('Position',[XYBOXpr2(1),XYBOXpr2(2),XYBOXpr2(3),XYBOXpr2(4)])
rectangle('Position',[XYBOXp1(1),XYBOXp1(2),XYBOXp1(3),XYBOXp1(4)])
rectangle('Position',[XYBOXp2(1),XYBOXp2(2),XYBOXp2(3),XYBOXp2(4)])
hold on;
%---


end



%-------------------------------------------%
% ONEDOTPLOT2
%-------------------------------------------%
function [] = ONEDOTPLOT2(stepN,Nsteps,xyG2,xyl2,XWIDE,YHIGH,...
				XYBOXpr1,XYBOXpr2,XYBOXp1,XYBOXp2)
%-------------------------------------------%
xlim = [0 XWIDE]; ylim = [-YHIGH YHIGH];
%---
 

LLg = 3;

if stepN > LLg
xp = xyl2(1,(stepN-LLg):stepN);
yp = xyl2(2,(stepN-LLg):stepN);
else
xp = xyl2(1,stepN);
yp = xyl2(2,stepN);
end


figure(1)
subplot(5,5,[3 25]), plot(xp, yp);
axis([xlim, ylim]);
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(get(gca,'XLabel'),'String','GluR2:Blue \bullet \bullet GluR1:Red')
rectangle('Position',[XYBOXpr1(1),XYBOXpr1(2),XYBOXpr1(3),XYBOXpr1(4)])
rectangle('Position',[XYBOXpr2(1),XYBOXpr2(2),XYBOXpr2(3),XYBOXpr2(4)])
rectangle('Position',[XYBOXp1(1),XYBOXp1(2),XYBOXp1(3),XYBOXp1(4)])
rectangle('Position',[XYBOXp2(1),XYBOXp2(2),XYBOXp2(3),XYBOXp2(4)])
hold on;
%---

% if mod(Nt,10)==0
% axis vis3d
% camorbit(10,0,'camera')
% drawnow
% end
 
end
















%-------------------------------------------%
% RECTANGLE OVERLAY
%-------------------------------------------%
function fig1box(fig1)
daspect([1,1,1])

% Create rectangle
annotation(fig1,'rectangle',...
	[0.656717687074824 0.254761904761905 0.0468537414966026 0.0539450816823097],...
	'LineWidth',1,...
	'FaceColor','flat',...
	'DisplayName','PSD4');

% Create rectangle
annotation(fig1,'rectangle',...
	[0.639285714285713 0.233333333333333 0.0803571428571431 0.099669794925269],...
	'LineWidth',1,...
	'FaceColor','flat',...
	'DisplayName','PSD1');

% Create rectangle
annotation(fig1,'rectangle',...
	[0.658928571428568 0.721758776503306 0.039285714285716 0.0544316996871702],...
	'LineWidth',1,...
	'FaceColor','flat',...
	'DisplayName','PSDALL');

% Create rectangle
annotation(fig1,'rectangle',...
	[0.637499999999998 0.693187347931874 0.0821428571428572 0.106812652068126],...
	'LineWidth',1,...
	'FaceColor','flat',...
	'DisplayName','PSD2');



end

function fig3box(fig3)

% Create rectangle
annotation(fig3,'rectangle',...
	[0.272959183673469 0.814285714285714 0.0431122448979596 0.056152241918665],...
	'LineWidth',1,...
	'FaceAlpha',0.2,...
	'FaceColor','flat',...
	'DisplayName','PSD1');

% Create rectangle
annotation(fig3,'rectangle',...
	[0.680357142857143 0.188095238095238 0.0875 0.116336461591934],...
	'LineWidth',1,...
	'FaceAlpha',0.2,...
	'FaceColor','flat',...
	'DisplayName','PSD2');

% Create rectangle
annotation(fig3,'rectangle',...
	[0.680357142857143 0.693187347931874 0.0875 0.118717413972888],...
	'LineWidth',1,...
	'FaceAlpha',0.2,...
	'FaceColor','flat',...
	'DisplayName','PSDALL');

% Create rectangle
annotation(fig3,'rectangle',...
	[0.272789115646258 0.521428571428572 0.0450680272108844 0.0563260340632587],...
	'LineWidth',1,...
	'FaceAlpha',0.2,...
	'FaceColor','flat',...
	'DisplayName','PSD4');


end


%-------------------------------------------%
% LIVE PLOT S1 S2 PSD
%-------------------------------------------%
function [] = PLOTS1S2(S1, S2)
%-------------------------------------------%
	figure(1);
	subplot(5,5,11), imagesc(S1);     % <--FIG--##
    colormap('jet')
    title('2D Particle Map');
    title('PSD1');
	figure(1);
    subplot(5,5,12), imagesc(S2); % <--FIG--##
    title('PSD2');
end

%-------------------------------------------%
% clusterplots
%-------------------------------------------%
function [] = clusterplots(S1, S2)
%-------------------------------------------%


S0 = 1;
Pkex0 = 1;
Pen0 = 1;
hk = 1;
SG1oc = 1;
S = 1;

%================================================%
% figure(1);
% subplot(5,5,11), imagesc(S1);
% colormap('bone')
% title('2D Particle Map');
% title('PSD1');
% figure(1);
% subplot(5,5,12), imagesc(S2);
% title('PSD2');
%================================================%
% figure(8) 
% set(gcf,'OuterPosition',[500 55 1200 800])
figure(1);
subplot(5,5,11), imagesc(S1);
%================================================%

% TOP LEFT (S0)
% subplot('Position',[.02 .54 .3 .45]),...
subplot(5,5,1),
imagesc(S0);
colormap('jet');
cbr = colorbar('location','South');
set(cbr(1), 'XColor', 'r');
text(1,1.5,'S0','FontSize',22,'Color',[.7 1 .7]);


% TOP MID (Pkex0)
% subplot('Position',[.35 .54 .3 .45]),...
subplot(5,5,2),
imagesc(Pkex0);
caxis(caxis);
cbr = colorbar('location','South');
set(cbr(1), 'XColor', 'r');
text(1,1.5,'Pkex','FontSize',22,'Color',[.7 1 .7]);


% TOP RIGHT (Pen0)
% subplot('Position',[.68 .54 .3 .45]),...
subplot(5,5,6),
imagesc(Pen0);
colormap('jet');
cbr = colorbar('location','South');
set(cbr(1), 'XColor', 'r');
text(1,1.5,'Pen0','FontSize',22,'Color',[.7 1 .7]);


%=======================================%


% BOTTOM LEFT (hk)
% subplot('Position',[.02 .05 .3 .45]),...
subplot(5,5,7),
imagesc(hk);
caxis(caxis);
cbr = colorbar('location','South');
set(cbr(1), 'XColor', 'g');
text(1,1.5,'hk','FontSize',22,'Color',[.7 1 .7]);


% BOTTOM MID (SG1oc)
% subplot('Position',[.35 .05 .3 .45]),...
subplot(5,5,16),
scfg1 = imagesc(SG1oc); hold on; 
colormap('jet');
cbr = colorbar('location','South');
set(cbr(1), 'XColor', 'g');
scfg2 = imagesc(Lhk);
alpha(.75);
hold off
text(1,1.5,'SG1oc + Lhk','FontSize',22,'Color',[.7 1 .7]);


% BOTTOM RIGHT (S)
% subplot('Position',[.68 .05 .3 .45]),...
subplot(5,5,17),
imagesc(S);
caxis(caxis);
cbr = colorbar('location','South');
set(cbr(1), 'XColor', 'g');
text(1,1.5,'S','FontSize',22,'Color',[.7 1 .7]);



%=======================================%
%{
%================================================%
% % TOP LEFT (S0)
% subplot('Position',[.02 .54 .3 .45]),...
% imagesc(S0);
% colormap('jet');
% cbr = colorbar('location','South');
% set(cbr(1), 'XColor', 'r');
% % xt = (get(gca,'XTick'))-.5;
% % set(gca,'XTickLabel', sprintf('%.0f|',xt))
% % yt = (get(gca,'YTick'))-.5;
% % set(gca,'YTickLabel', sprintf('%.0f|',yt))
% % haxes=axis;
% % ylim([haxes(3)+.5 haxes(4)-.5 ]);
% % xlim([haxes(1)+.5 haxes(2)-.5 ]);
% % grid on
% text(1,1.5,'S0','FontSize',22,'Color',[.7 1 .7]);
% text(XYLTpr1(1),XYLTpr1(2),...
% strcat(num2str(XYLTpr1(1)), '\bullet',num2str(XYLTpr1(2)),'\rightarrow'),...
% 'FontSize',12,'HorizontalAlignment','right');
% alpha('color'); alphamap('rampdown'); alphamap('increase',.1)
%=========================================================%
%%
%---------------------------%
% if stepN>500; keyboard; end
% keyboard
%---------------------------%
%======================================================================%
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf, ['outputfigs/CLUSTREC.png']);
%======================================================================%
%}
%=======================================%


end



%-------------------------------------------%
% GAUSSIAN TETHERING COLORMAP
%-------------------------------------------%
function [] = SLOTMAPFUN(G1STBASE,G1RTBASE,G1STLTP,G1RTLTP,G1BSMu,G1LSMu,...
G2STBASE,G2RTBASE,G2STLTP,G2RTLTP,G2BSMu,G2LSMu,GT,GTab)

saps = 0:4;
tails = 0:.5:4;
offsap = G1BSMu;

%param = [saps+tails(1) saps+tails(2) saps+tails(3) saps+tails(4) saps+tails(5)]-6

param = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G1BSMu;

cdfMx = 1-cdf('norm',0,param,1);
cdfMx2 = [linspace(0,1,numel(tails)); cdfMx];

cdfMxpad = padarray(cdfMx,[1 1],0,'pre');
cdfMxpad(:,1)=[0 saps]';
cdfMxpad(1,:)=[0 tails];

ProbCDF = cdfMx2(1,:);

	
%=========================================================%
% FIGURE SETUP
%---------------------------------------------------------%
fig99 = figure(99);
figure(fig99)
set(99,'Units','pixels') 
scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/1.8  scnsize(4)/7  scnsize(3)/2.3  scnsize(4)/1.2];
% OuterPosition(left, bottom, width, height)
set(fig99,'OuterPosition',pos1)
fig99 = figure(99);
figure(fig99)
set(gcf,'Color',[.9,.9,.9])
%=========================================================%
yt = [-1 saps];
xt = tails; 

subplot(5,5,[16 22]), subplot('Position',[.075 .052 .7 .3]),...
imagesc(cdfMx2)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
set(gca,'YTickLabel', sprintf('%.1f|',yt))
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(.8,1.8,strcat('\color[rgb]{.7 .7 .5}Gaussian \phi = ',...
	num2str(offsap)),'FontSize',18,'HorizontalAlignment','left');
text(6,1,strcat('REFERENCE ROW'),'FontSize',14.5,'HorizontalAlignment','center');
% text(0,1,strcat('P(on)= ', num2str(ProbCDF(:)')),'FontSize',14.5,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
colorbar('location','North')
text(.5,1,strcat('-------------------------------------------------------'),...
	'FontSize',28,'HorizontalAlignment','left',...
	'BackgroundColor',[1 1 1]);
hold on

G1STBL = plot((G1STBASE*2+1),2:6,...
	'MarkerFaceColor',[.9 .3 .5],'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',8,'Marker','d');
hold on
G1STLL = plot((G1STLTP*2+1),2:6,'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',24,'Marker','+');
hold on
G2STBL = plot((G2STBASE*2+1),2:6,...
	'MarkerFaceColor',[.2 .7 .8],'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',8,'Marker','o');
hold on
G2STLL = plot((G2STLTP*2+1),2:6,'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',18,'Marker','x');
hold on
STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');







%=========================================================%
% CDF SETUP		% P(G|S) = P(x??):{x?N(?+?,1)}
%---------------------------------------------------------%
saps = 0:4; tails = 0:.5:4;

G1BSMuTSmx = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G1BSMu;

G1LSMuTSmx = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G1LSMu;
	 
G2BSMuTSmx = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G2BSMu;

G2LSMuTSmx = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G2LSMu;
	 
%--------------------------------------%
% G1BSMu - G1 BASE 
%--------------------------------------%
G1BSMuCDFmx = 1-cdf('norm',0,G1BSMuTSmx,1);
G1BSMuCDF = [linspace(0,1,numel(tails)); G1BSMuCDFmx];

G1BSMuCDFmxpad = padarray(G1BSMuCDFmx,[1 1],0,'pre');
G1BSMuCDFmxpad(:,1)=[0 saps]';
G1BSMuCDFmxpad(1,:)=[0 tails];
G1BSMuCDFp = G1BSMuCDF(1,:);
%--------------------------------------%
%--------------------------------------%
% G1LSMu - G1 LTP 
%--------------------------------------%
G1LSMuCDFmx = 1-cdf('norm',0,G1LSMuTSmx,1);
G1LSMuCDF = [linspace(0,1,numel(tails)); G1LSMuCDFmx];

G1LSMuCDFmxpad = padarray(G1LSMuCDFmx,[1 1],0,'pre');
G1LSMuCDFmxpad(:,1)=[0 saps]';
G1LSMuCDFmxpad(1,:)=[0 tails];
G1LSMuCDFp = G1LSMuCDF(1,:);
%--------------------------------------%

%--------------------------------------%
% G2BSMu - G2 BASE 
%--------------------------------------%
G2BSMuCDFmx = 1-cdf('norm',0,G2BSMuTSmx,1);
G2BSMuCDF = [linspace(0,1,numel(tails)); G2BSMuCDFmx];

G2BSMuCDFmxpad = padarray(G2BSMuCDFmx,[1 1],0,'pre');
G2BSMuCDFmxpad(:,1)=[0 saps]';
G2BSMuCDFmxpad(1,:)=[0 tails];
G2BSMuCDFp = G2BSMuCDF(1,:);
%--------------------------------------%
%--------------------------------------%
% G2BSMu - G2 LTP 
%--------------------------------------%
G2LSMuCDFmx = 1-cdf('norm',0,G2LSMuTSmx,1);
G2LSMuCDF = [linspace(0,1,numel(tails)); G2LSMuCDFmx];

G2LSMuCDFmxpad = padarray(G2LSMuCDFmx,[1 1],0,'pre');
G2LSMuCDFmxpad(:,1)=[0 saps]';
G2LSMuCDFmxpad(1,:)=[0 tails];
G2LSMuCDFp = G2LSMuCDF(1,:);
%--------------------------------------%


%=========================================================%
%					FIGURE SETUP
%---------------------------------------------------------%
fig97 = figure(97);figure(fig97)
set(97,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig97,'OuterPosition',pos1)
fig97 = figure(97);
figure(fig97)
set(gcf,'Color',[.9,.9,.9])
%=========================================================%
yt = [-1 saps];
xt = tails;
% G1Bcolor=[.9 .3 .5];
% G2Bcolor=[.2 .7 .8];
G1Bcolor=[.9 .1 .6];
G2Bcolor=[.2 .7 .8];


%======================================%
% G1BS Mu COLORMAP FIGURE
%--------------------------------------%
subplot('Position',[.06 .07 .44 .40]),...
imagesc(G1BSMuCDF)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
%set(gca,'YTick',0:4)
%set(gca,'YTickLabel', sprintf('%.1f|',yt),'FontSize',14)
set(gca,'XTickLabel', sprintf('%.1f|',xt),'FontSize',14)
set(gca,'YTickLabel', ['N' sprintf('%.0f|',yt(2)) sprintf('%.1f|',yt(2:end))],'FontSize',14)
set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(1,1,strcat('\color[rgb]{.7 .8 .9}G1 basal \phi = ',...
	num2str(G1BSMu)),'FontSize',16,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
%colorbar('location','North')
hold on
%-------------------------%
G1STBL = plot((G1STBASE*2+1),2:6,...
	'MarkerFaceColor',G1Bcolor,'MarkerEdgeColor',G1Bcolor,'MarkerSize',10,'Marker','o');
hold on
% G1STLL = plot((G1STLTP*2+1),2:6,'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',24,'Marker','+');
% hold on
% G2STBL = plot((G2STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.2 .7 .8],'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',8,'Marker','o');
% hold on
% G2STLL = plot((G2STLTP*2+1),2:6,'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',18,'Marker','x');
% hold on
% STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
% 	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
% 	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
% set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');
% SAPMAPTITLE=title(char('Gaussian: ', num2str(offsap)));
% set(SAPMAPTITLE, 'FontSize', 14);
%======================================%





%======================================%
% G1LS Mu COLORMAP FIGURE
%--------------------------------------%
subplot('Position',[.55 .07 .44 .40]),...
imagesc(G1LSMuCDF)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
%set(gca,'YTick',0:4)
% set(gca,'YTickLabel', sprintf('%.1f|',yt),'FontSize',14)
set(gca,'XTickLabel', sprintf('%.1f|',xt),'FontSize',14)
set(gca,'YTickLabel', ['N' sprintf('%.0f|',yt(2)) sprintf('%.1f|',yt(2:end))],'FontSize',14)
% set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
% set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(1,1,strcat('\color[rgb]{.7 .8 .9}G1 LTP \phi = ',...
	num2str(G1LSMu)),'FontSize',16,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
%colorbar('location','North')
hold on
%-------------------------%
% G1STBL = plot((G1STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.9 .3 .5],'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',8,'Marker','d');
% hold on
G1STLL = plot((G1STLTP*2+1),2:6,...
	'MarkerFaceColor',G1Bcolor,'MarkerEdgeColor',G1Bcolor,'MarkerSize',10,'Marker','o');
hold on
% G2STBL = plot((G2STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.2 .7 .8],'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',8,'Marker','o');
% hold on
% G2STLL = plot((G2STLTP*2+1),2:6,'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',18,'Marker','x');
% hold on
% STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
% 	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
% 	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
% set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');
% SAPMAPTITLE=title(char('Gaussian: ', num2str(offsap)));
% set(SAPMAPTITLE, 'FontSize', 14);
%======================================%





%======================================%
% G2BS Mu COLORMAP FIGURE
%--------------------------------------%
subplot('Position',[.06 .57 .44 .40]),...
imagesc(G2BSMuCDF)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
%set(gca,'YTick',0:4)
% set(gca,'YTickLabel', sprintf('%.1f|',yt),'FontSize',14)
set(gca,'XTickLabel', sprintf('%.1f|',xt),'FontSize',14)
set(gca,'YTickLabel', ['N' sprintf('%.0f|',yt(2)) sprintf('%.1f|',yt(2:end))],'FontSize',14)
% set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
% set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(1,1,strcat('\color[rgb]{.7 .8 .9}G2 basal \phi = ',...
	num2str(G2BSMu)),'FontSize',16,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
%colorbar('location','North')
hold on
%-------------------------%
% G1STBL = plot((G1STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.9 .3 .5],'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',8,'Marker','d');
% hold on
% G1STLL = plot((G1STLTP*2+1),2:6,'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',24,'Marker','+');
% hold on
G2STBL = plot((G2STBASE*2+1),2:6,...
	'MarkerFaceColor',G2Bcolor,'MarkerEdgeColor',G2Bcolor,'MarkerSize',10,'Marker','o');
hold on
% G2STLL = plot((G2STLTP*2+1),2:6,'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',18,'Marker','x');
% hold on
% STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
% 	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
% 	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
% set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');
% SAPMAPTITLE=title(char('Gaussian: ', num2str(offsap)));
% set(SAPMAPTITLE, 'FontSize', 14);
%======================================%





%======================================%
% G2LS Mu COLORMAP FIGURE
%--------------------------------------%
subplot('Position',[.55 .57 .44 .40]),...
imagesc(G2LSMuCDF)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
%set(gca,'YTick',0:4)
% set(gca,'YTickLabel', sprintf('%.1f|',yt),'FontSize',14)
set(gca,'XTickLabel', sprintf('%.1f|',xt),'FontSize',14)
set(gca,'YTickLabel', ['N' sprintf('%.0f|',yt(2)) sprintf('%.1f|',yt(2:end))],'FontSize',14)
% set(gca,'YTickLabel', '-')
% set(gca,'XTickLabel', '-')
%axis off

% set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
% set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(1,1,strcat('\color[rgb]{.7 .8 .9}G2 LTP \phi = ',...
	num2str(G2LSMu)),'FontSize',16,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
%colorbar('location','North')
hold on
%-------------------------%
% G1STBL = plot((G1STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.9 .3 .5],'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',8,'Marker','d');
% hold on
% G1STLL = plot((G1STLTP*2+1),2:6,'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',24,'Marker','+');
% hold on
% G2STBL = plot((G2STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.2 .7 .8],'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',8,'Marker','o');
% hold on
G2STLL = plot((G2STLTP*2+1),2:6,...
	'MarkerFaceColor',G2Bcolor,'MarkerEdgeColor',G2Bcolor,'MarkerSize',10,'Marker','d');
hold on
% STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
% 	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
% 	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
% set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');
% SAPMAPTITLE=title(char('Gaussian: ', num2str(offsap)));
% set(SAPMAPTITLE, 'FontSize', 14);
%======================================%


end


%-------------------------------------------%
% GAUSSIAN RECRUITMENT COLORMAP
%-------------------------------------------%
function [] = SAPRECFUN(G1STBASE,G1RTBASE,G1STLTP,G1RTLTP,G1BSMu,G1LSMu,...
G2STBASE,G2RTBASE,G2STLTP,G2RTLTP,G2BSMu,G2LSMu,...
S1,S1sz, sap, stepN,GT,GTab,...
G1P1SAPM, G1SP1, LTP1onG1, LTP1offG1, LTP2onG1, LTP2offG1)
 
%---------------------------------------------------------%
% INITIALIZE CONSTANTS
%---------------------------------------------------------%
% sap1 = sap(1); % SAPdotsPSD1
% sap2 = sap(2); % SAPdotsPSD2
% sap3 = sap(3); % SAPbetaPSD1
% sap4 = sap(4); % SAPtauPSD1
% sap5 = sap(5); % SAPL1PSD1
% sap6 = sap(6); % SAPbetaPSD2
% sap7 = sap(7); % SAPtauPSD2
% sap8 = sap(8); % SAPL1PSD2
% sap9 = sap(9); % SAPmuPSD1
% sap10 = sap(10); % SAPmuPSD2
% sap11 = sap(11); % SAPdTPSD1
% sap12 = sap(12); % SAPdTPSD2
% sap13 = sap(13); % SAPrhoPSD1
% sap14 = sap(14); % SAPrhoPSD2
% sap15 = sap(15); % SAPrPSD1
% sap16 = sap(16); % SAPrPSD2
% sap17 = sap(17); % doDynamicLeP1
% sap18 = sap(18); % doDynamicLeP2
% sap19 = sap(19); % doLTPS1
% sap20 = sap(20); % doLTPS2
 
 
 


%%
%---------------------------------------------------------%
%			SAP CLUSTER MATRIX PROBABILITIES
%---------------------------------------------------------%
S1Le = sap(5);
doDynamicLe = sap(17);
if doDynamicLe
L1x = linspace(1.0,3.2);
Tsz = numel(S1);			% Total matrix space
Ssz = numel(find(S1>0));	% Filled matrix space
Psz = round(Ssz/Tsz*100)+1;	% Percent filled space
S1Le = L1x(Psz);
end
%---------------------------------------------------------%
%	TEST MASK AGAINST CLUSTER - GET CONVOLUTION MX
%---------------------------------------------------------%

S = S1;
S(8,7) = 0; S(10,10) = 0; S(5,12) = 0; S(3,14) = 1;


%-------------------------------%
hkMask=[0 1 0; 1 0 1; 0 1 0];
dT = sap(11);
Le = S1Le;		% hi = shrink
B = sap(3);		% hi = shrink
rho = sap(13);	% hi = grow
r = sap(15);	% hi = grow
mu = sap(9);	% hi = shrink
%-------------------------------%
Pmx = rand(size(S));
Soc = (S>0);
Snoc = (S<1);
Sno = ~Soc;
%---
hk = convn(Soc,hkMask,'same');
Lhk = hk-Le;
%---
Pex = 1 ./ (1+exp((-1*B)*Lhk));
Pkex = (1-Soc) .* ( rho * r * dT * Pex );
Pen = (Soc) .* mu * dT;
%---
S0 = S;
Pex0 = Pex;
Pkex0 = Pkex;
Pen0 = Pen;

%=========================================================%
% figure(9);
% set(gcf,'OuterPosition',[50 50 800 400])
% subplot('Position',[.05 .06 .45 .90]),...
% imagesc(PhkSnoc);
% colormap('bone')
% colorbar('location','North')
% subplot('Position',[.50 .06 .45 .90]),...
% imagesc(Phk);
% caxis(caxis);
% cbr = colorbar('location','North');
% set(cbr(1), 'XColor', 'r');
%=========================================================%

stepN=2;
%---------------------------------------------%
% GluA RT SAP recruitment (basal and LTP)
%---------------------------------------------%
G1RT = G1RTBASE;
if (stepN >= LTP1onG1 && stepN <= LTP1offG1); G1RT=G1RTLTP; end;
if (stepN >= LTP2onG1 && stepN <= LTP2offG1); G1RT=G1RTLTP; end;

G1P1SAPM = [1:3;4:6;7:9];
G1P1SAPM(1,2) = 18;
G1P1SAPM(2,2) = 12;
G1P1SAPM(3,2) = 28;
% close all;
%---------------------------------------------------------%
%       Get Matrix Locations of SAPs near AMPARs
%---------------------------------------------------------%
G1P1poly0 = G1P1SAPM(:,2);
[~, ~, G1P1poly1] = find(G1P1poly0);
G1P1poly2 = G1P1poly1;

for nx = 1:(numel(G1P1poly1))
    if G1P1poly1(nx)<=8 && G1P1poly1(nx)>0
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*0)-1;
    end
    if G1P1poly1(nx)<=16 && G1P1poly1(nx)>8
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*1)-1;
    end
    if G1P1poly1(nx)<=24 && G1P1poly1(nx)>16
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*2)-1;
    end
    if G1P1poly1(nx)<=32 && G1P1poly1(nx)>24
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*3)-1;
    end
    if G1P1poly1(nx)<=40 && G1P1poly1(nx)>32
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*4)-1;
    end
    if G1P1poly1(nx)<=48 && G1P1poly1(nx)>40
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*5)-1;
    end
    if G1P1poly1(nx)<=56 && G1P1poly1(nx)>48
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*6)-1;
    end
    if G1P1poly1(nx)<=64 && G1P1poly1(nx)>56
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*7)-1;
    end
end % for nx
%---------------------------------------------%

SG1oc = (S>2) .* 0.0;
G1P1n = numel(G1P1poly2);

%-------------------------------%
for polyid = 1:G1P1n
%-------------------------------%
G1P1IDo = G1P1poly2(polyid);

G1P1IDb = G1P1IDo+1;
G1P1IDr = G1P1IDo+16;
G1P1IDd = G1P1IDo+17;

if G1P1IDo>=16
	G1P1IDb = G1P1IDo+1;
	G1P1IDr = G1P1IDo+16+(round(rand)*-32); 
	G1P1IDd = G1P1IDr+1;
	% 50:50 WHETHER EAST OR WEST SAP
	% ALWAYS EAST SAP IF COMMENTED
	%G1P1IDd = G1P1IDo+17+(round(rand)*-34);
end
SG1oc(G1P1IDo)=1;
SG1oc(G1P1IDb)=1;
SG1oc(G1P1IDr)=1;
SG1oc(G1P1IDd)=1;
%-------------------------------%
end % for polyid = 1:G1P1n
%-------------------------------%

% Sno		SP non-occupied lattice locations
% SG1oc		G1 occupied lattice locations
% G1RT		G1 recruitment tail constant
% Pkex		SP P(exo) lattice

Gex = Sno .* SG1oc .* G1RT + Pkex;
Gp = 1 ./ (1+exp((-1*B)*Lhk)) + .0001;
Gpex = Gex .* Gp;

% Phkv	= 1-(1./(hk+G1RT+(1/(1+Le))));
% Gpex = Sno .* SG1oc .* Phkv;

%---
% Pkexo = Sno .* G1P1ex;
% Pkexo = Pkexo + Pkex0;

%---
Pkex = Gpex;
Pen = Pen0;
%---

%---------------------------------------------%
Sen = (Pen > Pmx);
Sex = (Pkex > Pmx);
S = (Soc-Sen) + Sex;
%---------------------------------------------%
%=========================================================%
figure(8) 
set(gcf,'OuterPosition',[500 55 1200 800])

%================================================%



% TOP LEFT (S0)
subplot('Position',[.02 .54 .3 .45]),...
imagesc(S0);
colormap('jet');
cbr = colorbar('location','South');
set(cbr(1), 'XColor', 'r');
text(1,1.5,'S0','FontSize',22,'Color',[.7 1 .7]);


% TOP MID (Pkex)
subplot('Position',[.35 .54 .3 .45]),...
imagesc(Pkex0);
caxis(caxis);
cbr = colorbar('location','South');
set(cbr(1), 'XColor', 'r');
text(1,1.5,'Pkex','FontSize',22,'Color',[.7 1 .7]);


% TOP RIGHT (Pen0)
subplot('Position',[.68 .54 .3 .45]),...
imagesc(Pen0);
colormap('jet');
cbr = colorbar('location','South');
set(cbr(1), 'XColor', 'r');
text(1,1.5,'Pen0','FontSize',22,'Color',[.7 1 .7]);


%=======================================%

hk = Sno .* hk;

% BOTTOM LEFT (hk)
subplot('Position',[.02 .05 .3 .45]),...
imagesc(hk);
caxis(caxis);
% cbr = colorbar('location','South');
% set(cbr(1), 'XColor', 'g');
% text(1,1.5,'hk','FontSize',22,'Color',[.7 1 .7]);


% BOTTOM MID (SG1oc)
subplot('Position',[.35 .05 .3 .45]),...
scfg1 = imagesc(SG1oc); hold on; 
colormap('jet');
cbr = colorbar('location','South');
set(cbr(1), 'XColor', 'g');
scfg2 = imagesc(Lhk);
alpha(.75);
hold off
text(1,1.5,'SG1oc + Lhk','FontSize',22,'Color',[.7 1 .7]);


% BOTTOM RIGHT (S)
subplot('Position',[.68 .05 .3 .45]),...
imagesc(S);
caxis(caxis);
cbr = colorbar('location','South');
set(cbr(1), 'XColor', 'g');
text(1,1.5,'S','FontSize',22,'Color',[.7 1 .7]);



%================================================%
% % TOP LEFT (S0)
% subplot('Position',[.02 .54 .3 .45]),...
% imagesc(S0);
% colormap('jet');
% cbr = colorbar('location','South');
% set(cbr(1), 'XColor', 'r');
% % xt = (get(gca,'XTick'))-.5;
% % set(gca,'XTickLabel', sprintf('%.0f|',xt))
% % yt = (get(gca,'YTick'))-.5;
% % set(gca,'YTickLabel', sprintf('%.0f|',yt))
% % haxes=axis;
% % ylim([haxes(3)+.5 haxes(4)-.5 ]);
% % xlim([haxes(1)+.5 haxes(2)-.5 ]);
% % grid on
% text(1,1.5,'S0','FontSize',22,'Color',[.7 1 .7]);
% text(XYLTpr1(1),XYLTpr1(2),...
% strcat(num2str(XYLTpr1(1)), '\bullet',num2str(XYLTpr1(2)),'\rightarrow'),...
% 'FontSize',12,'HorizontalAlignment','right');
% alpha('color'); alphamap('rampdown'); alphamap('increase',.1)
%=========================================================%
%%
%---------------------------%
% if stepN>500; keyboard; end

%---------------------------%
%======================================================================%
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, ['outputfigs/CLUSTREC.png']);
%======================================================================%
%%



%=========================================================%
figure(9);
set(gcf,'OuterPosition',[50 50 800 400])
subplot('Position',[.05 .06 .45 .90]),...
imagesc(PhkSnoc);
% colormap('gray');
colormap('jet');
colorbar('location','North')
subplot('Position',[.50 .06 .45 .90]),...
imagesc(Phk);
caxis(caxis);
cbr = colorbar('location','North');
set(cbr(1), 'XColor', 'r');
%=========================================================%












%=========================================================%
%=========================================================%
saps = 0:4;
tails = 0:.5:4;
offsap = G1BSMu;

param = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G1BSMu;

cdfMx = 1-cdf('norm',0,param,1);
cdfMx2 = [linspace(0,1,numel(tails)); cdfMx];

cdfMxpad = padarray(cdfMx,[1 1],0,'pre');
cdfMxpad(:,1)=[0 saps]';
cdfMxpad(1,:)=[0 tails];

ProbCDF = cdfMx2(1,:);

	
%=========================================================%
% FIGURE SETUP
%---------------------------------------------------------%
fig99 = figure(99);
figure(fig99)
set(99,'Units','pixels') 
scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/1.8  scnsize(4)/7  scnsize(3)/2.3  scnsize(4)/1.2];
% OuterPosition(left, bottom, width, height)
set(fig99,'OuterPosition',pos1)
fig99 = figure(99);
figure(fig99)
set(gcf,'Color',[.9,.9,.9])
%=========================================================%
yt = [-1 saps];
xt = tails; 

subplot(5,5,[16 22]), subplot('Position',[.075 .052 .7 .3]),...
imagesc(cdfMx2)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
set(gca,'YTickLabel', sprintf('%.1f|',yt))
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(.8,1.8,strcat('\color[rgb]{.7 .7 .5}Gaussian \phi = ',...
	num2str(offsap)),'FontSize',18,'HorizontalAlignment','left');
text(6,1,strcat('REFERENCE ROW'),'FontSize',14.5,'HorizontalAlignment','center');
% text(0,1,strcat('P(on)= ', num2str(ProbCDF(:)')),'FontSize',14.5,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
colorbar('location','North')
text(.5,1,strcat('-------------------------------------------------------'),...
	'FontSize',28,'HorizontalAlignment','left',...
	'BackgroundColor',[1 1 1]);
hold on

G1STBL = plot((G1STBASE*2+1),2:6,...
	'MarkerFaceColor',[.9 .3 .5],'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',8,'Marker','d');
hold on
G1STLL = plot((G1STLTP*2+1),2:6,'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',24,'Marker','+');
hold on
G2STBL = plot((G2STBASE*2+1),2:6,...
	'MarkerFaceColor',[.2 .7 .8],'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',8,'Marker','o');
hold on
G2STLL = plot((G2STLTP*2+1),2:6,'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',18,'Marker','x');
hold on
STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');







%=========================================================%
% CDF SETUP		% P(G|S) = P(x??):{x?N(?+?,1)}
%---------------------------------------------------------%
saps = 0:4; tails = 0:.5:4;

G1BSMuTSmx = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G1BSMu;

G1LSMuTSmx = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G1LSMu;
	 
G2BSMuTSmx = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G2BSMu;

G2LSMuTSmx = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G2LSMu;
	 
%--------------------------------------%
% G1BSMu - G1 BASE 
%--------------------------------------%
G1BSMuCDFmx = 1-cdf('norm',0,G1BSMuTSmx,1);
G1BSMuCDF = [linspace(0,1,numel(tails)); G1BSMuCDFmx];

G1BSMuCDFmxpad = padarray(G1BSMuCDFmx,[1 1],0,'pre');
G1BSMuCDFmxpad(:,1)=[0 saps]';
G1BSMuCDFmxpad(1,:)=[0 tails];
G1BSMuCDFp = G1BSMuCDF(1,:);
%--------------------------------------%
%--------------------------------------%
% G1LSMu - G1 LTP 
%--------------------------------------%
G1LSMuCDFmx = 1-cdf('norm',0,G1LSMuTSmx,1);
G1LSMuCDF = [linspace(0,1,numel(tails)); G1LSMuCDFmx];

G1LSMuCDFmxpad = padarray(G1LSMuCDFmx,[1 1],0,'pre');
G1LSMuCDFmxpad(:,1)=[0 saps]';
G1LSMuCDFmxpad(1,:)=[0 tails];
G1LSMuCDFp = G1LSMuCDF(1,:);
%--------------------------------------%

%--------------------------------------%
% G2BSMu - G2 BASE 
%--------------------------------------%
G2BSMuCDFmx = 1-cdf('norm',0,G2BSMuTSmx,1);
G2BSMuCDF = [linspace(0,1,numel(tails)); G2BSMuCDFmx];

G2BSMuCDFmxpad = padarray(G2BSMuCDFmx,[1 1],0,'pre');
G2BSMuCDFmxpad(:,1)=[0 saps]';
G2BSMuCDFmxpad(1,:)=[0 tails];
G2BSMuCDFp = G2BSMuCDF(1,:);
%--------------------------------------%
%--------------------------------------%
% G2BSMu - G2 LTP 
%--------------------------------------%
G2LSMuCDFmx = 1-cdf('norm',0,G2LSMuTSmx,1);
G2LSMuCDF = [linspace(0,1,numel(tails)); G2LSMuCDFmx];

G2LSMuCDFmxpad = padarray(G2LSMuCDFmx,[1 1],0,'pre');
G2LSMuCDFmxpad(:,1)=[0 saps]';
G2LSMuCDFmxpad(1,:)=[0 tails];
G2LSMuCDFp = G2LSMuCDF(1,:);
%--------------------------------------%


%=========================================================%
%					FIGURE SETUP
%---------------------------------------------------------%
fig97 = figure(97);figure(fig97)
set(97,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig97,'OuterPosition',pos1)
fig97 = figure(97);
figure(fig97)
set(gcf,'Color',[.9,.9,.9])
%=========================================================%
yt = [-1 saps];
xt = tails;
% G1Bcolor=[.9 .3 .5];
% G2Bcolor=[.2 .7 .8];
G1Bcolor=[.9 .1 .6];
G2Bcolor=[.2 .7 .8];


%======================================%
% G1BS Mu COLORMAP FIGURE
%--------------------------------------%
subplot('Position',[.06 .07 .44 .40]),...
imagesc(G1BSMuCDF)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
%set(gca,'YTick',0:4)
%set(gca,'YTickLabel', sprintf('%.1f|',yt),'FontSize',14)
set(gca,'XTickLabel', sprintf('%.1f|',xt),'FontSize',14)
set(gca,'YTickLabel', ['N' sprintf('%.0f|',yt(2)) sprintf('%.1f|',yt(2:end))],'FontSize',14)
set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(1,1,strcat('\color[rgb]{.7 .8 .9}G1 basal \phi = ',...
	num2str(G1BSMu)),'FontSize',16,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
%colorbar('location','North')
hold on
%-------------------------%
G1STBL = plot((G1STBASE*2+1),2:6,...
	'MarkerFaceColor',G1Bcolor,'MarkerEdgeColor',G1Bcolor,'MarkerSize',10,'Marker','o');
hold on
% G1STLL = plot((G1STLTP*2+1),2:6,'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',24,'Marker','+');
% hold on
% G2STBL = plot((G2STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.2 .7 .8],'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',8,'Marker','o');
% hold on
% G2STLL = plot((G2STLTP*2+1),2:6,'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',18,'Marker','x');
% hold on
% STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
% 	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
% 	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
% set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');
% SAPMAPTITLE=title(char('Gaussian: ', num2str(offsap)));
% set(SAPMAPTITLE, 'FontSize', 14);
%======================================%





%======================================%
% G1LS Mu COLORMAP FIGURE
%--------------------------------------%
subplot('Position',[.55 .07 .44 .40]),...
imagesc(G1LSMuCDF)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
%set(gca,'YTick',0:4)
% set(gca,'YTickLabel', sprintf('%.1f|',yt),'FontSize',14)
set(gca,'XTickLabel', sprintf('%.1f|',xt),'FontSize',14)
set(gca,'YTickLabel', ['N' sprintf('%.0f|',yt(2)) sprintf('%.1f|',yt(2:end))],'FontSize',14)
% set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
% set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(1,1,strcat('\color[rgb]{.7 .8 .9}G1 LTP \phi = ',...
	num2str(G1LSMu)),'FontSize',16,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
%colorbar('location','North')
hold on
%-------------------------%
% G1STBL = plot((G1STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.9 .3 .5],'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',8,'Marker','d');
% hold on
G1STLL = plot((G1STLTP*2+1),2:6,...
	'MarkerFaceColor',G1Bcolor,'MarkerEdgeColor',G1Bcolor,'MarkerSize',10,'Marker','o');
hold on
% G2STBL = plot((G2STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.2 .7 .8],'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',8,'Marker','o');
% hold on
% G2STLL = plot((G2STLTP*2+1),2:6,'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',18,'Marker','x');
% hold on
% STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
% 	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
% 	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
% set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');
% SAPMAPTITLE=title(char('Gaussian: ', num2str(offsap)));
% set(SAPMAPTITLE, 'FontSize', 14);
%======================================%





%======================================%
% G2BS Mu COLORMAP FIGURE
%--------------------------------------%
subplot('Position',[.06 .57 .44 .40]),...
imagesc(G2BSMuCDF)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
%set(gca,'YTick',0:4)
% set(gca,'YTickLabel', sprintf('%.1f|',yt),'FontSize',14)
set(gca,'XTickLabel', sprintf('%.1f|',xt),'FontSize',14)
set(gca,'YTickLabel', ['N' sprintf('%.0f|',yt(2)) sprintf('%.1f|',yt(2:end))],'FontSize',14)
% set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
% set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(1,1,strcat('\color[rgb]{.7 .8 .9}G2 basal \phi = ',...
	num2str(G2BSMu)),'FontSize',16,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
%colorbar('location','North')
hold on
%-------------------------%
% G1STBL = plot((G1STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.9 .3 .5],'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',8,'Marker','d');
% hold on
% G1STLL = plot((G1STLTP*2+1),2:6,'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',24,'Marker','+');
% hold on
G2STBL = plot((G2STBASE*2+1),2:6,...
	'MarkerFaceColor',G2Bcolor,'MarkerEdgeColor',G2Bcolor,'MarkerSize',10,'Marker','o');
hold on
% G2STLL = plot((G2STLTP*2+1),2:6,'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',18,'Marker','x');
% hold on
% STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
% 	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
% 	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
% set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');
% SAPMAPTITLE=title(char('Gaussian: ', num2str(offsap)));
% set(SAPMAPTITLE, 'FontSize', 14);
%======================================%





%======================================%
% G2LS Mu COLORMAP FIGURE
%--------------------------------------%
subplot('Position',[.55 .57 .44 .40]),...
imagesc(G2LSMuCDF)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
%set(gca,'YTick',0:4)
% set(gca,'YTickLabel', sprintf('%.1f|',yt),'FontSize',14)
set(gca,'XTickLabel', sprintf('%.1f|',xt),'FontSize',14)
set(gca,'YTickLabel', ['N' sprintf('%.0f|',yt(2)) sprintf('%.1f|',yt(2:end))],'FontSize',14)
% set(gca,'YTickLabel', '-')
% set(gca,'XTickLabel', '-')
%axis off

% set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
% set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(1,1,strcat('\color[rgb]{.7 .8 .9}G2 LTP \phi = ',...
	num2str(G2LSMu)),'FontSize',16,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
%colorbar('location','North')
hold on
%-------------------------%
% G1STBL = plot((G1STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.9 .3 .5],'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',8,'Marker','d');
% hold on
% G1STLL = plot((G1STLTP*2+1),2:6,'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',24,'Marker','+');
% hold on
% G2STBL = plot((G2STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.2 .7 .8],'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',8,'Marker','o');
% hold on
G2STLL = plot((G2STLTP*2+1),2:6,...
	'MarkerFaceColor',G2Bcolor,'MarkerEdgeColor',G2Bcolor,'MarkerSize',10,'Marker','d');
hold on
% STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
% 	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
% 	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
% set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');
% SAPMAPTITLE=title(char('Gaussian: ', num2str(offsap)));
% set(SAPMAPTITLE, 'FontSize', 14);
%======================================%


end



%===================================%
% HOMEOSTATIC
%===================================%
function [S1beta S2beta S1mu_array S2mu_array S1L1_array S2L1_array]...
	= HOMEOSTATIC(stepN,HShi,HSlow,TTDots,doCalcium,SPYNT,PSD1CaT,...
	S1num_epochs,S2num_epochs,S1mu,S2mu)

	if mod(stepN, 5) == 0
		PSDtot = SPYNT;
		PSDtot2 = SPYNT;
		if doCalcium
		PSDtot2 = SPYNT-(fix(sqrt(PSD1CaT)));		
		end
		%{
		% LINHi = linspace(60,80,(TTDots+1));
		% LINLo = linspace(25,45,(TTDots+1));
		% when L=1.1  tau=1.1  beta=60-80   higer=slower growth
        % when L=2.1  tau=2.1  beta=25-45   lower=slower shrink
		%}
		LINHi = linspace(60,80,(HSlow+2));
		LINLo = linspace(25,45,(TTDots-HShi+2));
		
		if PSDtot < HSlow
			Ltau = 1.1;
			beta1 = LINHi(PSDtot+1);
			S1mew=1/Ltau;
		elseif PSDtot > HShi
			Ltau = 2.1;
			beta1 = LINLo(PSDtot-HShi+1);
			S1mew=1/Ltau;
		else
			Ltau = 1.2;
			beta1 = 50;
			S1mew=S1mu;
		end
		
		
		if PSDtot2 < HSlow
			Ltau2 = 1.1;
			beta2 = LINHi(PSDtot+1);
			S2mew=1/Ltau2;
		elseif PSDtot2 > HShi
			Ltau2 = 2.1;
			beta2 = LINLo(PSDtot-HShi+1);
			S2mew=1/Ltau2;
		else
			Ltau2 = 1.2;
			beta2 = 50;
			S2mew=S2mu;
		end
		S1beta=beta1;   
        S2beta=beta2;
        S1mu_array=S1mew*ones(1,S1num_epochs);
        S2mu_array=S2mew*ones(1,S2num_epochs);
        S1L1_array=Ltau*ones(1,S1num_epochs);
        S2L1_array=Ltau2*ones(1,S2num_epochs);
	end


end
%===================================%


%===================================%
%		ZGEN FUNCTION
%===================================%
function [G1Z G2Z] = ZGEN(stepN,G1Z,G2Z,...
    G2N, xyG2,G1N, xyG1,...
	XYLTpr1,XYRTpr1,XYLBpr1,XYRBpr1,XYLTpr2,XYRTpr2,XYLBpr2,XYRBpr2,...
	XYLTp1,XYRTp1,XYLBp1,XYRBp1,XYLTp2,XYRTp2,XYLBp2,XYRBp2)



%--GLUR2--%
for j = 1:G2N 
%Inside PSD1
if (xyG2(1,j)>=XYLBp1(1) && xyG2(1,j)<=XYRTp1(1)) &&...
   (xyG2(2,j)>=XYLBp1(2) && xyG2(2,j)<=XYRTp1(2))
 
    G2Z(:,j) = 4;
 
%Inside PSD2
elseif  (xyG2(1,j)>=XYLBp2(1) && xyG2(1,j)<=XYRTp2(1)) &&...
        (xyG2(2,j)>=XYLBp2(2) && xyG2(2,j)<=XYRTp2(2))
    G2Z(:,j) = 4;
 
%Inside PSA1
elseif  (xyG2(1,j)>=XYLBpr1(1) && xyG2(1,j)<=XYRTpr1(1)) &&...
        (xyG2(2,j)>=XYLBpr1(2) && xyG2(2,j)<=XYRTpr1(2)) &&...
        ((xyG2(1,j)<=XYLBp1(1) || xyG2(1,j)>=XYRTp1(1))) &&...
        ((xyG2(2,j)<=XYLBp1(2) || xyG2(2,j)>=XYRTp1(2)))
    
    G2Z(:,j) = 2+(.1 * randi(15,1));
 
%Inside PSA2
elseif  (xyG2(1,j)>=XYLBpr2(1) && xyG2(1,j)<=XYRTpr2(1)) &&...
        (xyG2(2,j)>=XYLBpr2(2) && xyG2(2,j)<=XYRTpr2(2)) &&...
        ((xyG2(1,j)<=XYLBp2(1) || xyG2(1,j)>=XYRTp2(1))) &&...
        ((xyG2(2,j)<=XYLBp2(2) || xyG2(2,j)>=XYRTp2(2)))
    
    G2Z(:,j) = 2+(.1 * randi(15,1));
 
else
    G2Z(:,j) = 1;
end
 
end
 
%--GLUR1--%
for j = 1:G1N 
%Inside PSD1
if (xyG1(1,j)>=XYLBp1(1) && xyG1(1,j)<=XYRTp1(1)) &&...
   (xyG1(2,j)>=XYLBp1(2) && xyG1(2,j)<=XYRTp1(2))
 
    G1Z(:,j) = 4;
 
%Inside PSD2
elseif  (xyG1(1,j)>=XYLBp2(1) && xyG1(1,j)<=XYRTp2(1)) &&...
        (xyG1(2,j)>=XYLBp2(2) && xyG1(2,j)<=XYRTp2(2))
    G1Z(:,j) = 4;
 
%Inside PSA1
elseif  (xyG1(1,j)>=XYLBpr1(1) && xyG1(1,j)<=XYRTpr1(1)) &&...
        (xyG1(2,j)>=XYLBpr1(2) && xyG1(2,j)<=XYRTpr1(2)) &&...
        ((xyG1(1,j)<=XYLBp1(1) || xyG1(1,j)>=XYRTp1(1))) &&...
        ((xyG1(2,j)<=XYLBp1(2) || xyG1(2,j)>=XYRTp1(2)))
    
    G1Z(:,j) = 2+(.1 * randi(15,1));
 
%Inside PSA2
elseif  (xyG1(1,j)>=XYLBpr2(1) && xyG1(1,j)<=XYRTpr2(1)) &&...
        (xyG1(2,j)>=XYLBpr2(2) && xyG1(2,j)<=XYRTpr2(2)) &&...
        ((xyG1(1,j)<=XYLBp2(1) || xyG1(1,j)>=XYRTp2(1))) &&...
        ((xyG1(2,j)<=XYLBp2(2) || xyG1(2,j)>=XYRTp2(2)))
    
    G1Z(:,j) = 2+(.1 * randi(15,1));
 
else
    G1Z(:,j) = 1;
end
 
end
 

%{
%-------------------%
for j = 1:G2N
%-------------------%

		if xyG2(1,j)>=XcPr1lft && xyG2(1,j)<=XcPr1rit &&...
          xyG2(2,j)>=YrPr1top && xyG2(2,j)<=YrPr1bot &&...
		  xyG2(1,j)>=(XcPr1lft+PSDSZE(2,1)) && xyG2(1,j)<=(XcPr1rit-PSDSZE(2,1)) &&...
          xyG2(2,j)>=(YrPr1top+PSDSZE(2,1)) && xyG2(2,j)<=(YrPr1bot-PSDSZE(2,1))
		
			G2Z(:,j) = 3;
	
		elseif xyG2(1,j)>=XcPr1lft && xyG2(1,j)<=XcPr1rit &&...
          xyG2(2,j)>=YrPr2bot && xyG2(2,j)<=YrPr2top &&...
		  xyG2(1,j)>=(XcPr1lft+PSDSZE(2,2)) && xyG2(1,j)<=(XcPr1rit-PSDSZE(2,2)) &&...
          xyG2(2,j)>=(YrPr2bot+PSDSZE(2,2)) && xyG2(2,j)<=(YrPr2top-PSDSZE(2,2))

			G2Z(:,j) = 3;
			
			
		elseif xyG2(1,j)>=XcPr1lft && xyG2(1,j)<=XcPr1rit &&...
          xyG2(2,j)>=YrPr1top && xyG2(2,j)<=YrPr1bot &&...
          ((xyG2(1,j)<=(XcPr1lft+PSDSZE(2,1)) || xyG2(1,j)>=(XcPr1rit-PSDSZE(2,1))) ||...
          (xyG2(2,j)<=(YrPr1top+PSDSZE(2,1)) || xyG2(2,j)>=(YrPr1bot-PSDSZE(2,1))))
      
            G2Z(:,j) = 1+(.1 * randi(15,1));
            
        elseif xyG2(1,j)>=XcPr1lft && xyG2(1,j)<=XcPr1rit &&...
          xyG2(2,j)>=YrPr2bot && xyG2(2,j)<=YrPr2top &&...
          ((xyG2(1,j)<=(XcPr1lft+PSDSZE(2,2)) || xyG2(1,j)>=(XcPr1rit-PSDSZE(2,2))) ||...
          (xyG2(2,j)<=(YrPr2bot+PSDSZE(2,2)) || xyG2(2,j)>=(YrPr2top-PSDSZE(2,2))))
	  
            G2Z(:,j) = 1+(.1 * randi(15,1));
			
			
		else
		   G2Z(:,j) = 1;
		end % if
%-------------------%
end  %j = 1:G2N
%-------------------%



%-------------------%
for j = 1:G1N
%-------------------%
 
        if xyG1(1,j)>=XcPr1lft && xyG1(1,j)<=XcPr1rit &&...
          xyG1(2,j)>=YrPr1top && xyG1(2,j)<=YrPr1bot &&...
          xyG1(1,j)>=(XcPr1lft+PSDSZE(2,1)) && xyG1(1,j)<=(XcPr1rit-PSDSZE(2,1)) &&...
          xyG1(2,j)>=(YrPr1top+PSDSZE(2,1)) && xyG1(2,j)<=(YrPr1bot-PSDSZE(2,1))
        
            G1Z(:,j) = 3;
    
        elseif xyG1(1,j)>=XcPr1lft && xyG1(1,j)<=XcPr1rit &&...
          xyG1(2,j)>=YrPr2bot && xyG1(2,j)<=YrPr2top &&...
          xyG1(1,j)>=(XcPr1lft+PSDSZE(2,2)) && xyG1(1,j)<=(XcPr1rit-PSDSZE(2,2)) &&...
          xyG1(2,j)>=(YrPr2bot+PSDSZE(2,2)) && xyG1(2,j)<=(YrPr2top-PSDSZE(2,2))
 
            G1Z(:,j) = 3;
			
			
		elseif xyG1(1,j)>=XcPr1lft && xyG1(1,j)<=XcPr1rit &&...
          xyG1(2,j)>=YrPr1top && xyG1(2,j)<=YrPr1bot &&...
          ((xyG1(1,j)<=(XcPr1lft+PSDSZE(2,1)) || xyG1(1,j)>=(XcPr1rit-PSDSZE(2,1))) ||...
          (xyG1(2,j)<=(YrPr1top+PSDSZE(2,1)) || xyG1(2,j)>=(YrPr1bot-PSDSZE(2,1))))
      
            G1Z(:,j) = 1+(.1 * randi(15,1));
            
        elseif xyG1(1,j)>=XcPr1lft && xyG1(1,j)<=XcPr1rit &&...
          xyG1(2,j)>=YrPr2bot && xyG1(2,j)<=YrPr2top &&...
          ((xyG1(1,j)<=(XcPr1lft+PSDSZE(2,2)) || xyG1(1,j)>=(XcPr1rit-PSDSZE(2,2))) ||...
          (xyG1(2,j)<=(YrPr2bot+PSDSZE(2,2)) || xyG1(2,j)>=(YrPr2top-PSDSZE(2,2))))
      
            G1Z(:,j) = 1+(.1 * randi(15,1));
			
            
        else
           G1Z(:,j) = 1;
        end % if
%-------------------%
end  %j = 1:G1N
%-------------------%
%}
	
%================%
% if stepN == 100; keyboard; end
%================%
end
%===================================%


%===================================%
%		inboxfun
%===================================%
% Tests whether particles are in a box polygon
% and returns a logical vector
%===================================%
function [inbox] = inboxfunGENERAL(LB,RT,xyl)

if LB(1)>RT(1)
	LBt=LB;
	RTt=RT;
	LB(1)=RTt(1);
	RT(1)=LBt(1);
end
if LB(2)>RT(2)
	LBt=LB;
	RTt=RT;
	LB(2)=RTt(2);
	RT(2)=LBt(2);
end


sz = numel(xyl(1,:));
inbox = zeros(1,sz);

for s = 1:sz
	
	if xyl(1,s) > LB(1) && xyl(1,s) < RT(1) &&...
	   xyl(2,s) > LB(2) && xyl(2,s) < RT(2)
		
		inbox(s) = 1;
	end
   
end

inbox = inbox>0;

end
%===================================%



%-------------------------------------------%
% GLUR-SAPPOLYGON FUNCTIONS
%-------------------------------------------%
%{

%-------------------------------------------%
% GLUR-SAPPOLYGON FUNCTIONS
%-------------------------------------------%
% These SAP polygon functions take a SAPFmx that contains the
% number of SAPs in each SLOT across the PSD
% It cycles through the polygons in each PSD testing whether
% a given ampar is located in that region, and adds the current
% SAP value to a Mx with the same vector ID as the ampar
% it then returns this (2 x Ndot) Mx of values ranging from 0 to 4
%-------------------%
%---G1SAPPOLYGON---%
function [GluR1SP1 GluR1SP2 G1SAPLOC1 G1SAPLOC2]...
= G1SAPPOLYGON(stepN,...
PSD1Sxvecs,PSD1Syvecs,PSD2Sxvecs,PSD2Syvecs,...
xyG1,SAPFmx1,SAPFmx2)

    
GluR1Sd1 = ones(size(xyG1));
GluR1Sd2 = GluR1Sd1;

GluR1SLOCMXa = zeros(1,size(GluR1Sd1,2));
GluR1SLOCMXb = GluR1SLOCMXa;

[m,p,n] = size(PSD1Sxvecs);
	% currentcell = zeros(1,size(xyG1,2));
	currentcell = zeros(1,p);
	

for m = 1:p
	currentcell = currentcell+1;
    
PSD1SxvNOW = PSD1Sxvecs{m};
PSD1SyvNOW = PSD1Syvecs{m};
PSD2SxvNOW = PSD2Sxvecs{m};
PSD2SyvNOW = PSD2Syvecs{m};

G1INPSD1 = inpolygon(xyG1(1,:)',xyG1(2,:)',PSD1SxvNOW,PSD1SyvNOW);
G1INPSD2 = inpolygon(xyG1(1,:)',xyG1(2,:)',PSD2SxvNOW,PSD2SyvNOW);
 
GluR1Sd1(:,G1INPSD1) = GluR1Sd1(:,G1INPSD1)+(SAPFmx1(m));
GluR1Sd2(:,G1INPSD2) = GluR1Sd2(:,G1INPSD2)+(SAPFmx2(m));

GluR1SLOCMXa(1,G1INPSD1) = GluR1SLOCMXa(:,G1INPSD1)+(currentcell(m));
GluR1SLOCMXb(1,G1INPSD2) = GluR1SLOCMXb(:,G1INPSD2)+(currentcell(m));


end
 
GluR1SdP1 = GluR1Sd1-1;
GluR1SdP2 = GluR1Sd2-1;

GluR1SP1 = GluR1SdP1(1,:);
GluR1SP2 = GluR1SdP2(2,:);


G1SAPLOCP1 = [GluR1SP1;GluR1SLOCMXa];
G1SAPLOCP2 = [GluR1SP2;GluR1SLOCMXb];



%----
GluR1SP1 = GluR1SdP1(1,:);
GluR1SP2 = GluR1SdP2(2,:);
G1SAPLOC1 = G1SAPLOCP1(2,:);
G1SAPLOC2 = G1SAPLOCP2(2,:);




% if stepN == 201; keyboard; end
end

%---G2SAPPOLYGON---%
function [GluR2SP1 GluR2SP2 G2SAPLOC1 G2SAPLOC2]...
= G2SAPPOLYGON(stepN,...
PSD1Sxvecs,PSD1Syvecs,PSD2Sxvecs,PSD2Syvecs,...
xyG2,SAPFmx1,SAPFmx2)
 
    
GluR2Sd1 = ones(size(xyG2));
GluR2Sd2 = GluR2Sd1;
 
GluR2SLOCMXa = zeros(1,size(GluR2Sd1,2));
GluR2SLOCMXb = GluR2SLOCMXa;
 
[m,p,n] = size(PSD1Sxvecs);
    % currentcell = zeros(1,size(xyG2,2));
	currentcell = zeros(1,p);
	
for m = 1:p
    currentcell = currentcell+1;
    
PSD1SxvNOW = PSD1Sxvecs{m};
PSD1SyvNOW = PSD1Syvecs{m};
PSD2SxvNOW = PSD2Sxvecs{m};
PSD2SyvNOW = PSD2Syvecs{m};
 
G2INPSD1 = inpolygon(xyG2(1,:)',xyG2(2,:)',PSD1SxvNOW,PSD1SyvNOW);
G2INPSD2 = inpolygon(xyG2(1,:)',xyG2(2,:)',PSD2SxvNOW,PSD2SyvNOW);
 
GluR2Sd1(:,G2INPSD1) = GluR2Sd1(:,G2INPSD1)+(SAPFmx1(m));
GluR2Sd2(:,G2INPSD2) = GluR2Sd2(:,G2INPSD2)+(SAPFmx2(m));
 
GluR2SLOCMXa(1,G2INPSD1) = GluR2SLOCMXa(:,G2INPSD1)+(currentcell(m));
GluR2SLOCMXb(1,G2INPSD2) = GluR2SLOCMXb(:,G2INPSD2)+(currentcell(m));
 
 
end
 
GluR2SdP1 = GluR2Sd1-1;
GluR2SdP2 = GluR2Sd2-1;
 
GluR2SP1 = GluR2SdP1(1,:);
GluR2SP2 = GluR2SdP2(2,:);
 
 
G2SAPLOCP1 = [GluR2SP1;GluR2SLOCMXa];
G2SAPLOCP2 = [GluR2SP2;GluR2SLOCMXb];
 
 
 
%----
GluR2SP1 = GluR2SdP1(1,:);
GluR2SP2 = GluR2SdP2(2,:);
G2SAPLOC1 = G2SAPLOCP1(2,:);
G2SAPLOC2 = G2SAPLOCP2(2,:);
 
 
 
 
% if stepN == 100; keyboard; end
end
%===================================%


%}
%===================================%

%-------------------------------------------%
%			  AMPARS IN SLOTS
%-------------------------------------------%
%{

%-------------------------------------------%
%			  AMPARS IN SLOTS
%-------------------------------------------%
% GET THE SLOT DATA FOR AMPARS (FOR EACH PSD REGION SEPARATELY)
% THEN GET THE DWELL TIME VALUES AND THE MX IDs
% ROLL THE DICE FOR EACH GLUR PARTICLE IN A SLOT REGION
% UPDATE THE DIFFUSION RATES FOR THOSE GLUR RECEPTORS THAT MADE IT INTO A SLOT
%-------------------%
%----SLOTS GluR1----%
function [xysG1 GluR1_TdwellPSD GluR1_TdwellPERI GluR1_TdwellSPYN G1FSLOTS]...
	= G1SLOTS(stepN,G1N,xysG1,...
	GluR1_TdwellPSD,GluR1_TdwellPERI,GluR1_TdwellSPYN,...
	GluR1SP1,GluR1SP2,SAPFmx1,SAPFmx2,G1SAPLOC1, G1SAPLOC2)

%===========================%
% Slots
%---------------------------%   
PSD1slotN = numel(SAPFmx1);
PSD2slotN = numel(SAPFmx2);
%===========================%
%{

	[rPSD1,cDOT1,vPOLY1] = find(G1SAPLOC1);
	vrcPOLYDOT1 = [vPOLY1;rPSD1;cDOT1]';
	[rPSD2,cDOT2,vPOLY2] = find(G1SAPLOC2);
	vrcPOLYDOT2 = [vPOLY2;rPSD2;cDOT2]';
	
	[rPSD1IDpoly] = find(vrcPOLYDOT1(:,2)>0);	% PSD1 DOT-ID DATA
    [rPSD2IDpoly] = find(vrcPOLYDOT2(:,2)>0);	% PSD2 DOT-ID DATA
    vrcTIMEpoly1 = vrcPOLYDOT1(rPSD1IDpoly,:);	% PSD1 DWELL DATA ONLY
    vrcTIMEpoly2 = vrcPOLYDOT2(rPSD2IDpoly,:);	% PSD2 DWELL DATA ONLY
    vrcTIME1polys = sortrows(vrcTIMEpoly1,-1);	% PSD1: SORT DWELL DATA LONG TO SHORT
    vrcTIME2polys = sortrows(vrcTIMEpoly2,-1);	% PSD2: SORT DWELL DATA LONG TO SHORT
	
	if size(vrcTIME1polys,1) < PSD1slotN     % IF ALL SLOTS NOT FILLED MAKE BLANKS
        vrcTIME1polys = vertcat(vrcTIME1polys,zeros(PSD1slotN,3));
    end
    if size(vrcTIME2polys,1) < PSD2slotN
        vrcTIME2polys = vertcat(vrcTIME2polys,zeros(PSD2slotN,3));
	end 
	

% GET THE SLOT DATA FOR ONLY FOR GLUR1 AMPARS (FOR EACH PSD REGION SEPARATELY)
% THESE LABELS ARE A LITTLE MISLEADING BECAUSE HERE WE ARE ONLY
% INTERESTED IN MAKING SURE THE SLOT MX ID VALUES ARE RETAINED 
% TO MATCH THE DWELL TIME MX VALUES SORTED BELOW 

	[rPSD1,cDOT1,vTIME1] = find(GluR1SP1);
	vrcTIMEPSDDOT1 = [vTIME1;rPSD1;cDOT1]';
	[rPSD2,cDOT2,vTIME2] = find(GluR1SP2);
	vrcTIMEPSDDOT2 = [vTIME2;rPSD2;cDOT2]';
	
	[rPSD1ID] = find(vrcTIMEPSDDOT1(:,2)>0);	% PSD1 DOT-ID DATA
    [rPSD2ID] = find(vrcTIMEPSDDOT2(:,2)>0);	% PSD2 DOT-ID DATA
    vrcTIME1 = vrcTIMEPSDDOT1(rPSD1ID,:);		% PSD1 DWELL DATA ONLY
    vrcTIME2 = vrcTIMEPSDDOT2(rPSD2ID,:);		% PSD2 DWELL DATA ONLY
    vrcTIME1s = sortrows(vrcTIME1,-1);			% PSD1: SORT DWELL DATA LONG TO SHORT
    vrcTIME2s = sortrows(vrcTIME2,-1);			% PSD2: SORT DWELL DATA LONG TO SHORT
	
	if size(vrcTIME1s,1) < PSD1slotN     % IF ALL SLOTS NOT FILLED MAKE BLANKS
        vrcTIME1s = vertcat(vrcTIME1s,zeros(PSD1slotN,3));
    end
    if size(vrcTIME2s,1) < PSD2slotN
        vrcTIME2s = vertcat(vrcTIME2s,zeros(PSD2slotN,3));
	end 




	% (SLPSD1LOCkOFF: DOT)
	% SLPSD1LOCkOFF = SLPSD1((SLPSD1(:,1)>(2)),3);
    % SLPSD2LOCkOFF = SLPSD2((SLPSD2(:,1)>(2)),3);
	
	%}
	



% NOW GET ME THE DWELL TIME VALUES AND THE MX ID
% LOCATIONS OF THOSE GLUR1 RECEPTORS IN PSD1 OR PSD2

	[r,c,v] = find(GluR1_TdwellSPYN);	% WHICH PSD, WHICH AMPAR, HOW LONG
    rcv = [v r c];                      % ONLY AMPARS WITH DWELL TIMES (TIME,PSD,DOT)
    [r1] = find(rcv(:,2)<2);            % PSD1 DATA
    [r2] = find(rcv(:,2)>1);            % PSD2 DATA
    p1r = rcv(r1,:);                    % PSD1 DWELL DATA ONLY
    p2r = rcv(r2,:);                    % PSD2 DWELL DATA ONLY
    p1rs = sortrows(p1r,-1);            % PSD1: SORT DWELL DATA LONG TO SHORT
    p2rs = sortrows(p2r,-1);            % PSD2: SORT DWELL DATA LONG TO SHORT
    
    if size(p1rs,1) < PSD1slotN     % IF ALL SLOTS NOT FILLED MAKE BLANKS
        p1rs = vertcat(p1rs,zeros(PSD1slotN,3));
    end
    if size(p2rs,1) < PSD2slotN
        p2rs = vertcat(p2rs,zeros(PSD2slotN,3));
	end 
    
	% (SLPSD: TIME,PSD,DOT)
    SLPSD1 = p1rs(1:PSD1slotN,:);   %PSD1: GET THE N(slots) LONGEST DWELLING AMPARS
    SLPSD2 = p2rs(1:PSD2slotN,:);   %PSD2: GET THE N(slots) LONGEST DWELLING AMPARS
	

% LETS SEE SOME IDs FELLAS
	
	% (SLPSD1LOCkON: DOT)
	SLPSD1LOCkON = SLPSD1((SLPSD1(:,1)>0),3);
    SLPSD2LOCkON = SLPSD2((SLPSD2(:,1)>0),3);

	
	G1SAPLOCD1 = [G1SAPLOC1]';
	G1SAPLOCD2 = [G1SAPLOC2]';
	GluR1DSP1 = [GluR1SP1]';
	GluR1DSP2 = [GluR1SP2]';
	
GluR1DOTSP1 = [GluR1DSP1(SLPSD1LOCkON,:) SLPSD1LOCkON G1SAPLOCD1(SLPSD1LOCkON,:)];
GluR1DOTSP2 = [GluR1DSP2(SLPSD2LOCkON,:) SLPSD2LOCkON G1SAPLOCD2(SLPSD2LOCkON,:)];

% THESE ARE THE SLOTTED DOTS
	[uniV1 uniR1 uniC1] = unique(GluR1DOTSP1(:,3),'rows','first');
	GluR1DOTSP1 = GluR1DOTSP1(uniR1,:);
	GluR1DOTSP1 = GluR1DOTSP1(:,1:2);
	
	[uniV2 uniR2 uniC2] = unique(GluR1DOTSP2(:,3),'rows','first');
	GluR1DOTSP2 = GluR1DOTSP2(uniR2,:);
	GluR1DOTSP2 = GluR1DOTSP2(:,1:2);
	

	
% HERE IS WHERE WE ROLL THE DICE FOR EACH GLUR1 PARTICLE IN A SLOT REGION
% TO DETERMINE WHETHER IT BINDS TO THE SLOT

	GR1PROBPR1 = randn(numel(GluR1DOTSP1(:,1)),1)-4;
	GR1PROBP1 = [GluR1DOTSP1(:,1) GluR1DOTSP1(:,2) (GluR1DOTSP1(:,1)+GR1PROBPR1)];
	GR1PROBPRF1 = GR1PROBP1((GR1PROBP1(:,3)>0),2);
	
	GR1PROBPR2 = randn(numel(GluR1DOTSP2(:,1)),1)-4;
	GR1PROBP2 = [GluR1DOTSP2(:,1) GluR1DOTSP2(:,2) (GluR1DOTSP2(:,1)+GR1PROBPR2)];
	GR1PROBPRF2 = GR1PROBP2((GR1PROBP2(:,3)>0),2);
	

% NOW UPDATE THE DIFFUSION RATES FOR THOSE GLUR1 RECEPTORS
% THAT MADE IT INTO A SLOT

	xysG1(:,GR1PROBPRF1) = xysG1(:,GR1PROBPRF1)*(.001);
    xysG1(:,GR1PROBPRF2) = xysG1(:,GR1PROBPRF2)*(.001);
    

    
    % HOW MANY SLOTS ARE NOW FILLED?
    G1PSD1FSLOTS = numel(GR1PROBPRF1);
    G1PSD2FSLOTS = numel(GR1PROBPRF2);

	


G1FSLOTS = [G1PSD1FSLOTS G1PSD2FSLOTS 0 0];
%================%

% if stepN == 401; keyboard; end

end

%----SLOTS GluR2----%
function [xysG2 GluR2_TdwellPSD GluR2_TdwellPERI GluR2_TdwellSPYN G2FSLOTS]... 
	= G2SLOTS(stepN,G2N,xysG2,...
    GluR2_TdwellPSD,GluR2_TdwellPERI,GluR2_TdwellSPYN,...
    GluR2SP1,GluR2SP2,SAPFmx1,SAPFmx2,G2SAPLOC1, G2SAPLOC2)

%===========================%   
PSD1slotN = numel(SAPFmx1);
PSD2slotN = numel(SAPFmx2);
%===========================%
%{
 
    [rPSD1,cDOT1,vPOLY1] = find(G2SAPLOC1);
    vrcPOLYDOT1 = [vPOLY1;rPSD1;cDOT1]';
    [rPSD2,cDOT2,vPOLY2] = find(G2SAPLOC2);
    vrcPOLYDOT2 = [vPOLY2;rPSD2;cDOT2]';
    
    [rPSD1IDpoly] = find(vrcPOLYDOT1(:,2)>0);   % PSD1 DOT-ID DATA
    [rPSD2IDpoly] = find(vrcPOLYDOT2(:,2)>0);   % PSD2 DOT-ID DATA
    vrcTIMEpoly1 = vrcPOLYDOT1(rPSD1IDpoly,:);  % PSD1 DWELL DATA ONLY
    vrcTIMEpoly2 = vrcPOLYDOT2(rPSD2IDpoly,:);  % PSD2 DWELL DATA ONLY
    vrcTIME1polys = sortrows(vrcTIMEpoly1,-1);  % PSD1: SORT DWELL DATA LONG TO SHORT
    vrcTIME2polys = sortrows(vrcTIMEpoly2,-1);  % PSD2: SORT DWELL DATA LONG TO SHORT
    
    if size(vrcTIME1polys,1) < PSD1slotN     % IF ALL SLOTS NOT FILLED MAKE BLANKS
        vrcTIME1polys = vertcat(vrcTIME1polys,zeros(PSD1slotN,3));
    end
    if size(vrcTIME2polys,1) < PSD2slotN
        vrcTIME2polys = vertcat(vrcTIME2polys,zeros(PSD2slotN,3));
    end 
    
 
% GET THE SLOT DATA FOR ONLY FOR GLUR2 AMPARS (FOR EACH PSD REGION SEPARATELY)
% THESE LABELS ARE A LITTLE MISLEADING BECAUSE HERE WE ARE ONLY
% INTERESTED IN MAKING SURE THE SLOT MX ID VALUES ARE RETAINED 
% TO MATCH THE DWELL TIME MX VALUES SORTED BELOW 
 
    [rPSD1,cDOT1,vTIME1] = find(GluR2SP1);
    vrcTIMEPSDDOT1 = [vTIME1;rPSD1;cDOT1]';
    [rPSD2,cDOT2,vTIME2] = find(GluR2SP2);
    vrcTIMEPSDDOT2 = [vTIME2;rPSD2;cDOT2]';
    
    [rPSD1ID] = find(vrcTIMEPSDDOT1(:,2)>0);    % PSD1 DOT-ID DATA
    [rPSD2ID] = find(vrcTIMEPSDDOT2(:,2)>0);    % PSD2 DOT-ID DATA
    vrcTIME1 = vrcTIMEPSDDOT1(rPSD1ID,:);       % PSD1 DWELL DATA ONLY
    vrcTIME2 = vrcTIMEPSDDOT2(rPSD2ID,:);       % PSD2 DWELL DATA ONLY
    vrcTIME1s = sortrows(vrcTIME1,-1);          % PSD1: SORT DWELL DATA LONG TO SHORT
    vrcTIME2s = sortrows(vrcTIME2,-1);          % PSD2: SORT DWELL DATA LONG TO SHORT
    
    if size(vrcTIME1s,1) < PSD1slotN     % IF ALL SLOTS NOT FILLED MAKE BLANKS
        vrcTIME1s = vertcat(vrcTIME1s,zeros(PSD1slotN,3));
    end
    if size(vrcTIME2s,1) < PSD2slotN
        vrcTIME2s = vertcat(vrcTIME2s,zeros(PSD2slotN,3));
    end 
 
 
 
 
    % (SLPSD1LOCkOFF: DOT)
    % SLPSD1LOCkOFF = SLPSD1((SLPSD1(:,1)>(2)),3);
    % SLPSD2LOCkOFF = SLPSD2((SLPSD2(:,1)>(2)),3);
    
    %}
    

% GET ME THE DWELL TIME VALUES AND THE MX ID
% LOCATIONS OF THOSE GLUR2 RECEPTORS IN PSD1 OR PSD2
 
    [r,c,v] = find(GluR2_TdwellSPYN);    % WHICH PSD, WHICH AMPAR, HOW LONG
    rcv = [v r c];                      % ONLY AMPARS WITH DWELL TIMES (TIME,PSD,DOT)
    [r1] = find(rcv(:,2)<2);            % PSD1 DATA
    [r2] = find(rcv(:,2)>1);            % PSD2 DATA
    p1r = rcv(r1,:);                    % PSD1 DWELL DATA ONLY
    p2r = rcv(r2,:);                    % PSD2 DWELL DATA ONLY
    p1rs = sortrows(p1r,-1);            % PSD1: SORT DWELL DATA LONG TO SHORT
    p2rs = sortrows(p2r,-1);            % PSD2: SORT DWELL DATA LONG TO SHORT
    
    if size(p1rs,1) < PSD1slotN     % IF ALL SLOTS NOT FILLED MAKE BLANKS
        p1rs = vertcat(p1rs,zeros(PSD1slotN,3));
    end
    if size(p2rs,1) < PSD2slotN
        p2rs = vertcat(p2rs,zeros(PSD2slotN,3));
    end 
    
    % (SLPSD: TIME,PSD,DOT)
    SLPSD1 = p1rs(1:PSD1slotN,:);   %PSD1: GET THE N(slots) LONGEST DWELLING AMPARS
    SLPSD2 = p2rs(1:PSD2slotN,:);   %PSD2: GET THE N(slots) LONGEST DWELLING AMPARS
    
    
% LETS SEE SOME IDs FELLAS
    
    % (SLPSD1LOCkON: DOT)
    SLPSD1LOCkON = SLPSD1((SLPSD1(:,1)>0),3);
    SLPSD2LOCkON = SLPSD2((SLPSD2(:,1)>0),3);
 
    
    G2SAPLOCD1 = [G2SAPLOC1]';
    G2SAPLOCD2 = [G2SAPLOC2]';
    GluR2DSP1 = [GluR2SP1]';
    GluR2DSP2 = [GluR2SP2]';
    
GluR2DOTSP1 = [GluR2DSP1(SLPSD1LOCkON,:) SLPSD1LOCkON G2SAPLOCD1(SLPSD1LOCkON,:)];
GluR2DOTSP2 = [GluR2DSP2(SLPSD2LOCkON,:) SLPSD2LOCkON G2SAPLOCD2(SLPSD2LOCkON,:)];
 
% THESE ARE THE SLOTTED DOTS
    [uniV1 uniR2 uniC1] = unique(GluR2DOTSP1(:,3),'rows','first');
    GluR2DOTSP1 = GluR2DOTSP1(uniR2,:);
    GluR2DOTSP1 = GluR2DOTSP1(:,1:2);
    
    [uniV2 uniR2 uniC2] = unique(GluR2DOTSP2(:,3),'rows','first');
    GluR2DOTSP2 = GluR2DOTSP2(uniR2,:);
    GluR2DOTSP2 = GluR2DOTSP2(:,1:2);
    
 
    
% HERE IS WHERE WE ROLL THE DICE FOR EACH GLUR2 PARTICLE IN A SLOT REGION
% TO DETERMINE WHETHER IT BINDS TO THE SLOT
 
    GR2PROBPR2 = randn(numel(GluR2DOTSP1(:,1)),1)-2;
    GR2PROBP1 = [GluR2DOTSP1(:,1) GluR2DOTSP1(:,2) (GluR2DOTSP1(:,1)+GR2PROBPR2)];
    GR2PROBPRF1 = GR2PROBP1((GR2PROBP1(:,3)>0),2);
    
    GR2PROBPR2 = randn(numel(GluR2DOTSP2(:,1)),1)-2;
    GR2PROBP2 = [GluR2DOTSP2(:,1) GluR2DOTSP2(:,2) (GluR2DOTSP2(:,1)+GR2PROBPR2)];
    GR2PROBPRF2 = GR2PROBP2((GR2PROBP2(:,3)>0),2);
    
 
% NOW UPDATE THE DIFFUSION RATES FOR THOSE GLUR2 RECEPTORS
% THAT MADE IT INTO A SLOT
 
    xysG2(:,GR2PROBPRF1) = xysG2(:,GR2PROBPRF1)*(.001);
    xysG2(:,GR2PROBPRF2) = xysG2(:,GR2PROBPRF2)*(.001);
    
 
    
    % HOW MANY SLOTS ARE NOW FILLED?
    G2PSD1FSLOTS = numel(GR2PROBPRF1);
    G2PSD2FSLOTS = numel(GR2PROBPRF2);
 
   


 

G2FSLOTS = [G2PSD1FSLOTS G2PSD2FSLOTS 0 0];
%================%
% if stepN == 401; keyboard; end
end
%===================================%

%}
%===================================%

%===================================%
%----PERI SLOTS GluR1----%
%{
%===================================%
%----PERI SLOTS GluR1----%
function [xysG1 GluR1_TdwellPSD GluR1_TdwellPERI G1FSLOTS...
	] = G1PERISLOTS(stepN,G1N,xysG1,...
	GluR1_TdwellPSD,GluR1_TdwellPERI,...
	ko11,ko12,ko15,ko16,...
	doKo5,S1sum,S2sum,SAP5)


useSlotsPERIGR1 = doKo5;
%===================%
if useSlotsPERIGR1
%-------------------%
KonSpi1PERIGR1 = ko11;
KoffSpi1PERIGR1 = ko12;
KonSpi2PERIGR1 = ko15;	
KoffSpi2PERIGR1 = ko16; 

%===========================%
% SETUP
%---------------------------%   
PERI1slotN = fix(S1sum/5*SAP5(2));
PERI2slotN = fix(S2sum/5*SAP5(2));
%===========================%
%{

	[rPSD1,cDOT1,vPOLY1] = find(G1SAPLOC1);
	vrcPOLYDOT1 = [vPOLY1;rPSD1;cDOT1]';
	[rPSD2,cDOT2,vPOLY2] = find(G1SAPLOC2);
	vrcPOLYDOT2 = [vPOLY2;rPSD2;cDOT2]';
	
	[rPSD1IDpoly] = find(vrcPOLYDOT1(:,2)>0);	% PSD1 DOT-ID DATA
    [rPSD2IDpoly] = find(vrcPOLYDOT2(:,2)>0);	% PSD2 DOT-ID DATA
    vrcTIMEpoly1 = vrcPOLYDOT1(rPSD1IDpoly,:);	% PSD1 DWELL DATA ONLY
    vrcTIMEpoly2 = vrcPOLYDOT2(rPSD2IDpoly,:);	% PSD2 DWELL DATA ONLY
    vrcTIME1polys = sortrows(vrcTIMEpoly1,-1);	% PSD1: SORT DWELL DATA LONG TO SHORT
    vrcTIME2polys = sortrows(vrcTIMEpoly2,-1);	% PSD2: SORT DWELL DATA LONG TO SHORT
	
	if size(vrcTIME1polys,1) < PSD1slotN     % IF ALL SLOTS NOT FILLED MAKE BLANKS
        vrcTIME1polys = vertcat(vrcTIME1polys,zeros(PSD1slotN,3));
    end
    if size(vrcTIME2polys,1) < PSD2slotN
        vrcTIME2polys = vertcat(vrcTIME2polys,zeros(PSD2slotN,3));
	end 
	

% GET THE SLOT DATA FOR ONLY FOR GLUR1 AMPARS (FOR EACH PSD REGION SEPARATELY)
% THESE LABELS ARE A LITTLE MISLEADING BECAUSE HERE WE ARE ONLY
% INTERESTED IN MAKING SURE THE SLOT MX ID VALUES ARE RETAINED 
% TO MATCH THE DWELL TIME MX VALUES SORTED BELOW 

	[rPSD1,cDOT1,vTIME1] = find(GluR1SP1);
	vrcTIMEPSDDOT1 = [vTIME1;rPSD1;cDOT1]';
	[rPSD2,cDOT2,vTIME2] = find(GluR1SP2);
	vrcTIMEPSDDOT2 = [vTIME2;rPSD2;cDOT2]';
	
	[rPSD1ID] = find(vrcTIMEPSDDOT1(:,2)>0);	% PSD1 DOT-ID DATA
    [rPSD2ID] = find(vrcTIMEPSDDOT2(:,2)>0);	% PSD2 DOT-ID DATA
    vrcTIME1 = vrcTIMEPSDDOT1(rPSD1ID,:);		% PSD1 DWELL DATA ONLY
    vrcTIME2 = vrcTIMEPSDDOT2(rPSD2ID,:);		% PSD2 DWELL DATA ONLY
    vrcTIME1s = sortrows(vrcTIME1,-1);			% PSD1: SORT DWELL DATA LONG TO SHORT
    vrcTIME2s = sortrows(vrcTIME2,-1);			% PSD2: SORT DWELL DATA LONG TO SHORT
	
	if size(vrcTIME1s,1) < PSD1slotN     % IF ALL SLOTS NOT FILLED MAKE BLANKS
        vrcTIME1s = vertcat(vrcTIME1s,zeros(PSD1slotN,3));
    end
    if size(vrcTIME2s,1) < PSD2slotN
        vrcTIME2s = vertcat(vrcTIME2s,zeros(PSD2slotN,3));
	end 




	% (SLPSD1LOCkOFF: DOT)
	% SLPSD1LOCkOFF = SLPSD1((SLPSD1(:,1)>(2)),3);
    % SLPSD2LOCkOFF = SLPSD2((SLPSD2(:,1)>(2)),3);
	
	%}
%{
%===================%
if useSlotsPSDGR1
%-------------------% 

% NOW GET ME THE DWELL TIME VALUES AND THE MX ID
% LOCATIONS OF THOSE GLUR1 RECEPTORS IN PSD1 OR PSD2

	[r,c,v] = find(GluR1_TdwellPSD);	% WHICH PSD, WHICH AMPAR, HOW LONG
    rcv = [v r c];                      % ONLY AMPARS WITH DWELL TIMES (TIME,PSD,DOT)
    [r1] = find(rcv(:,2)<2);            % PSD1 DATA
    [r2] = find(rcv(:,2)>1);            % PSD2 DATA
    p1r = rcv(r1,:);                    % PSD1 DWELL DATA ONLY
    p2r = rcv(r2,:);                    % PSD2 DWELL DATA ONLY
    p1rs = sortrows(p1r,-1);            % PSD1: SORT DWELL DATA LONG TO SHORT
    p2rs = sortrows(p2r,-1);            % PSD2: SORT DWELL DATA LONG TO SHORT
    
    if size(p1rs,1) < PSD1slotN     % IF ALL SLOTS NOT FILLED MAKE BLANKS
        p1rs = vertcat(p1rs,zeros(PSD1slotN,3));
    end
    if size(p2rs,1) < PSD2slotN
        p2rs = vertcat(p2rs,zeros(PSD2slotN,3));
	end 
    
	% (SLPSD: TIME,PSD,DOT)
    SLPSD1 = p1rs(1:PSD1slotN,:);   %PSD1: GET THE N(slots) LONGEST DWELLING AMPARS
    SLPSD2 = p2rs(1:PSD2slotN,:);   %PSD2: GET THE N(slots) LONGEST DWELLING AMPARS
	
	
% LETS SEE SOME IDs FELLAS
	
	% (SLPSD1LOCkON: DOT)
	SLPSD1LOCkON = SLPSD1((SLPSD1(:,1)>0),3);
    SLPSD2LOCkON = SLPSD2((SLPSD2(:,1)>0),3);

	
	G1SAPLOCD1 = [G1SAPLOC1]';
	G1SAPLOCD2 = [G1SAPLOC2]';
	GluR1DSP1 = [GluR1SP1]';
	GluR1DSP2 = [GluR1SP2]';
	
GluR1DOTSP1 = [GluR1DSP1(SLPSD1LOCkON,:) SLPSD1LOCkON G1SAPLOCD1(SLPSD1LOCkON,:)];
GluR1DOTSP2 = [GluR1DSP2(SLPSD2LOCkON,:) SLPSD2LOCkON G1SAPLOCD2(SLPSD2LOCkON,:)];

% THESE ARE THE SLOTTED DOTS
	[uniV1 uniR1 uniC1] = unique(GluR1DOTSP1(:,3),'rows','first');
	GluR1DOTSP1 = GluR1DOTSP1(uniR1,:);
	GluR1DOTSP1 = GluR1DOTSP1(:,1:2);
	
	[uniV2 uniR2 uniC2] = unique(GluR1DOTSP2(:,3),'rows','first');
	GluR1DOTSP2 = GluR1DOTSP2(uniR2,:);
	GluR1DOTSP2 = GluR1DOTSP2(:,1:2);
	

	
% HERE IS WHERE WE ROLL THE DICE FOR EACH GLUR1 PARTICLE IN A SLOT REGION
% TO DETERMINE WHETHER IT BINDS TO THE SLOT

	GR1PROBPR1 = randn(numel(GluR1DOTSP1(:,1)),1)-4;
	GR1PROBP1 = [GluR1DOTSP1(:,1) GluR1DOTSP1(:,2) (GluR1DOTSP1(:,1)+GR1PROBPR1)];
	GR1PROBPRF1 = GR1PROBP1((GR1PROBP1(:,3)>0),2);
	
	GR1PROBPR2 = randn(numel(GluR1DOTSP2(:,1)),1)-4;
	GR1PROBP2 = [GluR1DOTSP2(:,1) GluR1DOTSP2(:,2) (GluR1DOTSP2(:,1)+GR1PROBPR2)];
	GR1PROBPRF2 = GR1PROBP2((GR1PROBP2(:,3)>0),2);
	

% NOW UPDATE THE DIFFUSION RATES FOR THOSE GLUR1 RECEPTORS
% THAT MADE IT INTO A SLOT

	xysG1(:,GR1PROBPRF1) = xysG1(:,GR1PROBPRF1)*(.001);
    xysG1(:,GR1PROBPRF2) = xysG1(:,GR1PROBPRF2)*(.001);
    

    
    % HOW MANY SLOTS ARE NOW FILLED?
    G1PSD1FSLOTS = numel(GR1PROBPRF1);
    G1PSD2FSLOTS = numel(GR1PROBPRF2);

	
% if stepN == 100; keyboard; end
%-------------------%   
end %useSlotsPSDGR1 
%===================%
%}
%===========================%


[rpx,cpx,vpx] = find(GluR1_TdwellPERI);
	rcv = [vpx rpx cpx]; % VALUE ROW COLUMN
	[r1] = find(rcv(:,2)<2);
	[r2] = find(rcv(:,2)>1);
	p1r = rcv(r1,:);
	p2r = rcv(r2,:);
	p1rs = sortrows(p1r,-1);
	p2rs = sortrows(p2r,-1);
	if size(p1rs,1) < PERI1slotN
		p1rs = vertcat(p1rs,zeros(PERI1slotN,3));
	end
	if size(p2rs,1) < PERI2slotN
		p2rs = vertcat(p2rs,zeros(PERI2slotN,3));
	end 



	SLPERI1 = p1rs(1:PERI1slotN,:); %SLOT LOCATION PERI1
	SLPERI2 = p2rs(1:PERI2slotN,:); %SLOT LOCATION PERI2

	SLPERI1LOCkon = SLPERI1((SLPERI1(:,1)>KonSpi1PERIGR1),3);
	SLPERI2LOCkon = SLPERI2((SLPERI2(:,1)>KonSpi2PERIGR1),3);
        
    SLPERI1LOCkoff = SLPERI1((SLPERI1(:,1)>(KoffSpi1PERIGR1+KonSpi1PERIGR1)),3);
    SLPERI2LOCkoff = SLPERI2((SLPERI2(:,1)>(KoffSpi2PERIGR1+KonSpi2PERIGR1)),3);
    
    xysG1(:,SLPERI1LOCkon) = xysG1(:,SLPERI1LOCkon)*(.001);
    xysG1(:,SLPERI2LOCkon) = xysG1(:,SLPERI2LOCkon)*(.001);
    xysG1(:,SLPERI1LOCkoff) = xysG1(:,SLPERI1LOCkoff)/(.001);
    xysG1(:,SLPERI2LOCkoff) = xysG1(:,SLPERI2LOCkoff)/(.001);

	G1PERI1FSLOTS = numel(SLPERI1LOCkon)-numel(SLPERI1LOCkoff);
	G1PERI2FSLOTS = numel(SLPERI2LOCkon)-numel(SLPERI2LOCkoff);


%================%
G1FSLOTS = [0 0 G1PERI1FSLOTS G1PERI2FSLOTS];
%================%
%-------------------%   
end %useSlotsPSDGR1
%===================%
end

%----PERI SLOTS GluR2----%
function [xysG2 GluR2_TdwellPSD GluR2_TdwellPERI G2FSLOTS...
    ] = G2PERISLOTS(stepN,G2N,xysG2,...
    GluR2_TdwellPSD,GluR2_TdwellPERI,...
    ko11,ko12,ko15,ko16,...
    doKo6,S1sum,S2sum,SAP5)
 
 
useSlotsPERIGR2 = doKo6;
%===================%
if useSlotsPERIGR2
%-------------------%
KonSpi1PERIGR2 = ko11;
KoffSpi1PERIGR2 = ko12;
KonSpi2PERIGR2 = ko15;
KoffSpi2PERIGR2 = ko16; 

%===========================%
% SETUP
%---------------------------%   
PERI1slotN = fix(S1sum/5*SAP5(2));
PERI2slotN = fix(S2sum/5*SAP5(2));
%===========================%  
%{
[rPSD1,cDOT1,vPOLY1] = find(G2SAPLOC1);
    vrcPOLYDOT1 = [vPOLY1;rPSD1;cDOT1]';
    [rPSD2,cDOT2,vPOLY2] = find(G2SAPLOC2);
    vrcPOLYDOT2 = [vPOLY2;rPSD2;cDOT2]';
    
    [rPSD1IDpoly] = find(vrcPOLYDOT1(:,2)>0);   % PSD1 DOT-ID DATA
    [rPSD2IDpoly] = find(vrcPOLYDOT2(:,2)>0);   % PSD2 DOT-ID DATA
    vrcTIMEpoly1 = vrcPOLYDOT1(rPSD1IDpoly,:);  % PSD1 DWELL DATA ONLY
    vrcTIMEpoly2 = vrcPOLYDOT2(rPSD2IDpoly,:);  % PSD2 DWELL DATA ONLY
    vrcTIME1polys = sortrows(vrcTIMEpoly1,-1);  % PSD1: SORT DWELL DATA LONG TO SHORT
    vrcTIME2polys = sortrows(vrcTIMEpoly2,-1);  % PSD2: SORT DWELL DATA LONG TO SHORT
    
    if size(vrcTIME1polys,1) < PSD1slotN     % IF ALL SLOTS NOT FILLED MAKE BLANKS
        vrcTIME1polys = vertcat(vrcTIME1polys,zeros(PSD1slotN,3));
    end
    if size(vrcTIME2polys,1) < PSD2slotN
        vrcTIME2polys = vertcat(vrcTIME2polys,zeros(PSD2slotN,3));
    end 
    
 
% GET THE SLOT DATA FOR ONLY FOR GLUR2 AMPARS (FOR EACH PSD REGION SEPARATELY)
% THESE LABELS ARE A LITTLE MISLEADING BECAUSE HERE WE ARE ONLY
% INTERESTED IN MAKING SURE THE SLOT MX ID VALUES ARE RETAINED 
% TO MATCH THE DWELL TIME MX VALUES SORTED BELOW 
 
    [rPSD1,cDOT1,vTIME1] = find(GluR2SP1);
    vrcTIMEPSDDOT1 = [vTIME1;rPSD1;cDOT1]';
    [rPSD2,cDOT2,vTIME2] = find(GluR2SP2);
    vrcTIMEPSDDOT2 = [vTIME2;rPSD2;cDOT2]';
    
    [rPSD1ID] = find(vrcTIMEPSDDOT1(:,2)>0);    % PSD1 DOT-ID DATA
    [rPSD2ID] = find(vrcTIMEPSDDOT2(:,2)>0);    % PSD2 DOT-ID DATA
    vrcTIME1 = vrcTIMEPSDDOT1(rPSD1ID,:);       % PSD1 DWELL DATA ONLY
    vrcTIME2 = vrcTIMEPSDDOT2(rPSD2ID,:);       % PSD2 DWELL DATA ONLY
    vrcTIME1s = sortrows(vrcTIME1,-1);          % PSD1: SORT DWELL DATA LONG TO SHORT
    vrcTIME2s = sortrows(vrcTIME2,-1);          % PSD2: SORT DWELL DATA LONG TO SHORT
    
    if size(vrcTIME1s,1) < PSD1slotN     % IF ALL SLOTS NOT FILLED MAKE BLANKS
        vrcTIME1s = vertcat(vrcTIME1s,zeros(PSD1slotN,3));
    end
    if size(vrcTIME2s,1) < PSD2slotN
        vrcTIME2s = vertcat(vrcTIME2s,zeros(PSD2slotN,3));
    end 
 
 
 
 
    % (SLPSD1LOCkOFF: DOT)
    % SLPSD1LOCkOFF = SLPSD1((SLPSD1(:,1)>(2)),3);
    % SLPSD2LOCkOFF = SLPSD2((SLPSD2(:,1)>(2)),3);
    
    %}
%{
%===================%
if useSlotsPSDGR2
%-------------------% 

% NOW GET ME THE DWELL TIME VALUES AND THE MX ID
% LOCATIONS OF THOSE GLUR2 RECEPTORS IN PSD1 OR PSD2
 
    [r,c,v] = find(GluR2_TdwellPSD);    % WHICH PSD, WHICH AMPAR, HOW LONG
    rcv = [v r c];                      % ONLY AMPARS WITH DWELL TIMES (TIME,PSD,DOT)
    [r1] = find(rcv(:,2)<2);            % PSD1 DATA
    [r2] = find(rcv(:,2)>1);            % PSD2 DATA
    p1r = rcv(r1,:);                    % PSD1 DWELL DATA ONLY
    p2r = rcv(r2,:);                    % PSD2 DWELL DATA ONLY
    p1rs = sortrows(p1r,-1);            % PSD1: SORT DWELL DATA LONG TO SHORT
    p2rs = sortrows(p2r,-1);            % PSD2: SORT DWELL DATA LONG TO SHORT
    
    if size(p1rs,1) < PSD1slotN     % IF ALL SLOTS NOT FILLED MAKE BLANKS
        p1rs = vertcat(p1rs,zeros(PSD1slotN,3));
    end
    if size(p2rs,1) < PSD2slotN
        p2rs = vertcat(p2rs,zeros(PSD2slotN,3));
    end 
    
    % (SLPSD: TIME,PSD,DOT)
    SLPSD1 = p1rs(1:PSD1slotN,:);   %PSD1: GET THE N(slots) LONGEST DWELLING AMPARS
    SLPSD2 = p2rs(1:PSD2slotN,:);   %PSD2: GET THE N(slots) LONGEST DWELLING AMPARS
    
    
% LETS SEE SOME IDs FELLAS
    
    % (SLPSD1LOCkON: DOT)
    SLPSD1LOCkON = SLPSD1((SLPSD1(:,1)>0),3);
    SLPSD2LOCkON = SLPSD2((SLPSD2(:,1)>0),3);
 
    
    G2SAPLOCD1 = [G2SAPLOC1]';
    G2SAPLOCD2 = [G2SAPLOC2]';
    GluR2DSP1 = [GluR2SP1]';
    GluR2DSP2 = [GluR2SP2]';
    
GluR2DOTSP1 = [GluR2DSP1(SLPSD1LOCkON,:) SLPSD1LOCkON G2SAPLOCD1(SLPSD1LOCkON,:)];
GluR2DOTSP2 = [GluR2DSP2(SLPSD2LOCkON,:) SLPSD2LOCkON G2SAPLOCD2(SLPSD2LOCkON,:)];
 
% THESE ARE THE SLOTTED DOTS
    [uniV1 uniR2 uniC1] = unique(GluR2DOTSP1(:,3),'rows','first');
    GluR2DOTSP1 = GluR2DOTSP1(uniR2,:);
    GluR2DOTSP1 = GluR2DOTSP1(:,1:2);
    
    [uniV2 uniR2 uniC2] = unique(GluR2DOTSP2(:,3),'rows','first');
    GluR2DOTSP2 = GluR2DOTSP2(uniR2,:);
    GluR2DOTSP2 = GluR2DOTSP2(:,1:2);
    
 
    
% HERE IS WHERE WE ROLL THE DICE FOR EACH GLUR2 PARTICLE IN A SLOT REGION
% TO DETERMINE WHETHER IT BINDS TO THE SLOT
 
    GR2PROBPR2 = randn(numel(GluR2DOTSP1(:,1)),1)-2;
    GR2PROBP1 = [GluR2DOTSP1(:,1) GluR2DOTSP1(:,2) (GluR2DOTSP1(:,1)+GR2PROBPR2)];
    GR2PROBPRF1 = GR2PROBP1((GR2PROBP1(:,3)>0),2);
    
    GR2PROBPR2 = randn(numel(GluR2DOTSP2(:,1)),1)-2;
    GR2PROBP2 = [GluR2DOTSP2(:,1) GluR2DOTSP2(:,2) (GluR2DOTSP2(:,1)+GR2PROBPR2)];
    GR2PROBPRF2 = GR2PROBP2((GR2PROBP2(:,3)>0),2);
    
 
% NOW UPDATE THE DIFFUSION RATES FOR THOSE GLUR2 RECEPTORS
% THAT MADE IT INTO A SLOT
 
    xysG2(:,GR2PROBPRF1) = xysG2(:,GR2PROBPRF1)*(.001);
    xysG2(:,GR2PROBPRF2) = xysG2(:,GR2PROBPRF2)*(.001);
    
 
    
    % HOW MANY SLOTS ARE NOW FILLED?
    G2PSD1FSLOTS = numel(GR2PROBPRF1);
    G2PSD2FSLOTS = numel(GR2PROBPRF2);
 
   

%-------------------%   
end %useSlotsPSDGR2 
%===================%

%}
%===========================%

[rpx,cpx,vpx] = find(GluR2_TdwellPERI);
	rcv = [vpx rpx cpx]; % VALUE ROW COLUMN
	[r1] = find(rcv(:,2)<2);
	[r2] = find(rcv(:,2)>1);
	p1r = rcv(r1,:);
	p2r = rcv(r2,:);
	p1rs = sortrows(p1r,-1);
	p2rs = sortrows(p2r,-1);
	if size(p1rs,1) < PERI1slotN
		p1rs = vertcat(p1rs,zeros(PERI1slotN,3));
	end
	if size(p2rs,1) < PERI2slotN
		p2rs = vertcat(p2rs,zeros(PERI2slotN,3));
	end 



	SLPERI1 = p1rs(1:PERI1slotN,:); %SLOT LOCATION PERI1
	SLPERI2 = p2rs(1:PERI2slotN,:); %SLOT LOCATION PERI2

	SLPERI1LOCkon = SLPERI1((SLPERI1(:,1)>KonSpi1PERIGR2),3);
	SLPERI2LOCkon = SLPERI2((SLPERI2(:,1)>KonSpi2PERIGR2),3);
        
    SLPERI1LOCkoff = SLPERI1((SLPERI1(:,1)>(KoffSpi1PERIGR2+KonSpi1PERIGR2)),3);
    SLPERI2LOCkoff = SLPERI2((SLPERI2(:,1)>(KoffSpi2PERIGR2+KonSpi2PERIGR2)),3);
    
    xysG2(:,SLPERI1LOCkon) = xysG2(:,SLPERI1LOCkon)*(.001);
    xysG2(:,SLPERI2LOCkon) = xysG2(:,SLPERI2LOCkon)*(.001);
    xysG2(:,SLPERI1LOCkoff) = xysG2(:,SLPERI1LOCkoff)/(.001);
    xysG2(:,SLPERI2LOCkoff) = xysG2(:,SLPERI2LOCkoff)/(.001);
 
    G2PERI1FSLOTS = numel(SLPERI1LOCkon)-numel(SLPERI1LOCkoff);
    G2PERI2FSLOTS = numel(SLPERI2LOCkon)-numel(SLPERI2LOCkoff);

 
%===================%
G2FSLOTS = [0 0 G2PERI1FSLOTS G2PERI2FSLOTS];
%===================%   
end %useSlotsPSDGR2
%===================%
end
%===================================%

%}
%===================================%

%===================================%
%----OLD SLOTS GluR1----%
%{
%===================================%
%----OLD SLOTS GluR1----%
function [xysG1 GluR1_TdwellPSD GluR1_TdwellPERI G1FSLOTS...
	] = G1SLOTSOLD(stepN,G1N,xysG1,...
	D,t,GluR1_TdwellPSD,GluR1_TdwellPERI,...
	ko1,ko2,ko3,ko4,ko5,ko6,ko7,ko8,ko9,ko10,ko11,ko12,ko13,ko14,ko15,ko16,...
	doKo1,doKo2,doKo3,doKo4,doKo5,doKo6,S1sum,S2sum,SAP5,...
	GluR1SP1,GluR1SP2,SAPFmx1,SAPFmx2,G1SAPLOC1, G1SAPLOC2)




%===========================%
% Slots
%---------------------------%   
PSD1slotN = fix(S1sum/5*SAP5(1));
PSD2slotN = fix(S2sum/5*SAP5(1));
PERI1slotN = fix(S1sum/5*SAP5(2));
PERI2slotN = fix(S2sum/5*SAP5(2));


PSD1slotN = numel(SAPFmx1);
PSD2slotN = numel(SAPFmx2);


G1PSD1FSLOTS = 0;
G1PSD2FSLOTS = 0;
G1PERI1FSLOTS = 0;
G1PERI2FSLOTS = 0;

useGluR1slots = doKo1;	useSlotsPSDGR1 = doKo3;	useSlotsPERIGR1 = doKo5;
KonSpi1PSDGR1 = ko9;	KoffSpi1PSDGR1 = ko10;	KonSpi1PERIGR1 = ko11;
KoffSpi1PERIGR1 = ko12;	KonSpi2PSDGR1 = ko13;	KoffSpi2PSDGR1 = ko14;
KonSpi2PERIGR1 = ko15;	KoffSpi2PERIGR1 = ko16; 
%===========================%

%===================%
if useSlotsPSDGR1
%-------------------%   

%{

	[rPSD1,cDOT1,vPOLY1] = find(G1SAPLOC1);
	vrcPOLYDOT1 = [vPOLY1;rPSD1;cDOT1]';
	[rPSD2,cDOT2,vPOLY2] = find(G1SAPLOC2);
	vrcPOLYDOT2 = [vPOLY2;rPSD2;cDOT2]';
	
	[rPSD1IDpoly] = find(vrcPOLYDOT1(:,2)>0);	% PSD1 DOT-ID DATA
    [rPSD2IDpoly] = find(vrcPOLYDOT2(:,2)>0);	% PSD2 DOT-ID DATA
    vrcTIMEpoly1 = vrcPOLYDOT1(rPSD1IDpoly,:);	% PSD1 DWELL DATA ONLY
    vrcTIMEpoly2 = vrcPOLYDOT2(rPSD2IDpoly,:);	% PSD2 DWELL DATA ONLY
    vrcTIME1polys = sortrows(vrcTIMEpoly1,-1);	% PSD1: SORT DWELL DATA LONG TO SHORT
    vrcTIME2polys = sortrows(vrcTIMEpoly2,-1);	% PSD2: SORT DWELL DATA LONG TO SHORT
	
	if size(vrcTIME1polys,1) < PSD1slotN     % IF ALL SLOTS NOT FILLED MAKE BLANKS
        vrcTIME1polys = vertcat(vrcTIME1polys,zeros(PSD1slotN,3));
    end
    if size(vrcTIME2polys,1) < PSD2slotN
        vrcTIME2polys = vertcat(vrcTIME2polys,zeros(PSD2slotN,3));
	end 
	

% GET THE SLOT DATA FOR ONLY FOR GLUR1 AMPARS (FOR EACH PSD REGION SEPARATELY)
% THESE LABELS ARE A LITTLE MISLEADING BECAUSE HERE WE ARE ONLY
% INTERESTED IN MAKING SURE THE SLOT MX ID VALUES ARE RETAINED 
% TO MATCH THE DWELL TIME MX VALUES SORTED BELOW 

	[rPSD1,cDOT1,vTIME1] = find(GluR1SP1);
	vrcTIMEPSDDOT1 = [vTIME1;rPSD1;cDOT1]';
	[rPSD2,cDOT2,vTIME2] = find(GluR1SP2);
	vrcTIMEPSDDOT2 = [vTIME2;rPSD2;cDOT2]';
	
	[rPSD1ID] = find(vrcTIMEPSDDOT1(:,2)>0);	% PSD1 DOT-ID DATA
    [rPSD2ID] = find(vrcTIMEPSDDOT2(:,2)>0);	% PSD2 DOT-ID DATA
    vrcTIME1 = vrcTIMEPSDDOT1(rPSD1ID,:);		% PSD1 DWELL DATA ONLY
    vrcTIME2 = vrcTIMEPSDDOT2(rPSD2ID,:);		% PSD2 DWELL DATA ONLY
    vrcTIME1s = sortrows(vrcTIME1,-1);			% PSD1: SORT DWELL DATA LONG TO SHORT
    vrcTIME2s = sortrows(vrcTIME2,-1);			% PSD2: SORT DWELL DATA LONG TO SHORT
	
	if size(vrcTIME1s,1) < PSD1slotN     % IF ALL SLOTS NOT FILLED MAKE BLANKS
        vrcTIME1s = vertcat(vrcTIME1s,zeros(PSD1slotN,3));
    end
    if size(vrcTIME2s,1) < PSD2slotN
        vrcTIME2s = vertcat(vrcTIME2s,zeros(PSD2slotN,3));
	end 




	% (SLPSD1LOCkOFF: DOT)
	% SLPSD1LOCkOFF = SLPSD1((SLPSD1(:,1)>(2)),3);
    % SLPSD2LOCkOFF = SLPSD2((SLPSD2(:,1)>(2)),3);
	
	%}
	
% NOW GET ME THE DWELL TIME VALUES AND THE MX ID
% LOCATIONS OF THOSE GLUR1 RECEPTORS IN PSD1 OR PSD2

	[r,c,v] = find(GluR1_TdwellPSD);	% WHICH PSD, WHICH AMPAR, HOW LONG
    rcv = [v r c];                      % ONLY AMPARS WITH DWELL TIMES (TIME,PSD,DOT)
    [r1] = find(rcv(:,2)<2);            % PSD1 DATA
    [r2] = find(rcv(:,2)>1);            % PSD2 DATA
    p1r = rcv(r1,:);                    % PSD1 DWELL DATA ONLY
    p2r = rcv(r2,:);                    % PSD2 DWELL DATA ONLY
    p1rs = sortrows(p1r,-1);            % PSD1: SORT DWELL DATA LONG TO SHORT
    p2rs = sortrows(p2r,-1);            % PSD2: SORT DWELL DATA LONG TO SHORT
    
    if size(p1rs,1) < PSD1slotN     % IF ALL SLOTS NOT FILLED MAKE BLANKS
        p1rs = vertcat(p1rs,zeros(PSD1slotN,3));
    end
    if size(p2rs,1) < PSD2slotN
        p2rs = vertcat(p2rs,zeros(PSD2slotN,3));
	end 
    
	% (SLPSD: TIME,PSD,DOT)
    SLPSD1 = p1rs(1:PSD1slotN,:);   %PSD1: GET THE N(slots) LONGEST DWELLING AMPARS
    SLPSD2 = p2rs(1:PSD2slotN,:);   %PSD2: GET THE N(slots) LONGEST DWELLING AMPARS
	
	
% LETS SEE SOME IDs FELLAS
	
	% (SLPSD1LOCkON: DOT)
	SLPSD1LOCkON = SLPSD1((SLPSD1(:,1)>0),3);
    SLPSD2LOCkON = SLPSD2((SLPSD2(:,1)>0),3);

	
	G1SAPLOCD1 = [G1SAPLOC1]';
	G1SAPLOCD2 = [G1SAPLOC2]';
	GluR1DSP1 = [GluR1SP1]';
	GluR1DSP2 = [GluR1SP2]';
	
GluR1DOTSP1 = [GluR1DSP1(SLPSD1LOCkON,:) SLPSD1LOCkON G1SAPLOCD1(SLPSD1LOCkON,:)];
GluR1DOTSP2 = [GluR1DSP2(SLPSD2LOCkON,:) SLPSD2LOCkON G1SAPLOCD2(SLPSD2LOCkON,:)];

% THESE ARE THE SLOTTED DOTS
	[uniV1 uniR1 uniC1] = unique(GluR1DOTSP1(:,3),'rows','first');
	GluR1DOTSP1 = GluR1DOTSP1(uniR1,:);
	GluR1DOTSP1 = GluR1DOTSP1(:,1:2);
	
	[uniV2 uniR2 uniC2] = unique(GluR1DOTSP2(:,3),'rows','first');
	GluR1DOTSP2 = GluR1DOTSP2(uniR2,:);
	GluR1DOTSP2 = GluR1DOTSP2(:,1:2);
	

	
% HERE IS WHERE WE ROLL THE DICE FOR EACH GLUR1 PARTICLE IN A SLOT REGION
% TO DETERMINE WHETHER IT BINDS TO THE SLOT

	GR1PROBPR1 = randn(numel(GluR1DOTSP1(:,1)),1)-4;
	GR1PROBP1 = [GluR1DOTSP1(:,1) GluR1DOTSP1(:,2) (GluR1DOTSP1(:,1)+GR1PROBPR1)];
	GR1PROBPRF1 = GR1PROBP1((GR1PROBP1(:,3)>0),2);
	
	GR1PROBPR2 = randn(numel(GluR1DOTSP2(:,1)),1)-4;
	GR1PROBP2 = [GluR1DOTSP2(:,1) GluR1DOTSP2(:,2) (GluR1DOTSP2(:,1)+GR1PROBPR2)];
	GR1PROBPRF2 = GR1PROBP2((GR1PROBP2(:,3)>0),2);
	

% NOW UPDATE THE DIFFUSION RATES FOR THOSE GLUR1 RECEPTORS
% THAT MADE IT INTO A SLOT

	xysG1(:,GR1PROBPRF1) = xysG1(:,GR1PROBPRF1)*(.001);
    xysG1(:,GR1PROBPRF2) = xysG1(:,GR1PROBPRF2)*(.001);
    

    
    % HOW MANY SLOTS ARE NOW FILLED?
    G1PSD1FSLOTS = numel(GR1PROBPRF1);
    G1PSD2FSLOTS = numel(GR1PROBPRF2);

	
% if stepN == 200; keyboard; end
%-------------------%   
end %useSlotsPSDGR1 
%===================%
 



%===================%
if useSlotsPERIGR1
%-------------------% 
    
        [rpx,cpx,vpx] = find(GluR1_TdwellPERI);
        rcv = [vpx rpx cpx]; % VALUE ROW COLUMN
        [r1] = find(rcv(:,2)<2);
        [r2] = find(rcv(:,2)>1);
        p1r = rcv(r1,:);
        p2r = rcv(r2,:);
        p1rs = sortrows(p1r,-1);
        p2rs = sortrows(p2r,-1);
        if size(p1rs,1) < PERI1slotN
            p1rs = vertcat(p1rs,zeros(PERI1slotN,3));
        end
        if size(p2rs,1) < PERI2slotN
            p2rs = vertcat(p2rs,zeros(PERI2slotN,3));
        end 
            
         
        
        SLPERI1 = p1rs(1:PERI1slotN,:); %SLOT LOCATION PERI1
        SLPERI2 = p2rs(1:PERI2slotN,:); %SLOT LOCATION PERI2
        
        SLPERI1LOCkon = SLPERI1((SLPERI1(:,1)>KonSpi1PERIGR1),3);
        SLPERI2LOCkon = SLPERI2((SLPERI2(:,1)>KonSpi2PERIGR1),3);
        
    SLPERI1LOCkoff = SLPERI1((SLPERI1(:,1)>(KoffSpi1PERIGR1+KonSpi1PERIGR1)),3);
    SLPERI2LOCkoff = SLPERI2((SLPERI2(:,1)>(KoffSpi2PERIGR1+KonSpi2PERIGR1)),3);
    
    xysG1(:,SLPERI1LOCkon) = xysG1(:,SLPERI1LOCkon)*(.001);
    xysG1(:,SLPERI2LOCkon) = xysG1(:,SLPERI2LOCkon)*(.001);
    xysG1(:,SLPERI1LOCkoff) = xysG1(:,SLPERI1LOCkoff)/(.001);
    xysG1(:,SLPERI2LOCkoff) = xysG1(:,SLPERI2LOCkoff)/(.001);

	G1PERI1FSLOTS = numel(SLPERI1LOCkon)-numel(SLPERI1LOCkoff);
	G1PERI2FSLOTS = numel(SLPERI2LOCkon)-numel(SLPERI2LOCkoff);
%-------------------%   
end %useSlotsPSDGR1
%===================%
 

G1FSLOTS = [G1PSD1FSLOTS G1PSD2FSLOTS G1PERI1FSLOTS G1PERI2FSLOTS];
%================%
end

%----OLD SLOTS GluR2----%
function [xysG2 GluR2_TdwellPSD GluR2_TdwellPERI G2FSLOTS...
    ] = G2SLOTSOLD(stepN,G2N,xysG2,...
    D,t,GluR2_TdwellPSD,GluR2_TdwellPERI,...
    ko1,ko2,ko3,ko4,ko5,ko6,ko7,ko8,ko9,ko10,ko11,ko12,ko13,ko14,ko15,ko16,...
    doKo1,doKo2,doKo3,doKo4,doKo5,doKo6,S1sum,S2sum,SAP5,...
    GluR2SP1,GluR2SP2,SAPFmx1,SAPFmx2,G2SAPLOC1, G2SAPLOC2)
 
 
 
 
%===========================%
% Slots
%---------------------------%   
PSD1slotN = fix(S1sum/5*SAP5(1));
PSD2slotN = fix(S2sum/5*SAP5(1));
PERI1slotN = fix(S1sum/5*SAP5(2));
PERI2slotN = fix(S2sum/5*SAP5(2));
 
 
PSD1slotN = numel(SAPFmx1);
PSD2slotN = numel(SAPFmx2);
 
 
G2PSD1FSLOTS = 0;
G2PSD2FSLOTS = 0;
G2PERI1FSLOTS = 0;
G2PERI2FSLOTS = 0;
 
useGluR2slots = doKo2;  useSlotsPSDGR2 = doKo4; useSlotsPERIGR2 = doKo5;
KonSpi1PSDGR2 = ko9;    KoffSpi1PSDGR2 = ko10;  KonSpi1PERIGR2 = ko11;
KoffSpi1PERIGR2 = ko12; KonSpi2PSDGR2 = ko13;   KoffSpi2PSDGR2 = ko14;
KonSpi2PERIGR2 = ko15;  KoffSpi2PERIGR2 = ko16; 
%===========================%
 
%===================%
if useSlotsPSDGR2
%-------------------%   
 
%{
 
    [rPSD1,cDOT1,vPOLY1] = find(G2SAPLOC1);
    vrcPOLYDOT1 = [vPOLY1;rPSD1;cDOT1]';
    [rPSD2,cDOT2,vPOLY2] = find(G2SAPLOC2);
    vrcPOLYDOT2 = [vPOLY2;rPSD2;cDOT2]';
    
    [rPSD1IDpoly] = find(vrcPOLYDOT1(:,2)>0);   % PSD1 DOT-ID DATA
    [rPSD2IDpoly] = find(vrcPOLYDOT2(:,2)>0);   % PSD2 DOT-ID DATA
    vrcTIMEpoly1 = vrcPOLYDOT1(rPSD1IDpoly,:);  % PSD1 DWELL DATA ONLY
    vrcTIMEpoly2 = vrcPOLYDOT2(rPSD2IDpoly,:);  % PSD2 DWELL DATA ONLY
    vrcTIME1polys = sortrows(vrcTIMEpoly1,-1);  % PSD1: SORT DWELL DATA LONG TO SHORT
    vrcTIME2polys = sortrows(vrcTIMEpoly2,-1);  % PSD2: SORT DWELL DATA LONG TO SHORT
    
    if size(vrcTIME1polys,1) < PSD1slotN     % IF ALL SLOTS NOT FILLED MAKE BLANKS
        vrcTIME1polys = vertcat(vrcTIME1polys,zeros(PSD1slotN,3));
    end
    if size(vrcTIME2polys,1) < PSD2slotN
        vrcTIME2polys = vertcat(vrcTIME2polys,zeros(PSD2slotN,3));
    end 
    
 
% GET THE SLOT DATA FOR ONLY FOR GLUR2 AMPARS (FOR EACH PSD REGION SEPARATELY)
% THESE LABELS ARE A LITTLE MISLEADING BECAUSE HERE WE ARE ONLY
% INTERESTED IN MAKING SURE THE SLOT MX ID VALUES ARE RETAINED 
% TO MATCH THE DWELL TIME MX VALUES SORTED BELOW 
 
    [rPSD1,cDOT1,vTIME1] = find(GluR2SP1);
    vrcTIMEPSDDOT1 = [vTIME1;rPSD1;cDOT1]';
    [rPSD2,cDOT2,vTIME2] = find(GluR2SP2);
    vrcTIMEPSDDOT2 = [vTIME2;rPSD2;cDOT2]';
    
    [rPSD1ID] = find(vrcTIMEPSDDOT1(:,2)>0);    % PSD1 DOT-ID DATA
    [rPSD2ID] = find(vrcTIMEPSDDOT2(:,2)>0);    % PSD2 DOT-ID DATA
    vrcTIME1 = vrcTIMEPSDDOT1(rPSD1ID,:);       % PSD1 DWELL DATA ONLY
    vrcTIME2 = vrcTIMEPSDDOT2(rPSD2ID,:);       % PSD2 DWELL DATA ONLY
    vrcTIME1s = sortrows(vrcTIME1,-1);          % PSD1: SORT DWELL DATA LONG TO SHORT
    vrcTIME2s = sortrows(vrcTIME2,-1);          % PSD2: SORT DWELL DATA LONG TO SHORT
    
    if size(vrcTIME1s,1) < PSD1slotN     % IF ALL SLOTS NOT FILLED MAKE BLANKS
        vrcTIME1s = vertcat(vrcTIME1s,zeros(PSD1slotN,3));
    end
    if size(vrcTIME2s,1) < PSD2slotN
        vrcTIME2s = vertcat(vrcTIME2s,zeros(PSD2slotN,3));
    end 
 
 
 
 
    % (SLPSD1LOCkOFF: DOT)
    % SLPSD1LOCkOFF = SLPSD1((SLPSD1(:,1)>(2)),3);
    % SLPSD2LOCkOFF = SLPSD2((SLPSD2(:,1)>(2)),3);
    
    %}
    
% NOW GET ME THE DWELL TIME VALUES AND THE MX ID
% LOCATIONS OF THOSE GLUR2 RECEPTORS IN PSD1 OR PSD2
 
    [r,c,v] = find(GluR2_TdwellPSD);    % WHICH PSD, WHICH AMPAR, HOW LONG
    rcv = [v r c];                      % ONLY AMPARS WITH DWELL TIMES (TIME,PSD,DOT)
    [r1] = find(rcv(:,2)<2);            % PSD1 DATA
    [r2] = find(rcv(:,2)>1);            % PSD2 DATA
    p1r = rcv(r1,:);                    % PSD1 DWELL DATA ONLY
    p2r = rcv(r2,:);                    % PSD2 DWELL DATA ONLY
    p1rs = sortrows(p1r,-1);            % PSD1: SORT DWELL DATA LONG TO SHORT
    p2rs = sortrows(p2r,-1);            % PSD2: SORT DWELL DATA LONG TO SHORT
    
    if size(p1rs,1) < PSD1slotN     % IF ALL SLOTS NOT FILLED MAKE BLANKS
        p1rs = vertcat(p1rs,zeros(PSD1slotN,3));
    end
    if size(p2rs,1) < PSD2slotN
        p2rs = vertcat(p2rs,zeros(PSD2slotN,3));
    end 
    
    % (SLPSD: TIME,PSD,DOT)
    SLPSD1 = p1rs(1:PSD1slotN,:);   %PSD1: GET THE N(slots) LONGEST DWELLING AMPARS
    SLPSD2 = p2rs(1:PSD2slotN,:);   %PSD2: GET THE N(slots) LONGEST DWELLING AMPARS
    
    
% LETS SEE SOME IDs FELLAS
    
    % (SLPSD1LOCkON: DOT)
    SLPSD1LOCkON = SLPSD1((SLPSD1(:,1)>0),3);
    SLPSD2LOCkON = SLPSD2((SLPSD2(:,1)>0),3);
 
    
    G2SAPLOCD1 = [G2SAPLOC1]';
    G2SAPLOCD2 = [G2SAPLOC2]';
    GluR2DSP1 = [GluR2SP1]';
    GluR2DSP2 = [GluR2SP2]';
    
GluR2DOTSP1 = [GluR2DSP1(SLPSD1LOCkON,:) SLPSD1LOCkON G2SAPLOCD1(SLPSD1LOCkON,:)];
GluR2DOTSP2 = [GluR2DSP2(SLPSD2LOCkON,:) SLPSD2LOCkON G2SAPLOCD2(SLPSD2LOCkON,:)];
 
% THESE ARE THE SLOTTED DOTS
    [uniV1 uniR2 uniC1] = unique(GluR2DOTSP1(:,3),'rows','first');
    GluR2DOTSP1 = GluR2DOTSP1(uniR2,:);
    GluR2DOTSP1 = GluR2DOTSP1(:,1:2);
    
    [uniV2 uniR2 uniC2] = unique(GluR2DOTSP2(:,3),'rows','first');
    GluR2DOTSP2 = GluR2DOTSP2(uniR2,:);
    GluR2DOTSP2 = GluR2DOTSP2(:,1:2);
    
 
    
% HERE IS WHERE WE ROLL THE DICE FOR EACH GLUR2 PARTICLE IN A SLOT REGION
% TO DETERMINE WHETHER IT BINDS TO THE SLOT
 
    GR2PROBPR2 = randn(numel(GluR2DOTSP1(:,1)),1)-2;
    GR2PROBP1 = [GluR2DOTSP1(:,1) GluR2DOTSP1(:,2) (GluR2DOTSP1(:,1)+GR2PROBPR2)];
    GR2PROBPRF1 = GR2PROBP1((GR2PROBP1(:,3)>0),2);
    
    GR2PROBPR2 = randn(numel(GluR2DOTSP2(:,1)),1)-2;
    GR2PROBP2 = [GluR2DOTSP2(:,1) GluR2DOTSP2(:,2) (GluR2DOTSP2(:,1)+GR2PROBPR2)];
    GR2PROBPRF2 = GR2PROBP2((GR2PROBP2(:,3)>0),2);
    
 
% NOW UPDATE THE DIFFUSION RATES FOR THOSE GLUR2 RECEPTORS
% THAT MADE IT INTO A SLOT
 
    xysG2(:,GR2PROBPRF1) = xysG2(:,GR2PROBPRF1)*(.001);
    xysG2(:,GR2PROBPRF2) = xysG2(:,GR2PROBPRF2)*(.001);
    
 
    
    % HOW MANY SLOTS ARE NOW FILLED?
    G2PSD1FSLOTS = numel(GR2PROBPRF1);
    G2PSD2FSLOTS = numel(GR2PROBPRF2);
 
   

% if stepN == 100; keyboard; end
%-------------------%   
end %useSlotsPSDGR2 
%===================%
% if stepN == 100; keyboard; end
 
 
 
%===================%
if useSlotsPERIGR2
%-------------------% 
    
        [rpx,cpx,vpx] = find(GluR2_TdwellPERI);
        rcv = [vpx rpx cpx]; % VALUE ROW COLUMN
        [r1] = find(rcv(:,2)<2);
        [r2] = find(rcv(:,2)>1);
        p1r = rcv(r1,:);
        p2r = rcv(r2,:);
        p1rs = sortrows(p1r,-1);
        p2rs = sortrows(p2r,-1);
        if size(p1rs,1) < PERI1slotN
            p1rs = vertcat(p1rs,zeros(PERI1slotN,3));
        end
        if size(p2rs,1) < PERI2slotN
            p2rs = vertcat(p2rs,zeros(PERI2slotN,3));
        end 
            
         
        
        SLPERI1 = p1rs(1:PERI1slotN,:); %SLOT LOCATION PERI1
        SLPERI2 = p2rs(1:PERI2slotN,:); %SLOT LOCATION PERI2
        
        SLPERI1LOCkon = SLPERI1((SLPERI1(:,1)>KonSpi1PERIGR2),3);
        SLPERI2LOCkon = SLPERI2((SLPERI2(:,1)>KonSpi2PERIGR2),3);
        
    SLPERI1LOCkoff = SLPERI1((SLPERI1(:,1)>(KoffSpi1PERIGR2+KonSpi1PERIGR2)),3);
    SLPERI2LOCkoff = SLPERI2((SLPERI2(:,1)>(KoffSpi2PERIGR2+KonSpi2PERIGR2)),3);
    
    xysG2(:,SLPERI1LOCkon) = xysG2(:,SLPERI1LOCkon)*(.001);
    xysG2(:,SLPERI2LOCkon) = xysG2(:,SLPERI2LOCkon)*(.001);
    xysG2(:,SLPERI1LOCkoff) = xysG2(:,SLPERI1LOCkoff)/(.001);
    xysG2(:,SLPERI2LOCkoff) = xysG2(:,SLPERI2LOCkoff)/(.001);
 
    G2PERI1FSLOTS = numel(SLPERI1LOCkon)-numel(SLPERI1LOCkoff);
    G2PERI2FSLOTS = numel(SLPERI2LOCkon)-numel(SLPERI2LOCkoff);
%-------------------%   
end %useSlotsPSDGR2
%===================%
 
 
G2FSLOTS = [G2PSD1FSLOTS G2PSD2FSLOTS G2PERI1FSLOTS G2PERI2FSLOTS];
%================%
end
%===================================%
%}
%===================================%


%===================================%
%			SAPS & SLOTS
%===================================%
% PRE-LOOP FUNCTIONS
%---
% SAPSLOTSETUP
% SAPMAP
% SAPMASKREPORT
% SAPPREALLOCATE
%--------------------
%{


%===================================%
% SAPSLOTSETUP
%===================================%
% Create initial SAP S-Clulster Matrix
% Take PSD Size (4) and multiply by 2 (8)
% Create a (8x8) Matrix of 1s 
% Then pad this Mx on all sides by Os
%-------------------%
% Establish the initial scalars, vectors, and matrices
% for the SAP cluster turnnover functions
%-------------------%
function [Nsteps,S1,S2,S1sum,S2sum,...
SC1beta,SC1mu,SC1r,SC1ro,SC1tau,...
SC2beta,SC2mu,SC2r,SC2ro,SC2tau...
] = SAPSLOTSETUP(Nsteps,...
sap,um,dot,SAPPADPSD1,SAPPADPSD2)

%{
SCdeltaT = 0.01;	% Shouval [.01]		Brad [.01]
SCbeta = 60;		% Shouval [60]		Brad [50]
SCtau = 1.0;		% Shouval [1.0]		Brad [1.8]
SCmu = 1/SCtau;		% Shouval [1/tau]	Brad [1/tau]
SCL = 1.5;			% Shouval [1.5]		Brad [1.2]
SCL2 = 1.1;			% Shouval [0.9]		Brad [1.1]
SCr = 10;			% Shouval [10]		Brad [15]
SCro = 0.95;		% Shouval [.95]		Brad [.90]
SCszi = 7;			% Shouval [8]		Brad [7]
SCszb = 17;			% Shouval [17]		Brad [17]

SCdeltaT = 0.01;	% Shouval [.01]		Brad [.01]
SCbeta = 50;		% Shouval [60]		Brad [50]
SCtau = 1.8;		% Shouval [1.0]		Brad [1.8]
SCmu = 1/SCtau;		% Shouval [1/tau]	Brad [1/tau]
SCL = 1.2;			% Shouval [1.5]		Brad [1.2]
SCL2 = 1.1;			% Shouval [0.9]		Brad [1.1]
SCr = 15;			% Shouval [10]		Brad [15]
SCro = 0.90;		% Shouval [.95]		Brad [.90]
SCszi = 7;			% Shouval [8]		Brad [7]
SCszb = 17;			% Shouval [17]		Brad [17]

%}

% TIME & SCALE
SAPSTEP = sap(11);
SAPSTEP1 = sap(11);
SAPSTEP2 = sap(12);

Steps = dot(3);			% Steps
TimeStep = dot(4)/1000; % TimeStep
Scale = dot(5);			% Scale
time=Nsteps;


% DENDRITIC FIELD
um1 = um(1); % denWidthX
um2 = um(2); % denHeightY
um3 = um(3); % PSD1um
um4 = um(4); % PSD2um
um5 = um(5); % PERI1um
um6 = um(6); % PERI2um
fsizeX = round(um1/Scale);
fsizeY = round(um2/Scale);
PSD1size = round(um3/Scale);
PSD2size = round(um4/Scale);
periPSD1size = round(um5/Scale);
periPSD2size = round(um6/Scale);

PSD1SAPF = PSD1size*2;
PSD2SAPF = PSD2size*2;
PSA1SAPF = periPSD1size*2;
PSA2SAPF = periPSD2size*2;


if ~SAPPADPSD1
SAPMXPSD1 = ones(PSD1SAPF);
S1=SAPMXPSD1;
end

if ~SAPPADPSD2
SAPMXPSD2 = ones(PSD2SAPF);
S2=SAPMXPSD2;
end

PADSAP1 = SAPPADPSD1*2;
PADSAP2 = SAPPADPSD2*2;

if SAPPADPSD1
SAPMXPSD1=padarray(ones(PSD1SAPF),[PADSAP1 PADSAP1], 0);
S1=SAPMXPSD1;
end

if SAPPADPSD2
SAPMXPSD2=padarray(ones(PSD2SAPF),[PADSAP2 PADSAP2], 0);
S2=SAPMXPSD2;
end


%===========================================%
sap1 = sap(1); % SAPdotsPSD1
sap2 = sap(2); % SAPdotsPSD2
sap3 = sap(3); % SAPbetaPSD1
sap4 = sap(4); % SAPtauPSD1
sap5 = sap(5); % SAPL1PSD1
sap6 = sap(6); % SAPbetaPSD2
sap7 = sap(7); % SAPtauPSD2
sap8 = sap(8); % SAPL1PSD2
sap9 = sap(9); % SAPmuPSD1
sap10 = sap(10); % SAPmuPSD2
sap11 = sap(11); % SAPdTPSD1
sap12 = sap(12); % SAPdTPSD2
sap13 = sap(13); % SAPrhoPSD1
sap14 = sap(14); % SAPrhoPSD2
sap15 = sap(15); % SAPrPSD1
sap16 = sap(16); % SAPrPSD2
sap17 = sap(17); % doDynamicLeP1
sap18 = sap(18); % doDynamicLeP2


h_mask=[0 1 0; 1 0 1; 0 1 0];


SC1szi = sap(1);		SC2szi = sap(2);
SC1beta = sap(3);		SC2beta = sap(6);	
SC1tau = sap(4);		SC2tau = sap(7);		
SC1mu = sap(9);			SC2mu = sap(10);	
SC1L = sap(5);			SC2L = sap(8);
SC1r = sap(15);			SC2r = sap(16);			
SC1ro = sap(13);		SC2ro = sap(14);

SC1L2 = 1.1;			SC2L2 = 1.1;	
SC1szb = 17;			SC2szb = 17;

SC1deltaT = TimeStep*SAPSTEP;	SC2deltaT = TimeStep*SAPSTEP;

%===========================================%
% CLUSTER SUMS FOR PSD POTENTIATION LEVEL
%-------------------------------------------%
S1sum = sum(S1(:));
S2sum = sum(S2(:));
%===========================================%


end
%===================================%


%===================================%
% SAPMAP
%===================================%
% Create two Cells that contain the X and Y coordinates respectively
% for polygon boxes distributed throughout the PSD
% Return Cell Mx to function caller
%-------------------%
function [PSD1Sxvecs PSD1Syvecs PSD2Sxvecs PSD2Syvecs] = SAPMAP(...
	um,dot,...
	SPYN1xv,SPYN1yv,SPYN2xv,SPYN2yv,...
	PSD1xv,PSD1yv,PSD2xv,PSD2yv,...
	PERI1xv,PERI1yv,PERI2xv,PERI2yv,SAPPADPSD1,SAPPADPSD2,SS)


%{
Scale = dot(5); % Scale
% DENDRITIC FIELD
um1 = um(1); % denWidthX
um2 = um(2); % denHeightY
um3 = um(3); % PSD1um
um4 = um(4); % PSD2um
um5 = um(5); % PERI1um
um6 = um(6); % PERI2um
PSD1sz = round(um3/Scale);
PSD2sz = round(um4/Scale);
periPSD1sz = round(um5/Scale);
periPSD2sz = round(um6/Scale);

PSD1XL = PSD1xv(1);
PSD1YL = PSD1yv(1);

PSD2XL = PSD2xv(1);
PSD2YL = PSD2yv(1);

% Get top left corner coordinates for SYN
PSD1XL = PSD1XL-SAPPADPSD1;
PSD1YL = PSD1YL+SAPPADPSD1;
PSD2XL = PSD2XL-SAPPADPSD2;
PSD2YL = PSD2YL+SAPPADPSD2;

% PSD1size = PSD1sz + (SAPPADPSD1*2);
% PSD2size = PSD2sz + (SAPPADPSD2*2);

% PSD1size = PSD1sz + (periPSD1sz*2);
% PSD2size = PSD2sz + (periPSD2sz*2);
%}

PSD1XL = PSD1xv(1);
PSD1YL = PSD1yv(1);

PSD2XL = PSD2xv(1);
PSD2YL = PSD2yv(1);

% Get top left corner coordinates for SYN
PSD1XL = PSD1XL-SAPPADPSD1;
PSD1YL = PSD1YL+SAPPADPSD1;
PSD2XL = PSD2XL-SAPPADPSD2;
PSD2YL = PSD2YL+SAPPADPSD2;

PSD1size = SS.S1um;
PSD2size = SS.S2um;


% Create two Cells that contain the X and Y coordinates respectively
% for the small polygons that are uniformly distributed throughout the PSD
% and return those Cells outside of the main loop
for m = 0:(PSD1size - 1)
for n = 0:(PSD1size - 1)
p = 1+n+(m*PSD1size);
PSD1Sxv = [(PSD1XL+m) (PSD1XL+(1+m)) (PSD1XL+(1+m)) (PSD1XL+m) (PSD1XL+m)]';
PSD1Syv = [(PSD1YL-n) (PSD1YL-n) (PSD1YL-(1+n)) (PSD1YL-(1+n)) (PSD1YL-n)]';
PSD1Sxvecs{p} = PSD1Sxv;
PSD1Syvecs{p} = PSD1Syv;
end
end


for m = 0:(PSD2size - 1)
for n = 0:(PSD2size - 1)
p = 1+n+(m*PSD2size);
PSD2Sxv = [(PSD2XL+m) (PSD2XL+(1+m)) (PSD2XL+(1+m)) (PSD2XL+m) (PSD2XL+m)]';
PSD2Syv = [(PSD2YL-n) (PSD2YL-n) (PSD2YL-(1+n)) (PSD2YL-(1+n)) (PSD2YL-n)]';
PSD2Sxvecs{p} = PSD2Sxv;
PSD2Syvecs{p} = PSD2Syv;
end
end


end
%===================================%


%===================================%
% SAPMASKREPORT
%===================================%
% Function takes SAP S-Cluster Matrix (S1)
% uses mask to count SAPs in each 1x1 PolyBox region
% returns condensed Field-Mx 1/4th of full S-Cluster field, 
% holds values (0-4) of number of SAPs in each PolyBox
%-------------------%
function [SAPFmx1 SAPFmx2 SAPmx1 SAPmx2] = SAPMASKREPORT(Nsteps,S1, S2,um,dot)

Scale = dot(5); % Scale
um3 = um(3); % PSD1um
um4 = um(4); % PSD2um
PSD1size = round(um3/Scale);
PSD2size = round(um4/Scale);
sap_mask=[1 1; 1 1];

[SS1PSD1 SS2PSD1] = size(S1);
[SS1PSD2 SS2PSD2] = size(S2);

SSPSD1 = round(SS1PSD1/2);
SSPSD2 = round(SS1PSD2/2);


% This constructs SAPFmx like: 
% SAPFmx = [(1,1) (1,3) (1,5)... (3,1) (3,3)....]

	SAPmx1_occ=(S1>0);  
	SAPmx1=convn(SAPmx1_occ,sap_mask,'same');
	SAPFmx1=zeros(SSPSD1);
	for n = 0:(SSPSD1-1)
	for m = 0:(SSPSD1-1)
	SAPFmx1(m+1,n+1) = SAPmx1(1+(m*2),1+(n*2));
	end
	end

	SAPmx2_occ=(S2>0);  
	SAPmx2=convn(SAPmx2_occ,sap_mask,'same');
	SAPFmx2=zeros(SSPSD2);
	for n = 0:(SSPSD2-1)
	for m = 0:(SSPSD2-1)
	SAPFmx2(m+1,n+1) = SAPmx2(1+(m*2),1+(n*2));
	end
	end

end
%===================================%


%===================================%
%		SAPPREALLOCATE
%===================================%
function [sap_mask SSPSD1 SSPSD2 GluR1SdP1 GluR1SdP2...
    GluR1SLOCMXa GluR1SLOCMXb g1polyN1 g2polyN1...
    GluR2SdP1 GluR2SdP2 GluR2SLOCMXa GluR2SLOCMXb] = SAPPREALLOCATE(...
    S1,S2,um,dot,runSAPPSD1,runSAPPSD2,PSD1Sxvecs,...
    xyG1,xyG2)


%-------------------
% SAPMASKREPORT
%-------------------
%Scale = dot(5); % Scale
%um3 = um(3); % PSD1um
%um4 = um(4); % PSD2um
%PSD1size = round(um3/Scale);
%PSD2size = round(um4/Scale);

sap_mask=[1 1; 1 1];
% [SS1PSD1 SS2PSD1] = size(S1);
% [SS1PSD2 SS2PSD2] = size(S2);
[SS1PSD1, ~] = size(S1);
[SS1PSD2, ~] = size(S2);
SSPSD1 = round(SS1PSD1/2); %!! SSPSD1 MAY NOT NEED ROUNDING
SSPSD2 = round(SS1PSD2/2); %!! SSPSD2 MAY NOT NEED ROUNDING



%-------------------
% G1SAPPOLYGON
%-------------------
GluR1SdP1 = zeros(size(xyG1));
GluR1SdP2 = GluR1SdP1;
 
GluR1SLOCMXa = zeros(1,size(GluR1SdP1,2));
GluR1SLOCMXb = GluR1SLOCMXa;

[m,g1polyN1,n] = size(PSD1Sxvecs);
    % currentcell = zeros(1,size(xyG1,2));


%-------------------
% G2SAPPOLYGON
%-------------------
GluR2SdP1 = zeros(size(xyG2));
GluR2SdP2 = GluR2SdP1;
 
GluR2SLOCMXa = zeros(1,size(GluR2SdP1,2));
GluR2SLOCMXb = GluR2SLOCMXa;
 
[m,g2polyN1,n] = size(PSD1Sxvecs);
    % currentcell = zeros(1,size(xyG2,2));



end

%}
%--------------------
% LOOPED FUNCTIONS
%---
% SUPERSLOT()
% S1_MainClusterFun()
% S2_MainClusterFun()
% MaskFun()
% inboxfun()
%--------------------


%----------------------------------%
% S1_MainClusterFun Old
%{
%----------------------------------%
% S1_MainClusterFun
%----------------------------------%
function [S1 G1SP1, SMx] = S1_MainClusterFun(S1Ph1, S1Ph2, S1Ph3, S1Ph4, S1Ph5,...
			S1, S1sz, sap, stepN, G1SP1,GT,GTab,...
			ActinTips,ActUpdate,TipCellA,AMask,hkMask,ACTINp,...
			G1xy,G2xy,SMx,MXSZ,G1V1,G1V2)

%=========================================================%
%	TEST MASK AGAINST CLUSTER - GET CONVOLUTION MX
%=========================================================%
% Pen = (S) * mu * dT
% Pex = 1/(1+exp(-B*hk))
% Pkex = (1-S) ( rho * r * dT * Pex)
%{
% B:	slope of Pex function						[hi = shrink]
% mu:	internalization rate						[hi = shrink]
% Le:	lattice energy repulsion constant			[hi = shrink]
% rho:	probability an endo pool sap is available	[hi = grow]
% r:	transition rate from endo to empty site		[hi = grow]
% hk:	mask energy field (matrix convelution)
% dT:	time-step interval
% S:	Surface SAP matrix

% S=padarray(ones(7),[4 4], 0);
% hkMask=[0 1 0; 1 0 1; 0 1 0];
% B = 80;			% hi = shrink
% Le = 1.21;		% hi = shrink
% rho = 0.9;		% hi = grow
% r = 10;			% hi = grow
% mu = .15;		% hi = shrink

% Pmx = rand(size(S));
% Soc = (S>0);
% Snoc = (S<1);
% 
% hk = convn(Soc,hkMask,'same');
% Lhk = hk-Le;
% 
% Pen = (Soc) .* mu * dT;
% Pex = 1 ./ (1+exp((-1*B)*Lhk));
% Pkex = (1-Soc) .* ( rho * r * dT * Pex );
% 
% Sen = (Pen>Pmx);
% Sex = (Pkex>Pmx);
% 
% S = (Soc-Sen) + Sex;
%}
%{
hkMask=[0 1 0; 1 0 1; 0 1 0];
S = S1;
dT = sap(11);

Le = S1Le;		% hi = shrink
B = sap(3);		% hi = shrink
rho = sap(13);	% hi = grow
r = sap(15);	% hi = grow
mu = sap(9);	% hi = shrink
%-------------------------------%
Pmx = rand(size(S));
Soc = (S>0);
Snoc = (S<1);
Sno = ~Soc;
%---
hk = convn(Soc,hkMask,'same');
Lhk = hk-Le;
%---
Pex = 1 ./ (1+exp((-1*B)*Lhk));
Pkex = (1-Soc) .* ( rho * r * dT * Pex );
Pkex0 = Pkex;
Pen = (Soc) .* mu * dT;
%---
S0 = S;
Pex0 = Pex;
Pkex0 = Pkex;
Pen0 = Pen;

%=========================================================%
% hk:	cluster force (2D field mask:Mx convelution)
% Le:	repulsion constant (lattice repulsion force)	[hi = shrink]
% Lhk:	(cluster force) - (repulsion constant)
% B:	slope of Pex(Lhk) function						[hi = shrink]
% rho:	probability an endo pool sap is available		[hi = grow]
% r:	transition rate from endo to empty site			[hi = grow]
% dT:	time-step interval
% mu:	internalization rate							[hi = shrink]

% S:	Surface SAP matrix (Soc:occupied | Snoc:not occupied)
%=========================================================%
%-------------------------------%
G1oc = (S>2) .* 0.0;
SG1oc = Pkex;
%-------------------------------%

%=========================================================%
doLTPS1 = sap(19);
if doLTPS1
%=========================================================%

%=========================================================%
% GluA RT SAP recruitment tails coefficent (basal and LTP)
%=========================================================%
G1RT = G1RTBASE;
if stepN >= LTP1onG1 && stepN <= LTP1offG1; G1RT=G1RTLTP; end
if stepN >= LTP2onG1 && stepN <= LTP2offG1; G1RT=G1RTLTP; end

%=========================================================%
%		Get Matrix Locations of SAPs near AMPARs
%=========================================================%
G1P1poly0 = G1P1SAPM(:,2);
[~, ~, G1P1poly1] = find(G1P1poly0);
G1P1poly2 = G1P1poly1;
%=============================================%
for nx = 1:(numel(G1P1poly1))
%=============================================%
    if G1P1poly1(nx)<=8 && G1P1poly1(nx)>0
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*0)-1;
    end
    if G1P1poly1(nx)<=16 && G1P1poly1(nx)>8
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*1)-1;
    end
    if G1P1poly1(nx)<=24 && G1P1poly1(nx)>16
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*2)-1;
    end
    if G1P1poly1(nx)<=32 && G1P1poly1(nx)>24
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*3)-1;
    end
    if G1P1poly1(nx)<=40 && G1P1poly1(nx)>32
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*4)-1;
    end
    if G1P1poly1(nx)<=48 && G1P1poly1(nx)>40
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*5)-1;
    end
    if G1P1poly1(nx)<=56 && G1P1poly1(nx)>48
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*6)-1;
    end
    if G1P1poly1(nx)<=64 && G1P1poly1(nx)>56
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*7)-1;
	end
%=============================================%
end % for nx
%=============================================%


%=============================================%
if G1RTBASE > 0
%=============================================%

%-------------------------------%
G1oc = (S>2) .* 0.0;
SG1oc = Pkex;
G1P1n = numel(G1P1poly2);
%-------------------------------%

%=============================================%
for polyid = 1:G1P1n
%=============================================%
G1P1IDo = G1P1poly2(polyid);

G1P1IDb = G1P1IDo+1;
G1P1IDr = G1P1IDo+16;
G1P1IDd = G1P1IDo+17;

if G1P1IDo>=16
	G1P1IDb = G1P1IDo+1;
	G1P1IDr = G1P1IDo+16+(round(rand)*-32); 
	G1P1IDd = G1P1IDr+1;
	% 50:50 WHETHER EAST OR WEST SAP
	% ALWAYS EAST SAP IF COMMENTED
	%G1P1IDd = G1P1IDo+17+(round(rand)*-34);
end
%---
G1oc(G1P1IDo)=1;
G1oc(G1P1IDb)=1;
G1oc(G1P1IDr)=1;
G1oc(G1P1IDd)=1;
%---
SG1oc(G1P1IDo)=1;
SG1oc(G1P1IDb)=1;
SG1oc(G1P1IDr)=1;
SG1oc(G1P1IDd)=1;
%=============================================%
end % for polyid = 1:G1P1n
%=============================================%
Gex = Sno .* SG1oc .* G1RT + Pkex;
Gp = 1 ./ (1+exp((-1*1)*hk)) + .001;
Gpex = Gex .* Gp;
%---
Pkex = Gpex;
%=============================================%
end % if G1RTBASE > 0
%=============================================%
end % if doLTPS1
%=========================================================%


%=========================================================%
%		FINAL CLUSTER PARAMETERS
%---------------------------------------------------------%
Sen = (Pen>Pmx);
Sex = (Pkex>Pmx);
S = (Soc-Sen) + Sex;
S1=S;

% S=padarray(ones(7),[4 4], 0);
%=========================================================%
%}
%---------------------------------------------%


S=S1;
dT = sap(11);
% hkMask=[0 1 0; 1 0 1; 0 1 0];


Lon = sap(21);	% On Energy (lower = more on events)
Bon = sap(22);	% On Neighbor-Independant Rate (new growth) (lower = more on)
Ron = sap(23);	% On Neighbor-Dependant Rate (cluster fill-in) (higher = more on)

Loff = sap(24);	% Off Energy (higher = more off events)
Boff = sap(25);	% Off Neighbor-Inependant Rate (uniform off)  (lower = more off)
Roff = sap(26);	% Off Neighbor-Dependant Rate (edge off) (higher = more off)

% %------------
% if mod(stepN,ActUpdate) == 0
% TipCellA = TipCellA+1;
% ACTINp = ActinTips{TipCellA};
% ACTINp = convn(ACTINp,AMask,'same');
% end
% %------------

%------------
Pmx = rand(size(S));
Soc = (S>0);
Sno = ~Soc;
%---
hk = convn(ACTINp,hkMask,'same');
%------------
Pon = 1 ./ (1+exp((hk-Lon).*(-Bon)));
Pkon = Sno .* ( Ron * dT * Pon );
%---
Son = (Pkon>Pmx);
%------------
Poff = 1 ./ (1+exp(((-hk)+Loff).*(-Boff)));
Pkoff = Soc .* ( Roff * dT * Poff );
%---
Soff = (Pkoff>Pmx);
%------------


%-------------------------------%
G1oc = zeros(S1sz);
%-------------------------------%

% %=========================================================%
% doLTPS1 = sap(19);
% if doLTPS1
% if GTon > 0
% 	
% G1oc(G1P2poly2)=1;
% %---------------------------------------------%
% GRhk = convn(G1oc,GhkMask,'same');
% GRk=(GRhk.*GTon);
% GSk=(GRhk.*GToff);
% %-------------------------------%
% Gon = Pkon+(GRk.*(Pkon + LTPv));
% Goff = Pkoff+(GSk.*Pkoff);
% %-------------------------------%
% Son = (Gon>Pmx);
% Soff = (Goff>Pmx);
% %-------------------------------%
% end % if G1RTBASE > 0
% %-------------------------------%
% end % if doLTPS1
% %=========================================================%


%---------------------------------------------%
%			FINAL CLUSTER PARAMETERS
%---------------------------------------------%
S = (Soc-Soff) + Son;
S1=S;
%---------------------------------------------%



%=============================================%
%				PLOT CLUSTER
%---------------------------------------------%
if mod(stepN,100) == 0
set(S1Ph1,'CData',S1);
drawnow
end
%---------------------------------------------%
%			ADDITIONAL CLUSTER PLOTS
%---------------------------------------------%
doS1PLOTS=1;
if doS1PLOTS
if mod(stepN,200) == 0
%---------------------------------------------%
figure(1);
%---

set(S1Ph2,'CData',hk);
drawnow

set(S1Ph3,'CData',Pkon);
drawnow

set(S1Ph4,'CData',Sno);
drawnow

% set(S1Ph5,'CData',GRhk);
% drawnow

%{
%=========================================================%
% if stepN == 2000; keyboard; end
% cbr = colorbar('location','South');
% set(cbr(1), 'XColor', 'r');
% text(1,1.5,'S0','FontSize',22,'Color',[.7 1 .7]);
% caxis(caxis);
%=========================================================%

%---
% LEFT POS-2 (hk)
%---
% subplot('Position',[.02 .54 .3 .45]),...
subplot(5,5,6),
imagesc(hk)
title('hk');

%---
% LEFT POS-3 (Pkon)
%---
% subplot('Position',[.35 .54 .3 .45]),...
subplot(5,5,11),
imagesc(Pkon);
title('Pkon');

%---
% LEFT POS-4 (Sno)
%---
% subplot('Position',[.68 .54 .3 .45]),...
subplot(5,5,16),
imagesc(Sno);
title('Sno');


%---
% LEFT POS-5 (G1oc)
%---
% subplot('Position',[.02 .05 .3 .45]),...
subplot(5,5,21),
imagesc(G1oc);
title('TETHERED G1');

%---
% LEFT POS-6 (Sex)
%---
% subplot('Position',[.35 .05 .3 .45]),...
% subplot(5,5,16),
% imagesc(Sex);
% title('Sex');

%---
% LEFT POS-4 (Sen)
%---
% subplot('Position',[.68 .54 .3 .45]),...
% subplot(5,5,16),
% imagesc(Sen);
% title('Sen');

%}
%---------------------------------------------%
end % if mod(stepN,100) == 0
end % if doS1PLOTS
%=========================================================%
end


%}
%----------------------------------%
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

