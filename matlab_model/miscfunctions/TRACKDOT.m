function [varargout] = DOTBOX2(dot,dr,um,sap,hr,ko,doUse,doRun,doKo,box,slt,stky)
format compact;
format short;
close all;
SNSZ = get(0,'ScreenSize');

doIz = 0;

doprofile = 0;
if doprofile
profile on;
end
%=========================================================%
%					INPUTS FROM GUI
%---------------------------------------------------------%

%{.
if 1-exist('dotbox','var')
dotbox=[30 200 5 0 0.1 1000 50 50 3 2.4 0.1 3 2.4 0.1 0.2 0.2 1 1 2 1000 2 1000 7 7 1 1 1 1];
end
%}


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
 
um1 = um(1); % denWidthX
um2 = um(2); % denHeightY
um3 = um(3); % PSD1um
um4 = um(4); % PSD2um
um5 = um(5); % PERI1um
um6 = um(6); % PERI2um
 
sap1 = sap(1); % SAPdotsPSD1
sap2 = sap(2); % SAPdotsPSD2
sap3 = sap(3); % SAPbetaPSD1
sap4 = sap(4); % SAPtauPSD1
sap5 = sap(5); % SAPL1PSD1
sap6 = sap(6); % SAPbetaPSD2
sap7 = sap(7); % SAPtauPSD2
sap8 = sap(8); % SAPL1PSD2

slt1 = slt(1); % G1PSDslotN
slt2 = slt(2); % G1PERIslotN
slt3 = slt(3); % G2PSDslotN
slt4 = slt(4); % G2PERIslotN

hr1 = hr(1); % HomeostaticLo
hr2 = hr(2); % HomeostaticHi
 
ko1 = ko(1); % KonSpi1PSDGR2
ko2 = ko(2); % KoffSpi1PSDGR2
ko9 = ko(9); % KonSpi1PSDGR1
ko10 = ko(10); % KoffSpi1PSDGR1

doUse1 = doUse(1); % useGluR1
doUse2 = doUse(2); % useGluR2
doUse3 = doUse(3); % runSAPPSD1
doUse4 = doUse(4); % runSAPPSD2
 
doRun1 = doRun(1); % run2Dplot
doRun2 = doRun(2); % run3Dplot
doRun3 = doRun(3); % runMSDtest
doRun4 = doRun(4); % runTraceSingleDot
doRun5 = doRun(5); % runMSDpopup      % POPUP
doRun6 = doRun(6); % runHomeostatic
doRun7 = doRun(7); % runPoissonsBox
 
doKo1 = doKo(1); % useGluR1slots
doKo2 = doKo(2); % useGluR2slots
doKo3 = doKo(3); % useSlotsPSDGR1
doKo4 = doKo(4); % useSlotsPSDGR2
doKo5 = doKo(5); % useSlotsPERIGR1
doKo6 = doKo(6); % useSlotsPERIGR2

box1 = box(1);	% GraphTime
box2 = box(2);	% AllowedTime
box3 = box(3);	% NumLoops
box4 = box(4);	% SteadyState


G1STBASE = stky(1);
G1RTBASE = stky(2);
G1STLTP = stky(3);
G1RTLTP = stky(4);
G1BSMu = stky(5);
G1LSMu = stky(6);

G2STBASE = stky(7);
G2RTBASE = stky(8);
G2STLTP = stky(9);
G2RTLTP = stky(10);
G2BSMu = stky(11);
G2LSMu = stky(12);

SAPPADPSD1 = stky(21);
SAPPADPSD2 = stky(22);


LTP1onG1 = stky(13);
LTP1offG1 = stky(15);
LTP2onG1 = stky(14);
LTP2offG1 = stky(16);

LTP1onG2 = stky(17);
LTP1offG2 = stky(19);
LTP2onG2 = stky(18);
LTP2offG2 = stky(20);


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

%%
%=========================================================%
%                  FUNCTION SWITCHES
%---------------------------------------------------------%
% DO PLOTS
doMainPlot = doRun1;
do3DPLOT = doRun2;

% PARTICLES?
Nsteps = dot3;
doGluR1 = doUse1;
doGluR2 = doUse2;

GluR1Ndots = dot1;		% GluR1 Particles
GluR2Ndots = dot2;		% GluR2 Particles

if ~doGluR1
	GluR1Ndots = 0;
end
if ~doGluR2
	GluR2Ndots = 0;
end

% DO PLOT ONE DOT?
doONEDOTPLOT = doRun4;
if doONEDOTPLOT
	GluR2Ndots = 1;
	clear GluR2xyds 
	clear GluR2xyl
	GluR2xyds = zeros(2,GluR2Ndots);
	GluR2xyl = ones(2,GluR2Ndots);
	xyl2 = GluR2xyl;			% save previous location
	Nsteps = 300;
	doMainPlot=0;
end

% DO TRACK MSD?
trackMSD = doRun3;
if trackMSD
	GluR2Ndots = 100;
	GluR1Ndots = 0;
	Nsteps = 100;
end

doHomeostatic = 0;

%=========================================================%
loops = dot6;
for re = 1:loops
%=========================================================%



%{
if 1-exist('dot1','var')
dot1 = 1;
end
%}
%---------------------------------------------------------%
%					INPUTS FROM GUI
%=========================================================%




%%
%=========================================================%
%				AMPAR DOT PARTICLE MATRIX
%---------------------------------------------------------%
%GLUR1
GluR1xyds = ones(2,GluR1Ndots);
GluR1xyl = ones(2,GluR1Ndots);
G1Z = ones(1,GluR1Ndots);

% GLUR2
GluR2xyds = ones(2,GluR2Ndots);
GluR2xyl = ones(2,GluR2Ndots);
G2Z = ones(1,GluR2Ndots);

TTDots = GluR2Ndots+GluR1Ndots;


%%
%-------------##########################------------------%
%{

If I want to generate a bunch of random steps with diffusion rate of:
D = 0.3
then I would first determine the stdev of the D distribution:
k = sqrt(d*D*t)
and then generate a bunch of random XY steps from that distribution:
Sxy = (k * randn(2,GluR2Ndots))
then if I wanted to increase or decrease that initial step size
I would need to determine the scalar value:
DSc = D/Dn   (DSc = scalar, Dn = D new)
then I can plug this DSc scalar into a function that will provide
a coefficient for adjusting each linear X and Y step:
Ls = 1/sqrt(DSc)
Now, if I multiply each original Sxy step by Ls, i will achieve Dn

Scale = dot5;				% scale of model
t = dot4/1000;				% time step
d = 2;                      % dimensions
D = esDGR1/Scale*t;			% Diffusion Rate (D = L² / 2d*t)
Ds = PSDDGR1/Scale;			% Diffusion Rate PSD
Dr = D/Ds;					% Ratio of D:Ds (1/Ls)^2;
Dn = D/Dr;					% new D after scaling L
k = sqrt(d*D);	            % stdev of D's step size distribution
MSD = 2*d*D;                % mean squared displacement
L = sqrt(2*d*D);            % average diagonal (2D) step size
Lx = L/sqrt(2);             % average linear (1D) step size
Ls = 1/sqrt(Dr);			% scales Lx values for Dn

% PSD DIFFUSION RATES
Dr_PSD1 = dr1/dr4;			% PSD-1 D Scalar Ratio of D:Ds
Dr_PSD2 = dr1/dr4;			% PSD-2 D Scalar Ratio of D:Ds
Dn_PSD1 = dr4/Scale;		% PSD-1 D value after scaling L
Dn_PSD2 = dr4/Scale;		% PSD-2 D value after scaling L
PSD1 = 1/sqrt(Dr_PSD1);		% PSD-1 D Scalar Function, scales Lx values for Dn
PSD2 = 1/sqrt(Dr_PSD2);		% PSD-2 D Scalar Function, scales Lx values for Dn

esDGR1=dr1;
spineDGR1=dr2;
PERIDGR1=dr3;
PSDDGR1=dr4;
esDGR2=dr5;
spineDGR2=dr6;
PERIDGR2=dr7;
PSDDGR2=dr8;
%}
%=========================================================%
%               STARTING PARAMETERS
%---------------------------------------------------------%
% BASE DIFFUSION RATES EQUATIONS
Scale = dot5;				% scale of model
t = dot4/1000;				% time step
d = 2;                      % dimensions
D = esDGR1*t/Scale;			% Diffusion Rate ES (D = L² / 2d*t)
Dp = PSDDGR1*t/Scale;		% Diffusion Rate PSD
Dr = D/Dp;					% Ratio of D:Ds (1/Ls)^2;
Dn = D/Dr;					% new D after scaling L
k = sqrt(d*D);	            % stdev of D's step size distribution
MSD = 2*d*D;                % mean squared displacement
L = sqrt(2*d*D);            % average diagonal (2D) step size
Lx = L/sqrt(2);             % average linear (1D) step size
Ls = 1/sqrt(Dr);			% scales Lx values for Dn

% GLUR DIFFUSION RATES
Ds = esDGR1*t/Scale;		% ExtraSynaptic Model-Scaled Diffusion Rate 
GR1Ds = PSDDGR1*t/Scale;	% GluR1 (PSD) Model-Scaled Diffusion Rate
GR2Ds = PSDDGR2*t/Scale;	% GluR2 (PSD) Model-Scaled Diffusion Rate
Dr_GR1 = Ds/GR1Ds;			% GluR1 (PSD) Ratio Ds:GR1Ds
Dr_GR2 = Ds/GR2Ds;			% GluR2 (PSD) Ratio Ds:GR2Ds
Dn_GR1 = Ds/Dr_GR1;			% GluR1 (PSD) D value after scaling L
Dn_GR2 = Ds/Dr_GR2;			% GluR2 (PSD) D value after scaling L
LsGR1 = 1/sqrt(Dr_GR1);		% GluR1 (PSD) D Lx-Scalar (Ls scales Lx so Dn_GR1 = Ds/Dr_GR1)
LsGR2 = 1/sqrt(Dr_GR2);		% GluR2 (PSD) D Lx-Scalar (Ls scales Lx so Dn_GR2 = Ds/Dr_GR2)

% PSD DIFFUSION RATES (DEFAULT: GLUR1-BASED)
Ds = esDGR1*t/Scale;		% ExtraSynaptic Model-Scaled Diffusion Rate 
PSD1Ds = PSDDGR1*t/Scale;	% (PSD1) Model-Scaled Diffusion Rate
PSD2Ds = PSDDGR1*t/Scale;	% (PSD2) Model-Scaled Diffusion Rate
Dr_PSD1 = Ds/PSD1Ds;		% (PSD1) Ratio Ds:PSD1Ds
Dr_PSD2 = Ds/PSD2Ds;		% (PSD2) Ratio Ds:PSD2Ds
Dn_PSD1 = Ds/Dr_PSD1;		% (PSD1) D value after scaling L
Dn_PSD2 = Ds/Dr_PSD2;		% (PSD2) D value after scaling L
PSD1 = 1/sqrt(Dr_PSD1);		% (PSD1) D Lx-Scalar (Ls scales Lx so Dn_PSD1 = Ds:PSD1Ds)
PSD2 = 1/sqrt(Dr_PSD2);		% (PSD2) D Lx-Scalar (Ls scales Lx so Dn_PSD2 = Ds:PSD2Ds)


% GLUR1 DIFFUSION RATES
DGR1esm = esDGR1*t/Scale;		% GluR1 (ESM) Model-Scaled Diffusion Rate
DGR1spy = spineDGR1*t/Scale;	% GluR1 (ESM) Model-Scaled Diffusion Rate
DGR1psa = PERIDGR1*t/Scale;		% GluR1 (PSA) Model-Scaled Diffusion Rate
DGR1psd = PSDDGR1*t/Scale;		% GluR1 (PSD) Model-Scaled Diffusion Rate
kGR1 = sqrt(d*DGR1esm);			% stdev of D's step size distribution
DrGR1psa = DGR1esm/DGR1psa;		% GluR1 (PSA) Ratio DGR1esm:DGR1psa
DrGR1psd = DGR1esm/DGR1psd;		% GluR1 (PSD) Ratio DGR1esm:DGR1psd
LsGR1psa = 1/sqrt(DrGR1psa);	% GluR1 (PSA) D Lx-Scalar
LsGR1psd = 1/sqrt(DrGR1psd);	% GluR1 (PSA) D Lx-Scalar

% GLUR2 DIFFUSION RATES
DGR2esm = esDGR2*t/Scale;		% GluR2 (ESM) Model-Scaled Diffusion Rate
DGR2spy = spineDGR2*t/Scale;	% GluR2 (ESM) Model-Scaled Diffusion Rate
DGR2psa = PERIDGR2*t/Scale;		% GluR2 (PSA) Model-Scaled Diffusion Rate
DGR2psd = PSDDGR2*t/Scale;		% GluR2 (PSD) Model-Scaled Diffusion Rate
kGR2 = sqrt(d*DGR2esm);			% stdev of D's step size distribution
DrGR2psa = DGR2esm/DGR2psa;		% GluR2 (PSA) Ratio DGR2esm:DGR2psa
DrGR2psd = DGR2esm/DGR2psd;		% GluR2 (PSD) Ratio DGR2esm:DGR2psd
LsGR2psa = 1/sqrt(DrGR2psa);	% GluR2 (PSA) D Lx-Scalar
LsGR2psd = 1/sqrt(DrGR2psd);	% GluR2 (PSA) D Lx-Scalar



%=========================================================%
% Synaptic Dwell-Time Values
%=========================================================%
GluR2_TdwellPSD = zeros(2,GluR2Ndots);
GluR1_TdwellPSD = zeros(2,GluR1Ndots);
GluR2_TdwellPERI = zeros(2,GluR2Ndots);
GluR1_TdwellPERI = zeros(2,GluR1Ndots);
GluR2_TdwellSPYN = zeros(2,GluR2Ndots);
GluR1_TdwellSPYN = zeros(2,GluR1Ndots);


%%
%=========================================================%
%			PARTICLE BOX EXIT SIMULATION (POISSONS BOX)
%---------------------------------------------------------%
if doRun7

pbox = [box1 box2 box3 box4 dot5 dot4 dot1 dot2...
    DGR1spy kGR1 DGR1psd DGR2spy kGR2 DGR2psd LsGR1psd LsGR2psd...
	doKo1 doKo2 ko9 ko10 ko1 ko2...
	sap1 sap2 slt1 slt2 slt3 slt4];

varargin = PoissonsBox(pbox);
GluR1exT30 = [varargin{1,1}]' 
GluR2exT30 = [varargin{1,2}]'
for k=1:2, varargout(k) = {eval(['GluR' int2str(k) 'exT30'])}; end
return
end




%=========================================================%
MSDstring = ['D  Dn_PSD1  Dn_PSD2  d  t  k  L  Lx'];
MSDdata = [D Dn_PSD1 Dn_PSD2 d t k L Lx];
%=========================================================%





%%
%=========================================================%
% PSD & FIELD MATRIX SETUP
%---------------------------------------------------------%
% DENDRITIC FIELD
fsizeX = round(um1/Scale);
fsizeY = round(um2/Scale);
PSD1size = round(um3/Scale);
PSD2size = round(um4/Scale);
periPSD1size = round(um5/Scale);
periPSD2size = round(um6/Scale);

PSD1CNTR = (PSD1size+(periPSD1size*2))/2;
PSD2CNTR = (PSD2size+(periPSD2size*2))/2;


%---##########################---%
%		FIELD MAP FUNCTION
%---##########################---%
[Yrows Xcols YrPr2bot XcPr2rit YrPr2top XcPr2lft...
YrPr1bot XcPr1rit YrPr1top XcPr1lft PSDfield...
PSDSZE Zfield XWIDE YHIGH...
PSD1WH PSD2WH PERI1WH PERI2WH SPYN1WH SPYN2WH...
XYLTpr1 XYRTpr1 XYLBpr1 XYRBpr1 XYLTpr2 XYRTpr2 XYLBpr2 XYRBpr2...
XYLTp1 XYRTp1 XYLBp1 XYRBp1 XYLTp2 XYRTp2 XYLBp2 XYRBp2...
XYBOXpr1 XYBOXpr2 XYBOXp1 XYBOXp2...
SPYN1xv SPYN1yv SPYN2xv SPYN2yv...
PSD1xv PSD1yv PSD2xv PSD2yv...
PERI1xv PERI1yv PERI2xv PERI2yv...
fPSD1 fPSD2...
] = FieldFun(fsizeX, fsizeY,PSD1size, PSD2size, periPSD1size, periPSD2size);
%---##########################---%


%---##########################---%
%		SLOT GAUSSIAN COLORMAP
%---##########################---%
SLOTMAPFUN(G1STBASE,G1RTBASE,G1STLTP,G1RTLTP,G1BSMu,G1LSMu,...
		   G2STBASE,G2RTBASE,G2STLTP,G2RTLTP,G2BSMu,G2LSMu);
%---##########################---%


%%
%-------------##########################------------------%
%              CLUSTER MODEL PARAMETERS
%=========================================================%
% SLOTS PER 5 SAPS
G1PSDSAP5 = slt1;	% G1PSDslotN
G1PERISAP5 = slt2; % G1PERIslotN
G2PSDSAP5 = slt3;	% G2PSDslotN
G2PERISAP5 = slt4; % G2PERIslotN
G1FSLOTS=[0 0 0 0];
G2FSLOTS=[0 0 0 0];
G1SAPLOC1 = zeros(1,GluR1Ndots);
G1SAPLOC2 = zeros(1,GluR1Ndots);
G2SAPLOC1 = zeros(1,GluR2Ndots);
G2SAPLOC2 = zeros(1,GluR2Ndots);

% SLOTS PER 5 SAPS
SAP5 = [G1PSDSAP5 G1PERISAP5 G2PSDSAP5 G2PERISAP5];

% DO SAP CLUSTERS?
runSAPPSD1 = doUse3;	
runSAPPSD2 = doUse4;
if ~runSAPPSD1;sap1=0;end;
if ~runSAPPSD2;sap2=0;end;
S1 = sap1^2; S1sum = sum(S1(:)); S1sumOrig = sum(S1(:)); 
S2 = sap2^2; S2sum = sum(S2(:)); S2sumOrig = sum(S2(:));
if ~runSAPPSD1 || ~runSAPPSD2
	doClusters = 0;
else 
	doClusters = 1;
end





%=============================================%
if doClusters
%---------------------------------------------%
%-------------------%
% SAPSLOTSETUP
%-------------------%
[Nsteps,h_mask,...
S1pop2_num,S1r,SC1L,...
S1,S1ro,SC1L2,S1L1,S1ro2,SC1beta,S1L1_array,S1L2,S1sum,SC1mu,...
S1L2_array,SC1r,S1beta,S1tau,SC1ro,S1delta_t,S1mem_num,S1mu,SC1tau...
S1mu_array,S1num_epochs,...
S2pop2_num,S2r,SC2L,...
S2,S2ro,SC2L2,S2L1,S2ro2,SC2beta,S2L1_array,S2L2,S2sum,SC2mu,...
S2L2_array,SC2r,S2beta,S2tau,SC2ro,S2delta_t,S2mem_num,S2mu,SC2tau...
S2mu_array,S2num_epochs] = SAPSLOTSETUP(Nsteps,...
sap1, sap2, sap3, sap6, sap4, sap7, sap5, sap8,runSAPPSD1,runSAPPSD2,...
um,Zfield,fPSD1,fPSD2,dot,SAPPADPSD1,SAPPADPSD2);
%-------------------%
% SAPMAP
%-------------------%
[PSD1Sxvecs PSD1Syvecs PSD2Sxvecs PSD2Syvecs] = SAPMAP(um,dot,...
SPYN1xv,SPYN1yv,SPYN2xv,SPYN2yv,...
PSD1xv,PSD1yv,PSD2xv,PSD2yv,...
PERI1xv,PERI1yv,PERI2xv,PERI2yv,SAPPADPSD1,SAPPADPSD2);
%-------------------%
% SAPMASKREPORT
%-------------------%
[SAPFmx1 SAPFmx2 SAPmx1 SAPmx2] = SAPMASKREPORT(Nsteps,...
	S1, S2,um,dot,runSAPPSD1,runSAPPSD2);
%-------------------%
% SAPPREALLOCATE
%-------------------%
[sap_mask SSPSD1 SSPSD2 GluR1SdP1 GluR1SdP2...
GluR1SLOCMXa GluR1SLOCMXb g1polyN1 g2polyN1...
GluR2SdP1 GluR2SdP2 GluR2SLOCMXa GluR2SLOCMXb] = SAPPREALLOCATE(...
S1,S2,um,dot,runSAPPSD1,runSAPPSD2,PSD1Sxvecs,...
GluR1xyl,GluR2xyl);

%-------------------%
G1SP1 = zeros(4,100);
G1SP2 = zeros(4,100);
G2SP1 = zeros(4,100);
G2SP2 = zeros(4,100);
G1LTParray = ones(1,numel(S1L2_array));
G2LTParray = ones(1,numel(S1L2_array));
%---------------------------------------------%
end % if runSAPPSD1 || runSAPPSD2
%=============================================%



%{
if runSAPPSD1 || runSAPPSD2
[Nsteps,S1pop2_num,S2beta,S2tau,SC2beta,SCr,initial_up,sz1_bound,...
S,S1r,S2delta_t,SC1L,SC2deltaT,SCro,initial_up1,sz2,...
S1,S1ro,S2mem_num,SC1L2,SC2mu,SCszb,initial_up2,sz2_bound,...
S1L1,S1ro2,S2mu,SC1beta,SC2r,SCszi,sap1,sz_bound,...
S1L1_array,S1s,S2mu_array,SC1deltaT,SC2ro,SCtau,sap2,time,...
S1L2,S1sum,S2num_epochs,SC1mu,SC2szb,center,sap3,...
S1L2_array,S1sumOrig,S2pop2_num,SC1r,SC2szi,center1,sap4,...
S1beta,S1tau,S2r,SC1ro,SC2tau,center2,sap5,...
S1delta_t,S2,S2ro,SC1szb,SCL,h_mask,sap6,...
S1mem_num,S2L1,S2ro2,SC1szi,SCL2,initial_down,sap7,...
S1mu,S2L1_array,S2s,SC1tau,SCbeta,initial_down1,sap8,...
S1mu_array,S2L2,S2sum,SC2L,SCdeltaT,initial_down2,sz,...
S1num_epochs,S2L2_array,S2sumOrig,SC2L2,SCmu,initial_sz,sz1] = SAPSLOTS(Nsteps,...
sap1, sap2, sap3, sap6, sap4, sap7, sap5, sap8,runSAPPSD1,runSAPPSD2,...
um,Zfield,pfPSD1,pfPSD2,dot);
end
	
	
	
	
if runSAPPSD1 || runSAPPSD2
[Nsteps,...
S1pop2_num,S2beta,S2tau,SC2beta,...
S1r,S2delta_t,SC1L,...
S1,S1ro,S2mem_num,SC1L2,SC2mu,...
S1L1,S1ro2,S2mu,SC1beta,SC2r,...
S1L1_array,S2mu_array,SC2ro,...
S1L2,S1sum,S2num_epochs,SC1mu,...
S1L2_array,S2pop2_num,SC1r,...
S1beta,S1tau,S2r,SC1ro,SC2tau,...
S1delta_t,S2,S2ro,h_mask,...
S1mem_num,S2L1,S2ro2,...
S1mu,S2L1_array,SC1tau,...
S1mu_array,S2L2,S2sum,SC2L,...
S1num_epochs,S2L2_array,SC2L2] = SAPSLOTS(Nsteps,...
sap1, sap2, sap3, sap6, sap4, sap7, sap5, sap8,runSAPPSD1,runSAPPSD2,...
um,Zfield,pfPSD1,pfPSD2,dot);
end
	
%}

%=========================================================%


%%
%===========================================%
% TIME AND PARTICLE COUNT VARIABLES
%-------------------------------------------%
ESNT = 0; SPY1N = 0; SPY2N = 0; SPYNT = 0;
SPY1N = 0; SPY2N = 0; SPYNT = SPY1N+SPY2N;
PSD1CaT = 0; PSD2CaT = 0;
PSD1VAL = SPYNT-PSD1CaT; PSD2VAL = SPYNT-PSD2CaT;


DATARATE = 10;
AveOver = 10;
SaveSteps = Nsteps/DATARATE;
dataset = zeros(SaveSteps,10);
SAPdata = zeros(SaveSteps,2);
Ddata = zeros(SaveSteps,2);
AMPARdata = zeros(SaveSteps,2);
GluRdata = zeros(SaveSteps,8);
G1SLOTSDATA = zeros(SaveSteps,4);
G2SLOTSDATA = zeros(SaveSteps,4);

%===========================================%

%===========================================%
% PREP GRAPHICS
%-------------------------------------------%
%{
% fig1 = figure(1);
% fig1box(fig1);
% datafig = figure(88);
% figure('Position',[1 SCREENSIZE(4)/8 SCREENSIZE(3)/8 SCREENSIZE(4)/8])
% set(fig1,'Renderer','painters')
% set(fig1,'DoubleBuffer','off')
% set(fig1,'WindowStyle','modal')
%}

figure(1)
fig1 = figure(1);
set(gcf,'OuterPosition',[300,100,900,800])




if do3DPLOT==0
figure(1);
subplot(5,5,[2 7]), map=imagesc(PSDfield);
maph = get(map,'child');
set(maph,'facea',.3);
set(gca, 'nextplot', 'add')
hold on;
end
%===========================================%


%=============%%%%%%%%%%%%%%%%%%%%%%%%%%%%=============%
%%					START Main Loop
%=============%%%%%%%%%%%%%%%%%%%%%%%%%%%%=============% 
stepN = 1;
for Nt = 1:Nsteps 
nn = Nt;
% if stepN == 201; keyboard; end
%===========================================%
    
  
    %-------------------------------%
    %   STEP DIRECTION & SIZE
    %-------------------------------%
	if doGluR1
	GluR1xyds = STEPxyds(GluR1Ndots, kGR1);
	end
	if doGluR2
    GluR2xyds = STEPxyds(GluR2Ndots, kGR2);
	end

	
	%-------------------------------%
    %     Homeostatic Functions
    %-------------------------------%
	if doHomeostatic
	if mod(stepN, 5) == 0
	[S1beta S2beta S1mu_array S2mu_array S1L1_array S2L1_array]...
	= HOMEOSTATIC(stepN,HShi,HSlow,TTDots,doCalcium,SPYNT,PSD1CaT,...
	S1num_epochs,S2num_epochs,S1mu,S2mu);
	end
	end % end doHomeostatic
	
	
	%-------------------------------%
    %     FRAP Functions
    %-------------------------------%
	%{
	doFRAP = 1;
	if doFRAP
	if stepN == 1
		GluR2xyl = ones(2,GluR2Ndots).*rand(2,GluR2Ndots); 
		GluR2xyl(1,:) = XYLBp1(1)+GluR2xyl(1,:);
		GluR2xyl(2,:) = XYLBp1(2)+GluR2xyl(2,:);
	end
	if stepN < 200
		LsGR2psd = .2;
		GluR1xyds = GluR1xyds.*.0001;
	end
	% if stepN == 5; keyboard; end
	end
	%}
	
	
	%-------------------------------%
    %     SUPERSLOT FUNCTION
    %-------------------------------%
	if doClusters
	[GluR1xyl GluR2xyl GluR1xyds GluR2xyds...
	SAPFmx1 SAPFmx2 SAPmx1 SAPmx2 G1FSLOTS G2FSLOTS...
	GluR1_TdwellSPYN GluR1_TdwellPSD GluR1_TdwellPERI...
	GluR2_TdwellSPYN GluR2_TdwellPSD GluR2_TdwellPERI...
	G1INSPYN1 G1INSPYN2 G1INPSD1 G1INPSD2 G1INPERI1 G1INPERI2...
	G2INSPYN1 G2INSPYN2 G2INPSD1 G2INPSD2 G2INPERI1 G2INPERI2...
	G1P1SAPM G1P2SAPM G2P1SAPM G2P2SAPM]...
	= SUPERSLOT(stepN,t,GluR1xyl,GluR2xyl,GluR1xyds,GluR2xyds,...
	GluR1_TdwellPSD,GluR1_TdwellPERI,GluR1_TdwellSPYN,LsGR1psa,LsGR1psd,...
	GluR2_TdwellPSD,GluR2_TdwellPERI,GluR2_TdwellSPYN,LsGR2psa,LsGR2psd,...	
	S1,S2,runSAPPSD1,runSAPPSD2,...
	PSD1Sxvecs,PSD1Syvecs,PSD2Sxvecs,PSD2Syvecs,...
	SAPFmx1,SAPFmx2,...
	sap_mask,SSPSD1,SSPSD2,GluR1SdP1,GluR1SdP2,...
	GluR1SLOCMXa,GluR1SLOCMXb,g1polyN1,g2polyN1,...
	GluR2SdP1,GluR2SdP2,GluR2SLOCMXa,GluR2SLOCMXb,...
	LTP1onG1,LTP1offG1,LTP2onG1,LTP2offG1,...
	LTP1onG2,LTP1offG2,LTP2onG2,LTP2offG2,...
	XYLBp1,XYRTp1,XYLBp2,XYRTp2,...
	G1STBASE,G1RTBASE,G1STLTP,G1RTLTP,G1BSMu,G1LSMu,...
	G2STBASE,G2RTBASE,G2STLTP,G2RTLTP,G2BSMu,G2LSMu);
	end
	%===================================%
	% if stepN == 401; keyboard; end
	
	

    %===================================%
    %     SAP CLUSTERS
    %-----------------------------------%
	if doClusters	% ?? doClusters ??
	%-----------------------------------%
	%----S1_MainClusterFun----%
	if runSAPPSD1
	[S1 G1SP1] = S1_MainClusterFun(...
	h_mask, S1L1, S1beta,... 
    S1mu, S1delta_t, S1ro, S1r, S1, nn,...
	G1P1SAPM,G1P2SAPM,G2P1SAPM,G2P2SAPM,...
	G1SP1,G2SP1,G1LTParray,G2LTParray,...
	G1STBASE,G1RTBASE,G1STLTP,G1RTLTP,G1BSMu,G1LSMu,...
	G2STBASE,G2RTBASE,G2STLTP,G2RTLTP,G2BSMu,G2LSMu,...
	LTP1onG1,LTP1offG1,LTP2onG1,LTP2offG1);
	end
	%----S2_MainClusterFun----%
	if runSAPPSD2
	S2 = S2_MainClusterFun(S2pop2_num, h_mask, S2L1_array, S2L2_array, S2beta,... 
	S2mu_array, S2delta_t, S2ro, S2r, S2ro2, S2, nn,...
	G1SAPLOC1,G1SAPLOC2,G2SAPLOC1,G2SAPLOC2,...
	SAPFmx1,SAPFmx2,SAPmx1,SAPmx2,...
	G1P1SAPM,G1P2SAPM,G2P1SAPM,G2P2SAPM);
	end
	
	%----S1sum & S2sum----%
	S1sum = sum(S1(:)); S2sum = sum(S2(:));
	%----PLOTS1S2----%
	if doMainPlot; if mod(stepN, 10) == 0;PLOTS1S2(S1, S2);end;end
	if ~doMainPlot; if mod(stepN, 100) == 0;PLOTS1S2(S1, S2);end;end
	
	%{
	%----SAP MASKS----%
	[SAPFmx1 SAPFmx2 SAPmx1 SAPmx2] = SAPMASKREPORT(Nsteps,...
	S1, S2,um,dot,runSAPPSD1,runSAPPSD2);
	%}
	%-----------------------------------%
	end				% ?? doClusters ??
	%===================================%

	% if stepN == 201; keyboard; end
	
	%{
	%=============================================%
    %     MOVE GLURs SLOTS & DIFFUSION RATE MODS
    %---------------------------------------------%
	
	%=============================================%
				%----PSD SLOTS----%
		%===================================%
		
	%=========================%
	if doKo(1)	% GLUR-1
	%-------------------%
	
	%----G1SAPPOLYGON----%
	[GluR1SP1 GluR1SP2 G1SAPLOC1 G1SAPLOC2]...
	= G1SAPPOLYGON(stepN,...
	PSD1Sxvecs,PSD1Syvecs,PSD2Sxvecs,PSD2Syvecs,...
	GluR1xyl,SAPFmx1,SAPFmx2);

	% useSlotsPSDGR1 = doKo3;
	%----G1SLOTS----%
	[GluR1xyds GluR1_TdwellPSD GluR1_TdwellPERI GluR1_TdwellSPYN G1FSLOTS]...
	= G1SLOTS(stepN,GluR1Ndots,GluR1xyds,...
	GluR1_TdwellPSD,GluR1_TdwellPERI,GluR1_TdwellSPYN,...
	GluR1SP1,GluR1SP2,SAPFmx1,SAPFmx2,G1SAPLOC1, G1SAPLOC2);

	%-------------------%
	end			% GLUR-1
	%=========================%
	
	%=========================%
	if doKo(2)	% GLUR-2
	%-------------------%
	
	%----G2SAPPOLYGON----%
	[GluR2SP1 GluR2SP2 G2SAPLOC1 G2SAPLOC2]...
	= G2SAPPOLYGON(stepN,...
	PSD1Sxvecs,PSD1Syvecs,PSD2Sxvecs,PSD2Syvecs,...
	GluR2xyl,SAPFmx1,SAPFmx2);

	% useSlotsPSDGR2 = doKo4;
	%----G2SLOTS----%
	[GluR2xyds GluR2_TdwellPSD GluR2_TdwellPERI GluR2_TdwellSPYN G2FSLOTS]... 
	= G2SLOTS(stepN,GluR2Ndots,GluR2xyds,...
    GluR2_TdwellPSD,GluR2_TdwellPERI,GluR2_TdwellSPYN,...
    GluR2SP1,GluR2SP2,SAPFmx1,SAPFmx2,G2SAPLOC1, G2SAPLOC2);

	%-------------------%
	end			% GLUR-2
	%=========================%
	% if stepN == 200; keyboard; end
	
	
	%=============================================%
				%----PERI SLOTS----%
		%===================================%
	if doKo5	% USE PERI SLOTS GR1
	%----PERI SLOTS GluR1----%
	[GluR1xyds GluR1_TdwellPSD GluR1_TdwellPERI G1FSLOTS...
	] = G1PERISLOTS(stepN,GluR1Ndots,GluR1xyds,...
	GluR1_TdwellPSD,GluR1_TdwellPERI,...
	ko11,ko12,ko15,ko16,...
	doKo5,S1sum,S2sum,SAP5);
	end
	if doKo6	% USE PERI SLOTS GR2
	%----PERI SLOTS GluR2----%
	[GluR2xyds GluR2_TdwellPSD GluR2_TdwellPERI G2FSLOTS...
    ] = G2PERISLOTS(stepN,GluR2Ndots,GluR2xyds,...
    GluR2_TdwellPSD,GluR2_TdwellPERI,...
    ko11,ko12,ko15,ko16,...
    doKo6,S1sum,S2sum,SAP5);
	end
	%=============================================%
	%}
	%{
		%===================================%
    %     GLUR+SAP >> SLOTS
    %-----------------------------------%
	if doKo(1)	% ?? useGluR1slots ??
	%-----------------------------------%
	%----G1SAPPOLYGON----%
	[GluR1SdP1 GluR1SdP2 G1SAPLOCP1 G1SAPLOCP2]...
	= G1SAPPOLYGON(stepN,...
	PSD1Sxvecs,PSD1Syvecs,PSD2Sxvecs,PSD2Syvecs,...
	GluR1xyl,SAPFmx1,SAPFmx2);
	%----G1SLOTS----%
	[GluR1xyds GluR1_TdwellPSD GluR1_TdwellPERI G1FSLOTS...
	] = G1SLOTS(stepN,GluR1Ndots,GluR1xyds,...
	D,t,GluR1_TdwellPSD,GluR1_TdwellPERI,...
	ko1,ko2,ko3,ko4,ko5,ko6,ko7,ko8,ko9,ko10,ko11,ko12,ko13,ko14,ko15,ko16,...
	doKo1,doKo2,doKo3,doKo4,doKo5,doKo6,S1sum,S2sum,SAP5,...
	GluR1SdP1,GluR1SdP2,SAPFmx1,SAPFmx2,...
	G1SAPLOCP1, G1SAPLOCP2);
	%-----------------------------------%
	end			% ?? useGluR1slots ??
	%===================================%
	if doKo(2)	% ?? useGluR2slots ??
	%-----------------------------------%
	%----G2SAPPOLYGON----%
	[GluR2SdP1 GluR2SdP2]...
	= G2SAPPOLYGON(stepN,...
	PSD1Sxvecs,PSD1Syvecs,PSD2Sxvecs,PSD2Syvecs,...
	GluR2xyl,SAPFmx1,SAPFmx2);
	%----G2SLOTS----%
	[GluR2xyds GluR2_TdwellPSD GluR2_TdwellPERI G2FSLOTS...
    ] = G2SLOTS(stepN,GluR2Ndots,GluR2xyds,...
    D,t,GluR2_TdwellPSD,GluR2_TdwellPERI,...
    ko1,ko2,ko3,ko4,ko5,ko6,ko7,ko8,ko9,ko10,ko11,ko12,ko13,ko14,ko15,ko16,...
    doKo1,doKo2,doKo3,doKo4,doKo5,doKo6,S1sum,S2sum,SAP5,...
    GluR2SdP1,GluR2SdP2,SAPFmx1,SAPFmx2);
	%-----------------------------------%
	end			% ?? useGluR2slots ??
	%===================================%
	%}
	
	
	
	%=============================================%
    %		MOVE PARTICLES MAIN FUNCTION
	%=============================================%
	[GluR1xyds GluR1xyl GluR2xyds GluR2xyl]...
    = MOVEGLUR(stepN,XWIDE,YHIGH,...
	GluR1Ndots,GluR1xyds,GluR1xyl,...
    GluR2Ndots,GluR2xyds,GluR2xyl);
	
	%{
				%----MOVE GluR1----%
		%===================================%
	[GluR1xyds GluR1xyl GluR1_TdwellPSD GluR1_TdwellPERI GluR1_TdwellSPYN...
	G1INSPYN1 G1INSPYN2 G1INPSD1 G1INPSD2 G1INPERI1 G1INPERI2]...
    = MOVEGLUR1(stepN,GluR1Ndots,GluR1xyds,GluR1xyl, PSD1, PSD2,...
    D,t,GluR1_TdwellPSD,GluR1_TdwellPERI,GluR1_TdwellSPYN,XWIDE,YHIGH,...
	SPYN1xv,SPYN1yv,SPYN2xv,SPYN2yv,...
	PSD1xv,PSD1yv,PSD2xv,PSD2yv,...
	PERI1xv,PERI1yv,PERI2xv,PERI2yv,...
	LsGR1psa,LsGR1psd,stky);
	%=============================================%
				%----MOVE GluR2----%
		%===================================%
	[GluR2xyds GluR2xyl GluR2_TdwellPSD GluR2_TdwellPERI GluR2_TdwellSPYN...
	G2INSPYN1 G2INSPYN2 G2INPSD1 G2INPSD2 G2INPERI1 G2INPERI2]...
	= MOVEGLUR2(stepN,GluR2Ndots,GluR2xyds,GluR2xyl, PSD1, PSD2,...
    D,t,GluR2_TdwellPSD,GluR2_TdwellPERI,GluR2_TdwellSPYN,XWIDE,YHIGH,...
	SPYN1xv,SPYN1yv,SPYN2xv,SPYN2yv,...
	PSD1xv,PSD1yv,PSD2xv,PSD2yv,...
	PERI1xv,PERI1yv,PERI2xv,PERI2yv,...
	LsGR2psa,LsGR2psd,stky);
	%=============================================%
	%}
	
	
	
	
	
	%-------------------------------%
    %   3D VISUALS
    %-------------------------------%
	if do3DPLOT
	[G1Z G2Z] = ZGEN(stepN,G1Z,G2Z,...
    GluR2Ndots, GluR2xyl,GluR1Ndots, GluR1xyl,...
	XYLTpr1,XYRTpr1,XYLBpr1,XYRBpr1,XYLTpr2,XYRTpr2,XYLBpr2,XYRBpr2,...
	XYLTp1,XYRTp1,XYLBp1,XYRBp1,XYLTp2,XYRTp2,XYLBp2,XYRBp2);

	PLOT3DS(stepN,GluR2xyl,GluR1xyl,XWIDE,YHIGH,G1Z,G2Z,Zfield);
	end
	

    %-------------------------------%
    %   MAIN PLOT XYL SIMULATION
    %-------------------------------%
	if doMainPlot
	MAINPLOT(stepN,GluR2xyl,GluR1xyl,...
	XWIDE,YHIGH,XYBOXpr1,XYBOXpr2,XYBOXp1,XYBOXp2);
	end

	if ~doMainPlot && ~doONEDOTPLOT
	if mod(stepN, 200) == 0
	MAINPLOT(stepN,GluR2xyl,GluR1xyl,...
	XWIDE,YHIGH,XYBOXpr1,XYBOXpr2,XYBOXp1,XYBOXp2);
	end
	end
    
	

	%-------------------------------%
    %   ONE DOT PLOT GRAPHICS
    %-------------------------------%
	if doONEDOTPLOT
	doMainPlot=0;
	xyl2(:,stepN) = GluR2xyl(:,1);
	ONEDOTPLOT2(stepN,Nsteps,GluR2xyl,xyl2,XWIDE,YHIGH,...
				XYBOXpr1,XYBOXpr2,XYBOXp1,XYBOXp2)
	end
	
	
  %---------------------------------------------%
  %       START PARTICLE COUNTERS
  %---------------------------------------------%
  if mod(stepN, DATARATE) == 0
	  
	% if stepN == 201; keyboard; end
    %-------------------------------%
    %       PARTICLE COUNTERS
    %-------------------------------%
	[GR1SPY1N GR1SPY2N GR1PER1N GR1PER2N GR1PSD1N GR1PSD2N...
	GR2SPY1N GR2SPY2N GR2PER1N GR2PER2N GR2PSD1N GR2PSD2N...
	GR1PERN GR1PSDN GR2PERN GR2PSDN GR1SPYN GR2SPYN...
	SPY1N SPY2N SPYNT ESNT EST]...
	= dotCount(stepN,GluR2Ndots, GluR1Ndots,...
	G2INSPYN1,G2INSPYN2,G2INPSD1,G2INPSD2,G2INPERI1,G2INPERI2,...
	G1INSPYN1,G1INSPYN2,G1INPSD1,G1INPSD2,G1INPERI1,G1INPERI2);
	
	%-------------------------------%
    %   GRAPH PARTICLE COUNTS
    %-------------------------------%
	if doMainPlot
	GRAPHCOUNTS(stepN,GluR2Ndots,doGluR1,doGluR2,doRun1,doRun2,doRun6,...
	HSlow, HShi, G1FSLOTS,G2FSLOTS,SAP5,S1sum,S2sum,doKo,...
	EST,SPY1N,SPY2N);
	end
	
	%-------------------------------%
    % SAVE DATA AND STREAM OUTPUT
    %-------------------------------%
       
       data = [stepN/Nsteps ESNT SPY1N SPY2N SPYNT;... 
			   (Nsteps-stepN) D/10 Dn_PSD1/10 Dn_PSD2/10 k];
       % data1(stepN/DATARATE,:) = data(1,:);
       % data2(stepN/DATARATE,:) = data(2,:);
	   % dataset = cat(2, data1, data2);
       dataset(stepN/DATARATE,:) = [data(1,:) data(2,:)];
	   SAPdata(stepN/DATARATE,:) = [S1sum S2sum];
	   Ddata(stepN/DATARATE,:) = [Dn_PSD1/10 Dn_PSD2/10];
	   AMPARdata(stepN/DATARATE,:) = [GR1SPYN GR2SPYN];
	   GluRdata(stepN/DATARATE,:) = [GR1PSDN GR1PERN GR2PSDN GR2PERN...
								GR1SPY1N GR1SPY2N GR2SPY1N GR2SPY2N];
	   G1SLOTSDATA(stepN/DATARATE,:) = G1FSLOTS;
	   G2SLOTSDATA(stepN/DATARATE,:) = G2FSLOTS; 
	   
	   if mod(stepN, DATARATE*10) == 0
	   head1 = ('Time%2go    CntDwn');
	   disp(head1);
	   disp([data(1,1) data(2,1)])
	   end
       
       snapnow;
	end
 %---------------------------------------------%
 %       END PARTICLE COUNTERS
 %---------------------------------------------%
	
	

%=============%%%%%%%%%%%%%%%%%%%%%%%%%%%%=============%
% if stepN == 400; keyboard; end 
stepN = stepN+1;
nn = stepN;
end
%=============%%%%%%%%%%%%%%%%%%%%%%%%%%%%=============%
%%					END Main Loop
%=============%%%%%%%%%%%%%%%%%%%%%%%%%%%%=============%

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
%=========================================================%
end % for re = 1:loops
%=========================================================%


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
%===========================================%
%			OUTPUT GRAPHS
%===========================================%
% filter(ones(1,10)/10,1,G1SLOTSDATA(:,1)) % running mean
% set(gcf,'OuterPosition',[100,100,1100,800])
% [left, bottom, width, height]
% SNSZ = [1, 1, 1680, 1028;]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			FIGURE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',...
	[SNSZ(3)/SNSZ(3)	SNSZ(4)/SNSZ(4)		SNSZ(3)/2	SNSZ(4)/2.3]);


%%%%%%%%%%%%%%%%%%
figure(2)
subplot(2,2,1), subplot('Position',[.05 .55 .4 .4]),...
	plot([dataset(:,3) dataset(:,4) dataset(:,5) dataset(:,2)])
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
F2AT = title('Receptors Per Region'); leg1=legend('SPINE1', 'SPINE2', 'SPINES', 'ES');
set(F2AT, 'FontSize', 14);
set(leg1,'Location','NorthWest');
%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%
figure(2)
subplot(2,2,3), subplot('Position',[.55 .75 .4 .2]),...
	plot([GluRdata(:,1) GluRdata(:,2)])
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
GluR1PSDvsPERI=title('GluR1 in PSD vs PERI'); leg1=legend('PSD', 'PERI');
set(GluR1PSDvsPERI, 'FontSize', 14);
set(leg1,'Location','NorthWest');
	
figure(2)
subplot(2,2,3), subplot('Position',[.55 .55 .4 .2]),...
	plot([GluRdata(:,3) GluRdata(:,4)])
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
GluR2PSDvsPERI=title('GluR2 PSD vs PERI'); leg1=legend('PSD', 'PERI');
set(GluR2PSDvsPERI, 'FontSize', 14);
set(leg1,'Location','NorthWest');
%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%
figure(2)
subplot(2,2,3), subplot('Position',[.05 .05 .4 .4]),...
	plot([GluRdata(:,5) GluRdata(:,6) GluRdata(:,7) GluRdata(:,8)])
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
F2CT = title('GluR Subtypes in SPINE1 vs SPINE2'); 
set(F2CT, 'FontSize', 14);
leg1=legend('SPINE1-GluR1', 'SPINE2-GluR1',	'SPINE1-GluR2', 'SPINE2-GluR2');
set(leg1,'Location','NorthWest');
%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%
figure(2)
subplot(2,2,4), subplot('Position',[.55 .05 .4 .4]),...
	plot([dataset(:,3) dataset(:,4) dataset(:,5)]);
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
F2DT=title('SPINE Totals'); leg2=legend('SPINE-1', 'SPINE-2', 'Total');
set(F2DT, 'FontSize', 14);
set(leg2,'Location','SouthWest');
%%%%%%%%%%%%%%%%%%





%{
% figure(2)
% subplot(2,2,2), subplot('Position',[.55 .55 .4 .4]),...
% 	plot([GluRdata(:,1) GluRdata(:,2) GluRdata(:,3) GluRdata(:,4)])
% xt = ((get(gca,'XTick'))*10);
% set(gca,'XTickLabel', sprintf('%.1f|',xt))
% title('GluR Subtypes in PSD vs periPSD'); leg1=legend('GR1PSDN', 'GR1PERN', 'GR2PSDN', 'GR2PERN');
% set(leg1,'Location','NorthWest');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			FIGURE 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set(gcf,'OuterPosition',[50,50,1100,800])
figure('Position',...
	[SNSZ(3)/SNSZ(3)	SNSZ(4)/2			SNSZ(3)/2	SNSZ(4)/2.3]);

%%%%%%%%%%%%%%%%%%
figure(3)
subplot(2,2,1), subplot('Position',[.05 .55 .4 .4]),...
	plot([AMPARdata(:,1) AMPARdata(:,2)]);
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
F3AT=title('AMPAR SUBTYPE TOTALS IN SPINES'); leg3=legend('GluR1', 'GluR2');
set(F3AT, 'FontSize', 14);
set(leg3,'Location','SouthWest');
%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%
figure(3)
subplot(2,2,3), subplot('Position',[.55 .75 .4 .2]),...
	plot([GluRdata(:,5) GluRdata(:,6)])
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'YLabel'),'String','GLUR1 PARTICLES')
F3B2T=title('GluR1');
set(F3B2T, 'FontSize', 14);
leg1=legend('SPINE1', 'SPINE2');
set(leg1,'Location','NorthWest');

figure(3)
subplot(2,2,3), subplot('Position',[.55 .55 .4 .2]),...
	plot([GluRdata(:,7) GluRdata(:,8)])
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'YLabel'),'String','GLUR2 PARTICLES')
F3B1T=title('GluR2');
set(F3B1T, 'FontSize', 14);
leg1=legend('SPINE1', 'SPINE2');
set(leg1,'Location','NorthWest');
	
%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%
figure(3)
subplot(2,2,3), subplot('Position',[.55 .25 .4 .2]),...
	plot([G1SLOTSDATA(:,1) G1SLOTSDATA(:,2) G1SLOTSDATA(:,3) G1SLOTSDATA(:,4)])
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'YLabel'),'String','FILLED SLOTS')
F3B1T=title('GluR1');
set(F3B1T, 'FontSize', 14);
leg1=legend('PSD-1', 'PSD-2','PERI-1','PERI-2');
set(leg1,'Location','NorthWest');
		


figure(3)
subplot(2,2,3), subplot('Position',[.55 .05 .4 .2]),...
	plot([G2SLOTSDATA(:,1) G2SLOTSDATA(:,2) G2SLOTSDATA(:,3) G2SLOTSDATA(:,4)])
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','FILLED SLOTS')
F3B2T=title('GluR2');
set(F3B2T, 'FontSize', 14);
leg1=legend('PSD-1', 'PSD-2','PERI-1','PERI-2');
set(leg1,'Location','NorthWest');
%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%
figure(3)
subplot(2,2,3), subplot('Position',[.05 .05 .4 .4]),...
	plot([SAPdata(:,1) SAPdata(:,2)]);
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','SAP PARTICLES')
F3CT=title('SAP Cluster Particles'); leg4=legend('SAP1', 'SAP2');
set(F3CT, 'FontSize', 14);
set(leg4,'Location','SouthWest');
%%%%%%%%%%%%%%%%%%




%{
figure(3)
subplot(2,2,4), subplot('Position',[.55 .05 .4 .4]),...
	plot([(Ddata(:,1)) (Ddata(:,2))]);
yt = get(gca,'YTick');
set(gca,'YTickLabel', sprintf('%.4f|',yt))
set(get(gca,'YLabel'),'String','scaled diffusion rate')
xt = (get(gca,'XTick'))*10;
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time-Step')
F3DT=title('PSD DIFFUSION RATE'); leg4=legend('D-SPINE1', 'D-SPINE2');
set(F3DT, 'FontSize', 16);
set(leg4,'Location','SouthWest');
%}
%-------------##########################------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			FIGURE 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',...
	[SNSZ(3)/SNSZ(3)	SNSZ(4)/SNSZ(4)		SNSZ(3)/2	SNSZ(4)/2.3]);


%%%%%%%%%%%%%%%%%%
figure(4)
subplot(2,2,1), subplot('Position',[.05 .55 .4 .4]),...
	plot([DATAdataset(:,3) DATAdataset(:,4) DATAdataset(:,5) DATAdataset(:,2)])
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
F2AT = title('Receptors Per Region'); leg1=legend('SPINE1', 'SPINE2', 'SPINES', 'ES');
set(F2AT, 'FontSize', 14);
set(leg1,'Location','NorthWest');
%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%
figure(4)
subplot(2,2,3), subplot('Position',[.55 .75 .4 .2]),...
	plot([DATAGluRdata(:,1) DATAGluRdata(:,2)])
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
GluR1PSDvsPERI=title('GluR1 in PSD vs PERI'); leg1=legend('PSD', 'PERI');
set(GluR1PSDvsPERI, 'FontSize', 14);
set(leg1,'Location','NorthWest');
	
figure(4)
subplot(2,2,3), subplot('Position',[.55 .55 .4 .2]),...
	plot([DATAGluRdata(:,3) DATAGluRdata(:,4)])
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
GluR2PSDvsPERI=title('GluR2 PSD vs PERI'); leg1=legend('PSD', 'PERI');
set(GluR2PSDvsPERI, 'FontSize', 14);
set(leg1,'Location','NorthWest');
%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%
figure(4)
subplot(2,2,3), subplot('Position',[.05 .05 .4 .4]),...
	plot([DATAGluRdata(:,5) DATAGluRdata(:,6) DATAGluRdata(:,7) DATAGluRdata(:,8)])
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
F2CT = title('GluR Subtypes in SPINE1 vs SPINE2'); 
set(F2CT, 'FontSize', 14);
leg1=legend('SPINE1-GluR1', 'SPINE2-GluR1',	'SPINE1-GluR2', 'SPINE2-GluR2');
set(leg1,'Location','NorthWest');
%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%
figure(4)
subplot(2,2,4), subplot('Position',[.55 .05 .4 .4]),...
	plot([DATAdataset(:,3) DATAdataset(:,4) DATAdataset(:,5)]);
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
F2DT=title('SPINE Totals'); leg2=legend('SPINE-1', 'SPINE-2', 'Total');
set(F2DT, 'FontSize', 14);
set(leg2,'Location','SouthWest');
%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			FIGURE 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set(gcf,'OuterPosition',[50,50,1100,800])
figure('Position',...
	[SNSZ(3)/SNSZ(3)	SNSZ(4)/2			SNSZ(3)/2	SNSZ(4)/2.3]);

%%%%%%%%%%%%%%%%%%
figure(5)
subplot(2,2,1), subplot('Position',[.05 .55 .4 .4]),...
	plot([DATAAMPARdata(:,1) DATAAMPARdata(:,2)]);
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
F3AT=title('AMPAR SUBTYPE TOTALS IN SPINES'); leg3=legend('GluR1', 'GluR2');
set(F3AT, 'FontSize', 14);
set(leg3,'Location','SouthWest');
%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%
figure(5)
subplot(2,2,3), subplot('Position',[.55 .75 .4 .2]),...
	plot([DATAGluRdata(:,5) DATAGluRdata(:,6)])
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'YLabel'),'String','GLUR1 PARTICLES')
F3B2T=title('GluR1');
set(F3B2T, 'FontSize', 14);
leg1=legend('SPINE1', 'SPINE2');
set(leg1,'Location','NorthWest');

figure(5)
subplot(2,2,3), subplot('Position',[.55 .55 .4 .2]),...
	plot([DATAGluRdata(:,7) DATAGluRdata(:,8)])
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'YLabel'),'String','GLUR2 PARTICLES')
F3B1T=title('GluR2');
set(F3B1T, 'FontSize', 14);
leg1=legend('SPINE1', 'SPINE2');
set(leg1,'Location','NorthWest');
	
%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%
figure(5)
subplot(2,2,3), subplot('Position',[.55 .25 .4 .2]),...
	plot([DATAG1SLOTSDATA(:,1) DATAG1SLOTSDATA(:,2) DATAG1SLOTSDATA(:,3) DATAG1SLOTSDATA(:,4)])
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'YLabel'),'String','FILLED SLOTS')
F3B1T=title('GluR1');
set(F3B1T, 'FontSize', 14);
leg1=legend('PSD-1', 'PSD-2','PERI-1','PERI-2');
set(leg1,'Location','NorthWest');
		


figure(5)
subplot(2,2,3), subplot('Position',[.55 .05 .4 .2]),...
	plot([DATAG2SLOTSDATA(:,1) DATAG2SLOTSDATA(:,2) DATAG2SLOTSDATA(:,3) DATAG2SLOTSDATA(:,4)])
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','FILLED SLOTS')
F3B2T=title('GluR2');
set(F3B2T, 'FontSize', 14);
leg1=legend('PSD-1', 'PSD-2','PERI-1','PERI-2');
set(leg1,'Location','NorthWest');
%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%
figure(5)
subplot(2,2,3), subplot('Position',[.05 .05 .4 .4]),...
	plot([DATASAPdata(:,1) DATASAPdata(:,2)]);
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','SAP PARTICLES')
F3CT=title('SAP Cluster Particles'); leg4=legend('SAP1', 'SAP2');
set(F3CT, 'FontSize', 14);
set(leg4,'Location','SouthWest');
%%%%%%%%%%%%%%%%%%
%===========================================%



%%
%===========================================%
%   IF TRACK MSD
%-------------------------------%
% TRACK MSD AND STEP SIZE MEANS (0:no 1:yes)
trackMSD = doRun(3);
if trackMSD

dT = t;
GluR2Ndots = 100;
Nsteps = GluR2Ndots;
tracks = cell(GluR2Ndots, 1);
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
    GluR2xyds = STEPxyds(GluR2Ndots, k);
	%-------------------------------%
    %     DO MANUAL STEP SIZES
    %-------------------------------%
	if MANstepsize
		GluR2xyds = CIRCSTEPxyds(GluR2Ndots, k);
	end % xyd = DIRxyd(xyd);
	
		
	[GluR2xyl GluR2xyds] = MSDAMPARSTEP(GluR2Ndots, GluR2xyds, GluR2xyl,... 
	testPSD1MSD, testPSD2MSD, testESMSD, PSD1, PSD2, D, Dn_PSD1, Dn_PSD2);


    [tracks] = MSDfun(stepN, Nsteps, tracks, GluR2xyds);

    
    if trackStepMeans
     if mod(stepN, 10) == 0
        StepMeans = meanStep(GluR2xyds, k, D);
        head = ['MeanAbsStep CalcHalfNorm MeanPythStep CalcPythStep CalcHalfPythStep'];
       disp(head)
       disp(StepMeans)
     end
    end


 stepN = stepN+1;
end
output1 = tracks;
[rawMSDout] = MSDfunction(tracks, GluR2Ndots, Nsteps, D, Dn_PSD1, Dn_PSD2, dT, k, MSDtest);


%===========================================================%
%               LIVE PARTICLE DIFFUSION
%-----------------------------------------------------------%
for Nt = 1:Nsteps

    GluR2xyds = STEPxyds2(GluR2Ndots, k);
	[GluR2xyl] = AMPARSTEP2(GluR2Ndots, GluR2xyds, GluR2xyl);
	MAINPLOT2(GluR2xyl, lims);

end
%===========================================================%


%===========================================================%
%               MSD RANDOM STEPS ANALYSIS
%-----------------------------------------------------------%
tracks = cell(GluR2Ndots, 1);

stepN = 1;
for Nt = 1:Nsteps 

	GluR2xyds = STEPxyds2(GluR2Ndots, k);
	[GluR2xyl] = AMPARSTEP2(GluR2Ndots, GluR2xyds, GluR2xyl);
    [tracks] = MSDfun2(stepN, Nsteps, tracks, GluR2xyds);

stepN = stepN+1;
end
[ESMSDout] = MSDfunction2(tracks,GluR2Ndots,Nsteps,D,Dn,L,dT,k,Scale,MSDtest);

%===========================================================%
%               MSD UNIFORM STEPS ANALYSIS
%-----------------------------------------------------------%
stepN = 1;
for t = 1:Nsteps 

	GluR2xyds = UNIFORMSTEPS2(GluR2Ndots, Lx);
	[GluR2xyl GluR2xyds] = SCALEUNIFORMSTEPS2(GluR2Ndots, GluR2xyds, GluR2xyl, Ls, MSDtest);
    [tracks] = MSDfun2(stepN, Nsteps, tracks, GluR2xyds);

stepN = stepN+1;
end
[PSDMSDout] = MSDfunction2(tracks,GluR2Ndots,Nsteps,D,Dn,L,dT,k,Scale,MSDtest);

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



%----------------------------------------%
%   Do Izhikevich Spiking
%----------------------------------------%
if doIz
spikes=DATAdataset(:,3:4);
IzhFun(spikes);
end



%%
%----------------------------------------%
%   OUTPUT RETURNS FOR TOP LEVEL FUNCTION
%----------------------------------------%

% varargout = {GluRdata;SAPdata;Ddata;AMPARdata};
varargout = {DATAdataset;DATASAPdata;DATADdata;DATAAMPARdata;DATAGluRdata;...
	DATAG1SLOTSDATA;DATAG2SLOTSDATA};
assignin('base', 'varargout', varargout)

%{
time1 = datestr(clock, 0);
time2 = num2str(time1(1,[1:2 16:17 19:20]));
time3 = strcat(['BradsModelData' time2 '.txt']);
time4 = strcat(['BradsModelData' time2 '.txt']);
XLSdata1 = dataset({varargout{1}});
export(XLSdata1,'file',eval('time3'));
%}

if doprofile
profile viewer;
end
%-------------############################------------------%
end %		  ##    END MAIN FUNCTION   ##
%-------------############################------------------% 














%%						 SUBROUTINES
%=================########################=================%
%				  ##		INDEX		##
%=================########################=================% 
%				1. BROWNIAN MOTION FUNCTIONS
%				2. SAPS & SLOTS
%				3. PARTICLE COUNTERS
%				4. DENDRITIC FIELD MAPPING TOOLS
%				5. PLOTTING AND LIVE SIMULATION
%				6. MSD BROWNIAN MOTION ANALYSIS TOOLS
%=================############################=================%



%%			  BROWNIAN MOTION FUNCTIONS
%-------------##########################------------------%
%			 PARTICLE MOVEMENT FUNCTIONS
%-------------##########################------------------%


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
function [GluR1xyds GluR1xyl GluR2xyds GluR2xyl]...
    = MOVEGLUR(stepN,XWIDE,YHIGH,...
	GluR1Ndots,GluR1xyds,GluR1xyl,...
    GluR2Ndots,GluR2xyds,GluR2xyl)

%===================================%
GluR1xyl = GluR1xyl+GluR1xyds;
 
for j = 1:GluR1Ndots 
    if GluR1xyl(1,j)>(XWIDE) || GluR1xyl(1,j)<(0)
            GluR1xyl(1,j) = uint8(sign(GluR1xyl(1,j)))*(XWIDE);
    elseif GluR1xyl(2,j)>(YHIGH) || GluR1xyl(2,j)<(-YHIGH)
            GluR1xyl(2,j) = sign(GluR1xyl(2,j))*(YHIGH);
    end    
        
end
%===================================%
GluR2xyl = GluR2xyl+GluR2xyds;

for j = 1:GluR2Ndots 
	if GluR2xyl(1,j)>(XWIDE) || GluR2xyl(1,j)<(0)
            GluR2xyl(1,j) = uint8(sign(GluR2xyl(1,j)))*(XWIDE);
	elseif GluR2xyl(2,j)>(YHIGH) || GluR2xyl(2,j)<(-YHIGH)
            GluR2xyl(2,j) = sign(GluR2xyl(2,j))*(YHIGH);
	end	   
        
end
%===================================%
% if stepN == 500; keyboard; end
end
%===================================%



%%					 SAPS & SLOTS
%-------------##########################------------------%
%			SAPS SLOTS & S-CLUSTER FUNCTIONS
%-------------##########################------------------%

%-----------------------------------------%
%			PRE-LOOP SETUP
%-----------------------------------------%

%===================================%
% SAPSLOTSETUP
%===================================%
% Establish the initial scalars, vectors, and matrices
% for the SAP cluster turnnover functions
%-------------------%
function [Nsteps,h_mask,...
S1pop2_num,S1r,SC1L,...
S1,S1ro,SC1L2,S1L1,S1ro2,SC1beta,S1L1_array,S1L2,S1sum,SC1mu,...
S1L2_array,SC1r,S1beta,S1tau,SC1ro,S1delta_t,S1mem_num,S1mu,SC1tau...
S1mu_array,S1num_epochs,...
S2pop2_num,S2r,SC2L,...
S2,S2ro,SC2L2,S2L1,S2ro2,SC2beta,S2L1_array,S2L2,S2sum,SC2mu,...
S2L2_array,SC2r,S2beta,S2tau,SC2ro,S2delta_t,S2mem_num,S2mu,SC2tau...
S2mu_array,S2num_epochs] = SAPSLOTSETUP(Nsteps,...
sap1, sap2, sap3, sap6, sap4, sap7, sap5, sap8,runSAPPSD1,runSAPPSD2,...
um,Zfield,fPSD1,fPSD2,dot,SAPPADPSD1,SAPPADPSD2)

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



SAPSTEP = .1;
Steps = dot(3); % Steps
TimeStep = dot(4)/1000; % TimeStep
Scale = dot(5); % Scale
% Scale = Scale/10;

time=Nsteps;
h_mask=[0 1 0; 1 0 1; 0 1 0];

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




%{
if SAPPADPSD1
SAPMXPSD1=padarray(ones(PSD1SAPF),[PSA1SAPF PSA1SAPF], 0);
S1=SAPMXPSD1;
end

if SAPPADPSD2
SAPMXPSD2=padarray(ones(PSD2SAPF),[PSA2SAPF PSA2SAPF], 0);
S2=SAPMXPSD2;
end
%}


%===========================================%
SC1szi = sap1;			SC2szi = sap2;
SC1beta = sap3;			SC2beta = sap6;		
SC1tau = sap4;			SC2tau = sap7;		
SC1mu = 1/SC1tau;		SC2mu = 1/SC2tau;	
SC1L = sap5;			SC2L = sap8;
SC1L2 = 1.1;			SC2L2 = 1.1;		
SC1r = 15;				SC2r = 15;			
SC1ro = 0.90;			SC2ro = 0.90;			
SC1szb = 17;			SC2szb = 17;

SC1deltaT = TimeStep*SAPSTEP;	SC2deltaT = TimeStep*SAPSTEP;

%===========================================%
% if runSAPPSD1 || runSAPPSD2 %   S1 PSD1
%===========================================%

%==================%
%	PSD1 S1 SAP
%==================%

S1delta_t=SC1deltaT; % smaller yields slower changes [0.1]
S1num_epochs=floor(time/S1delta_t)+1; % 1001
S1ro2=ones(1,S1num_epochs)*0.00;
S1ro2(floor(15.5/S1delta_t):floor(16.5/S1delta_t))=0.02;

S1beta=SC1beta;	% slope of on-rate 'h' probability insertion [60]
				% beta>100 low P() of psudopodia growth; beta<20 high P()

S1tau=SC1tau;	% off rate  random degredation [1] lower=fast degrade
S1mu=SC1mu;		% on rate   random degredation [1]
 
S1r=SC1r;		% receptor transition rate into unoccupied site [10]
S1ro=SC1ro;		% cellular concentration [.95]
 
S1L1=SC1L;		% repulsive lattice constant [1.2 - 2]
				% [3] MAKES IT SHRINK [1.2] MAKES IT GROW
S1L2=SC1L2;		% repulsive lattice constant [0.9] active during LTP window
 
S1L1_array=S1L1*ones(1,S1num_epochs);
S1L2_array=S1L2*ones(1,S1num_epochs);
S1mu_array=S1mu*ones(1,S1num_epochs);

S1s=S1;
S1mem_num=zeros(1,S1num_epochs);
S1mem_num(1)=sum(S1(:));
S1pop2_num=zeros(1,S1num_epochs);


%==================%
%	PSD2 S2 SAP
%==================%

S2delta_t=SC2deltaT; % smaller yields slower changes [0.1]
S2num_epochs=floor(time/S2delta_t)+1; % 1001
S2ro2=ones(1,S2num_epochs)*0.00;
S2ro2(floor(15.5/S2delta_t):floor(16.5/S2delta_t))=0.02;

S2beta=SC2beta;	% slope of on-rate 'h' probability insertion [60]
				% beta>100 low P() of psudopodia growth; beta<20 high P()

S2tau=SC2tau;	% off rate  random degredation [1] lower=fast degrade
S2mu=SC2mu;		% on rate	random degredation [1]

S2r=SC2r;		% receptor transition rate into unoccupied site [10]
S2ro=SC2ro;		% cellular concentration [.95]

S2L1=SC2L;		% repulsive lattice constant [1.2 - 2]
				% [3] MAKES IT SHRINK [1.2] MAKES IT GROW
S2L2=SC2L2;		% repulsive lattice constant [0.9]

S2L1_array=S2L1*ones(1,S2num_epochs);
S2L2_array=S2L2*ones(1,S2num_epochs);
S2mu_array=S2mu*ones(1,S2num_epochs);
S2s = S2;

S2mem_num=zeros(1,S2num_epochs);
S2mem_num(1)=sum(S2(:));
S2pop2_num=zeros(1,S2num_epochs);

%{
sz2_bound=SC2szb;
sz2=sz2_bound-2;
S2=zeros(sz2+2);
center2=floor(sz2_bound/2)+1;
initial_down2=center2-floor(SC2szi/2);
initial_up2=initial_down2+SC2szi-1;
S2(initial_down2:initial_up2,initial_down2:initial_up2)=1;

S2delta_t=SC2deltaT; % smaller yields slower changes [0.1]
S2num_epochs=floor(time/S2delta_t)+1; % 1001
S2ro2=ones(1,S2num_epochs)*0.00;
S2ro2(floor(15.5/S2delta_t):floor(16.5/S2delta_t))=0.02;

S2beta=SC2beta;		% slope of on-rate 'h' probability insertion [60]
				% beta>100 low P() of psudopodia growth; beta<20 high P()

S2tau=SC2tau;		% off rate  random degredation [1] lower=fast degrade
S2mu=SC2mu;	% on rate	random degredation [1]

S2r=SC2r;			% receptor transition rate into unoccupied site [10]
S2ro=SC2ro;		% cellular concentration [.95]

% THIS MAKES THE CLUSTERS GROW OR SHRINK
S2L1=SC2L;		% repulsive lattice constant [1.2 - 2]
				% [3] MAKES IT SHRINK [1.2] MAKES IT GROW
S2L2=SC2L2;		% repulsive lattice constant [0.9]

S2L1_array=S2L1*ones(1,S2num_epochs);
S2L2_array=S2L2*ones(1,S2num_epochs);
S2mu_array=S2mu*ones(1,S2num_epochs);
S2s = S2;

S2mem_num=zeros(1,S2num_epochs);
S2mem_num(1)=sum(S2(:));
S2pop2_num=zeros(1,S2num_epochs);
%}
%===========================================%
% end % if runSAPPSD2 %   S2 PSD2
%===========================================%

% if ~runSAPPSD1
% S1=0;
% end
% 
% if ~runSAPPSD2
% S2=0;
% end

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
% for the small polygons that are uniformly distributed throughout the PSD
% and return those Cells outside of the main loop
%-------------------%
function [PSD1Sxvecs PSD1Syvecs PSD2Sxvecs PSD2Syvecs] = SAPMAP(...
	um,dot,...
	SPYN1xv,SPYN1yv,SPYN2xv,SPYN2yv,...
	PSD1xv,PSD1yv,PSD2xv,PSD2yv,...
	PERI1xv,PERI1yv,PERI2xv,PERI2yv,SAPPADPSD1,SAPPADPSD2)



Scale = dot(5); % Scale
% DENDRITIC FIELD
um1 = um(1); % denWidthX
um2 = um(2); % denHeightY
um3 = um(3); % PSD1um
um4 = um(4); % PSD2um
um5 = um(5); % PERI1um
um6 = um(6); % PERI2um
PSD1size = round(um3/Scale);
PSD2size = round(um4/Scale);
periPSD1size = round(um5/Scale);
periPSD2size = round(um6/Scale);

PSD1SAPF = PSD1size*2;
PSD2SAPF = PSD2size*2;
PSA1SAPF = periPSD1size*2;
PSA2SAPF = periPSD2size*2;

PSD1XL = PSD1xv(1);
PSD1YL = PSD1yv(1);

PSD2XL = PSD2xv(1);
PSD2YL = PSD2yv(1);

PSD1XL = PSD1XL-SAPPADPSD1;
PSD1YL = PSD1YL+SAPPADPSD1;
PSD2XL = PSD2XL-SAPPADPSD2;
PSD2YL = PSD2YL+SAPPADPSD2;

PSD1size = PSD1size + (SAPPADPSD1*2);
PSD2size = PSD2size + (SAPPADPSD2*2);



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

%{
for m = 0:(PSD1size-1)
for n = 0:(PSD1size-1)
p = 1+n+(m*PSD1size);
PSD1Sxv = [(PSD1XL+m) (PSD1XL+(1+m)) (PSD1XL+(1+m)) (PSD1XL+m) (PSD1XL+m)]';
PSD1Syv = [(PSD1YL-n) (PSD1YL-n) (PSD1YL-(1+n)) (PSD1YL-(1+n)) (PSD1YL-n)]';
PSD1Sxvecs{p} = PSD1Sxv;
PSD1Syvecs{p} = PSD1Syv;
end
end




for m = 0:(PSD2size-1)
for n = 0:(PSD2size-1)
p = 1+n+(m*PSD2size);
PSD2Sxv = [(PSD2XL+m) (PSD2XL+(1+m)) (PSD2XL+(1+m)) (PSD2XL+m) (PSD2XL+m)]';
PSD2Syv = [(PSD2YL-n) (PSD2YL-n) (PSD2YL-(1+n)) (PSD2YL-(1+n)) (PSD2YL-n)]';
PSD2Sxvecs{p} = PSD2Sxv;
PSD2Syvecs{p} = PSD2Syv;
end
end
%}



end
%===================================%


%===================================%
% SAPMASKREPORT
%===================================%
% This function takes the S1 & S2 Mx for SAP clusters
% and uses a mask to tally the number of SAPs in each cluster
% and returns a condensed Field Mx 1/4th the size, which holds values
% from 0 to 4 depending on the number of SAPs in a defined ROI
%-------------------%
function [SAPFmx1 SAPFmx2 SAPmx1 SAPmx2] = SAPMASKREPORT(Nsteps,...
	S1, S2,um,dot,runSAPPSD1,runSAPPSD2)

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



if runSAPPSD1
	SAPmx1_occ=(S1>0);  
	SAPmx1=convn(SAPmx1_occ,sap_mask,'same');
	SAPFmx1=zeros(SSPSD1);
	for n = 0:(SSPSD1-1)
	for m = 0:(SSPSD1-1)
	SAPFmx1(m+1,n+1) = SAPmx1(1+(m*2),1+(n*2));
	end
	end
end

if runSAPPSD2
	SAPmx2_occ=(S2>0);  
	SAPmx2=convn(SAPmx2_occ,sap_mask,'same');
	SAPFmx2=zeros(SSPSD2);
	for n = 0:(SSPSD2-1)
	for m = 0:(SSPSD2-1)
	SAPFmx2(m+1,n+1) = SAPmx2(1+(m*2),1+(n*2));
	end
	end
end

if ~runSAPPSD1
	SAPmx1_occ=(S1>0);  
	SAPmx1=convn(SAPmx1_occ,sap_mask,'same');
	SAPFmx1=zeros(SSPSD1);
	SAPFmx1 = SAPFmx1.*0;
end

if ~runSAPPSD2
	SAPmx2_occ=(S2>0);  
	SAPmx2=convn(SAPmx2_occ,sap_mask,'same');
	SAPFmx2=zeros(SSPSD2);
	SAPFmx2 = SAPFmx2.*0;
end


%{

if runSAPPSD1
	SAPmx1_occ=(S1>0);  
	SAPmx1=convn(SAPmx1_occ,sap_mask,'same');
	SAPFmx1=zeros(PSD1size);
	for n = 0:(PSD1size-1)
	for m = 0:(PSD1size-1)
	SAPFmx1(m+1,n+1) = SAPmx1(1+(m*2),1+(n*2));
	end
	end
end

if runSAPPSD2
	SAPmx2_occ=(S2>0);  
	SAPmx2=convn(SAPmx2_occ,sap_mask,'same');
	SAPFmx2=zeros(PSD2size);
	for n = 0:(PSD2size-1)
	for m = 0:(PSD2size-1)
	SAPFmx2(m+1,n+1) = SAPmx2(1+(m*2),1+(n*2));
	end
	end
end




SAPSTEP = .1;
Steps = dot(3); % Steps
TimeStep = dot(4)/1000; % TimeStep

um1 = um(1); % denWidthX
um2 = um(2); % denHeightY

um5 = um(5); % PERI1um
um6 = um(6); % PERI2um
fsizeX = round(um1/Scale);
fsizeY = round(um2/Scale);

periPSD1size = round(um5/Scale);
periPSD2size = round(um6/Scale);

PSD1SAPF = PSD1size*2;
PSD2SAPF = PSD2size*2;
PSA1SAPF = periPSD1size*2;
PSA2SAPF = periPSD2size*2;


%=======================================%
%				GLUR1 PSD				%
%=======================================%
if doKo(1) %doKo(1); % useGluR1slots
	
GluR1Sd1 = ones(size(GluR1xyl));
GluR1Sd2 = GluR1Sd1;
[m,p,n] = size(PSD1Sxvecs);

for m = 1:p
	
PSD1SxvNOW = PSD1Sxvecs{m};
PSD1SyvNOW = PSD1Syvecs{m};
PSD2SxvNOW = PSD2Sxvecs{m};
PSD2SyvNOW = PSD1SyvNOW.*-1;

G1INPSD1 = inpolygon(GluR1xyl(1,:)',GluR1xyl(2,:)',PSD1SxvNOW,PSD1SyvNOW);
G1INPSD2 = inpolygon(GluR1xyl(1,:)',GluR1xyl(2,:)',PSD2SxvNOW,PSD2SyvNOW);

GluR1Sd1(:,G1INPSD1) = GluR1Sd1(:,G1INPSD1)+(SAPFmx1(m));
GluR1Sd2(:,G1INPSD2) = GluR1Sd2(:,G1INPSD2)+(SAPFmx2(m));

end

GluR1SdP1 = GluR1Sd1-1;
GluR1SdP2 = GluR1Sd2-1;

end %if %doKo(1); % useGluR1slots
%=======================================%

%=======================================%
%				GLUR2 PSD				%
%=======================================%
if doKo(2) %doKo(2); % useGluR2slots

GluR2Sd1 = ones(size(GluR2xyl));
GluR2Sd2 = GluR2Sd1;
[m,p,n] = size(PSD1Sxvecs);

for m = 1:p
	
PSD1SxvNOW = PSD1Sxvecs{m};
PSD1SyvNOW = PSD1Syvecs{m};
PSD2SxvNOW = PSD2Sxvecs{m};
PSD2SyvNOW = PSD1SyvNOW.*-1;

G2INPSD1 = inpolygon(GluR2xyl(1,:)',GluR2xyl(2,:)',PSD1SxvNOW,PSD1SyvNOW);
G2INPSD2 = inpolygon(GluR2xyl(1,:)',GluR2xyl(2,:)',PSD2SxvNOW,PSD2SyvNOW);

GluR2Sd1(:,G2INPSD1) = GluR2Sd1(:,G2INPSD1)+(SAPFmx1(m));
GluR2Sd2(:,G2INPSD2) = GluR2Sd2(:,G2INPSD2)+(SAPFmx2(m));

end

GluR2SdP1 = GluR2Sd1-1;
GluR2SdP2 = GluR2Sd2-1;

end %doKo(2); % useGluR2slots
%=======================================%

%}
end
%===================================%


%===================================%
%		SAPPREALLOCATE
%===================================%
function [sap_mask SSPSD1 SSPSD2 GluR1SdP1 GluR1SdP2...
    GluR1SLOCMXa GluR1SLOCMXb g1polyN1 g2polyN1...
    GluR2SdP1 GluR2SdP2 GluR2SLOCMXa GluR2SLOCMXb] = SAPPREALLOCATE(...
    S1,S2,um,dot,runSAPPSD1,runSAPPSD2,PSD1Sxvecs,...
    GluR1xyl,GluR2xyl)


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
GluR1SdP1 = zeros(size(GluR1xyl));
GluR1SdP2 = GluR1SdP1;
 
GluR1SLOCMXa = zeros(1,size(GluR1SdP1,2));
GluR1SLOCMXb = GluR1SLOCMXa;

[m,g1polyN1,n] = size(PSD1Sxvecs);
    % currentcell = zeros(1,size(GluR1xyl,2));


%-------------------
% G2SAPPOLYGON
%-------------------
GluR2SdP1 = zeros(size(GluR2xyl));
GluR2SdP2 = GluR2SdP1;
 
GluR2SLOCMXa = zeros(1,size(GluR2SdP1,2));
GluR2SLOCMXb = GluR2SLOCMXa;
 
[m,g2polyN1,n] = size(PSD1Sxvecs);
    % currentcell = zeros(1,size(GluR2xyl,2));



end


%-----------------------------------------%
%			LOOPED SAP FUNCTIONS
%-----------------------------------------%

%===================================%
%		SUPERSLOT
%===================================%
function [GluR1xyl GluR2xyl GluR1xyds GluR2xyds...
	SAPFmx1 SAPFmx2 SAPmx1 SAPmx2 G1FSLOTS G2FSLOTS...
	GluR1_TdwellSPYN GluR1_TdwellPSD GluR1_TdwellPERI...
	GluR2_TdwellSPYN GluR2_TdwellPSD GluR2_TdwellPERI...
	G1INSPYN1 G1INSPYN2 G1INPSD1 G1INPSD2 G1INPERI1 G1INPERI2...
	G2INSPYN1 G2INSPYN2 G2INPSD1 G2INPSD2 G2INPERI1 G2INPERI2...
	G1P1SAPM G1P2SAPM G2P1SAPM G2P2SAPM]...
	= SUPERSLOT(stepN,t,GluR1xyl,GluR2xyl,GluR1xyds,GluR2xyds,...
	GluR1_TdwellPSD,GluR1_TdwellPERI,GluR1_TdwellSPYN,LsGR1psa,LsGR1psd,...
	GluR2_TdwellPSD,GluR2_TdwellPERI,GluR2_TdwellSPYN,LsGR2psa,LsGR2psd,...	
	S1,S2,runSAPPSD1,runSAPPSD2,...
	PSD1Sxvecs,PSD1Syvecs,PSD2Sxvecs,PSD2Syvecs,...
	SAPFmx1,SAPFmx2,...
	sap_mask,SSPSD1,SSPSD2,GluR1SdP1,GluR1SdP2,...
	GluR1SLOCMXa,GluR1SLOCMXb,g1polyN1,g2polyN1,...
	GluR2SdP1,GluR2SdP2,GluR2SLOCMXa,GluR2SLOCMXb,...
	LTP1onG1,LTP1offG1,LTP2onG1,LTP2offG1,...
	LTP1onG2,LTP1offG2,LTP2onG2,LTP2offG2,...
	XYLBp1,XYRTp1,XYLBp2,XYRTp2,...
	G1STBASE,G1RTBASE,G1STLTP,G1RTLTP,G1BSMu,G1LSMu,...
	G2STBASE,G2RTBASE,G2STLTP,G2RTLTP,G2BSMu,G2LSMu)





%-------------------
% SAPMASKREPORT
%-------------------
%{
% SAPmx1 = 16x16 convolution Mx (X checks X,E,S,SE)
% NW N NE
%  W X E
% SW S SE

% SAPFmx1 = 8x8 condensed Mx (always 1/4 of SAPmx1)
% takes all odd row and column values
% keeps "odd" discards odd+1 odd+16 odd+17
% odd   odd+16
% odd+1 odd+17

dumMx = [1:16]'
for mx = 1:15
nx = mx+1;
dumMx(:,nx) = dumMx(:,mx)+16
end
dumMx2 = permute(dumMx,[2 1])
dumMx3 = reshape(dumMx2,16,[])

INPOLYGON REFERENCE
ODDID = ID+(ID-1)
ODDID+1
ODDID+16
ODDID+17


S1(2:5,2:5) = 1
if runSAPPSD1
	SAPmx1_occ=(S1>0);  
	SAPmx1=convn(SAPmx1_occ,sap_mask,'same');
	SAPFmx1=zeros(SSPSD1);
	for n = 0:(SSPSD1-1)
	for m = 0:(SSPSD1-1)
	SAPFmx1(m+1,n+1) = SAPmx1(1+(m*2),1+(n*2));
	end
	end
end

%}

if runSAPPSD1
	SAPmx1_occ=(S1>0);  
	SAPmx1=convn(SAPmx1_occ,sap_mask,'same');
	SAPFmx1=zeros(SSPSD1);
	for n = 0:(SSPSD1-1)
	for m = 0:(SSPSD1-1)
	SAPFmx1(m+1,n+1) = SAPmx1(1+(m*2),1+(n*2));
	end
	end
end

if runSAPPSD2
	SAPmx2_occ=(S2>0);  
	SAPmx2=convn(SAPmx2_occ,sap_mask,'same');
	SAPFmx2=zeros(SSPSD2);
	for n = 0:(SSPSD2-1)
	for m = 0:(SSPSD2-1)
	SAPFmx2(m+1,n+1) = SAPmx2(1+(m*2),1+(n*2));
	end
	end
end

if ~runSAPPSD1
	SAPmx1_occ=(S1>0);  
	SAPmx1=convn(SAPmx1_occ,sap_mask,'same');
	SAPFmx1=zeros(SSPSD1);
	SAPFmx1 = SAPFmx1.*0;
end

if ~runSAPPSD2
	SAPmx2_occ=(S2>0);  
	SAPmx2=convn(SAPmx2_occ,sap_mask,'same');
	SAPFmx2=zeros(SSPSD2);
	SAPFmx2 = SAPFmx2.*0;
end


%-------------------
% G1SAPPOLYGON
%-------------------
%{

SAPFmx1 & SAPFmx2 = 8x8 Mx of #SAPs in each 64 SLOTs of SPYN1 & SPYN2
this G1SAPPOLYGON subroutine takes an 1x64 Cell of 1x1 polygons

G1INPSD1 holds logic array testing if a GluR in a polygon (on each iteration)

GluR1SdP1 takes G1INPSD1 logic test and SAPFmx1 and stores the #SAP
	each GluR is currently surrounded by. cycles down each column

GluR1SLOCMXa stores the ID of the polygon each GluR is located in
The polygon check starts with the top left and cycles down each column
like so...

%----Test Polygon Check Order---%
fig2 = figure(2);
currentcell = zeros(1,g1polyN1);
for m = 1:g1polyN1
currentcell = currentcell+1;
PSD1SxvNOW = PSD1Sxvecs{m}; % CURRENT 1x1 POLYGON (NOT JUST PSD)
PSD1SyvNOW = PSD1Syvecs{m};
PSD2SxvNOW = PSD2Sxvecs{m};
PSD2SyvNOW = PSD2Syvecs{m};
G1INPSD1 = inpolygon(GluR1xyl(1,:)',GluR1xyl(2,:)',PSD1SxvNOW,PSD1SyvNOW);
G1INPSD2 = inpolygon(GluR1xyl(1,:)',GluR1xyl(2,:)',PSD2SxvNOW,PSD2SyvNOW);
GluR1SdP1(:,G1INPSD1) = GluR1SdP1(:,G1INPSD1)+(SAPFmx1(m));
GluR1SdP2(:,G1INPSD2) = GluR1SdP2(:,G1INPSD2)+(SAPFmx2(m));
GluR1SLOCMXa(1,G1INPSD1) = GluR1SLOCMXa(:,G1INPSD1)+(currentcell(m));
GluR1SLOCMXb(1,G1INPSD2) = GluR1SLOCMXb(:,G1INPSD2)+(currentcell(m));
figure(fig2)
axesm miller
plotm(PSD1SyvNOW,PSD1SxvNOW,'r')
% draw now
hold on
pause(.1)
end
%----Test Polygon Check Order---%




%----Test Polygon Check Order---%
fig2 = figure(2);
currentcell = zeros(1,g1polyN1);
for m = 1:g1polyN1
currentcell = currentcell+1;
PSD1SxvNOW = PSD1Sxvecs{m}; % CURRENT 1x1 POLYGON (NOT JUST PSD)
PSD1SyvNOW = PSD1Syvecs{m};
PSD2SxvNOW = PSD2Sxvecs{m};
PSD2SyvNOW = PSD2Syvecs{m};

XYLB1 = [PSD1SxvNOW(4) PSD1SyvNOW(4)];
XYRT1 = [PSD1SxvNOW(2) PSD1SyvNOW(2)];
XYLB2 = [PSD2SxvNOW(4) PSD2SyvNOW(4)];
XYRT2 = [PSD2SxvNOW(2) PSD2SyvNOW(2)];

G1INPSD1 = inboxfun(XYLB1,XYRT1,GluR1xyl);
G1INPSD2 = inboxfun(XYLB2,XYRT2,GluR1xyl);

GluR1SdP1(:,G1INPSD1) = GluR1SdP1(:,G1INPSD1)+(SAPFmx1(m));
GluR1SdP2(:,G1INPSD2) = GluR1SdP2(:,G1INPSD2)+(SAPFmx2(m));
GluR1SLOCMXa(1,G1INPSD1) = GluR1SLOCMXa(:,G1INPSD1)+(currentcell(m));
GluR1SLOCMXb(1,G1INPSD2) = GluR1SLOCMXb(:,G1INPSD2)+(currentcell(m));
figure(fig2)
axesm miller
plotm(PSD1SyvNOW,PSD1SxvNOW,'r')
% draw now
hold on
pause(.1)
end
%----Test Polygon Check Order---%
%}



currentcell = zeros(1,g1polyN1);
for m = 1:g1polyN1
currentcell = currentcell+1;

PSD1SxvNOW = PSD1Sxvecs{m}; % CURRENT 1x1 POLYGON (NOT JUST PSD)
PSD1SyvNOW = PSD1Syvecs{m};
PSD2SxvNOW = PSD2Sxvecs{m};
PSD2SyvNOW = PSD2Syvecs{m};

XYLB1 = [PSD1SxvNOW(4) PSD1SyvNOW(4)]; % MAKE SEPARATE LOOPS FOR PSD1/PSD2
XYRT1 = [PSD1SxvNOW(2) PSD1SyvNOW(2)];
XYLB2 = [PSD2SxvNOW(4) PSD2SyvNOW(4)];
XYRT2 = [PSD2SxvNOW(2) PSD2SyvNOW(2)];
G1IN1 = inboxfun(XYLB1,XYRT1,GluR1xyl)';
G1IN2 = inboxfun(XYLB2,XYRT2,GluR1xyl)';

% G1INPSD1 = inpolygon(GluR1xyl(1,:)',GluR1xyl(2,:)',PSD1SxvNOW,PSD1SyvNOW);
% G1INPSD2 = inpolygon(GluR1xyl(1,:)',GluR1xyl(2,:)',PSD2SxvNOW,PSD2SyvNOW);

GluR1SdP1(:,G1IN1) = GluR1SdP1(:,G1IN1)+(SAPFmx1(m));
GluR1SdP2(:,G1IN2) = GluR1SdP2(:,G1IN2)+(SAPFmx2(m));
GluR1SLOCMXa(1,G1IN1) = GluR1SLOCMXa(:,G1IN1)+(currentcell(m));
GluR1SLOCMXb(1,G1IN2) = GluR1SLOCMXb(:,G1IN2)+(currentcell(m));
end

GluR1INSPYN1 = GluR1SLOCMXa>0;
GluR1INSPYN2 = GluR1SLOCMXb>0;

GluR1SP1 = GluR1SdP1(1,:); %!! GluR1SdP1 doesn't need to by 2 rows
GluR1SP2 = GluR1SdP2(2,:);
G1SAPLOC1 = GluR1SLOCMXa;
G1SAPLOC2 = GluR1SLOCMXb;

% GLUR-ID , POLY ID , SAPS NEAR GLUR , GLUR DWELL TIME
G1P1SAPMX = [(1:numel(GluR1SP1))' G1SAPLOC1' GluR1SP1' GluR1_TdwellSPYN(1,:)'];
G1P2SAPMX = [([1:numel(GluR1SP2)]') G1SAPLOC2' GluR1SP2' GluR1_TdwellSPYN(2,:)'];
G1P1SAPMX = sortrows(G1P1SAPMX,-4);
G1P2SAPMX = sortrows(G1P2SAPMX,-4);


%-------------------
% G2SAPPOLYGON
%-------------------
currentcell = zeros(1,g2polyN1);	
for m = 1:g2polyN1
currentcell = currentcell+1;
    
PSD1SxvNOW = PSD1Sxvecs{m};
PSD1SyvNOW = PSD1Syvecs{m};
PSD2SxvNOW = PSD2Sxvecs{m};
PSD2SyvNOW = PSD2Syvecs{m};

XYLB1 = [PSD1SxvNOW(4) PSD1SyvNOW(4)];
XYRT1 = [PSD1SxvNOW(2) PSD1SyvNOW(2)];
XYLB2 = [PSD2SxvNOW(4) PSD2SyvNOW(4)];
XYRT2 = [PSD2SxvNOW(2) PSD2SyvNOW(2)];
G2IN1 = inboxfun(XYLB1,XYRT1,GluR2xyl)';
G2IN2 = inboxfun(XYLB2,XYRT2,GluR2xyl)';
 
% G2INPSD1 = inpolygon(GluR2xyl(1,:)',GluR2xyl(2,:)',PSD1SxvNOW,PSD1SyvNOW);
% G2INPSD2 = inpolygon(GluR2xyl(1,:)',GluR2xyl(2,:)',PSD2SxvNOW,PSD2SyvNOW);
 
GluR2SdP1(:,G2IN1) = GluR2SdP1(:,G2IN1)+(SAPFmx1(m));
GluR2SdP2(:,G2IN2) = GluR2SdP2(:,G2IN2)+(SAPFmx2(m));
 
GluR2SLOCMXa(1,G2IN1) = GluR2SLOCMXa(:,G2IN1)+(currentcell(m));
GluR2SLOCMXb(1,G2IN2) = GluR2SLOCMXb(:,G2IN2)+(currentcell(m));
end

GluR2INSPYN1 = GluR2SLOCMXa>0;
GluR2INSPYN2 = GluR2SLOCMXb>0;

GluR2SP1 = GluR2SdP1(1,:); %!! GluR2SdP1 doesn't need to by 2 rows
GluR2SP2 = GluR2SdP2(2,:);
G2SAPLOC1 = GluR2SLOCMXa;
G2SAPLOC2 = GluR2SLOCMXb;

% GLUR-ID , POLY ID , SAPS NEAR GLUR , GLUR DWELL TIME
G2P1SAPMX = [([1:numel(GluR2SP1)]') G2SAPLOC1' GluR2SP1' GluR2_TdwellSPYN(1,:)'];
G2P2SAPMX = [([1:numel(GluR2SP2)]') G2SAPLOC2' GluR2SP2' GluR2_TdwellSPYN(2,:)'];
G2P1SAPMX = sortrows(G2P1SAPMX,-4);
G2P2SAPMX = sortrows(G2P2SAPMX,-4);

%===================================%
%-------------------%
% GLUR SLOT TEST
%-------------------%
G1CDF = G1STBASE-G1BSMu;
G2CDF = G2STBASE-G2BSMu;
if stepN >= LTP1onG1 && stepN <= LTP1offG1; 
		G1CDF=G1STLTP-G1LSMu;
end
if stepN >= LTP2onG1 && stepN <= LTP2offG1; 
		G1CDF=G1STLTP-G1LSMu; 
end
if stepN >= LTP1onG2 && stepN <= LTP1offG2; 
		G2CDF=G2STLTP-G2LSMu;
end
if stepN >= LTP2onG2 && stepN <= LTP2offG2; 
		G2CDF=G2STLTP-G2LSMu; 
end

[~, uniG1P1SAPMX, ~] = unique(G1P1SAPMX(:,2),'rows','first');
[~, uniG1P2SAPMX, ~] = unique(G1P2SAPMX(:,2),'rows','first');
[~, uniG2P1SAPMX, ~] = unique(G2P1SAPMX(:,2),'rows','first');
[~, uniG2P2SAPMX, ~] = unique(G2P2SAPMX(:,2),'rows','first');
G1P1SAPMx = G1P1SAPMX(uniG1P1SAPMX,:);
G1P2SAPMx = G1P2SAPMX(uniG1P2SAPMX,:);
G2P1SAPMx = G2P1SAPMX(uniG2P1SAPMX,:);
G2P2SAPMx = G2P2SAPMX(uniG2P2SAPMX,:);

G1P1SAPMxP = G1P1SAPMx;
PtGR1P1 = randn(numel(G1P1SAPMx(:,1)),1)+(G1P1SAPMxP(:,3)+G1CDF);
G1P1SAPMxP(:,5) = PtGR1P1;
G1P1SAPMxPx = G1P1SAPMxP((G1P1SAPMxP(:,5)>0),1);

G1P2SAPMxP = G1P2SAPMx;
PtGR1P2 = randn(numel(G1P2SAPMx(:,1)),1)+(G1P2SAPMxP(:,3)+G1CDF);
G1P2SAPMxP(:,5) = PtGR1P2;
G1P2SAPMxPx = G1P2SAPMxP((G1P2SAPMxP(:,5)>0),1);

G2P1SAPMxP = G2P1SAPMx;
PtGR2P1 = randn(numel(G2P1SAPMx(:,1)),1)+(G2P1SAPMxP(:,3)+G2CDF);
G2P1SAPMxP(:,5) = PtGR2P1;
G2P1SAPMxPx = G2P1SAPMxP((G2P1SAPMxP(:,5)>0),1);

G2P2SAPMxP = G2P2SAPMx;
PtGR2P2 = randn(numel(G2P2SAPMx(:,1)),1)+(G2P2SAPMxP(:,3)+G2CDF);
G2P2SAPMxP(:,5) = PtGR2P2;
G2P2SAPMxPx = G2P2SAPMxP((G2P2SAPMxP(:,5)>0),1);


G1P1SAPM = sortrows(G1P1SAPMX,1);
G1P2SAPM = sortrows(G1P2SAPMX,1);
G2P1SAPM = sortrows(G2P1SAPMX,1);
G2P2SAPM = sortrows(G2P2SAPMX,1);


%===================================%
% if stepN == 300; keyboard; end
%===================================%
% FULL IN-SLOT DATA: GLUR-ID , POLY-ID , SAPS , DWELL
G1P1SAPM = G1P1SAPM(G1P1SAPMxPx,:);
G1P2SAPM = G1P2SAPM(G1P2SAPMxPx,:);
G2P1SAPM = G2P1SAPM(G2P1SAPMxPx,:);
G2P2SAPM = G2P2SAPM(G2P2SAPMxPx,:);
%===================================%


[G1P1dw,~,~] = find(G1P1SAPM(:,4) <= 3);
G1P1SAPM(G1P1dw,:) = [];
[G1P2dw,~,~] = find(G1P2SAPM(:,4) <= 3);
G1P2SAPM(G1P2dw,:) = [];
[G2P1dw,~,~] = find(G2P1SAPM(:,4) <= 3);
G2P1SAPM(G2P1dw,:) = [];
[G2P2dw,~,~] = find(G2P2SAPM(:,4) <= 3);
G2P2SAPM(G2P2dw,:) = [];


% GLUR-ID , POLY ID , SAPS NEAR GLUR , GLUR DWELL TIME
G1P1xxx = [G1P1SAPM [(11+zeros(numel(G1P1SAPM(:,1)),1))]];
G1P2xxx = [G1P2SAPM [(12+zeros(numel(G1P2SAPM(:,1)),1))]];
G2P1xxx = [G2P1SAPM [(21+zeros(numel(G2P1SAPM(:,1)),1))]];
G2P2xxx = [G2P2SAPM [(22+zeros(numel(G2P2SAPM(:,1)),1))]];

G1G2P1 = [G1P1xxx; G2P1xxx];
G1G2P2 = [G1P2xxx; G2P2xxx];

G1G2P1xMX = sortrows(G1G2P1,-4);
G1G2P2xMX = sortrows(G1G2P2,-4);
[~, uniG1G2P1, ~] = unique(G1G2P1xMX(:,2),'rows','first');
[~, uniG1G2P2, ~] = unique(G1G2P2xMX(:,2),'rows','first');
G1G2P1xxMX = G1G2P1(uniG1G2P1,:);
G1G2P2xxMX = G1G2P2(uniG1G2P2,:);

G1P1xxID = G1G2P1xxMX((G1G2P1xxMX(:,5)<15),1);
G2P1xxID = G1G2P1xxMX((G1G2P1xxMX(:,5)>15),1);
G1P2xxID = G1G2P2xxMX((G1G2P2xxMX(:,5)<15),1);
G2P2xxID = G1G2P2xxMX((G1G2P2xxMX(:,5)>15),1);


G1P1SLOTTED = G1P1xxID;
G1P2SLOTTED = G1P2xxID;
G2P1SLOTTED = G2P1xxID;
G2P2SLOTTED = G2P2xxID;


% G1P1SLOTTED = G1P1SAPM(:,1);
% G1P2SLOTTED = G1P2SAPM(:,1);
% G2P1SLOTTED = G2P1SAPM(:,1);
% G2P2SLOTTED = G2P2SAPM(:,1);




%===================================%
GluR1xyds(:,G1P1SLOTTED) = GluR1xyds(:,G1P1SLOTTED)*(.001);
GluR1xyds(:,G1P2SLOTTED) = GluR1xyds(:,G1P2SLOTTED)*(.001);
GluR2xyds(:,G2P1SLOTTED) = GluR2xyds(:,G2P1SLOTTED)*(.001);
GluR2xyds(:,G2P2SLOTTED) = GluR2xyds(:,G2P2SLOTTED)*(.001);
% HOW MANY SLOTS ARE NOW FILLED?
G1PSD1FSLOTS = numel(G1P1SLOTTED);
G1PSD2FSLOTS = numel(G1P2SLOTTED);
G1FSLOTS = [G1PSD1FSLOTS G1PSD2FSLOTS 0 0];
G2PSD1FSLOTS = numel(G2P1SLOTTED);
G2PSD2FSLOTS = numel(G2P2SLOTTED);
G2FSLOTS = [G2PSD1FSLOTS G2PSD2FSLOTS 0 0];




%===================================%
% DWELL TIME COUNTER
%===================================%
% IF INPOLYGON GETS TOO SLOW, THIS MAY REPLACE THE MOVEGLUR BELOW
% G1P1Tdwell = zeros(1,numel(GluR1xyl(1,:)))';
% G1P1INaPOLY = G1P1SAPMX(:,2)>0;
% G1P1Tdwell(G1P1INaPOLY) = G1P1Tdwell(G1P1INaPOLY)+t;
% G1P1Tdwell(1,~G1P1INaPOLY) = 0;


%-------------------
% MOVEGLUR1
%-------------------
%===================%
% SPINE %
%-------%
% G1INSPYN1 = inboxfun(XYLBpr1,XYRTpr1,GluR1xyl);
% G1INSPYN2 = inboxfun(XYLBpr2,XYRTpr2,GluR1xyl);

G1INSPYN1 = GluR1INSPYN1;
G1INSPYN2 = GluR1INSPYN2;

GluR1_TdwellSPYN(1,G1INSPYN1) = GluR1_TdwellSPYN(1,G1INSPYN1)+t;
GluR1_TdwellSPYN(1,~G1INSPYN1) = 0;
GluR1_TdwellSPYN(2,G1INSPYN2) = GluR1_TdwellSPYN(2,G1INSPYN2)+t;
GluR1_TdwellSPYN(2,~G1INSPYN2) = 0;

%===================%
% PSD   %
%-------%
G1INPSD1 = inboxfun(XYLBp1,XYRTp1,GluR1xyl);
G1INPSD2 = inboxfun(XYLBp2,XYRTp2,GluR1xyl);
 
GluR1xyds(:,G1INPSD1) = GluR1xyds(:,G1INPSD1)*(LsGR1psd);
GluR1xyds(:,G1INPSD2) = GluR1xyds(:,G1INPSD2)*(LsGR1psd);
 
GluR1_TdwellPSD(1,G1INPSD1) = GluR1_TdwellPSD(1,G1INPSD1)+t;
GluR1_TdwellPSD(1,~G1INPSD1) = 0;
GluR1_TdwellPSD(2,G1INPSD2) = GluR1_TdwellPSD(2,G1INPSD2)+t;
GluR1_TdwellPSD(2,~G1INPSD2) = 0;

%===================%
% PERI  %
%-------%
G1INPERI1 = G1INSPYN1 > G1INPSD1;
G1INPERI2 = G1INSPYN2 > G1INPSD2;

GluR1xyds(:,G1INPERI1) = GluR1xyds(:,G1INPERI1)*(LsGR1psa);
GluR1xyds(:,G1INPERI2) = GluR1xyds(:,G1INPERI2)*(LsGR1psa);

GluR1_TdwellPERI(1,G1INPERI1) = GluR1_TdwellPERI(1,G1INPERI1)+t;
GluR1_TdwellPERI(1,~G1INPERI1) = 0;
GluR1_TdwellPERI(2,G1INPERI2) = GluR1_TdwellPERI(2,G1INPERI2)+t;
GluR1_TdwellPERI(2,~G1INPERI2) = 0;
%===================%



%-------------------
% MOVEGLUR2
%-------------------
%===================%
% SPINE %
%-------%
% G2INSPYN1 = inboxfun(XYLBpr1,XYRTpr1,GluR2xyl);
% G2INSPYN2 = inboxfun(XYLBpr2,XYRTpr2,GluR2xyl);

G2INSPYN1 = GluR2INSPYN1;
G2INSPYN2 = GluR2INSPYN2;

GluR2_TdwellSPYN(1,G2INSPYN1) = GluR2_TdwellSPYN(1,G2INSPYN1)+t;
GluR2_TdwellSPYN(1,~G2INSPYN1) = 0;
GluR2_TdwellSPYN(2,G2INSPYN2) = GluR2_TdwellSPYN(2,G2INSPYN2)+t;
GluR2_TdwellSPYN(2,~G2INSPYN2) = 0;

%===================%
% PSD   %
%-------%
G2INPSD1 = inboxfun(XYLBp1,XYRTp1,GluR2xyl);
G2INPSD2 = inboxfun(XYLBp2,XYRTp2,GluR2xyl);

GluR2xyds(:,G2INPSD1) = GluR2xyds(:,G2INPSD1)*(LsGR2psd);
GluR2xyds(:,G2INPSD2) = GluR2xyds(:,G2INPSD2)*(LsGR2psd);

GluR2_TdwellPSD(1,G2INPSD1) = GluR2_TdwellPSD(1,G2INPSD1)+t;
GluR2_TdwellPSD(1,~G2INPSD1) = 0;
GluR2_TdwellPSD(2,G2INPSD2) = GluR2_TdwellPSD(2,G2INPSD2)+t;
GluR2_TdwellPSD(2,~G2INPSD2) = 0;

%===================%
% PERI  %
%-------%
G2INPERI1 = G2INSPYN1 > G2INPSD1;
G2INPERI2 = G2INSPYN2 > G2INPSD2;

GluR2xyds(:,G2INPERI1) = GluR2xyds(:,G2INPERI1)*(LsGR2psa);
GluR2xyds(:,G2INPERI2) = GluR2xyds(:,G2INPERI2)*(LsGR2psa);

GluR2_TdwellPERI(1,G2INPERI1) = GluR2_TdwellPERI(1,G2INPERI1)+t;
GluR2_TdwellPERI(1,~G2INPERI1) = 0;
GluR2_TdwellPERI(2,G2INPERI2) = GluR2_TdwellPERI(2,G2INPERI2)+t;
GluR2_TdwellPERI(2,~G2INPERI2) = 0;
%===================%
%===================================%
% DWELL TIME COUNTER
%===================================%



%{
%-------------------
% MOVEGLUR1
%-------------------
G1INSPYN1 = inpolygon(GluR1xyl(1,:)',GluR1xyl(2,:)',SPYN1xv,SPYN1yv);
G1INSPYN2 = inpolygon(GluR1xyl(1,:)',GluR1xyl(2,:)',SPYN2xv,SPYN2yv);

GluR1_TdwellSPYN(1,G1INSPYN1) = GluR1_TdwellSPYN(1,G1INSPYN1)+t;
GluR1_TdwellSPYN(1,~G1INSPYN1) = 0;
GluR1_TdwellSPYN(2,G1INSPYN2) = GluR1_TdwellSPYN(2,G1INSPYN2)+t;
GluR1_TdwellSPYN(2,~G1INSPYN2) = 0;

%===================%
% PSD   %
%-------%
G1INPSD1 = inpolygon(GluR1xyl(1,:)',GluR1xyl(2,:)',PSD1xv,PSD1yv);
G1INPSD2 = inpolygon(GluR1xyl(1,:)',GluR1xyl(2,:)',PSD2xv,PSD2yv);
 
GluR1xyds(:,G1INPSD1) = GluR1xyds(:,G1INPSD1)*(LsGR1psd);
GluR1xyds(:,G1INPSD2) = GluR1xyds(:,G1INPSD2)*(LsGR1psd);
 
GluR1_TdwellPSD(1,G1INPSD1) = GluR1_TdwellPSD(1,G1INPSD1)+t;
GluR1_TdwellPSD(1,~G1INPSD1) = 0;
GluR1_TdwellPSD(2,G1INPSD2) = GluR1_TdwellPSD(2,G1INPSD2)+t;
GluR1_TdwellPSD(2,~G1INPSD2) = 0;

%===================%
% PERI  %
%-------%
G1INPERI1 = inpolygon(GluR1xyl(1,:)',GluR1xyl(2,:)',PERI1xv,PERI1yv);
G1INPERI2 = inpolygon(GluR1xyl(1,:)',GluR1xyl(2,:)',PERI2xv,PERI2yv);
 
GluR1xyds(:,G1INPERI1) = GluR1xyds(:,G1INPERI1)*(LsGR1psa);
GluR1xyds(:,G1INPERI2) = GluR1xyds(:,G1INPERI2)*(LsGR1psa);

GluR1_TdwellPERI(1,G1INPERI1) = GluR1_TdwellPERI(1,G1INPERI1)+t;
GluR1_TdwellPERI(1,~G1INPERI1) = 0;
GluR1_TdwellPERI(2,G1INPERI2) = GluR1_TdwellPERI(2,G1INPERI2)+t;
GluR1_TdwellPERI(2,~G1INPERI2) = 0;
%===================%



%-------------------
% MOVEGLUR2
%-------------------
%===================%
% SPINE %
%-------%
G2INSPYN1 = inpolygon(GluR2xyl(1,:)',GluR2xyl(2,:)',SPYN1xv,SPYN1yv);
G2INSPYN2 = inpolygon(GluR2xyl(1,:)',GluR2xyl(2,:)',SPYN2xv,SPYN2yv);

GluR2_TdwellSPYN(1,G2INSPYN1) = GluR2_TdwellSPYN(1,G2INSPYN1)+t;
GluR2_TdwellSPYN(1,~G2INSPYN1) = 0;
GluR2_TdwellSPYN(2,G2INSPYN2) = GluR2_TdwellSPYN(2,G2INSPYN2)+t;
GluR2_TdwellSPYN(2,~G2INSPYN2) = 0;

%===================%
% PSD   %
%-------%
G2INPSD1 = inpolygon(GluR2xyl(1,:)',GluR2xyl(2,:)',PSD1xv,PSD1yv);
G2INPSD2 = inpolygon(GluR2xyl(1,:)',GluR2xyl(2,:)',PSD2xv,PSD2yv);

GluR2xyds(:,G2INPSD1) = GluR2xyds(:,G2INPSD1)*(LsGR2psd);
GluR2xyds(:,G2INPSD2) = GluR2xyds(:,G2INPSD2)*(LsGR2psd);

GluR2_TdwellPSD(1,G2INPSD1) = GluR2_TdwellPSD(1,G2INPSD1)+t;
GluR2_TdwellPSD(1,~G2INPSD1) = 0;
GluR2_TdwellPSD(2,G2INPSD2) = GluR2_TdwellPSD(2,G2INPSD2)+t;
GluR2_TdwellPSD(2,~G2INPSD2) = 0;

%===================%
% PERI  %
%-------%
G2INPERI1 = inpolygon(GluR2xyl(1,:)',GluR2xyl(2,:)',PERI1xv,PERI1yv);
G2INPERI2 = inpolygon(GluR2xyl(1,:)',GluR2xyl(2,:)',PERI2xv,PERI2yv);

GluR2xyds(:,G2INPERI1) = GluR2xyds(:,G2INPERI1)*(LsGR2psa);
GluR2xyds(:,G2INPERI2) = GluR2xyds(:,G2INPERI2)*(LsGR2psa);

GluR2_TdwellPERI(1,G2INPERI1) = GluR2_TdwellPERI(1,G2INPERI1)+t;
GluR2_TdwellPERI(1,~G2INPERI1) = 0;
GluR2_TdwellPERI(2,G2INPERI2) = GluR2_TdwellPERI(2,G2INPERI2)+t;
GluR2_TdwellPERI(2,~G2INPERI2) = 0;
%===================%
%===================================%
% DWELL TIME COUNTER
%===================================%
%}

% if stepN == 300; keyboard; end
end



%===================================%
% S_MainClusterFun
%===================================%
% For each PSD area independently update the current
% S1 & S2  SAP clusters and return this S-matrix 
% (these ioMx represent the current surface SAPs) 
% The raw matrices will be plotted using a function further below
%-------------------%
% S1_MainClusterFun
%-------------------%
function [S1 G1SP1] = S1_MainClusterFun(...
	h_mask, S1L1, S1beta,... 
    S1mu, S1delta_t, S1ro, S1r, S1, nn,...
	G1P1SAPM,G1P2SAPM,G2P1SAPM,G2P2SAPM,...
	G1SP1,G2SP1,G1LTParray,G2LTParray,...
	G1STBASE,G1RTBASE,G1STLTP,G1RTLTP,G1BSMu,G1LSMu,...
	G2STBASE,G2RTBASE,G2STLTP,G2RTLTP,G2BSMu,G2LSMu,...
	LTP1onG1,LTP1offG1,LTP2onG1,LTP2offG1)


%-------------------------------%
G1P1poly = G1P1SAPM(:,2);
[~, ~, G1P1poly] = find(G1P1poly);
for nx = 1:(numel(G1P1poly))
	if G1P1poly(nx)<=8 && G1P1poly(nx)>0
		G1P1poly(nx) = (G1P1poly(nx)*2)+(16*0)-1;
	end
	if G1P1poly(nx)<=16 && G1P1poly(nx)>8
		G1P1poly(nx) = (G1P1poly(nx)*2)+(16*1)-1;
	end
	if G1P1poly(nx)<=24 && G1P1poly(nx)>16
		G1P1poly(nx) = (G1P1poly(nx)*2)+(16*2)-1;
	end
	if G1P1poly(nx)<=32 && G1P1poly(nx)>24
		G1P1poly(nx) = (G1P1poly(nx)*2)+(16*3)-1;
	end
	if G1P1poly(nx)<=40 && G1P1poly(nx)>32
		G1P1poly(nx) = (G1P1poly(nx)*2)+(16*4)-1;
	end
	if G1P1poly(nx)<=48 && G1P1poly(nx)>40
		G1P1poly(nx) = (G1P1poly(nx)*2)+(16*5)-1;
	end
	if G1P1poly(nx)<=56 && G1P1poly(nx)>48
		G1P1poly(nx) = (G1P1poly(nx)*2)+(16*6)-1;
	end
	if G1P1poly(nx)<=64 && G1P1poly(nx)>56
		G1P1poly(nx) = (G1P1poly(nx)*2)+(16*7)-1;
	end
end




%-------------------------------%
L1x = linspace(1.0,3.2);
Tsz = numel(S1);			% Total matrix space
Ssz = numel(find(S1>0));	% Filled matrix space
Psz = round(Ssz/Tsz*100);	% Percent filled space
S1L1 = L1x(Psz);
%-------------------------------%
S1_occ=(S1>0);									% current surface SAPs
h=convn(S1_occ,h_mask,'same');					% energy field from mask
dE1=S1L1-h;										% energy diff for pop 1
P1=1./(1+exp(dE1*S1beta));						% cond. exocytosis for pop 1
Pen=S1_occ.*(1-exp(-S1mu*S1delta_t));			% large delta_t approx
P_rand=rand(size(S1));
S1_en=(Pen>P_rand);
P_ex1=(1-S1_occ).*(S1ro*S1r*S1delta_t*P1);
%-------------------------------%



G1RT = G1RTBASE;
if nn >= LTP1onG1 && nn <= LTP1offG1; 
		G1RT=G1RTLTP;
end
if nn >= LTP2onG1 && nn <= LTP2offG1; 
		G1RT=G1RTLTP; 
end


%------------------------%
if G1RTBASE > 0
%------------------------%
	
	G1P1ex = P_ex1;
	G1P1n = numel(G1P1poly);
	for polyid = 1:G1P1n
	G1P1IDo = G1P1poly(polyid);
	
	G1P1IDb = G1P1IDo+1;
	G1P1IDr = G1P1IDo+16;
	G1P1IDd = G1P1IDo+17;
	
	if G1P1IDo>=16
		G1P1IDb = G1P1IDo+1;
		G1P1IDr = G1P1IDo+16+(round(rand)*-32);
		G1P1IDd = G1P1IDr+1;
		% G1P1IDd = G1P1IDo+17+(round(rand)*-34);
	end
	
	G1P1ex(G1P1IDo)=G1RT;
	G1P1ex(G1P1IDb)=G1RT;
	G1P1ex(G1P1IDr)=G1RT;
	G1P1ex(G1P1IDd)=G1RT;
	end
	%-----
	S1_NONocc=(S1==0);
	G1P_ex1 = S1_NONocc.*G1P1ex;
	%-----
	P_ex1=G1P_ex1; % Enact LTP Manipulation
	%-----
%------------------------%
end % GluR1 G1RTBASE > 0
%------------------------%



%-------------------------------%
S1_ex=(P_ex1>P_rand);
S1=(S1_occ-S1_en).*S1+S1_ex;
%-------------------------------%

 
 
%if nn == 400; keyboard; end
end

%-------------------%
% S2_MainClusterFun
%-------------------%
function [S2] = S2_MainClusterFun(S2pop2_num, h_mask, S2L1_array, S2L2_array, S2beta,... 
    S2mu_array, S2delta_t, S2ro, S2r, S2ro2, S2, nn,...
	G1SAPLOC1,G1SAPLOC2,G2SAPLOC1,G2SAPLOC2,...
	SAPFmx1,SAPFmx2,SAPmx1,SAPmx2,...
	G1P1SAPM,G1P2SAPM,G2P1SAPM,G2P2SAPM)

  S2_occ=(S2>0);  
  h=convn(S2_occ,h_mask,'same');
  dE1=S2L1_array(nn)-h; % energy diff for pop #1
  dE2=S2L2_array(nn)-h; % energy diff for pop #2
  P1=1./(1+exp(dE1*S2beta)); % conditional exocytosis for pop 1
  P2=1./(1+exp(dE2*S2beta)); % conditional exocytosis for pop 2
  
  %Pen=S2_occ.*(mu_array(nn)*delta_t);
  Pen=S2_occ.*(1-exp(-S2mu_array(nn)*S2delta_t)); % large delta_t approx
  P_rand=rand(size(S2));
  S2_en=(Pen>P_rand);
  
  
  
  P_ex1=(1-S2_occ).*(S2ro*S2r*S2delta_t*P1);
  %P_ex1=(1-S_occ).*(1-exp(-ro*r*delta_t*P1)); % large delta_t approx
  S2_ex=(P_ex1>P_rand);
  
  S2_occ2=S2_occ+S2_ex;
  %P_ex2=(1-S2_occ2).*(ro2(nn)*r*delta_t*P2); % after normal exo and endo, check for pop2 exo
  P_ex2=(1-S2_occ2).*(1-exp(-S2ro2(nn)*S2r*S2delta_t*P2)); % after normal exo and endo, check for pop2 exo
  S2_ex2=2*(P_ex2>P_rand);
  S2=(S2_occ-S2_en).*S2+S2_ex+S2_ex2;
 
 

end
%===================================%




%{
function [S1 G1SP1] = S1_MainClusterFun(S1pop2_num, h_mask, S1L1_array, S1L2_array, S1beta,... 
    S1mu_array, S1delta_t, S1ro, S1r, S1ro2, S1, nn,...
	G1SAPLOC1,G1SAPLOC2,G2SAPLOC1,G2SAPLOC2,...
	SAPFmx1,SAPFmx2,SAPmx1,SAPmx2,...
	G1P1SAPM,G1P2SAPM,G2P1SAPM,G2P2SAPM,...
	G1SP1,G2SP1,G1LTParray,G2LTParray,...
	G1STBASE,G1RTBASE,G1STLTP,G1RTLTP,G1BSMu,G1LSMu,...
	G2STBASE,G2RTBASE,G2STLTP,G2RTLTP,G2BSMu,G2LSMu,...
	LTP1onG1,LTP1offG1,LTP2onG1,LTP2offG1)
%{
mxshape = [1:16]';
mxshape2 = mxshape;
for x=1:15
	y=x+1;
mxshape2(:,y) = mxshape2(:,x)+16;
end

mshape = [1:8]';
mshape2 = mshape;
for x=1:7
	y=x+1;
mshape2(:,y) = mshape2(:,x)+8;
end  
%}

%-------------------------------%
G1P1poly = G1P1SAPM(:,2);
[~, ~, G1P1poly] = find(G1P1poly);
for nx = 1:(numel(G1P1poly))
	if G1P1poly(nx)<=8 && G1P1poly(nx)>0
		G1P1poly(nx) = (G1P1poly(nx)*2)+(16*0)-1;
	end
	if G1P1poly(nx)<=16 && G1P1poly(nx)>8
		G1P1poly(nx) = (G1P1poly(nx)*2)+(16*1)-1;
	end
	if G1P1poly(nx)<=24 && G1P1poly(nx)>16
		G1P1poly(nx) = (G1P1poly(nx)*2)+(16*2)-1;
	end
	if G1P1poly(nx)<=32 && G1P1poly(nx)>24
		G1P1poly(nx) = (G1P1poly(nx)*2)+(16*3)-1;
	end
	if G1P1poly(nx)<=40 && G1P1poly(nx)>32
		G1P1poly(nx) = (G1P1poly(nx)*2)+(16*4)-1;
	end
	if G1P1poly(nx)<=48 && G1P1poly(nx)>40
		G1P1poly(nx) = (G1P1poly(nx)*2)+(16*5)-1;
	end
	if G1P1poly(nx)<=56 && G1P1poly(nx)>48
		G1P1poly(nx) = (G1P1poly(nx)*2)+(16*6)-1;
	end
	if G1P1poly(nx)<=64 && G1P1poly(nx)>56
		G1P1poly(nx) = (G1P1poly(nx)*2)+(16*7)-1;
	end
end


%-------------------------------%
% NOTES
%{
% To make SAP round or amorphous, set L = 2.05
and uncomment this section

if mod(nn,2)==0
h_mask = [1 0 1; 0 0 0; 1 0 1];
end

if mod(nn,4)==0
h_mask = [0 0 0; 1 0 1; 0 0 0];
end

if mod(nn,3)==0
S1L1_array(nn) = 3;
%S1beta = 60;
end

if mod(nn,4)==0
S1L1_array(nn) = 3;
%S1beta = 60;
end


%-------------------------------%
S1_occ=(S1>0);									% current surface SAPs
h=convn(S1_occ,h_mask,'same');					% energy field from mask

dE1=S1L1_array(nn)-h;							% energy diff for pop 1
dE2=S1L2_array(nn)-h;							% energy diff for pop 2
P1=1./(1+exp(dE1*S1beta));						% cond. exocytosis for pop 1
P2=1./(1+exp(dE2*S1beta));						% cond. exocytosis for pop 2
Pen=S1_occ.*(1-exp(-S1mu_array(nn)*S1delta_t)); % large delta_t approx
P_rand=rand(size(S1));
S1_en=(Pen>P_rand);
P_ex1=(1-S1_occ).*(S1ro*S1r*S1delta_t*P1);
%-------------------------------%
%}
%{
S1beta = 50

L1 -- Lattice repulsion constant

h -- the energy field from the mask; Mx of int vals from 0 - 4

dE1 -- Lattice repulsion constant minus the h energy field
dE1 = L1-h

> L1 is usually 1.2 and since the h-energy field ranges from
> values 0 to 4 here are some typical vals of dE1 when a GAP
> is surrounded by 0 1 2 3 or 4 SAPs

N surrounding SAPs:   0    1     2     3     4
L1(1.1):  dE1  =    [1.1  0.1  -0.9  -1.9  -2.9]
L1(1.2):  dE1  =    [1.2  0.2  -0.8  -1.8  -2.8]
L1(1.3):  dE1  =    [1.3  0.3  -0.7  -1.7  -2.7]

> as we will see, dE1 values closer to zero become more relevant
> because...

P1 -- conditional exocytosis probability Mx
P1 = 1 / (1  +  exp(dE1*S1beta))

> positive dE1 vals >  1  will never exo
> negative dE1 vals < -1  will always exo
> since...
> P1 results from:  1 / 1+exp()   
> and  exp()  inflates/deflates vals log-quickly  [exp(1)=2.7183 | exp(-1)=0.3679]
> even a dE1 val of .1 has a low P(exo)
> because typically...
> S1beta = 50
> thus:	 .1 * 50 = 5
> thus:  exp(5) = 150
> makes... P(exo) i.e. P1 = 1/(1+exp(dE1*S1beta))
> P1 = 1/151 = .006



G1LTParray
%}
%-------------------------------%

% %-------------------------------%
S1_occ=(S1>0);									% current surface SAPs
h=convn(S1_occ,h_mask,'same');					% energy field from mask
dE1=S1L1_array(nn)-h;							% energy diff for pop 1
dE2=S1L2_array(nn)-h;							% energy diff for pop 2
P1=1./(1+exp(dE1*S1beta));						% cond. exocytosis for pop 1
P2=1./(1+exp(dE2*S1beta));						% cond. exocytosis for pop 2
Pen=S1_occ.*(1-exp(-S1mu_array(nn)*S1delta_t)); % large delta_t approx
P_rand=rand(size(S1));
S1_en=(Pen>P_rand);
P_ex1=(1-S1_occ).*(S1ro*S1r*S1delta_t*P1);
% %-------------------------------%





%-------------------------------%
% Build LTP Manipulation
%-------------------------------%
%{
% % GluR1 G1RTBASE > 0
% if G1RTBASE > 0
% 	G1P1ex = P_ex1;
% 	G1P1n = numel(G1P1poly);
% 	for polyid = 1:G1P1n
% 	G1P1IDo = G1P1poly(polyid);
% 	G1P1IDb = G1P1IDo+1;
% 	G1P1IDr = G1P1IDo+16;
% 	G1P1IDd = G1P1IDo+17;
% 	G1P1ex(G1P1IDo)=G1RTBASE;
% 	G1P1ex(G1P1IDb)=G1RTBASE;
% 	G1P1ex(G1P1IDr)=G1RTBASE;
% 	G1P1ex(G1P1IDd)=G1RTBASE;
% 	end
% 	%-----
% 	S1_NONocc=(S1==0);
% 	G1P_ex1 = S1_NONocc.*G1P1ex;
% 	%-----
% 	P_ex1=G1P_ex1; % Enact LTP Manipulation
% 	%-----
% end
%}


G1RT = G1RTBASE;
if nn >= LTP1onG1 && nn <= LTP1offG1; 
		G1RT=G1RTLTP;
end
if nn >= LTP2onG1 && nn <= LTP2offG1; 
		G1RT=G1RTLTP; 
end


%------------------------%
if G1RTBASE > 0
%------------------------%
	
	G1P1ex = P_ex1;
	G1P1n = numel(G1P1poly);
	for polyid = 1:G1P1n
	G1P1IDo = G1P1poly(polyid);
	
	G1P1IDb = G1P1IDo+1;
	G1P1IDr = G1P1IDo+16;
	G1P1IDd = G1P1IDo+17;
	
	if G1P1IDo>=16
		G1P1IDb = G1P1IDo+1;
		G1P1IDr = G1P1IDo+16+(round(rand)*-32);
		G1P1IDd = G1P1IDr+1;
		% G1P1IDd = G1P1IDo+17+(round(rand)*-34);
	end
	
	G1P1ex(G1P1IDo)=G1RT;
	G1P1ex(G1P1IDb)=G1RT;
	G1P1ex(G1P1IDr)=G1RT;
	G1P1ex(G1P1IDd)=G1RT;
	end
	%-----
	S1_NONocc=(S1==0);
	G1P_ex1 = S1_NONocc.*G1P1ex;
	%-----
	P_ex1=G1P_ex1; % Enact LTP Manipulation
	%-----
%------------------------%
end % GluR1 G1RTBASE > 0
%------------------------%


%{
% GluR1 G1STAB == 2
if G1STAB == 2
	G1P1ex = P_ex1;
	G1P1n = numel(G1P1poly);
	for polyid = 1:G1P1n
	G1P1IDo = G1P1poly(polyid);
	G1P1IDb = G1P1IDo+1;
	G1P1IDr = G1P1IDo+16;
	G1P1IDd = G1P1IDo+17;
	G1P1ex(G1P1IDo)=.5;
	G1P1ex(G1P1IDb)=.5;
	G1P1ex(G1P1IDr)=.5;
	G1P1ex(G1P1IDd)=.5;
	end
	%-----
	S1_NONocc=(S1==0);
	G1P_ex1 = S1_NONocc.*G1P1ex;
	%-----
	P_ex1=G1P_ex1; % Enact LTP Manipulation
	%-----
end

% GluR1 G1STAB == 1
if G1STAB == 1
	G1P1ex = P_ex1;
	G1P1n = numel(G1P1poly);
	for polyid = 1:G1P1n
	G1P1IDo = G1P1poly(polyid);
	G1P1IDb = G1P1IDo+1;
	G1P1IDr = G1P1IDo+16;
	G1P1IDd = G1P1IDo+17;
	G1P1ex(G1P1IDo)=.3;
	G1P1ex(G1P1IDb)=.3;
	G1P1ex(G1P1IDr)=.3;
	G1P1ex(G1P1IDd)=.3;
	end
	%-----
	S1_NONocc=(S1==0);
	G1P_ex1 = S1_NONocc.*G1P1ex;
	%-----
	P_ex1=G1P_ex1; % Enact LTP Manipulation
	%-----
end
%}

%-------------------------------%
S1_ex=(P_ex1>P_rand);
S1_occ2=S1_occ+S1_ex;
P_ex2=(1-S1_occ2).*(1-exp(-S1ro2(nn)*S1r*S1delta_t*P2)); % after exo/endo, check pop2 exo
S1_ex2=2*(P_ex2>P_rand);
S1=(S1_occ-S1_en).*S1+S1_ex+S1_ex2;
%-------------------------------%

%{
S1_occ=(S1>0);  
% S1mem_num(nn)=sum(S1_occ(:));
% S1pop2_num(nn)=sum((S1(:)>1));
h=convn(S1_occ,h_mask,'same');
dE1=S1L1_array(nn)-h; % energy diff for pop #1
dE2=S1L2_array(nn)-h; % energy diff for pop #2
P1=1./(1+exp(dE1*S1beta)); % conditional exocytosis for pop 1
P2=1./(1+exp(dE2*S1beta)); % conditional exocytosis for pop 2

%Pen=S1_occ.*(mu_array(nn)*delta_t);
Pen=S1_occ.*(1-exp(-S1mu_array(nn)*S1delta_t)); % large delta_t approx
P_rand=rand(size(S1));
S1_en=(Pen>P_rand);
  
  

P_rand=G1P1ex;


P_ex1=(1-S1_occ).*(S1ro*S1r*S1delta_t*P1);
%P_ex1=(1-S1_occ).*(S1ro*S1r*S1delta_t*P1);
%P_ex1=(1-S_occ).*(1-exp(-ro*r*delta_t*P1)); % large delta_t approx
S1_ex=(P_ex1>P_rand);

S1_occ2=S1_occ+S1_ex;
%P_ex2=(1-S1_occ2).*(ro2(nn)*r*delta_t*P2); % after normal exo and endo, check for pop2 exo
P_ex2=(1-S1_occ2).*(1-exp(-S1ro2(nn)*S1r*S1delta_t*P2)); % after normal exo and endo, check for pop2 exo
S1_ex2=2*(P_ex2>P_rand);
S1=(S1_occ-S1_en).*S1+S1_ex+S1_ex2;
%}
 
 
%if nn == 400; keyboard; end
end
%}





%%				   SPECIALIZED TOOLS
%-------------##########################------------------%
%				   SPECIALIZED TOOLS
%-------------##########################------------------%


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
% IzhFun SPIKE GENERATOR
%===================================%
function [] = IzhFun(spikes)

global v;
izhikevich(spikes)

end
%===================================%


%===================================%
%		ZGEN FUNCTION
%===================================%
function [G1Z G2Z] = ZGEN(stepN,G1Z,G2Z,...
    GluR2Ndots, GluR2xyl,GluR1Ndots, GluR1xyl,...
	XYLTpr1,XYRTpr1,XYLBpr1,XYRBpr1,XYLTpr2,XYRTpr2,XYLBpr2,XYRBpr2,...
	XYLTp1,XYRTp1,XYLBp1,XYRBp1,XYLTp2,XYRTp2,XYLBp2,XYRBp2)



%--GLUR2--%
for j = 1:GluR2Ndots 
%Inside PSD1
if (GluR2xyl(1,j)>=XYLBp1(1) && GluR2xyl(1,j)<=XYRTp1(1)) &&...
   (GluR2xyl(2,j)>=XYLBp1(2) && GluR2xyl(2,j)<=XYRTp1(2))
 
    G2Z(:,j) = 4;
 
%Inside PSD2
elseif  (GluR2xyl(1,j)>=XYLBp2(1) && GluR2xyl(1,j)<=XYRTp2(1)) &&...
        (GluR2xyl(2,j)>=XYLBp2(2) && GluR2xyl(2,j)<=XYRTp2(2))
    G2Z(:,j) = 4;
 
%Inside PSA1
elseif  (GluR2xyl(1,j)>=XYLBpr1(1) && GluR2xyl(1,j)<=XYRTpr1(1)) &&...
        (GluR2xyl(2,j)>=XYLBpr1(2) && GluR2xyl(2,j)<=XYRTpr1(2)) &&...
        ((GluR2xyl(1,j)<=XYLBp1(1) || GluR2xyl(1,j)>=XYRTp1(1))) &&...
        ((GluR2xyl(2,j)<=XYLBp1(2) || GluR2xyl(2,j)>=XYRTp1(2)))
    
    G2Z(:,j) = 2+(.1 * randi(15,1));
 
%Inside PSA2
elseif  (GluR2xyl(1,j)>=XYLBpr2(1) && GluR2xyl(1,j)<=XYRTpr2(1)) &&...
        (GluR2xyl(2,j)>=XYLBpr2(2) && GluR2xyl(2,j)<=XYRTpr2(2)) &&...
        ((GluR2xyl(1,j)<=XYLBp2(1) || GluR2xyl(1,j)>=XYRTp2(1))) &&...
        ((GluR2xyl(2,j)<=XYLBp2(2) || GluR2xyl(2,j)>=XYRTp2(2)))
    
    G2Z(:,j) = 2+(.1 * randi(15,1));
 
else
    G2Z(:,j) = 1;
end
 
end
 
%--GLUR1--%
for j = 1:GluR1Ndots 
%Inside PSD1
if (GluR1xyl(1,j)>=XYLBp1(1) && GluR1xyl(1,j)<=XYRTp1(1)) &&...
   (GluR1xyl(2,j)>=XYLBp1(2) && GluR1xyl(2,j)<=XYRTp1(2))
 
    G1Z(:,j) = 4;
 
%Inside PSD2
elseif  (GluR1xyl(1,j)>=XYLBp2(1) && GluR1xyl(1,j)<=XYRTp2(1)) &&...
        (GluR1xyl(2,j)>=XYLBp2(2) && GluR1xyl(2,j)<=XYRTp2(2))
    G1Z(:,j) = 4;
 
%Inside PSA1
elseif  (GluR1xyl(1,j)>=XYLBpr1(1) && GluR1xyl(1,j)<=XYRTpr1(1)) &&...
        (GluR1xyl(2,j)>=XYLBpr1(2) && GluR1xyl(2,j)<=XYRTpr1(2)) &&...
        ((GluR1xyl(1,j)<=XYLBp1(1) || GluR1xyl(1,j)>=XYRTp1(1))) &&...
        ((GluR1xyl(2,j)<=XYLBp1(2) || GluR1xyl(2,j)>=XYRTp1(2)))
    
    G1Z(:,j) = 2+(.1 * randi(15,1));
 
%Inside PSA2
elseif  (GluR1xyl(1,j)>=XYLBpr2(1) && GluR1xyl(1,j)<=XYRTpr2(1)) &&...
        (GluR1xyl(2,j)>=XYLBpr2(2) && GluR1xyl(2,j)<=XYRTpr2(2)) &&...
        ((GluR1xyl(1,j)<=XYLBp2(1) || GluR1xyl(1,j)>=XYRTp2(1))) &&...
        ((GluR1xyl(2,j)<=XYLBp2(2) || GluR1xyl(2,j)>=XYRTp2(2)))
    
    G1Z(:,j) = 2+(.1 * randi(15,1));
 
else
    G1Z(:,j) = 1;
end
 
end
 

%{
%-------------------%
for j = 1:GluR2Ndots
%-------------------%

		if GluR2xyl(1,j)>=XcPr1lft && GluR2xyl(1,j)<=XcPr1rit &&...
          GluR2xyl(2,j)>=YrPr1top && GluR2xyl(2,j)<=YrPr1bot &&...
		  GluR2xyl(1,j)>=(XcPr1lft+PSDSZE(2,1)) && GluR2xyl(1,j)<=(XcPr1rit-PSDSZE(2,1)) &&...
          GluR2xyl(2,j)>=(YrPr1top+PSDSZE(2,1)) && GluR2xyl(2,j)<=(YrPr1bot-PSDSZE(2,1))
		
			G2Z(:,j) = 3;
	
		elseif GluR2xyl(1,j)>=XcPr1lft && GluR2xyl(1,j)<=XcPr1rit &&...
          GluR2xyl(2,j)>=YrPr2bot && GluR2xyl(2,j)<=YrPr2top &&...
		  GluR2xyl(1,j)>=(XcPr1lft+PSDSZE(2,2)) && GluR2xyl(1,j)<=(XcPr1rit-PSDSZE(2,2)) &&...
          GluR2xyl(2,j)>=(YrPr2bot+PSDSZE(2,2)) && GluR2xyl(2,j)<=(YrPr2top-PSDSZE(2,2))

			G2Z(:,j) = 3;
			
			
		elseif GluR2xyl(1,j)>=XcPr1lft && GluR2xyl(1,j)<=XcPr1rit &&...
          GluR2xyl(2,j)>=YrPr1top && GluR2xyl(2,j)<=YrPr1bot &&...
          ((GluR2xyl(1,j)<=(XcPr1lft+PSDSZE(2,1)) || GluR2xyl(1,j)>=(XcPr1rit-PSDSZE(2,1))) ||...
          (GluR2xyl(2,j)<=(YrPr1top+PSDSZE(2,1)) || GluR2xyl(2,j)>=(YrPr1bot-PSDSZE(2,1))))
      
            G2Z(:,j) = 1+(.1 * randi(15,1));
            
        elseif GluR2xyl(1,j)>=XcPr1lft && GluR2xyl(1,j)<=XcPr1rit &&...
          GluR2xyl(2,j)>=YrPr2bot && GluR2xyl(2,j)<=YrPr2top &&...
          ((GluR2xyl(1,j)<=(XcPr1lft+PSDSZE(2,2)) || GluR2xyl(1,j)>=(XcPr1rit-PSDSZE(2,2))) ||...
          (GluR2xyl(2,j)<=(YrPr2bot+PSDSZE(2,2)) || GluR2xyl(2,j)>=(YrPr2top-PSDSZE(2,2))))
	  
            G2Z(:,j) = 1+(.1 * randi(15,1));
			
			
		else
		   G2Z(:,j) = 1;
		end % if
%-------------------%
end  %j = 1:GluR2Ndots
%-------------------%



%-------------------%
for j = 1:GluR1Ndots
%-------------------%
 
        if GluR1xyl(1,j)>=XcPr1lft && GluR1xyl(1,j)<=XcPr1rit &&...
          GluR1xyl(2,j)>=YrPr1top && GluR1xyl(2,j)<=YrPr1bot &&...
          GluR1xyl(1,j)>=(XcPr1lft+PSDSZE(2,1)) && GluR1xyl(1,j)<=(XcPr1rit-PSDSZE(2,1)) &&...
          GluR1xyl(2,j)>=(YrPr1top+PSDSZE(2,1)) && GluR1xyl(2,j)<=(YrPr1bot-PSDSZE(2,1))
        
            G1Z(:,j) = 3;
    
        elseif GluR1xyl(1,j)>=XcPr1lft && GluR1xyl(1,j)<=XcPr1rit &&...
          GluR1xyl(2,j)>=YrPr2bot && GluR1xyl(2,j)<=YrPr2top &&...
          GluR1xyl(1,j)>=(XcPr1lft+PSDSZE(2,2)) && GluR1xyl(1,j)<=(XcPr1rit-PSDSZE(2,2)) &&...
          GluR1xyl(2,j)>=(YrPr2bot+PSDSZE(2,2)) && GluR1xyl(2,j)<=(YrPr2top-PSDSZE(2,2))
 
            G1Z(:,j) = 3;
			
			
		elseif GluR1xyl(1,j)>=XcPr1lft && GluR1xyl(1,j)<=XcPr1rit &&...
          GluR1xyl(2,j)>=YrPr1top && GluR1xyl(2,j)<=YrPr1bot &&...
          ((GluR1xyl(1,j)<=(XcPr1lft+PSDSZE(2,1)) || GluR1xyl(1,j)>=(XcPr1rit-PSDSZE(2,1))) ||...
          (GluR1xyl(2,j)<=(YrPr1top+PSDSZE(2,1)) || GluR1xyl(2,j)>=(YrPr1bot-PSDSZE(2,1))))
      
            G1Z(:,j) = 1+(.1 * randi(15,1));
            
        elseif GluR1xyl(1,j)>=XcPr1lft && GluR1xyl(1,j)<=XcPr1rit &&...
          GluR1xyl(2,j)>=YrPr2bot && GluR1xyl(2,j)<=YrPr2top &&...
          ((GluR1xyl(1,j)<=(XcPr1lft+PSDSZE(2,2)) || GluR1xyl(1,j)>=(XcPr1rit-PSDSZE(2,2))) ||...
          (GluR1xyl(2,j)<=(YrPr2bot+PSDSZE(2,2)) || GluR1xyl(2,j)>=(YrPr2top-PSDSZE(2,2))))
      
            G1Z(:,j) = 1+(.1 * randi(15,1));
			
            
        else
           G1Z(:,j) = 1;
        end % if
%-------------------%
end  %j = 1:GluR1Ndots
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








%%					PARTICLE COUNTERS
%-------------##########################------------------%
%					PARTICLE COUNTERS
%-------------##########################------------------%

%-------------------------------%
%       MAIN PARTICLE COUNTER
%-------------------------------%
function [GR1SPY1N GR1SPY2N GR1PER1N GR1PER2N GR1PSD1N GR1PSD2N...
		GR2SPY1N GR2SPY2N GR2PER1N GR2PER2N GR2PSD1N GR2PSD2N...
		GR1PERN GR1PSDN GR2PERN GR2PSDN GR1SPYN GR2SPYN...
		SPY1N SPY2N SPYNT ESNT EST]...
		= dotCount(stepN,GluR2Ndots, GluR1Ndots,...
		G2INSPYN1,G2INSPYN2,G2INPSD1,G2INPSD2,G2INPERI1,G2INPERI2,...
		G1INSPYN1,G1INSPYN2,G1INPSD1,G1INPSD2,G1INPERI1,G1INPERI2)
%================%
% if stepN == 600; keyboard; end
%================%


GR1SPY1N = 0; GR1SPY2N = 0; 
GR1PER1N = 0; GR1PER2N = 0; 
GR1PSD1N = 0; GR1PSD2N = 0; 
GR2SPY1N = 0; GR2SPY1N = 0; 
GR2PER1N = 0; GR2PER2N = 0; 
GR2PSD1N = 0; GR2PSD2N = 0; 
GR1PERN = 0; GR1PSDN = 0; 
GR2PERN = 0; GR2PSDN = 0; 
GR1SPYN = 0; GR2SPYN = 0; 
SPY1N = 0; SPY2N = 0;
SPYNT = 0; ESNT = 0; EST = 0;



%================%
GR1SPY1N = sum(G1INSPYN1);	% GR1PSD1
GR1SPY2N = sum(G1INSPYN2);	% GR1PSD2
GR1PER1N = sum(G1INPERI1);	% GR1PS1T
GR1PER2N = sum(G1INPERI2);	% GR1PS2T
GR1PSD1N = sum(G1INPSD1);	% GR1PSD1T
GR1PSD2N = sum(G1INPSD2);	% GR1PSD2T
%==========%
GR2SPY1N = sum(G2INSPYN1);	% GR2SPY1N
GR2SPY2N = sum(G2INSPYN2);	% GR2SPY2N
GR2PER1N = sum(G2INPERI1);	% GR2PER1N
GR2PER2N = sum(G2INPERI2);	% GR2PER2N
GR2PSD1N = sum(G2INPSD1);	% GR2PSD1N
GR2PSD2N = sum(G2INPSD2);	% GR2PSD2N
%================%

GR1PERN = GR1PER1N + GR1PER2N;	% GR1PERN
GR1PSDN = GR1PSD1N + GR1PSD2N;	% GR1PSDN
GR2PERN = GR2PER1N + GR2PER2N;	% GR2PERN
GR2PSDN = GR2PSD1N + GR2PSD2N;	% GR2PSDN

GR1SPYN = GR1SPY1N + GR1SPY2N;	% GR1SPYN
GR2SPYN = GR2SPY1N + GR2SPY2N;	% GR2SPYN
%GR1SPYN = GR1PERN + GR1PSDN;
%GR2SPYN = GR2PERN + GR2PSDN;

SPY1N = GR1SPY1N + GR2SPY1N;	% SPY1N SPY1N
SPY2N = GR1SPY2N + GR2SPY2N;	% SPY2N SPY2N
%SPY1N = GR1PER1N+GR1PSD1N+GR2PER1N+GR2PSD1N;
%SPY2N = GR1PER2N+GR1PSD2N+GR2PER2N+GR2PSD2N;

SPYNT = SPY1N + SPY2N;							% SPYNT SPYNT
ESNT = (GluR2Ndots+GluR1Ndots)-(SPY1N+SPY2N);	% ESNT
EST = [ESNT SPY1N SPY2N SPYNT];					% EST



end
%================================%



%-------------------------------%
%   GRAPH PARTICLE COUNTERS
%-------------------------------%
function [] = GRAPHCOUNTS(stepN,GluR2Ndots,doGluR1,doGluR2,doRun1,doRun2,doRun6,...
	HSlow, HShi, G1FSLOTS,G2FSLOTS,SAP5,S1sum,S2sum,doKo,...
	EST,SPY1N,SPY2N)



doPLOT2D = doRun1;
doPLOT3D = doRun2;
doHOMEO = doRun6;

if GluR2Ndots == 0
	GluR2Ndots = 50;
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
	set(gca,'Ylim',[0 GluR2Ndots],'XLim',[.5 4.5])
	title('AMPAR#:   EC    -   SPINE1   -     SPINE2  -    SPINES');
	end
end

if doPLOT3D == 1
	if mod(stepN, 100) == 0
	figure(1)
	subplot(5,5,[21 22]);
	bar(1:4, [EST'], .8, 'EdgeColor',[1 0.5 0.5]);
	set(gca,'Ylim',[0 GluR2Ndots],'XLim',[.5 4.5])
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





%%			4. DENDRITIC FIELD MAPPING TOOLS
%-------------##########################------------------%
%		FIELD MAP FOR DENDRITE AREA AND PSD AREAS
%-------------##########################------------------%

% FieldFun
function [Yrows Xcols YrPr2bot XcPr2rit YrPr2top XcPr2lft...
YrPr1bot XcPr1rit YrPr1top XcPr1lft PSDfield...
PSDSZE Zfield XWIDE YHIGH...
PSD1WH PSD2WH PERI1WH PERI2WH SPYN1WH SPYN2WH...
XYLTpr1 XYRTpr1 XYLBpr1 XYRBpr1 XYLTpr2 XYRTpr2 XYLBpr2 XYRBpr2...
XYLTp1 XYRTp1 XYLBp1 XYRBp1 XYLTp2 XYRTp2 XYLBp2 XYRBp2...
XYBOXpr1 XYBOXpr2 XYBOXp1 XYBOXp2...
SPYN1xv SPYN1yv SPYN2xv SPYN2yv...
PSD1xv PSD1yv PSD2xv PSD2yv...
PERI1xv PERI1yv PERI2xv PERI2yv...
fPSD1 fPSD2...
] = FieldFun(fsizeX, fsizeY,PSD1size, PSD2size, periPSD1size, periPSD2size)


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
PSDSZE = [PSD1size PSD2size; periPSD1size periPSD2size];
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
xlim = [0 XcolsWIDE];
ylim = [-YrowsHIGH YrowsHIGH];
%---
SNSZ = get(0,'ScreenSize');
fig99 = figure(99);

figure(fig99)
subplot(5,5,[3 20]), subplot('Position',[.45 .40 .5 .59]),...
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

figure(fig99)
subplot(5,5,[1 6]), imagesc(Zfield);     % <--FIG--##
colormap('bone')
title('PSD1');
figure(fig99)
subplot(5,5,[2 7]), imagesc(PSDfield); % <--FIG--##
title('PSD2');

set(gcf,'Position',...
	[SNSZ(3)/1.5	SNSZ(4)/5		SNSZ(3)/3	SNSZ(4)/1.5]);


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



%=================================%
%{
plot(SPYN1xv,SPYN1yv,GluR2xyl(1,G2INSPYN1),GluR2xyl(2,G2INSPYN1),...
    'r+',GluR2xyl(1,~G2INSPYN1),GluR2xyl(2,~G2INSPYN1),'bo');
hold on
plot(SPYN2xv,SPYN2yv,GluR2xyl(1,G2INSPYN2),GluR2xyl(2,G2INSPYN2),...
    'r+',GluR2xyl(1,~G2INSPYN2),GluR2xyl(2,~G2INSPYN2),'bo');
hold off

plot(PSD1xv,PSD1yv,GluR2xyl(1,G2INPSD1),GluR2xyl(2,G2INPSD1),...
	'r+',GluR2xyl(1,~G2INPSD1),GluR2xyl(2,~G2INPSD1),'bo');
hold on
plot(PSD2xv,PSD2yv,GluR2xyl(1,G2INPSD2),GluR2xyl(2,G2INPSD2),...
	'r+',GluR2xyl(1,~G2INPSD2),GluR2xyl(2,~G2INPSD2),'bo');
hold off

plot(PERI1xv,PERI1yv,GluR2xyl(1,G2INPERI1),GluR2xyl(2,G2INPERI1),...
    'r+',GluR2xyl(1,~G2INPERI1),GluR2xyl(2,~G2INPERI1),'bo');
hold on
plot(PERI2xv,PERI2yv,GluR2xyl(1,G2INPERI2),GluR2xyl(2,G2INPERI2),...
    'r+',GluR2xyl(1,~G2INPERI2),GluR2xyl(2,~G2INPERI2),'bo');
hold off
%}
%{
plot(SPYN1xv,SPYN1yv,GluR1xyl(1,G1INSPYN1),GluR1xyl(2,G1INSPYN1),...
    'r+',GluR1xyl(1,~G1INSPYN1),GluR1xyl(2,~G1INSPYN1),'bo');
hold on
plot(SPYN2xv,SPYN2yv,GluR1xyl(1,G1INSPYN2),GluR1xyl(2,G1INSPYN2),...
    'r+',GluR1xyl(1,~G1INSPYN2),GluR1xyl(2,~G1INSPYN2),'bo');
hold off
 
plot(PSD1xv,PSD1yv,GluR1xyl(1,G1INPSD1),GluR1xyl(2,G1INPSD1),...
    'r+',GluR1xyl(1,~G1INPSD1),GluR1xyl(2,~G1INPSD1),'bo');
hold on
plot(PSD2xv,PSD2yv,GluR1xyl(1,G1INPSD2),GluR1xyl(2,G1INPSD2),...
    'r+',GluR1xyl(1,~G1INPSD2),GluR1xyl(2,~G1INPSD2),'bo');
hold off
 
plot(PERI1xv,PERI1yv,GluR1xyl(1,G1INPERI1),GluR1xyl(2,G1INPERI1),...
    'r+',GluR1xyl(1,~G1INPERI1),GluR1xyl(2,~G1INPERI1),'bo');
hold on
plot(PERI2xv,PERI2yv,GluR1xyl(1,G1INPERI2),GluR1xyl(2,G1INPERI2),...
    'r+',GluR1xyl(1,~G1INPERI2),GluR1xyl(2,~G1INPERI2),'bo');
hold off
%}
%=================================%
end



%%			  PLOTTING AND LIVE SIMULATION
%-------------##########################------------------%
%				VISUALIZATION FUNCTIONS
%-------------##########################------------------%

%-------------------------------------------%
% PLOT Particle Motion
%-------------------------------------------%
function [] = MAINPLOT(stepN,GluR2xyl,GluR1xyl,...
	XWIDE,YHIGH,XYBOXpr1,XYBOXpr2,XYBOXp1,XYBOXp2)
%-------------------------------------------%
xlim = [0 XWIDE];
ylim = [-YHIGH YHIGH];
%---

figure(1)
subplot(5,5,[3 25]), 
scatter(GluR2xyl(1,:),GluR2xyl(2,:),5,[0 0 1]);
hold on;
subplot(5,5,[3 25]), 
scatter(GluR1xyl(1,:),GluR1xyl(2,:),5,[1 0 0]);
axis([xlim, ylim]);
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(get(gca,'XLabel'),'String','GluR2:Blue \bullet \bullet GluR1:Red')
rectangle('Position',[XYBOXpr1(1),XYBOXpr1(2),XYBOXpr1(3),XYBOXpr1(4)])
rectangle('Position',[XYBOXpr2(1),XYBOXpr2(2),XYBOXpr2(3),XYBOXpr2(4)])
rectangle('Position',[XYBOXp1(1),XYBOXp1(2),XYBOXp1(3),XYBOXp1(4)])
rectangle('Position',[XYBOXp2(1),XYBOXp2(2),XYBOXp2(3),XYBOXp2(4)])
hold off;



%=================================%
%           3D PLOT
%---------------------------------%
%{
%=================================%
%           3D PLOT
%---------------------------------%
%     MAIN DOTS
%----------------------%
if do3DPLOT
figure(1);
subplot(5,5,[1 7]), 
%gscatter(GluR2xyl(1,:),GluR2xyl(2,:)); view(20, 30);
scatter3(GluR2xyl(1,:),GluR2xyl(2,:),G2Z(:),'.'), view(20, 30)
axis normal;
grid off
axis([xlim, ylim, zlim]);
set(gca, 'Box', 'on');
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
%-------%
%}
%=================================%

end


%-------------------------------------------%
% PLOT3DS
%-------------------------------------------%
function [] = PLOT3DS(stepN,GluR2xyl,GluR1xyl,XWIDE,YHIGH,G1Z,G2Z,Zfield)
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
ZI = griddata(GluR2xyl(1,:),GluR2xyl(2,:),G2Z(:),XI,YI,'v4');
end

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





%======================%
% GLUR2 DOTS
%----------------------%
figure(1);
subplot('Position',[.02 .61 .4 .35]),...
scatter3(GluR2xyl(1,:),GluR2xyl(2,:),G2Z(:),'.','MarkerEdgeColor',[0 0 1]),...
view(70, 10)
hold on;
subplot('Position',[.02 .61 .4 .35]),...
scatter3(GluR1xyl(1,:),GluR1xyl(2,:),G1Z(:),'.','MarkerEdgeColor',[1 0 0]),...
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
scatter3(GluR2xyl(1,:),GluR2xyl(2,:),G2Z(:),'.'), view(60, 10)
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
scatter3(GluR1xyl(1,:),GluR1xyl(2,:),G1Z(:),'.'), view(60, 10)
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
function [] = ONEDOTPLOT(stepN,Nsteps,D, GluR2xyl, xyl2, XWIDE,YHIGH,...
				XYBOXpr1,XYBOXpr2,XYBOXp1,XYBOXp2)
%-------------------------------------------%
xlim = [0 XWIDE];
ylim = [-YHIGH YHIGH];
%---

xylcat = cat(2, xyl2,GluR2xyl)';
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
function [] = ONEDOTPLOT2(stepN,Nsteps,GluR2xyl,xyl2,XWIDE,YHIGH,...
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
    colormap('bone')
    title('2D Particle Map');
    title('PSD1');
	figure(1);
    subplot(5,5,12), imagesc(S2); % <--FIG--##
    title('PSD2');
end



%-------------------------------------------%
% GAUSSIAN SLOT COLORMAP FUN
%-------------------------------------------%
function [] = SLOTMAPFUN(G1STBASE,G1RTBASE,G1STLTP,G1RTLTP,G1BSMu,G1LSMu,...
G2STBASE,G2RTBASE,G2STLTP,G2RTLTP,G2BSMu,G2LSMu)
		 
saps = 0:4;
tails = 0:.5:4;
offsap = G1BSMu;

%param = [saps+tails(1) saps+tails(2) saps+tails(3) saps+tails(4) saps+tails(5)]-6

param = [tails+saps(1); tails+saps(2); tails+saps(3);...
	     tails+saps(4); tails+saps(5)]-offsap;

cdfMx = 1-cdf('norm',0,param,1);
cdfMx2 = [linspace(0,1,numel(tails)); cdfMx];

cdfMxpad = padarray(cdfMx,[1 1],0,'pre');
cdfMxpad(:,1)=[0 saps]';
cdfMxpad(1,:)=[0 tails]

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
%=========================================================%
yt = [-1 saps];
xt = tails;

subplot(5,5,[16 22]), subplot('Position',[.075 .045 .7 .3]),...
imagesc(cdfMx2)
%axis image
colormap('bone')
set(gca,'YTickLabel', sprintf('%.1f|',yt))
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'YLabel'),'String','CLOSE  SAPS')
set(get(gca,'XLabel'),'String','STICK PARAMETER')

text(1,0,...
strcat('Gaussian On µ : ', num2str(offsap)),'FontSize',14,'HorizontalAlignment','left');

text(0,1,...
strcat('P(on)= ', num2str(ProbCDF(:)')),'FontSize',14.5,'HorizontalAlignment','left');


% SAPMAPTITLE=title(char('Gaussian: ', num2str(offsap)));
% set(SAPMAPTITLE, 'FontSize', 14);
%=========================================================%
	
end












%%		  MSD BROWNIAN MOTION ANALYSIS TOOLS
%-------------##########################------------------%
%			 DIFFUSION ANALYSIS FUNCTIONS
%-------------##########################------------------%
%-------------------------------------------%
% MSD PARTICLE MOVEMENT GENERATION
%-------------------------------------------%
function [GluR2xyl GluR2xyds] = MSDAMPARSTEP(GluR2Ndots, GluR2xyds, GluR2xyl,... 
	testPSD1MSD, testPSD2MSD, testESMSD, PSD1, PSD2, D, Dn_PSD1, Dn_PSD2)
	
	for j = 1:GluR2Ndots
        
		if testPSD1MSD
           GluR2xyds(:,j) = GluR2xyds(:,j)*PSD1;
        elseif testPSD2MSD
            GluR2xyds(:,j) = GluR2xyds(:,j)*PSD2;
		elseif testESMSD
            GluR2xyds(:,j) = GluR2xyds(:,j);
		end
		
	   
        GluR2xyl(:,j) = GluR2xyl(:,j)+GluR2xyds(:,j);
	end

end

%-------------------------------------------%
% MEAN SQUARED DISPLACEMENT FUNCTION
%-------------------------------------------%
function [tracks] = MSDfun(stepN, Nsteps, tracks, GluR2xyds)
    time = (0:Nsteps-1)';
    xymsd = GluR2xyds';
    xymsd = cumsum(xymsd,1);
    tracks{stepN} = [time xymsd];

end


%-------------------------------------------%
% MSDfunction - FINAL TRACKS ANALYSIS
%-------------------------------------------%
function [MSDout] = MSDfunction(tracks, GluR2Ndots, Nsteps, D, Dn_PSD1, Dn_PSD2, dT, k, MSDtest)

SPACE_UNITS = 'µm';
TIME_UNITS = 's';
N_PARTICLES = GluR2Ndots;
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



%-------------############################------------------%
%			  ##   MSD2 SUBFUNCTIONS    ##
%-------------############################------------------% 

%-------------------------------------------%
% STEP SIZE GENERATOR
%-------------------------------------------%
function GluR2xyds = STEPxyds2(GluR2Ndots, k)

    GluR2xyds = (k * randn(2,GluR2Ndots));

end


%-------------------------------------------%
% MOVE PARTICLES MAIN FUNCTION
%-------------------------------------------%
function [GluR2xyl] = AMPARSTEP2(GluR2Ndots, GluR2xyds, GluR2xyl)
	
	for j = 1:GluR2Ndots
        GluR2xyl(:,j) = GluR2xyl(:,j)+GluR2xyds(:,j);
	end
	
end


%-------------------------------------------%
% LIVE DIFFUSION PLOT
%-------------------------------------------%
function [] = MAINPLOT2(GluR2xyl, lims)
%-------------------------------------------%
xlim = [-lims lims];
ylim = [-lims lims];
zlim = [-5 5];

%=================================%
%       MAIN 2D PLOT
%---------------------------------%
figure(1)
subplot(2,1,1), 
AMPARPlot = gscatter(GluR2xyl(1,:),GluR2xyl(2,:));
axis([xlim, ylim]);
set(AMPARPlot,'marker','.','markersize',[6],'color',[1 0 0])


%=================================%
%           3D PLOT
%---------------------------------%
figure(1);
subplot(2,1,2), 
gscatter(GluR2xyl(1,:),GluR2xyl(2,:)); view(20, 30);
axis normal;
grid off
axis([xlim, ylim, zlim]);
set(gca, 'Box', 'on');

end


%-------------------------------------------%
% MANUAL STEP SIZE FUNCTION
%-------------------------------------------%
function GluR2xyds = UNIFORMSTEPS2(GluR2Ndots, Lx)
%-------------------------------------------%

   Lx(1:2,1:GluR2Ndots) = Lx;
   xyd = randi([0 1],GluR2Ndots,2)';
   xyd(xyd == 0) = -1;
   GluR2xyds = (Lx.*xyd);
   
end


%-------------------------------------------%
% MSD SCALED STEPS FUNCTION
%-------------------------------------------%
function [GluR2xyl GluR2xyds] = SCALEUNIFORMSTEPS2(GluR2Ndots, GluR2xyds, GluR2xyl, Ls, MSDtest)


if MSDtest(1)
	Ls = 1;	
end
	
	for j = 1:GluR2Ndots
        GluR2xyds(:,j) = GluR2xyds(:,j)*Ls;
        GluR2xyl(:,j) = GluR2xyl(:,j)+GluR2xyds(:,j);
	end	
	
end


%-------------------------------------------%
% MSD TRACKS GENERATOR
%-------------------------------------------%
function [tracks] = MSDfun2(stepN, Nsteps, tracks, GluR2xyds)
    time = (0:Nsteps-1)';
    xymsd = GluR2xyds';
    xymsd = cumsum(xymsd,1);
    tracks{stepN} = [time xymsd];

end


%-------------------------------------------%
% MSD TRACKS ANALYSIS
%-------------------------------------------%
function [Dest] = MSDfunction2(tracks,GluR2Ndots,Nsteps,D,Dn,L,dT,k,Scale,MSDtest)


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
N_PARTICLES = GluR2Ndots;
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









%%					CODE GRAVE YARD
%-------------##########################------------------%
%			   UNUSED OR ARCHIVED CODE
%-------------##########################------------------%

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
GluR1xyl,SAPFmx1,SAPFmx2)

    
GluR1Sd1 = ones(size(GluR1xyl));
GluR1Sd2 = GluR1Sd1;

GluR1SLOCMXa = zeros(1,size(GluR1Sd1,2));
GluR1SLOCMXb = GluR1SLOCMXa;

[m,p,n] = size(PSD1Sxvecs);
	% currentcell = zeros(1,size(GluR1xyl,2));
	currentcell = zeros(1,p);
	

for m = 1:p
	currentcell = currentcell+1;
    
PSD1SxvNOW = PSD1Sxvecs{m};
PSD1SyvNOW = PSD1Syvecs{m};
PSD2SxvNOW = PSD2Sxvecs{m};
PSD2SyvNOW = PSD2Syvecs{m};

G1INPSD1 = inpolygon(GluR1xyl(1,:)',GluR1xyl(2,:)',PSD1SxvNOW,PSD1SyvNOW);
G1INPSD2 = inpolygon(GluR1xyl(1,:)',GluR1xyl(2,:)',PSD2SxvNOW,PSD2SyvNOW);
 
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
GluR2xyl,SAPFmx1,SAPFmx2)
 
    
GluR2Sd1 = ones(size(GluR2xyl));
GluR2Sd2 = GluR2Sd1;
 
GluR2SLOCMXa = zeros(1,size(GluR2Sd1,2));
GluR2SLOCMXb = GluR2SLOCMXa;
 
[m,p,n] = size(PSD1Sxvecs);
    % currentcell = zeros(1,size(GluR2xyl,2));
	currentcell = zeros(1,p);
	
for m = 1:p
    currentcell = currentcell+1;
    
PSD1SxvNOW = PSD1Sxvecs{m};
PSD1SyvNOW = PSD1Syvecs{m};
PSD2SxvNOW = PSD2Sxvecs{m};
PSD2SyvNOW = PSD2Syvecs{m};
 
G2INPSD1 = inpolygon(GluR2xyl(1,:)',GluR2xyl(2,:)',PSD1SxvNOW,PSD1SyvNOW);
G2INPSD2 = inpolygon(GluR2xyl(1,:)',GluR2xyl(2,:)',PSD2SxvNOW,PSD2SyvNOW);
 
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
function [GluR1xyds GluR1_TdwellPSD GluR1_TdwellPERI GluR1_TdwellSPYN G1FSLOTS]...
	= G1SLOTS(stepN,GluR1Ndots,GluR1xyds,...
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

	GluR1xyds(:,GR1PROBPRF1) = GluR1xyds(:,GR1PROBPRF1)*(.001);
    GluR1xyds(:,GR1PROBPRF2) = GluR1xyds(:,GR1PROBPRF2)*(.001);
    

    
    % HOW MANY SLOTS ARE NOW FILLED?
    G1PSD1FSLOTS = numel(GR1PROBPRF1);
    G1PSD2FSLOTS = numel(GR1PROBPRF2);

	


G1FSLOTS = [G1PSD1FSLOTS G1PSD2FSLOTS 0 0];
%================%

% if stepN == 401; keyboard; end

end

%----SLOTS GluR2----%
function [GluR2xyds GluR2_TdwellPSD GluR2_TdwellPERI GluR2_TdwellSPYN G2FSLOTS]... 
	= G2SLOTS(stepN,GluR2Ndots,GluR2xyds,...
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
 
    GluR2xyds(:,GR2PROBPRF1) = GluR2xyds(:,GR2PROBPRF1)*(.001);
    GluR2xyds(:,GR2PROBPRF2) = GluR2xyds(:,GR2PROBPRF2)*(.001);
    
 
    
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
function [GluR1xyds GluR1_TdwellPSD GluR1_TdwellPERI G1FSLOTS...
	] = G1PERISLOTS(stepN,GluR1Ndots,GluR1xyds,...
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

	GluR1xyds(:,GR1PROBPRF1) = GluR1xyds(:,GR1PROBPRF1)*(.001);
    GluR1xyds(:,GR1PROBPRF2) = GluR1xyds(:,GR1PROBPRF2)*(.001);
    

    
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
    
    GluR1xyds(:,SLPERI1LOCkon) = GluR1xyds(:,SLPERI1LOCkon)*(.001);
    GluR1xyds(:,SLPERI2LOCkon) = GluR1xyds(:,SLPERI2LOCkon)*(.001);
    GluR1xyds(:,SLPERI1LOCkoff) = GluR1xyds(:,SLPERI1LOCkoff)/(.001);
    GluR1xyds(:,SLPERI2LOCkoff) = GluR1xyds(:,SLPERI2LOCkoff)/(.001);

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
function [GluR2xyds GluR2_TdwellPSD GluR2_TdwellPERI G2FSLOTS...
    ] = G2PERISLOTS(stepN,GluR2Ndots,GluR2xyds,...
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
 
    GluR2xyds(:,GR2PROBPRF1) = GluR2xyds(:,GR2PROBPRF1)*(.001);
    GluR2xyds(:,GR2PROBPRF2) = GluR2xyds(:,GR2PROBPRF2)*(.001);
    
 
    
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
    
    GluR2xyds(:,SLPERI1LOCkon) = GluR2xyds(:,SLPERI1LOCkon)*(.001);
    GluR2xyds(:,SLPERI2LOCkon) = GluR2xyds(:,SLPERI2LOCkon)*(.001);
    GluR2xyds(:,SLPERI1LOCkoff) = GluR2xyds(:,SLPERI1LOCkoff)/(.001);
    GluR2xyds(:,SLPERI2LOCkoff) = GluR2xyds(:,SLPERI2LOCkoff)/(.001);
 
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
function [GluR1xyds GluR1_TdwellPSD GluR1_TdwellPERI G1FSLOTS...
	] = G1SLOTSOLD(stepN,GluR1Ndots,GluR1xyds,...
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

	GluR1xyds(:,GR1PROBPRF1) = GluR1xyds(:,GR1PROBPRF1)*(.001);
    GluR1xyds(:,GR1PROBPRF2) = GluR1xyds(:,GR1PROBPRF2)*(.001);
    

    
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
    
    GluR1xyds(:,SLPERI1LOCkon) = GluR1xyds(:,SLPERI1LOCkon)*(.001);
    GluR1xyds(:,SLPERI2LOCkon) = GluR1xyds(:,SLPERI2LOCkon)*(.001);
    GluR1xyds(:,SLPERI1LOCkoff) = GluR1xyds(:,SLPERI1LOCkoff)/(.001);
    GluR1xyds(:,SLPERI2LOCkoff) = GluR1xyds(:,SLPERI2LOCkoff)/(.001);

	G1PERI1FSLOTS = numel(SLPERI1LOCkon)-numel(SLPERI1LOCkoff);
	G1PERI2FSLOTS = numel(SLPERI2LOCkon)-numel(SLPERI2LOCkoff);
%-------------------%   
end %useSlotsPSDGR1
%===================%
 

G1FSLOTS = [G1PSD1FSLOTS G1PSD2FSLOTS G1PERI1FSLOTS G1PERI2FSLOTS];
%================%
end

%----OLD SLOTS GluR2----%
function [GluR2xyds GluR2_TdwellPSD GluR2_TdwellPERI G2FSLOTS...
    ] = G2SLOTSOLD(stepN,GluR2Ndots,GluR2xyds,...
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
 
    GluR2xyds(:,GR2PROBPRF1) = GluR2xyds(:,GR2PROBPRF1)*(.001);
    GluR2xyds(:,GR2PROBPRF2) = GluR2xyds(:,GR2PROBPRF2)*(.001);
    
 
    
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
    
    GluR2xyds(:,SLPERI1LOCkon) = GluR2xyds(:,SLPERI1LOCkon)*(.001);
    GluR2xyds(:,SLPERI2LOCkon) = GluR2xyds(:,SLPERI2LOCkon)*(.001);
    GluR2xyds(:,SLPERI1LOCkoff) = GluR2xyds(:,SLPERI1LOCkoff)/(.001);
    GluR2xyds(:,SLPERI2LOCkoff) = GluR2xyds(:,SLPERI2LOCkoff)/(.001);
 
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


