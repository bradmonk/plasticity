function [varargout] = MAINBOX(dot,dr,um,sap,hr,ko,doUse,doRun,doKo,box,slt,stky,GT,GTab)
clc; format compact; format short; close all; SNSZ = get(0,'ScreenSize');


%=========================================================%
%					INPUTS FROM GUI
%---------------------------------------------------------%
doIz = doRun(10);
doprofile = doRun(11);

if doprofile
profile on;
end

%{.
% IF NO GUI INPUTS
if nargin < 1 

dot = [150.00 50.00 3600.00 1000.00 0.10 1.00];
dr = [0.130 0.100 0.050 0.040 0.130 0.100 0.050 0.040];
um = [3.000 6.000 0.400 0.400 0.200 0.200];
sap = [7.000 7.000 50.000 1.800 1.300 50.000 1.800 1.300];
hr = [30 60];
ko = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
doUse = [1 1 1 1];
doRun = [0 0 0 0 1 0 0];
doKo = [1 1 1 1 0 0];
box = [30 200 5 0];
slt = [1 1 1 1];
stky = [1.500 0.200 4.000 0.400 6.900 5.800 4.000 0.000 2.000...
	0.000 5.800 5.800 1800.000 1.000 2100.000 1.000 1800.000...
	1.000 2100.000 1.000 2.000 2.000];

disp('Using internal (non-GUI) parameter constants:');
% sprintf('%.3f|',strcat('dot:', num2str(dot)))
disp('dot...');sprintf('%.3f|',dot)
disp('dr...');sprintf('%.3f|',dr)
disp('um...');sprintf('%.3f|',um)
disp('sap...');sprintf('%.3f|',sap)
disp('hr...');sprintf('%.0f|',hr)
disp('ko...');sprintf('%.0f|',ko)
disp('doUse...');sprintf('%.0f|',doUse)
disp('doRun...');sprintf('%.0f|',doRun)
disp('doKo...');sprintf('%.0f|',doKo)
disp('box...');sprintf('%.0f|',box)
disp('slt...');sprintf('%.0f|',slt)
disp('stky...');sprintf('%.3f|',stky)

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

sap21 = sap(21); % LonS1
sap22 = sap(22); % BonS1
sap23 = sap(23); % RonS1
sap24 = sap(24); % LoffS1
sap25 = sap(25); % BoffS1
sap26 = sap(26); % RoffS1
sap27 = sap(27); % LonS2
sap28 = sap(28); % BonS2
sap29 = sap(29); % RonS2
sap30 = sap(30); % LoffS2
sap31 = sap(31); % BoffS2
sap32 = sap(32); % RoffS2


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
doRun8 = doRun(8); % doFieldFig
doRun9 = doRun(9); % doSlotColormap

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


SAPPADPSD1 = stky(21);
SAPPADPSD2 = stky(22);


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

LTP1onG1 = stky(13);
LTP1offG1 = stky(15);
LTP2onG1 = stky(14);
LTP2offG1 = stky(16);

LTP1onG2 = stky(17);
LTP1offG2 = stky(19);
LTP2onG2 = stky(18);
LTP2offG2 = stky(20);


GT1on = GT(1);
GT1off = GT(2);
GT1LTPv = GT(3);
GT2on = GT(4);
GT2off = GT(5);
GT2LTPv = GT(6);
GT1masktab = GTab{1};
GT2masktab = GTab{2};





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
	% Nsteps = 300;
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
%			PARTICLE BOX EXIT SIMULATION (EXITBOX)
%---------------------------------------------------------%
if doRun7

pbox = [box1 box2 box3 box4 dot5 dot4 dot1 dot2...
    DGR1spy kGR1 DGR1psd DGR2spy kGR2 DGR2psd LsGR1psd LsGR2psd...
	doKo1 doKo2 ko9 ko10 ko1 ko2...
	sap1 sap2 slt1 slt2 slt3 slt4];

varargin = EXITBOX(pbox,um);
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
] = FieldFun(doRun,fsizeX, fsizeY,PSD1size, PSD2size, periPSD1size, periPSD2size);
%---##########################---%





%%
%-------------##########################------------------%
%              CLUSTER MODEL PARAMETERS
%=========================================================%
%==================================%
mask3x3 = [1 1 1; 1 1 1; 1 1 1];
oVx = round(.5:.5:8);
oVx = [oVx;oVx];
rVx = [oVx;oVx+8;oVx+16;oVx+24;oVx+32;oVx+40;oVx+48;oVx+56];
%==================================%
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
[Nsteps,S1,S2,S1sum,S2sum,...
SC1beta,SC1mu,SC1r,SC1ro,SC1tau,...
SC2beta,SC2mu,SC2r,SC2ro,SC2tau...
] = SAPSLOTSETUP(Nsteps,...
sap,um,dot,SAPPADPSD1,SAPPADPSD2);
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
S1sz = numel(S1(:,1));
S2sz = numel(S2(:,1));
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
% G2SP1 = zeros(4,100);
% G2SP2 = zeros(4,100);
% G1LTParray = ones(1,numel(S1L2_array));
% G2LTParray = ones(1,numel(S1L2_array));
%---------------------------------------------%
end % if runSAPPSD1 || runSAPPSD2
%=============================================%

%---##########################---%
%		SLOT GAUSSIAN COLORMAP
%---##########################---%
if doRun(9)
SLOTMAPFUN(G1STBASE,G1RTBASE,G1STLTP,G1RTLTP,G1BSMu,G1LSMu,...
G2STBASE,G2RTBASE,G2STLTP,G2RTLTP,G2BSMu,G2LSMu,GT,GTab);

stepN = 1;
G1P1SAPM = [1:3;4:6;7:9];
	   
SAPRECFUN(G1STBASE,G1RTBASE,G1STLTP,G1RTLTP,G1BSMu,G1LSMu,...
G2STBASE,G2RTBASE,G2STLTP,G2RTLTP,G2BSMu,G2LSMu,...
S1,S1sz, sap, stepN,GT,GTab,...
G1P1SAPM, G1SP1, LTP1onG1, LTP1offG1, LTP2onG1, LTP2offG1);

G1P1SAPM = [];
end
%---##########################---%



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


xlim = [0 XWIDE];
ylim = [-YHIGH YHIGH];

subplot(5,5,[3 25]), 
G2Ph1 = scatter(GluR2xyl(1,:),GluR2xyl(2,:),5,[0 0 1]);
hold on;
subplot(5,5,[3 25]), 
G1Ph1 = scatter(GluR1xyl(1,:),GluR1xyl(2,:),5,[1 0 0]);
axis([xlim, ylim]);
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(get(gca,'XLabel'),'String','GluR2:Blue \bullet \bullet GluR1:Red')
rectangle('Position',[XYBOXpr1(1),XYBOXpr1(2),XYBOXpr1(3),XYBOXpr1(4)])
rectangle('Position',[XYBOXpr2(1),XYBOXpr2(2),XYBOXpr2(3),XYBOXpr2(4)])
rectangle('Position',[XYBOXp1(1),XYBOXp1(2),XYBOXp1(3),XYBOXp1(4)])
rectangle('Position',[XYBOXp2(1),XYBOXp2(2),XYBOXp2(3),XYBOXp2(4)])
hold off;





%-- S1 Cluster Plots --
subplot(5,5,1); S1Ph1 = imagesc(S1);
colormap('hot');
title('SYNAPSE-1 CLUSTERS');

subplot(5,5,6),
S1Ph2 = imagesc(S1);
title('hk');

subplot(5,5,11),
S1Ph3 = imagesc(S1);
title('Pkon');

subplot(5,5,16),
S1Ph4 = imagesc(S1);
title('Sno');

subplot(5,5,21),
S1Ph5 = imagesc(S1);
title('TETHERED G1');

%-- S2 Cluster Plots --
subplot(5,5,2); 
S2Ph1 = imagesc(S2);
colormap('hot');
title('SYNAPSE-2 CLUSTERS');


subplot(5,5,7),
S2Ph2 = imagesc(S2);
title('hk');


subplot(5,5,12),
S2Ph3 = imagesc(S2);
title('Pkon');


subplot(5,5,17),
S2Ph4 = imagesc(S2);
title('Sno');


subplot(5,5,22),
S2Ph5 = imagesc(S2);
title('TETHERED G1');


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
	SAPFmx1 SAPFmx2 G1FSLOTS G2FSLOTS...
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
	mask3x3,rVx,SSPSD1,SSPSD2,GluR1SdP1,GluR1SdP2,...
	GluR1SLOCMXa,GluR1SLOCMXb,g1polyN1,g2polyN1,...
	GluR2SdP1,GluR2SdP2,GluR2SLOCMXa,GluR2SLOCMXb,...
	LTP1onG1,LTP1offG1,LTP2onG1,LTP2offG1,...
	LTP1onG2,LTP1offG2,LTP2onG2,LTP2offG2,...
	XYLBp1,XYRTp1,XYLBp2,XYRTp2,...
	G1STBASE,G1STLTP,G1BSMu,G1LSMu,...
	G2STBASE,G2STLTP,G2BSMu,G2LSMu);
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
	[S1 G1SP1] = S1_MainClusterFun(S1Ph1, S1Ph2, S1Ph3, S1Ph4, S1Ph5,...
	S1, S1sz, sap, stepN, G1P1SAPM, G1SP1,GT,GTab);
	end
	%----S2_MainClusterFun----%
	if runSAPPSD2
	[S2 G1SP2] = S2_MainClusterFun(S2Ph1, S2Ph2, S2Ph3, S2Ph4, S2Ph5,...
    S2, S2sz, sap, stepN, G1P2SAPM, G1SP2,GT,GTab);
	end
		
	%----S1sum & S2sum----%
	S1sum = sum(S1(:)); S2sum = sum(S2(:));
	%----PLOTS1S2----%
	% if doMainPlot; if mod(stepN, 10) == 0;PLOTS1S2(S1, S2);end;end
	% if ~doMainPlot; if mod(stepN, 100) == 0;PLOTS1S2(S1, S2);end;end
	%-----------------------------------%
	end				% ?? doClusters ??
	%===================================%

	% if stepN == 2000; keyboard; end
	
	
	%=============================================%
    %		MOVE PARTICLES MAIN FUNCTION
	%=============================================%
	[GluR1xyds GluR1xyl GluR2xyds GluR2xyl]...
    = MOVEGLUR(stepN,XWIDE,YHIGH,...
	GluR1Ndots,GluR1xyds,GluR1xyl,...
    GluR2Ndots,GluR2xyds,GluR2xyl);
	
	
	
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
	MAINPLOT(G1Ph1, G2Ph1,stepN,GluR2xyl,GluR1xyl,...
	XWIDE,YHIGH,XYBOXpr1,XYBOXpr2,XYBOXp1,XYBOXp2);
	end

	if ~doMainPlot && ~doONEDOTPLOT
	if mod(stepN, 200) == 0
	MAINPLOT(G1Ph1, G2Ph1,stepN,GluR2xyl,GluR1xyl,...
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
       
       %snapnow;
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

if doprofile
profile viewer;
end

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
BMData1= reDATAdatasetFull;
BMData2= reDATASAPdataFull;
BMData3= reDATADdataFull;
BMData4= reDATAAMPARdataFull;
BMData5= reDATAGluRdataFull;
BMData6= reDATAG1SLOTSDATAFull;
BMData7= reDATAG2SLOTSDATAFull;







%{
% UNCOMMENT TO REINITIALIZE DATA

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





%%
%---------------
if ~doRun(3);
%---------------

ChopStart = 1;
if ChopStart

% chop = 3;
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

end		% if ChopStart


% for p = 1:numel(reDATAdataset)
% reDATAdataset = smooth(reDATAdataset{p},5);
% end

%%
%======================================================================%
%					SAVE FIGURE 1 PLOTS
%======================================================================%
set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,'STARShiP2','png');
saveas(gcf, ['outputfigs/FIGURE1.png']);
%======================================================================%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%======================================================================%
%					FINAL OUTPUT (SPLINE) FIGURE 1 OF 2
%======================================================================%
%							FIGURE SETUP
%----------------------------------------------------------------------%
%----------------------------------------------------------------------%
fig21 = figure(21);
set(21,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig21,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
cOLOR = [c1; c2; c3; c4];


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
%======================================================================%



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
%======================================================================%
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
%======================================================================%



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
%======================================================================%





%===========================================================%
% FIG1 BOTTOM RIGHT: Synapse Particle Counts
%===========================================================%
sbpos = [.55 .09 .4 .38]; ptype = 4;
cOLOR = [c1; c2; c1; c2; c3; c2; c3; c4];
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
%======================================================================%


%======================================================================%
set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,'STARShiP1','png');
saveas(gcf, ['outputfigs/STARShiP1.png']);
% printpreview
%======================================================================%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%======================================================================%
%					FINAL OUTPUT (SPLINE) FIGURE 2 OF 2
%======================================================================%
%							FIGURE SETUP
%----------------------------------------------------------------------%
fig12 = figure(12);
set(12,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig12,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% WAS FIG1 BOTTOM LEFT: GluR Subtypes in Synapse 1v2
%===========================================================%
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
set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,'STARShiP2','png');
saveas(gcf, ['outputfigs/STARShiP2.png']);

%{
% FName = gcf;
% OName = 'STARShiP';
% save(['/Users/bradleymonk/Documents/MatLab/BradsModel/outputfigs/'...
% 	OName '.mat'],...
% 	'BMData1', 'BMData2', 'BMData3', 'BMData4', 'BMData5', 'BMData6',...
% 	'AveOver', 'DATARATE','t');
% saveas(FName, ['synth1/memo/' OName '.png']);

% AveOver = 10;
% DATARATE = 100;
% t = .1;
% BMData1 = reDATAdataset;
% BMData2 = reDATASAPdata;
% BMData3 = reDATADdata;
% BMData4 = reDATAAMPARdata;
% BMData5 = reDATAGluRdata;
% BMData6 = reDATAG1SLOTSDATA;
% BMData7 = reDATAG2SLOTSDATA;
% saveas(gcf,...
% '/Users/bradleymonk/Documents/MatLab/BradsModel/outputfigs/STARShiP2',...
% 'png');
% fname = 'D:\path1\path2';
% saveas(gca, fullfile(fname, filename), 'jpeg');
% set(gcf, 'PaperPositionMode', 'auto');
% print -depsc2 STARShiP3.eps
% print('-depsc2','-tiff','-r300','STARShiP3')
% fixPSlinestyle('STARShiP3.eps', 'STARShiP3.eps');
% saveas(gcf,'STARShiP3.tiff', 'tiffn')
%}
%======================================================================%


%---------------
end % if ~doRun(3);
%---------------

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%======================================================================%
%					FINAL OUTPUT (NONAVERAGED) FIGURE 1 OF 2
%======================================================================%
%======================================================================%
%{
%======================================================================%
%							FIGURE SETUP
%----------------------------------------------------------------------%
fig2 = figure(2);
set(2,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig2,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%======================================================================%
subplot(2,2,1), subplot('Position',[.055 .57 .4 .38]),...
RecepPerRegFIG1=plot([dataset(:,3) dataset(:,4) dataset(:,5) dataset(:,2)]);
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
F2AT = title('Receptors Per Region'); leg1=legend('Synapse 1', 'Synapse 2', 'Synapse Total', 'ES');
set(F2AT, 'FontSize', 14);
set(leg1,'Location','NorthWest');
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
MS1 = 5; 
MS2 = 5;
set(RecepPerRegFIG1(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(RecepPerRegFIG1(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS2,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(RecepPerRegFIG1(3),'LineStyle','-','Color', c3,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
set(RecepPerRegFIG1(4),'LineStyle','-','Color',c4,'LineWidth',1,...
'Marker','o','MarkerSize',MS2,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('Extrasynaptic vs Synaptic AMPARs');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%======================================================================%


%======================================================================%
subplot(2,2,2), subplot('Position',[.55 .77 .4 .18]),...
G1PSDvPSAFIG1=plot([GluRdata(:,1) GluRdata(:,2)]);
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
GluR1PSDvsPERI=title('GluR1 in PSD vs PSA'); leg1=legend('PSD', 'PSA');
set(GluR1PSDvsPERI, 'FontSize', 14);
set(leg1,'Location','NorthWest');
set(gca,'XTickLabel',{'-'})
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
MS1 = 5; 
MS2 = 5;
set(G1PSDvPSAFIG1(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(G1PSDvPSAFIG1(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
% set(XXXXXXXX(3),'LineStyle','-','Color', c3,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
% set(XXXXXXXX(4),'LineStyle','-','Color',c4,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
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
%===========================================================%
subplot('Position',[.55 .57 .4 .18]),...
G1PSDvPSAFIG2=plot([GluRdata(:,3) GluRdata(:,4)]);
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
GluR2PSDvsPERI=title('GluR2 PSD vs PSA'); leg1=legend('PSD', 'PSA');
set(GluR2PSDvsPERI, 'FontSize', 14);
set(leg1,'Location','NorthWest');
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
MS1 = 5; 
MS2 = 5;
set(G1PSDvPSAFIG2(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(G1PSDvPSAFIG2(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
% set(XXXXXXXX(3),'LineStyle','-','Color', c3,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
% set(XXXXXXXX(4),'LineStyle','-','Color',c4,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
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
%======================================================================%


%======================================================================%
subplot(2,2,3), subplot('Position',[.05 .09 .4 .38]),...
Spy1vSpy2FIG1=plot([GluRdata(:,5) GluRdata(:,6) GluRdata(:,7) GluRdata(:,8)]);
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
F2CT = title('GluR Subtypes in Synapse 1v2'); 
set(F2CT, 'FontSize', 14);
leg1=legend('G1/2 Synapse 1', 'G1/2 Synapse 2',	'G2/3 Synapse 1', 'G2/3 Synapse 2');
set(leg1,'Location','NorthWest');
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
MS1 = 5; 
MS2 = 5;
set(Spy1vSpy2FIG1(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(Spy1vSpy2FIG1(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(Spy1vSpy2FIG1(3),'LineStyle','-','Color', c3,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
set(Spy1vSpy2FIG1(4),'LineStyle','-','Color',c4,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('GluR Subtypes in Synapse 1v2');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%======================================================================%


%======================================================================%
subplot(2,2,4), subplot('Position',[.55 .09 .4 .38]),...
SpyTotalsFIG1=plot([dataset(:,3) dataset(:,4) dataset(:,5)]);
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
F2DT=title('Synapse Totals'); leg2=legend('Synapse 1', 'Synapse 2', 'Synapse Total');
set(F2DT, 'FontSize', 14);
set(leg2,'Location','SouthWest');
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
MS1 = 5; 
MS2 = 5;
set(SpyTotalsFIG1(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(SpyTotalsFIG1(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(SpyTotalsFIG1(3),'LineStyle','-','Color', c3,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
% set(XXXXXXXX(4),'LineStyle','-','Color',c4,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('Synapse Particle Counts');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%======================================================================%

%======================================================================%
set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 NonMeanOutput1.eps
fixPSlinestyle('NonMeanOutput1.eps', 'NonMeanOutput1.eps');
%======================================================================%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%======================================================================%
%					FINAL OUTPUT (NONAVERAGED) FIGURE 2 OF 2
%======================================================================%
%							FIGURE SETUP
%----------------------------------------------------------------------%
fig3 = figure(3);
set(3,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig3,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%======================================================================%
subplot(2,2,1), subplot('Position',[.05 .57 .4 .38]),...
AMPARTinSpyFIG1=plot([AMPARdata(:,1) AMPARdata(:,2)]);
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','Particles')
F3AT=title('AMPAR subtype totals in synapses'); 
leg3=legend('GluR1', 'GluR2');
set(F3AT, 'FontSize', 14);
set(leg3,'Location','NorthWest');
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
MS1 = 5; 
MS2 = 5;
set(AMPARTinSpyFIG1(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(AMPARTinSpyFIG1(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
% set(XXXXXXXX(3),'LineStyle','-','Color', c3,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
% set(XXXXXXXX(4),'LineStyle','-','Color',c4,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('GluR Subtypes in Synapses');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%======================================================================%


%======================================================================%
% subplot(2,2,2), 
subplot('Position',[.55 .77 .4 .18]),...
G1DotsInSPYFIG1=plot([GluRdata(:,5) GluRdata(:,6)]);
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'YLabel'),'String','GluR1 Particles')
F3B2T=title('GluR1');
set(F3B2T, 'FontSize', 14);
leg1=legend('Synapse 1', 'Synapse 2');
set(leg1,'Location','NorthWest');
set(gca,'XTickLabel',{'-'})
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
MS1 = 5; 
MS2 = 5;
set(G1DotsInSPYFIG1(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(G1DotsInSPYFIG1(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
% set(XXXXXXXX(3),'LineStyle','-','Color', c3,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
% set(XXXXXXXX(4),'LineStyle','-','Color',c4,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('GluR Subtypes in Synapse 1v2');
hXLabel = xlabel('-');
hYLabel = ylabel('GluR1 Particles');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%===========================================================%
subplot('Position',[.55 .57 .4 .18]),...
G2DotsInSPYFIG2=plot([GluRdata(:,7) GluRdata(:,8)]);
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'YLabel'),'String','GluR2 Particles')
F3B1T=title('GluR2 Particles');
set(F3B1T, 'FontSize', 14);
leg1=legend('Synapse 1', 'Synapse 2');
set(leg1,'Location','NorthWest');
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
MS1 = 5; 
MS2 = 5;
set(G2DotsInSPYFIG2(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(G2DotsInSPYFIG2(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
% set(XXXXXXXX(3),'LineStyle','-','Color', c3,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
% set(XXXXXXXX(4),'LineStyle','-','Color',c4,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('-');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('GluR2 Particles');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%======================================================================%


%======================================================================%
% subplot(2,2,3), 
subplot('Position',[.55 .29 .4 .18]),...
%G1FSLOTSFIG1=plot([G1SLOTSDATA(:,1) G1SLOTSDATA(:,2) G1SLOTSDATA(:,3) G1SLOTSDATA(:,4)]);
G1FSLOTSFIG1=plot([G1SLOTSDATA(:,1) G1SLOTSDATA(:,2)]);
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'YLabel'),'String','Occupied Slots')
F3B1T=title('Slots Occupied by GluR1');
set(F3B1T, 'FontSize', 14);
leg1=legend('PSD 1', 'PSD 2');
set(leg1,'Location','NorthWest');
set(gca,'XTickLabel',{'-'})
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
MS1 = 5; 
MS2 = 5;
set(G1FSLOTSFIG1(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(G1FSLOTSFIG1(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
% set(G1FSLOTSFIG1(3),'LineStyle','-','Color', c3,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
% set(G1FSLOTSFIG1(4),'LineStyle','-','Color',c4,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('Occupied Synaptic Slots');
hXLabel = xlabel('-');
hYLabel = ylabel('Particles');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%===========================================================%
sbpos = [.55 .09 .4 .18];
subplot('Position',[.55 .09 .4 .18]),...
G2FSLOTSFIG1=plot([G2SLOTSDATA(:,1) G2SLOTSDATA(:,2)]);
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','Occupied Slots')
F3B2T=title('Slots Occupied by GluR2');
set(F3B2T, 'FontSize', 14);
leg1=legend('PSD 1', 'PSD 2');
set(leg1,'Location','NorthWest');
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
MS1 = 5; 
MS2 = 5;
set(G2FSLOTSFIG1(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(G2FSLOTSFIG1(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
% set(G2FSLOTSFIG1(3),'LineStyle','-','Color', c3,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
% set(G2FSLOTSFIG1(4),'LineStyle','-','Color',c4,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('-');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);

%======================================================================%


%======================================================================%
subplot('Position',[.055 .09 .4 .38]),...
SAPCLPTFIG1=plot([SAPdata(:,1) SAPdata(:,2)]);
xt = (get(gca,'XTick'))*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','SAP Particles')
F3CT=title('Synaptic SAP Cluster Size'); leg4=legend('Synapse 1', 'Synapse 2');
set(F3CT, 'FontSize', 14);
set(leg4,'Location','NorthWest');
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
MS1 = 5; 
MS2 = 5;
set(SAPCLPTFIG1(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(SAPCLPTFIG1(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
% set(XXXXXXXX(3),'LineStyle','-','Color', c3,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
% set(XXXXXXXX(4),'LineStyle','-','Color',c4,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('Synaptic SAP Cluster Size');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%======================================================================%

%======================================================================%
set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 NonMeanOutput2.eps
fixPSlinestyle('NonMeanOutput2.eps', 'NonMeanOutput2.eps');
%======================================================================%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%======================================================================%
%					AVERAGED OUTPUTS FIGURE 1 OF 2
%======================================================================%
%							FIGURE SETUP
%----------------------------------------------------------------------%
fig4 = figure(4);
set(4,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig4,'OuterPosition',pos1)
set(gcf,'Color',[.95,.95,.95])
%----------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%======================================================================%
subplot(2,2,1),subplot('Position',[.055 .57 .4 .38]),...
RPRFIG1 = plot([DATAdataset(:,3) DATAdataset(:,4) DATAdataset(:,5) DATAdataset(:,2)]);
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
F2AT = title('Receptors Per Region'); leg1=legend('Synapse 1', 'Synapse 2', 'Synapse Total', 'ES');
set(F2AT, 'FontSize', 14);
set(leg1,'Location','NorthWest');
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
MS1 = 5; 
MS2 = 5;
set(RPRFIG1(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(RPRFIG1(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(RPRFIG1(3),'LineStyle','-','Color', c3,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
set(RPRFIG1(4),'LineStyle','-','Color',c4,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('Extrasynaptic vs Synaptic AMPARs');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%======================================================================%




%======================================================================%
% subplot(2,2,2), 
subplot('Position',[.55 .77 .4 .18]),...
G1PSvPAF1=plot([DATAGluRdata(:,1) DATAGluRdata(:,2)]);
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
GluR1PSDvsPERI=title('GluR1 in PSD vs PSA'); leg1=legend('PSD', 'PSA');
set(GluR1PSDvsPERI, 'FontSize', 14);
set(leg1,'Location','NorthWest');
set(gca,'XTickLabel',{'-'})
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
% MS1 = 6; 
% MS2 = 6;
set(G1PSvPAF1(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(G1PSvPAF1(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
% set(XXXXXXXX(3),'LineStyle','-','Color', c3,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
% set(XXXXXXXX(4),'LineStyle','-','Color',c4,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
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
%===========================================================%
subplot('Position',[.55 .57 .4 .18]),...
G2PSvPAF1=plot([DATAGluRdata(:,3) DATAGluRdata(:,4)]);
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
GluR2PSDvsPERI=title('GluR2 PSD vs PSA'); leg1=legend('PSD', 'PSA');
set(GluR2PSDvsPERI, 'FontSize', 14);
set(leg1,'Location','NorthWest');
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
% MS1 = 6; 
% MS2 = 6;
set(G2PSvPAF1(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(G2PSvPAF1(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
% set(XXXXXXXX(3),'LineStyle','-','Color', c3,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
% set(XXXXXXXX(4),'LineStyle','-','Color',c4,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
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
%======================================================================%


%======================================================================%
% subplot(2,2,3), subplot('Position',[.05 .077 .4 .38]),...
subplot('Position',[.05 .09 .4 .38]),...
GR1SpyF1=plot([DATAGluRdata(:,5) DATAGluRdata(:,6) DATAGluRdata(:,7) DATAGluRdata(:,8)]);
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','Particles')
F2CT = title('GluR Subtypes in Synapse 1v2'); 
set(F2CT, 'FontSize', 14);
leg1=legend('G1/2 Synapse 1', 'G1/2 Synapse 2',	'G2/3 Synapse 1', 'G2/3 Synapse 2');
set(leg1,'Location','NorthWest');
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
% MS1 = 6; 
% MS2 = 6;
set(GR1SpyF1(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(GR1SpyF1(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(GR1SpyF1(3),'LineStyle','-','Color', c3,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
set(GR1SpyF1(4),'LineStyle','-','Color',c4,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('GluR Subtypes in Synapse 1v2');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%======================================================================%


%======================================================================%
% subplot(2,2,4), subplot('Position',[.55 .077 .4 .38]),...
subplot('Position',[.55 .09 .4 .38]),...
SpyTF1=plot([DATAdataset(:,3) DATAdataset(:,4) DATAdataset(:,5)]);
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
F2DT=title('Synapse Totals'); leg2=legend('Synapse 1', 'Synapse 2', 'Synapse Total');
set(F2DT, 'FontSize', 14);
set(leg2,'Location','SouthWest');
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
% MS1 = 6; 
% MS2 = 6;
set(SpyTF1(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(SpyTF1(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(SpyTF1(3),'LineStyle','-','Color', c3,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
% set(XXXXXXXX(4),'LineStyle','-','Color',c4,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('Synapse Particle Counts');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%======================================================================%


%======================================================================%
set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 MeanOutput1.eps
fixPSlinestyle('MeanOutput1.eps', 'MeanOutput2.eps');
%======================================================================%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%======================================================================%
%					AVERAGED OUTPUTS FIGURE 2 OF 2
%======================================================================%
%							FIGURE SETUP
%----------------------------------------------------------------------%
fig5 = figure(5);
set(5,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig5,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%======================================================================%
subplot(2,2,1),subplot('Position',[.05 .57 .4 .38]),...
GRSubsSpinesFig1=plot([DATAAMPARdata(:,1) DATAAMPARdata(:,2)]);
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','PARTICLES')
F3AT=title('AMPAR SUBTYPE TOTALS IN SynapseS'); 
leg3=legend('GluR1/2', 'GluR2/3');
set(F3AT, 'FontSize', 14);
set(leg3,'Location','NorthWest');
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
MS1 = 5; 
MS2 = 5;
set(GRSubsSpinesFig1(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(GRSubsSpinesFig1(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
% set(XXXXXXXX(3),'LineStyle','-','Color', c3,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
% set(XXXXXXXX(4),'LineStyle','-','Color',c4,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('GluR Subtypes in Synapses');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%======================================================================%


%======================================================================%
subplot('Position',[.55 .765 .4 .18]),...
GluInSpineFig1=plot([DATAGluRdata(:,5) DATAGluRdata(:,6)]);
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'YLabel'),'String','GluR1/2 Particles')
F3B2T=title('GluR1');
set(F3B2T, 'FontSize', 14);
leg1=legend('Synapse 1', 'Synapse 2');
set(leg1,'Location','NorthWest');
set(gca,'XTickLabel',{'-'})
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
% MS1 = 6; 
% MS2 = 6;
set(GluInSpineFig1(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(GluInSpineFig1(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
% set(XXXXXXXX(3),'LineStyle','-','Color', c3,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
% set(XXXXXXXX(4),'LineStyle','-','Color',c4,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
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
%===========================================================%
subplot(2,2,3), subplot('Position',[.55 .57 .4 .18]),...
GluInSpineFig2=plot([DATAGluRdata(:,7) DATAGluRdata(:,8)]);
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'YLabel'),'String','GluR2/3 Particles')
F3B1T=title('GluR2/3');
set(F3B1T, 'FontSize', 14);
leg1=legend('Synapse 1', 'Synapse 2');
set(leg1,'Location','NorthWest');
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
% MS1 = 6; 
% MS2 = 6;
set(GluInSpineFig2(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(GluInSpineFig2(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
% set(XXXXXXXX(3),'LineStyle','-','Color', c3,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
% set(XXXXXXXX(4),'LineStyle','-','Color',c4,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
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
%======================================================================%





%======================================================================%
%{
%======================================================================%
subplot('Position',[.55 .265 .4 .19]),...
FilledSlotsFig1=plot([DATAG1SLOTSDATA(:,1) DATAG1SLOTSDATA(:,2) ...
					DATAG1SLOTSDATA(:,3) DATAG1SLOTSDATA(:,4)]);
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'YLabel'),'String','Occupied Slots')
F3B1T=title('Slots Occupied by GluR1');
set(F3B1T, 'FontSize', 14);
leg1=legend('PSD 1', 'PSD 2','PSA 1','PSA 2');
set(leg1,'Location','NorthWest');
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
% MS1 = 6; 
% MS2 = 6;
set(FilledSlotsFig1(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(FilledSlotsFig1(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(FilledSlotsFig1(3),'LineStyle','-','Color', c3,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
set(FilledSlotsFig1(4),'LineStyle','-','Color',c4,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('Occupied Slots');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%===========================================%
subplot('Position',[.55 .055 .4 .19]),...
FilledSlotsFig2=plot([DATAG2SLOTSDATA(:,1) DATAG2SLOTSDATA(:,2) ...
					DATAG2SLOTSDATA(:,3) DATAG2SLOTSDATA(:,4)]);
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','Occupied Slots')
F3B2T=title('Slots Occupied by GluR2');
set(F3B2T, 'FontSize', 14);
leg1=legend('PSD 1', 'PSD 2','PSA 1','PSA 2');
set(leg1,'Location','NorthWest');
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
% MS1 = 6; 
% MS2 = 6;
set(FilledSlotsFig2(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(FilledSlotsFig2(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(FilledSlotsFig2(3),'LineStyle','-','Color', c3,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
set(FilledSlotsFig2(4),'LineStyle','-','Color',c4,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('-');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%======================================================================%
%}
sbpos=[.55 .09 .4 .38];
subplot(2,2,4), subplot('Position',[.55 .09 .4 .38]),...
FilledSlotsF1=plot([DATAG1SLOTSDATA(:,1) DATAG1SLOTSDATA(:,2) ...
					DATAG2SLOTSDATA(:,1) DATAG2SLOTSDATA(:,2)]);
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'YLabel'),'String','Occupied Slots')
F3B1T=title('Slots Occupied by GluR');
set(F3B1T, 'FontSize', 14);
leg1=legend('G1/2 PSD 1', 'G1/2 PSD 2','G2/3 PSD 1', 'G2/3 PSD 2');
set(leg1,'Location','NorthWest');
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
% MS1 = 6; 
% MS2 = 6;
set(FilledSlotsF1(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(FilledSlotsF1(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(FilledSlotsF1(3),'LineStyle','-','Color', c3,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
set(FilledSlotsF1(4),'LineStyle','-','Color',c4,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('Occupied Slots');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%===========================================%













%======================================================================%
subplot(2,2,3), subplot('Position',[.055 .09 .4 .38]),...
SAPNdotsFig2 = plot([DATASAPdata(:,1) DATASAPdata(:,2)]);
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','SAP PARTICLES')
F3CT=title('Synaptic SAP Cluster Size'); leg4=legend('Synapse 1', 'Synapse 2');
set(F3CT, 'FontSize', 14);
set(leg4,'Location','NorthWest');
%------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
% MS1 = 6; 
% MS2 = 6;
set(SAPNdotsFig2(1),'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(SAPNdotsFig2(2),'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
% set(XXXXXXXX(3),'LineStyle','-','Color', c3,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
% set(XXXXXXXX(4),'LineStyle','-','Color',c4,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('Synaptic SAP Cluster Size');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Particles');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%======================================================================%


%======================================================================%
set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 MeanOutput2.eps
fixPSlinestyle('MeanOutput2.eps', 'MeanOutput2.eps');
%======================================================================%
%}
%======================================================================%
%======================================================================%
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

% Get top left corner coordinates for SYN
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
% Function takes SAP S-Cluster Matrix (S1)
% uses mask to count SAPs in each 1x1 PolyBox region
% returns condensed Field-Mx 1/4th of full S-Cluster field, 
% holds values (0-4) of number of SAPs in each PolyBox
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


% This constructs SAPFmx like: 
% SAPFmx = [(1,1) (1,3) (1,5)... (3,1) (3,3)....]
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
	SAPFmx1 SAPFmx2 G1FSLOTS G2FSLOTS...
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
	mask3x3,rVx,SSPSD1,SSPSD2,GluR1SdP1,GluR1SdP2,...
	GluR1SLOCMXa,GluR1SLOCMXb,g1polyN1,g2polyN1,...
	GluR2SdP1,GluR2SdP2,GluR2SLOCMXa,GluR2SLOCMXb,...
	LTP1onG1,LTP1offG1,LTP2onG1,LTP2offG1,...
	LTP1onG2,LTP1offG2,LTP2onG2,LTP2offG2,...
	XYLBp1,XYRTp1,XYLBp2,XYRTp2,...
	G1STBASE,G1STLTP,G1BSMu,G1LSMu,...
	G2STBASE,G2STLTP,G2BSMu,G2LSMu)

% THIS WILL ONLY WORK IF THE SYNAPSE IS A TOTAL SIZE OF 8X8
% IF THE TOTAL PSD+PSA IS NOT==8 THE ENTIRE SAP/SLOT SUBROUTINE 
% WILL MALFUNCTION!!!!!!



%-------------------
% SAPMASKREPORT
%-------------------
%==================================%
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
% THIS WILL ONLY WORK IF THE SAP CLUSTER FIELD IS TOTAL SIZE IS 8X8
% IF THE TOTAL PSD+PSA IS NOT==8 THE ENTIRE SAP/SLOT SUBROUTINE 
% WILL MALFUNCTION!!!!!!


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

%==================================%
%{
%==================================%

Smx = padarray(ones(8),[4 4],0);
fMxSz = size(Smx,1);
pMxSz = fMxSz/2;

mask3x3 = [1 1 1; 1 1 1; 1 1 1];
convnSMx = convn(Smx,mask3x3,'same');

rc = 0;
convnV = [];
for r = 1:2:fMxSz; for c = 1:2:fMxSz
rc = rc+1;
convnV(rc) = mean(mean(convnSMx(r:r+1,c:c+1)));
end; end;

convnRMx = round(reshape(convnV,pMxSz,[]) / 2);

gaussian = fspecial('gaussian');

%==================================%
%}
%==================================%
if runSAPPSD1
%------------
S1oc=(S1>0);
SAPmx1 = convn(S1oc,mask3x3,'same');
convnS1Mx = accumarray(rVx(:), SAPmx1(:));
SAPFmx1 = reshape(convnS1Mx,8,[]) ./8;
%------------
end
%==================================%
if runSAPPSD2
%------------
S2oc=(S2>0);
SAPmx2 = convn(S2oc,mask3x3,'same');
convnS2Mx = accumarray(rVx(:), SAPmx2(:));
SAPFmx2 = reshape(convnS2Mx,8,[]) ./8;
%------------
end
%==================================%



%-------------------%
% old conv
%{

sap_mask2 = [1 1 1; 1 1 1; 1 1 1];

if runSAPPSD1
	SAPmx1_occ=(S1>0);							% SAPmx1_occ 16x16 sap occupied lattice
	% SAPmx1=convn(SAPmx1_occ,sap_mask,'same');	% SAPmx1	 16x16 convolution
	SAPmx1 = round(convn(SAPmx1_occ,sap_mask2,'same')./2);
	SAPFmx1=zeros(SSPSD1);						% SAPFmx1	 8x8 0s preallocation
	for n = 0:(SSPSD1-1)
	for m = 0:(SSPSD1-1)
	SAPFmx1(m+1,n+1) = SAPmx1(1+(m*2),1+(n*2));	% SAPFmx1	 8x8 
	end
	end
end

if runSAPPSD2
	SAPmx2_occ=(S2>0);  
	% SAPmx2=convn(SAPmx2_occ,sap_mask,'same');
	SAPmx2 = round(convn(SAPmx2_occ,sap_mask2,'same')./2);
	SAPFmx2=zeros(SSPSD2);
	for n = 0:(SSPSD2-1)
	for m = 0:(SSPSD2-1)
	SAPFmx2(m+1,n+1) = SAPmx2(1+(m*2),1+(n*2));
	end
	end
end



%==================================%
mask3x3 = [1 1 1; 1 1 1; 1 1 1];
%==================================%
if runSAPPSD1
%------------
S1oc=(S1>0);
fMxS1z = size(S1oc,1);
pMxS1z = fMxS1z/2;
%----------------------------------%
convnS1Mx = convn(S1oc,mask3x3,'same');
loc = 0;
convnS1v = zeros(pMxS1z);
for r = 1:2:fMxS1z; for c = 1:2:fMxS1z
loc = loc+1;
convnS1v(loc) = mean(mean(convnS1Mx(r:r+1,c:c+1)));
end; end;

SAPFmx1 = round(convnS1v / 2) ;
SAPmx1 = convnS1Mx;
%------------
end;%if runSAPPSD1
%==================================%
%==================================%
if runSAPPSD2
%------------
S2oc=(S2>0);
fMxS2z = size(S2oc,1);
pMxS2z = fMxS2z/2;
%----------------------------------%
convnS2Mx = convn(S2oc,mask3x3,'same');
loc = 0;
convnS2v = zeros(pMxS2z);
for r = 1:2:fMxS2z; for c = 1:2:fMxS2z
loc = loc+1;
convnS2v(loc) = mean(mean(convnS2Mx(r:r+1,c:c+1)));
end; end;

SAPFmx2 = round(convnS2v / 2) ;
SAPmx2 = convnS2Mx;
%------------
end;%if runSAPPSD2
%==================================%








%==================================%
mask2x2 = [1 1; 1 1];
%==================================%
if runSAPPSD1
%------------
S1oc=(S1>0);
fMxS1z = size(S1oc,1);
pMxS1z = fMxS1z/2;
%----------------------------------%
loc = 0;
SAPmx1 = zeros(pMxS1z);
for r = 1:2:fMxS1z; for c = 1:2:fMxS1z
loc = loc+1;
SAPmx1(loc) = convn(S1oc(r:r+1,c:c+1),mask2x2,'valid');
end; end;

%------------
end;%if runSAPPSD1
%==================================%
%==================================%
if runSAPPSD2
%------------
S2oc=(S2>0);
fMxS2z = size(S2oc,1);
pMxS2z = fMxS2z/2;
%----------------------------------%
loc = 0;
SAPmx2 = zeros(pMxS2z);
for r = 1:2:fMxS2z; for c = 1:2:fMxS2z
loc = loc+1;
SAPmx2(loc) = convn(S2oc(r:r+1,c:c+1),mask2x2,'valid');
end; end;

%------------
end;%if runSAPPSD2
%==================================%






% if runSAPPSD1
% 	SAPmx1_occ=(S1>0);							% SAPmx1_occ 16x16 sap occupied lattice
% 	SAPmx1=convn(SAPmx1_occ,sap_mask,'same');	% SAPmx1	 16x16 convolution
% 	SAPFmx1=zeros(SSPSD1);						% SAPFmx1	 8x8 0s preallocation
% 	for n = 0:(SSPSD1-1)
% 	for m = 0:(SSPSD1-1)
% 	SAPFmx1(m+1,n+1) = SAPmx1(1+(m*2),1+(n*2));	% SAPFmx1	 8x8 
% 	end
% 	end
% end

% if runSAPPSD2
% 	SAPmx2_occ=(S2>0);  
% 	SAPmx2=convn(SAPmx2_occ,sap_mask,'same');
% 	SAPFmx2=zeros(SSPSD2);
% 	for n = 0:(SSPSD2-1)
% 	for m = 0:(SSPSD2-1)
% 	SAPFmx2(m+1,n+1) = SAPmx2(1+(m*2),1+(n*2));
% 	end
% 	end
% end

%}
%-------------------%

if ~runSAPPSD1
	SAPmx1_occ=(S1>0);  
	SAPmx1=round(convn(SAPmx1_occ,sap_mask2,'same')./2);
	SAPFmx1=zeros(SSPSD1);
	SAPFmx1 = SAPFmx1.*0;
end

if ~runSAPPSD2
	SAPmx2_occ=(S2>0);  
	SAPmx2=round(convn(SAPmx2_occ,sap_mask2,'same')./2);
	SAPFmx2=zeros(SSPSD2);
	SAPFmx2 = SAPFmx2.*0;
end



%-------------------
% G1SAPPOLYGON
%-------------------
%{

% SAPFmx1 & SAPFmx2 = 8x8 Mx of #SAPs in each 64 SLOTs of SPYN1 & SPYN2
% this G1SAPPOLYGON subroutine takes a 1x64 Cell of 1x1 polygons
% 
% G1INPSD1 holds logic array testing if a GluR in a polygon (on each iteration)
% 
% GluR1SdP1 takes G1INPSD1 logic test and SAPFmx1 and stores the #SAP
% 	each GluR is currently surrounded by. cycles down each column
% 
% GluR1SLOCMXa stores the ID of the polygon each GluR is located in
% The polygon check starts with the top left and cycles down each column
% like so...

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

GluR1SP1 = GluR1SdP1(1,:); %!! GluR1SdP1 doesn't need 2 rows
GluR1SP2 = GluR1SdP2(2,:);
G1SAPLOC1 = GluR1SLOCMXa;
G1SAPLOC2 = GluR1SLOCMXb;

% GLUR-ID , POLY ID , SAPS NEAR GLUR , GLUR DWELL TIME
G1P1SAPMX = [(1:numel(GluR1SP1))' G1SAPLOC1' GluR1SP1' GluR1_TdwellSPYN(1,:)'];
G1P2SAPMX = [([1:numel(GluR1SP2)]') G1SAPLOC2' GluR1SP2' GluR1_TdwellSPYN(2,:)'];
G1P1SAPMX = sortrows(G1P1SAPMX,-4);
G1P2SAPMX = sortrows(G1P2SAPMX,-4);

% if stepN>2000; keyboard; end;
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

% G1P1SAPMX: GLUR-ID , POLY ID , SAPS NEAR GLUR , GLUR DWELL TIME
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
% G1P1poly = G1P1SAPM(:,2);
% [~, ~, G1P1poly] = find(G1P1poly);
% if numel(G1P1poly)>1 && stepN>400; keyboard; end
%===================================%



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

%----------------------------------%
% S1_MainClusterFun
%----------------------------------%
function [S1 G1SP1] = S1_MainClusterFun(S1Ph1, S1Ph2, S1Ph3, S1Ph4, S1Ph5,...
	S1, S1sz, sap, stepN, G1P1SAPM, G1SP1,GT,GTab)

%=========================================================%
% INITIALIZE CONSTANTS
%=========================================================%
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

%-------------------------------%
GhkMask=GTab{1};
GTon = GT(1);
GToff = GT(2);
LTPv = GT(3);
%-------------------------------%

%=========================================================%
% COMPUTE DYNAMIC LATTICE REPULSION FORCE FROM SAP DENSITY
%=========================================================%
S1Le = sap(5);
doDynamicLe = sap(17);
if doDynamicLe
L1x = linspace(1.0,3.2);
Tsz = numel(S1);			% Total matrix space
Ssz = numel(find(S1>0));	% Filled matrix space
Psz = round(Ssz/Tsz*100)+1;	% Percent filled space
S1Le = L1x(Psz);
end

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
%---------------------------------------------%
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

S=S1;
dT = sap(11);
hkMask=[0 1 0; 1 0 1; 0 1 0];


Lon = sap(21);	% On Energy (lower = more on events)
Bon = sap(22);	% On Neighbor-Independant Rate (new growth) (lower = more on)
Ron = sap(23);	% On Neighbor-Dependant Rate (cluster fill-in) (higher = more on)

Loff = sap(24);	% Off Energy (higher = more off events)
Boff = sap(25);	% Off Neighbor-Inependant Rate (uniform off)  (lower = more off)
Roff = sap(26);	% Off Neighbor-Dependant Rate (edge off) (higher = more off)


Pmx = rand(size(S));
Soc = (S>0);
Sno = ~Soc;
hk = convn(Soc,hkMask,'same');

Lhon = (hk-Lon) .* (-Bon);
Pon = 1 ./ (1+exp(Lhon));
Pkon = Sno .* ( Ron * dT * Pon );
Son = (Pkon>Pmx);

Lhoff = ((-hk)+Loff) .* (-Boff);
Poff = 1 ./ (1+exp(Lhoff));
Pkoff = Soc .* ( Roff * dT * Poff );
Soff = (Pkoff>Pmx);


%-------------------------------%
G1oc = zeros(S1sz);
%-------------------------------%

%=========================================================%
doLTPS1 = sap(19);
if doLTPS1
%=========================================================%

%=========================================================%
%		Get Matrix Locations of SAPs near AMPARs
%=========================================================%
G1P1poly0 = G1P1SAPM(:,2);
[~, ~, G1P1poly1] = find(G1P1poly0);
G1P1poly2 = G1P1poly1;
%=============================================%
for nx = 1:(numel(G1P1poly1))
%---------------------------------------------%
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
	if G1P1poly1(nx)<=72 && G1P1poly1(nx)>64
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*8)-1;
	end
%---------------------------------------------%
end % for nx
%=============================================%


%=============================================%
if GTon > 0
%---------------------------------------------%

G1oc(G1P1poly2)=1;

%-------------------------------%
GRhk = convn(G1oc,GhkMask,'same');
GRk=(GRhk.*GTon);
GSk=(GRhk.*GToff);
%-------------------------------%
Gon = Pkon+(GRk.*(Pkon + LTPv));
Goff = Pkoff+(GSk.*Pkoff);
%-------------------------------%
Son = (Gon>Pmx);
Soff = (Goff>Pmx);
%-------------------------------%



%-------------------------------%
% if stepN == 2000; keyboard; end;
%{
%=========================================================%
% GluA RT SAP recruitment tails coefficent (basal and LTP)
%=========================================================%
% G1RT = G1RTBASE;
% if stepN >= LTP1onG1 && stepN <= LTP1offG1; G1RT=G1RTLTP; end
% if stepN >= LTP2onG1 && stepN <= LTP2offG1; G1RT=G1RTLTP; end
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


%=============================================%
for polyid = 1:G1P1n
%=============================================%
G1P1IDo = G1P1poly2(polyid);

% G1P1IDb = G1P1IDo+1;
% G1P1IDr = G1P1IDo+16;
% G1P1IDd = G1P1IDo+17;
% 
% if G1P1IDo>=16
% 	G1P1IDb = G1P1IDo+1;
% 	G1P1IDr = G1P1IDo+16+(round(rand)*-32); 
% 	G1P1IDd = G1P1IDr+1;
% 	% 50:50 WHETHER EAST OR WEST SAP
% 	% ALWAYS EAST SAP IF COMMENTED
% 	%G1P1IDd = G1P1IDo+17+(round(rand)*-34);
% end
%---
G1oc(G1P1IDo)=1;
% G1oc(G1P1IDb)=1;
% G1oc(G1P1IDr)=1;
% G1oc(G1P1IDd)=1;
%---
% SG1oc(G1P1IDo)=1;
% SG1oc(G1P1IDb)=1;
% SG1oc(G1P1IDr)=1;
% SG1oc(G1P1IDd)=1;
%=============================================%
end % for polyid = 1:G1P1n
%=============================================%
% Gex = Sno .* SG1oc .* G1RT + Pkon;
% Gp = 1 ./ (1+exp((-1*1)*hk)) + .001;
% Gpex = Gex .* Gp;
% %---
% 
% Son = (Gpex>Pmx);


%}
%-------------------------------%
end % if G1RTBASE > 0
%-------------------------------%
end % if doLTPS1
%=========================================================%


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

set(S1Ph5,'CData',GRhk);
drawnow

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

%----------------------------------%
% S2_MainClusterFun
%----------------------------------%
function [S2 G1SP2] = S2_MainClusterFun(S2Ph1, S2Ph2, S2Ph3, S2Ph4, S2Ph5,...
    S2, S2sz, sap, stepN, G1P2SAPM, G1SP2,GT,GTab)
 
%=========================================================%
% INITIALIZE CONSTANTS
%=========================================================%
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
% sap17 = sap(17); % doDynamicLeP2
% sap18 = sap(18); % doDynamicLeP2
% sap19 = sap(19); % doLTPS2
% sap20 = sap(20); % doLTPS2
GhkMask=GTab{2};
GTon = GT(4);
GToff = GT(5);
LTPv = GT(6);

%=========================================================%
% COMPUTE DYNAMIC LATTICE REPULSION FORCE FROM SAP DENSITY
%=========================================================%
S2Le = sap(8);
doDynamicLe = sap(18);
if doDynamicLe
L1x = linspace(1.0,3.2);
Tsz = numel(S2);            % Total matrix space
Ssz = numel(find(S2>0));    % Filled matrix space
Psz = round(Ssz/Tsz*100)+1; % Percent filled space
S2Le = L1x(Psz);
end
 
%=========================================================%
%   TEST MASK AGAINST CLUSTER - GET CONVOLUTION MX
%=========================================================%
S=S2;
dT = sap(12);
hkMask=[0 1 0; 1 0 1; 0 1 0];

Lon = sap(27);	% On Energy (lower = more on events)
Bon = sap(28);	% On Neighbor-Independant Rate (new growth) (lower = more on)
Ron = sap(29);	% On Neighbor-Dependant Rate (cluster fill-in) (higher = more on)

Loff = sap(30);	% Off Energy (higher = more off events)
Boff = sap(31);	% Off Neighbor-Inependant Rate (uniform off)  (lower = more off)
Roff = sap(32);	% Off Neighbor-Dependant Rate (edge off) (higher = more off)


Pmx = rand(size(S));
Soc = (S>0);
Sno = ~Soc;
hk = convn(Soc,hkMask,'same');

Lhon = (hk-Lon) .* (-Bon);
Pon = 1 ./ (1+exp(Lhon));
Pkon = Sno .* ( Ron * dT * Pon );
Son = (Pkon>Pmx);

Lhoff = ((-hk)+Loff) .* (-Boff);
Poff = 1 ./ (1+exp(Lhoff));
Pkoff = Soc .* ( Roff * dT * Poff );
Soff = (Pkoff>Pmx);

% S = (Soc-Soff) + Son;

%-------------------------------%
G1oc = zeros(S2sz);
%-------------------------------%

%=========================================================%
doLTPS2 = sap(20);
if doLTPS2 
%---------------------------------------------------------%
%       Get Matrix Locations of SAPs near AMPARs
%---------------------------------------------------------%
G1P2poly0 = G1P2SAPM(:,2);
[~, ~, G1P2poly1] = find(G1P2poly0);
G1P2poly2 = G1P2poly1;

%=============================================%
for nx = 1:(numel(G1P2poly1))
%=============================================%
    if G1P2poly1(nx)<=8 && G1P2poly1(nx)>0
        G1P2poly2(nx) = (G1P2poly1(nx)*2)+(S2sz*0)-1;
    end
    if G1P2poly1(nx)<=16 && G1P2poly1(nx)>8
        G1P2poly2(nx) = (G1P2poly1(nx)*2)+(S2sz*1)-1;
    end
    if G1P2poly1(nx)<=24 && G1P2poly1(nx)>16
        G1P2poly2(nx) = (G1P2poly1(nx)*2)+(S2sz*2)-1;
    end
    if G1P2poly1(nx)<=32 && G1P2poly1(nx)>24
        G1P2poly2(nx) = (G1P2poly1(nx)*2)+(S2sz*3)-1;
    end
    if G1P2poly1(nx)<=40 && G1P2poly1(nx)>32
        G1P2poly2(nx) = (G1P2poly1(nx)*2)+(S2sz*4)-1;
    end
    if G1P2poly1(nx)<=48 && G1P2poly1(nx)>40
        G1P2poly2(nx) = (G1P2poly1(nx)*2)+(S2sz*5)-1;
    end
    if G1P2poly1(nx)<=56 && G1P2poly1(nx)>48
        G1P2poly2(nx) = (G1P2poly1(nx)*2)+(S2sz*6)-1;
    end
    if G1P2poly1(nx)<=64 && G1P2poly1(nx)>56
        G1P2poly2(nx) = (G1P2poly1(nx)*2)+(S2sz*7)-1;
	end
%=============================================%
end % for nx
%=============================================%



%=============================================%
if GTon > 0
%---------------------------------------------%

G1oc(G1P2poly2)=1;

%-------------------------------%
GRhk = convn(G1oc,GhkMask,'same');
GRk=(GRhk.*GTon);
GSk=(GRhk.*GToff);
%-------------------------------%
Gon = Pkon+(GRk.*(Pkon + LTPv));
Goff = Pkoff+(GSk.*Pkoff);
%-------------------------------%
Son = (Gon>Pmx);
Soff = (Goff>Pmx);

%-------------------------------%
end % if G1RTBASE > 0
%-------------------------------%
end % if doLTPS1
%=========================================================%


%=========================================================%
%		FINAL CLUSTER PARAMETERS
%---------------------------------------------------------%
S = (Soc-Soff) + Son;
S2=S;
%=========================================================%



%=========================================================%
%					PLOT CLUSTER
%---------------------------------------------------------%
if mod(stepN,100) == 0
set(S2Ph1,'CData',S2);
drawnow
end

%=========================================================%
%			ADDITIONAL CLUSTER PLOTS
%---------------------------------------------------------%
doS1PLOTS=1;
if doS1PLOTS
if mod(stepN,200) == 0
figure(1);
%---

set(S2Ph2,'CData',hk);
drawnow

set(S2Ph3,'CData',Pkon);
drawnow

set(S2Ph4,'CData',Sno);
drawnow

set(S2Ph5,'CData',G1oc);
drawnow

end % if mod(stepN,100) == 0
end % if doS1PLOTS
%=========================================================%
%%




 
%---------------------------% 
% if nn == 400; keyboard; end
% if G1P2n>1 && nn>1000; keyboard; end
%---------------------------%
end

%===================================%





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

xylLB1 = xyl(1,:) > LB(1);
xylRT1 = xyl(1,:) < RT(1);
xylLB2 = xyl(2,:) > LB(2);
xylRT2 = xyl(2,:) < RT(2);

inbox = xylLB1 & xylRT1 & xylLB2 & xylRT2;

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
] = FieldFun(doRun,fsizeX, fsizeY,PSD1size, PSD2size, periPSD1size, periPSD2size)


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




%=========================================================%
if doRun(8)
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
function [] = MAINPLOT(G1Ph1, G2Ph1, stepN,GluR2xyl,GluR1xyl,...
	XWIDE,YHIGH,XYBOXpr1,XYBOXpr2,XYBOXp1,XYBOXp2)
%-------------------------------------------%

set(G2Ph1,'XData',GluR2xyl(1,:),'YData',GluR2xyl(2,:));
drawnow

set(G1Ph1,'XData',GluR1xyl(1,:),'YData',GluR1xyl(2,:));
drawnow


%{
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
%}


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

%{.
if mod(stepN, 20) == 0 && stepN>=100
figure(36);
ti = -XWIDE:2:XWIDE; 
[XI,YI] = meshgrid(ti,ti);
ZI = griddata(GluR2xyl(1,:),GluR2xyl(2,:),G2Z(:),XI,YI,'v4');
hold on,
grid off,
axis([xlim, ylim, zlim]),
set(gca, 'Box', 'on'),
set(gca,'xticklabel',[]),
set(gca,'yticklabel',[]),
set(gca,'zticklabel',[]),
mesh(XI,YI,ZI),view(70, 10)
scatter3(GluR2xyl(1,:),GluR2xyl(2,:),G2Z(:),'.','MarkerEdgeColor',[0 0 1]),...
view(70, 10)
hold on;
scatter3(GluR1xyl(1,:),GluR1xyl(2,:),G1Z(:),'.','MarkerEdgeColor',[1 0 0]),...
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


