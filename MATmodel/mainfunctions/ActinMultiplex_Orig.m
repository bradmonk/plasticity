function [varargout] = ActinMultiplex_Orig(LBR,TIME,SIZE,DOES,REVA,AMX,MSK,AMS)
%% ActinMultiplex

clc, close all; scsz = get(0,'ScreenSize');

if (~isdeployed); cd(fileparts(which(mfilename))); end


%% -- WELCOME TO THE ActinMultiplex TOOLBOX

    if exist('LBR','var') && nargin > 0
        fprintf('RUNNING ActinMultiplex TOOLBOX VERSION 1.3 \n');
    else
        fprintf('DEPLOYING GUI: ActinMultiplexGUI.m \n');
        ActinMultiplexGUI()
      return
    end


%%
%-------------------------------%
% Timing-Related Parameters
%-------------------------------%
dT = TIME(2);
datarate = TIME(3);
viewTime = TIME(4);
%-------------------------------%
% PSD Cluster Size
%-------------------------------%
PSDsz = SIZE(1);
PSAsz = SIZE(2);
SN0 = PSDsz^2;
SYNsz = PSDsz+(PSAsz*2);
S=padarray(ones(PSDsz),[PSAsz PSAsz], 0);
S0=S;


%-------------------------------%
% doItems
%-------------------------------%
doGaussianMask = 1;
doRev = REVA(1);
doPlot = DOES(1);
doFluorPlot=DOES(3);
doNoff = 0; 

%-------------------------------%
% Revolving Particle Entry
%-------------------------------%
MuT = LBR(6) * dT;
Rev = REVA(2);	
if doRev
dT = dT*Rev;
MuT = LBR(6)*dT*Rev;
end



%===============================================%
%				AMPAR STUFF
%===============================================%
%{
doAMPARs = GLU(1);
% amparate = GLU(3);
AMPARN = GLU(4);

GhkMask=GTab;
GTon = GT(1);
GToff = GT(2);
LTPv = GT(3);

if doAMPARs
Fh3 = figure(13);
set(Fh3,'OuterPosition',(scsz./[2e-1 .2 4 4]))
Ph3 = imagesc(S0);
end
%}
%------------
% if doAMPARs (goes in loop)
%{
%====================================%
if doAMPARs; %if mod(stepN, amparate) == 0;
%-------------------------------%
SG1oc = zeros(SYNsz);
%-------------------------------%
GRPOS=randi([1 (SYNsz*SYNsz)],1,AMPARN);
SG1oc(GRPOS)=1;
GRhk = convn(SG1oc,GhkMask,'same');
GRk=(GRhk.*GTon); GSk=(GRhk.*GToff);
%-------------------------------%
Gon = Pkon+(GRk.*(Pkon+LTPv));
Goff = Pkoff+(GSk.*Pkoff);

Son = (Gon>Pmx);
Soff = (Goff>Pmx);

if doPlot;	if mod(stepN, 500) == 0
set(Ph3,'CData',GRk);
drawnow
end;		end;
%-------------------------------%
end; %end;
%====================================%
%}
%------------
% MASKS
%{
% imagesc(ActinTips{140})
%-------------------------------%
% hkMask=[0 1 0; 1 0 1; 0 1 0];
AMask=ones(AMX{5});
hkMask=ones(AMX{6});
% ActinBranches = 50;
% ActinSlotsV = 2;
% ACTINpos = zeros(PSDsz);
% ACTINrand=randi([1 (PSDsz*PSDsz)],1,ActinBranches);
% ACTINpos(ACTINrand)=ActinSlotsV;
% ACTINp=padarray(ACTINpos,[PSAsz PSAsz], 0);
%-------------------------------%
%===============================================%
%}
%------------
%-------------------------------%


%===============================================%
%				RATE PARAMETERS
%===============================================%
Lon = LBR(1);
Loff = LBR(2);
Bon = LBR(3);
Boff = LBR(4);
Ron = LBR(5);
Roff = LBR(6);
%--------



%===============================================%
%	ACTIN MULTIPLEX TIP DATA
%===============================================%
global ATs;
global Ax;
global AMx;
global AFMx;

AMx = AMX;

loadActinTips = AMX{1};
generateActinTips = AMX{2};

%----------
if loadActinTips
savetipdir(AMX)
ATdat = which(AMX{3});
load(ATdat);
assignin('base', 'ATs', ATs)
assignin('base', 'Ax', Ax)

LivePlot(Ax{1},Ax{2},Ax{3},Ax{4},Ax{5},Ax{6},Ax{7})
end
%----------



%----------
if generateActinTips
ActSteps = AMX{4};
[ActinTips, AxLP, Ax, AFMx] = ActinMainStage2(ActSteps,AMX,AMS);

varargout={0};
return
%----------
if AMS{5} || ~AMX{25}
varargout={0};
return; 
end
%----------

assignin('base', 'ATs', ActinTips)
ATs = evalin('base', 'ATs');
assignin('base', 'Ax', AxLP)
Ax = evalin('base', 'Ax');

assignin('base', 'AMx', AMx)
AMx = evalin('base', 'AMx');

assignin('base', 'AFMx', AFMx)
AFMx = evalin('base', 'AFMx');


savetipdir(AMX,2)
%save('ATdata.mat', 'ATs')
end

if ~loadActinTips
	ATs = evalin('base', 'ATs');
	Ax = evalin('base', 'Ax');
	AMx = evalin('base', 'AMx');
end
%----------


%------------
% SHRINK MX
%{
% doActMxShrink = 1;
% SkTo = 150;
% if doActMxShrink
% 	ATn = ATs{1};
% 	Asz = size(ATn);
% 	Ash = Asz - SkTo;
% 	Adel = round(linspace(1,Asz(1),Ash(1)));
% 	for n = 1:numel(ATs)
% 		ATn = ATs{n};
% 		ATn(Adel,:) = [];
% 		ATn(:,Adel) = [];
% 		ATs(n) = {ATn};
% 	end
% end


% for n = 1:numel(ATs)	
% 	SkTo = 150;
% 	ATn = ATs{n};
% 	Asz = size(ATn,1);
% 	TipsN = nnz(ATn);
% 	Tdc = sum(ATn,1);
% 	Tdr = sum(ATn,2);
% 	Tdcv = Tdc > 0;
% 	Tdrv = Tdr > 0;
% 	TdCR = [Tdcv' Tdrv];
% 	TdCRv = sum(TdCR,2);
% 	Tdf = find(TdCRv);
% 	while Asz > SkTo
% 		ARc = randi([1,Asz],1);
% 		ARok = sum(Tdf == ARc);
% 		ARok = ARok<1;
% 		if ARok
% 			ATn(ARc,:) = [];
% 			ATn(:,ARc) = [];
% 		end
% 		Asz = size(ATn,1);
% 	end
% 	ATs(n) = {ATn};
% end
%}
%------------


%------------
ActinTips = ATs;
TipCellN = numel(ActinTips);
TipCellA = TipCellN - AMX{7};
NumTipCells = TipCellN - TipCellA;
%------------
TrimActMx = size(ATs{1},1)+1;
ACTINp = ActinTips{TipCellA};
ACTINp(:,TrimActMx:end) = [];
ACTINp(TrimActMx:end,:) = [];
%------------
%-------------------------------%
% Mask Setup
%-------------------------------%
[hkMask dT LBR] = MaskFun(S,dT,LBR,doGaussianMask,PSAsz,scsz,AMX,MSK);
%hkMask=ones(AMX{6});


AMask=ones(AMX{5});
Ahk = convn(ACTINp,AMask,'same');
S = (Ahk>0).*1.0;
%------------

%-------------------------
% DEDUCE NSteps
%-------------------------
ActUpdate = AMX{8};
NSteps = ActUpdate * NumTipCells - ActUpdate - 1;
%-------------------------


%================================================%
%				RUN CLUSTERING
%------------------------------------------------%

doClustering = AMX{28};
if ~doClustering
	varargout = {AMX,MSK,LBR};
	return
end





%================================================%
%				FIGURE SETUP
%------------------------------------------------%
if doPlot
%--------
sz = [2.5e-3 2.5e-3 1.6 1.2];
Fh1 = FigSetup(11,sz);
Ph1 = imagesc(S);
colormap('bone')
hold on;
%--------
[PSDY,PSDX] = find(ACTINp);
Ph2 = scatter(PSDX,PSDY,100, [.9 .08 .24],'filled');
hold off
%--------
end
%------------------------------------------------%

% RGB/255
% [255 20 147] ./ 255

%====================================================================%
%						MAIN CLUSTERING LOOP
%====================================================================%
for stepN = 1:NSteps
%--------------------------------------------------------------------%

%------------
if mod(stepN,ActUpdate) == 0
TipCellA = TipCellA+1;
ACTINp = ActinTips{TipCellA};
ACTINp = convn(ACTINp,AMask,'same');
end
%------------


%------------
Pmx = rand(size(S));
Soc = (S>0);
Sno = ~Soc;
%---
Ah = convn(ACTINp,hkMask,'same');
%------------
APon = 1 ./ (1+exp((Ah-Lon).*(-Bon)));
APkon = Sno .* ( Ron * dT * APon );
%---
Son = (APkon>Pmx);
%------------
APoff = 1 ./ (1+exp(((-Ah)+Loff).*(-Boff)));
APkoff = Soc .* ( Roff * dT * APoff );
%---
Soff = (APkoff>Pmx);
%------------

S = (Soc-Soff) + Son;

% if doPlot
if mod(stepN, viewTime) == 0
set(Ph1,'CData',S);
drawnow
%--------
[PSDY,PSDX] = find(ACTINp);
set(Ph2,'XData',PSDX,'YData',PSDY);
drawnow;
%--------
end
% end

%--------------------------------------------------------------------%
end % end main clustering loop
%====================================================================%



varargout = {AMX,MSK,LBR};
%====================================================================%
end % end main function
%####################################################################%





%####################################################################%
%
%							ACTIN MAIN STAGE
%
%####################################################################%
function [varargout] = ActinMainStage2(ActSteps,AMX,AMS,varargin)
scsz = get(0,'ScreenSize');

BTs = [];
AFMx = [];
Nsteps = ActSteps;
dT = AMX{47};

SaveTipsAfter = AMX{19};
SaveTipsRate = AMX{20};
SvTpMx = AMX{25};

doActCounts = AMX{35};
doDelOrig = AMX{39};
doDelOrigT = AMX{40};

%-- LIVE PLOTS --
doLiveTipPlot = AMS{1};
doLiveHullPlot = AMS{2};
LiveTipMod = AMS{3};
LiveHullMod = AMS{4};

LivePlotMod = AMS{6};
TriHullMod  = AMS{7};
doTM = 0;
%-------


runTest = AMS{5};
if runTest
    Nsteps = 30000;
    LivePlotMod = 1000;
    SaveTipsAfter = 26000;
    SaveTipsRate = 10;
    doLiveHullPlot = 0;
    doLiveTipPlot = 0;
end




%---------------------------------------------%
%				ANGLE SETUP
%---------------------------------------------%
% unitsratio('rad', 'deg')
d2r = 1/(180/pi);
Ov = [2 29 57 85 112 140 168 195 223 251 278 306 334];
Ov = [Ov -15 -40 -65 -100 -125 -150 -175 -205 -230 -260 -295 -320 -350];

%----------------------------------------
% Spine Dimensions
%----------------------------------------
SPYnXY = AMX{10}/2;	% Spy neck XY radius
SPYhXY = AMX{13}/2;	% Spy head XY radius
SPYhZN = AMX{11};	% Spy head north
SPYhZS = AMX{12};	% Spy head south


PSDproxy = AMX{15};
PSDproxyXY = AMX{27};
inPSD = SPYhZN - PSDproxy;

SPYH = [SPYhXY SPYhXY];
AcMx = zeros(SPYhXY*2/5,SPYhXY*2/5);

dims = [SPYnXY SPYhZN SPYhZS SPYhXY SPYhXY PSDproxy inPSD PSDproxyXY];
%----------------------------------------



%-----------------------------------------------------------------------------------------%
%                       INITIALIZE AND TAG STARTING FILAMENTS
                NStFils = AMX{42};      Actin = zeros(NStFils,20);  
%-----------------------------------------------------------------------------------------%
%                               Actin(Nfil,19)
% N  Xa  Xo  Xt  Ya  Yo  Yt  Za  Zo  Zt  MomID  ID  Fkd Born Died Lif MaxL MeanL Null Lgth
% 1  2   3   4   5   6   7   8   9   10  11     12  13  14   15   16  17   18    19   20
%-----------------------------------------------------------------------------------------%


% Starting Length Loc & Angles
StartMonos = AMX{41};
fXYo = AMX{45};
fZo = AMX{43};
fXYa = AMX{46};
fZa = AMX{44};

% Branching Angles
TPi = d2r*AMX{9};
PPi = 0;

%----------------------------------------
% MAKES STARTING FILAMENTS
Actin = MakeStartFils(Actin,NStFils,StartMonos,d2r,fZa,fXYo,SPYhZN,SPYhZS,SPYhXY);
%----------------------------------------

%----------------------------------------
TagN = numel(Actin(:,1));
Actin(:,12) = 1:TagN;
TagN = TagN+1;
%----------------------------------------


%% ---- ACTIN SIZES ----
% NOTES
%{
From: Andre Kamkin, Irina Kiseleva
Springer Science & Business Media, Nov 18, 2010 - Biochemistry - 395 pages
Chapter Authors: Luo and Robinson
2.2 Microstructures and Deformations of the Actin Cyctoskeleton
- Gactin is 5 nm in diameter
- Factin filaments are 8 nm wide with a left-handed helical morphology
- 13 actin monomers per pseudo-repeat
- 1 pseudo-repeat length of 37 nm 
- Alternatively the actin filament can be considered to have a right-handed 
helical structure with two strands slowly twisting around each other. 
Each actin monomer is rotated 166 degrees with respect to its nearest neighbors 
across the strand (Holmes 1990). Within the strand, subdomains 2 and 4 contact 
subdomains 1 and 3 in the next monomer in the strand, and each monomer reaches 
across to the other strand through a hydrophobic plug that links the two 
strands together. 

37 nm / 13 p = 2.85 nm/p
13 p / 37 nm = 0.35 p/nm

From: Lodish, Principles of Molecular Biology
Adapted from C. E. Schutt et al. 1993, Nature 365:810
- There are 2 strands with 14 units per strand, 
- so 28 monomers in each 360 degree turn,
- with a length of 72 nm per 360 degree turn

72 nm / 28 p = 2.57 nm/p
28 p / 72 nm = 0.39 p/nm

Factin particle size in filaments:
0.369 p/nm   (Factin particles per nm of filament)
2.71  nm/p    (nm of filament per Factin particle)

%}
%-------

p_per_nm = 0.369;    % (Factin particles per nm of filament)
nm_per_p = 2.71;     % (nm of filament per Factin particle)
Actin(:,20) = Actin(:,1) .* nm_per_p;   % Store filament length
%------------------------

%----------------------------------------
% TRIG: branch XYZ tip coordinates
Actin(:,4) = Actin(:,1) .* nm_per_p .* sin(Actin(:,8)) .* cos(Actin(:,2)) + Actin(:,3);
Actin(:,7) = Actin(:,1) .* nm_per_p .* sin(Actin(:,8)) .* sin(Actin(:,2)) + Actin(:,6);
Actin(:,10) = Actin(:,1) .* nm_per_p .* cos(Actin(:,8)) + Actin(:,9);
%----------------------------------------



%----------------------------------------
ACTs = Actin;
oActin = Actin;
DiedACTs = Actin(1,:);
%----------------------------------------






%% ---- Conversion Factors & VCP Function ----
% >> pNuM = VCP(vol,uM,pN)
%-------------------------
% NOTES
%{
 
% [pNuM] = VCP(vol,uM,pN);
% INPUTS
% vol: volume in um^3
% uM: concentration in uM
% pN: particle count
% 
% enter zero for the unknown value
% if pN is unknown enter... pN = VCP(.1,10,0)
% if uM is unknown enter... uM = VCP(.1,0,6e5)
 
%}
%-------

MOL = 6.022e23;     % 1 mol Avagadro's number
mol = 6e23;			% 1 mol rounded



%% ---- SPINE VOLUME ----
% NOTES
%{
%---------------------------------------------------
% A typical spine volume is 0.1 um^3 or 1e-16 L
% volume of cylinder
% V = pi * r^2 * h
%---------------------------------------------------


% MATH - Spine Volume

This spine has a neck volume equivalent to a cylendar of dimensions: 
50d x 200

And a head volume equivalent to a 3D disk with dimensions:
100d x 100


% Units
mol = 6e23;		% in N
mM = 1e-3;
uM = 1e-6;
upM = 1e3;
dnM = 1e-3;

% 1 cm^3 = 1 mL
% SI units conversion from m^3 to L are:
% 1 cm^3 = 1 mL
% To convert cubic volume to mL you may first need to convert
% by an order of magnitude; here's a reminder of cubic volume conversions:
% 1 m^3 = 1 cm^3 * .01^3;		% 1e-6
% 1 m^3 = 1 mm^3 * .001^3;		% 1e-9
% 1 cm^3 = 1 mm^3 * .1^3;		% 1e-3
% 1 cm^3 = 1 um^3 * .0001^3;	% 1e-12
% Thus, if dendritic spines have an average volume of 0.1 um^3
% that would be equivalent to X uL
%
% 0.1 um^3 * (1 cm^3 / 1e12 um^3) * (1 mL / 1 cm^3) * (1000 uL / 1 mL)
% 0.1*(1/1e12)*(1/1)*(1000/1)
% 0.1 um^3 = .1e-12 cm^3 = .1e-9 uL
%
% and 0.1 um^3 equivalent to X L
% 
% 0.1*(1/1e12)*(1/1)*(1/1000)
% 1e-16


% Spine volume (0.1 um^3) in L and uL
SpyV = 1e-16;	% in L
SpyVu = .1e-9;	% in uL

% Actin Polymerization (12 p/µM*s)
% Act_Na = 1e5;
% BeKanT = 1000;

% Actin Depolymerization (2 p/s)	
BeKdnT = 200;
DePSum = 0;

% Actin Poly Math
% Cytosolic concentration of actin in cells ranges from .1 to .5 mM
% Given a spine volume of 'SpyV' we can find how many actin monomers
% are in an average spine:
%
% .1 mM (mmol/L) * (1 mol / 1000 mmol) * SpyV = 1e-17 mol
% 1e-17 mol * (6e23 units / 1 mol) = 6e3 monomer units

molAct = AMX{16};
Act_N = molAct * (1/1000) * SpyV * 6e23; % 6e3 monomer units

% we can check our math starting with a set number of actin monomers
% and calculate the spine molarity (6e3 monomer units as an example):
% 
% 6e3 units/SpyV * (1 mol / 6e23 units) * (1000 mmol / 1 mol)
% 6e3/SpyV*(1/6e23) 

Act_uM = Act_N/SpyV * (1/6e23) * (1000/1);	% 1.6e-10 

% Act_uM = Act_N / SpyVu / mol;	% 1.6e-10 
% Act_N = Act_N / SpyVu / mol;
BeKa = 12 * Act_uM * dT;
BeKd = 2 * dT;

volume of cylinder
V = pi * r^2 * h

volume of a sphere
V = (4*pi*r^3)/3



%}
%-------


Vneck = pi * SPYnXY^2 * SPYhZS;
Vhead = pi * SPYhXY^2 * (SPYhZN-SPYhZS);
SpyV = (Vneck+Vhead) * 1e-24;	% nm^3 to L conversion (nm: 1e-24; um: 1e-15)
SpyV = 1e-16;




%% ---- ACTIN CONCENTRATIONS ----
% NOTES
%{

NA

%}
%-------

uM_Act = AMX{16};                           % Actin uM
GActinN0 = uM_Act / 1e6 * MOL * SpyV;       % t0 N Gactin monomer units (6e4)
Act_uM = GActinN0 / SpyV *(1/MOL)*1e6;      % Check uM_Act == Act_uM

Nfi = numel(Actin(:,1));                    % current number of branches
FActinN = sum(Actin(:,1));                  % current number FActins
GActinN = ceil(GActinN0) - FActinN;         % current number GActins


%% ---- ACTIN POLYMERIZATION & DEPOLYMERIZATION RATES ----
% NOTES
%{
%----------------------------------------------------
% Pollard:	
% Actin+	ATP      ADPP       ADP		| MEAN
% BeKa: 	11.6     3.4        2.9		| 6.0 p/µM*s
% BeKd:		1.4      0.2		5.4		| 2.3 p/s
% PeKa:		1.3      0.11       0.09	| 0.5 p/µM*s
% PeKd: 	0.8      0.02       0.25	| 0.4 p/s
%           BeKaATP  BeKaADPP   BeKaADP
%           BeKdATP  BeKdADPP   BeKdADP
%           PeKaATP  PeKaADPP   PeKaADP
%           PeKdATP  PeKdADPP   PeKdADP
%----------------------------------------------------
%}
%-------

BeKaATP = AMX{60};  BeKaADPP = AMX{61};  BeKaADP = AMX{62};
BeKdATP = AMX{63};  BeKdADPP = AMX{64};  BeKdADP = AMX{65};
PeKaATP = AMX{66};  PeKaADPP = AMX{67};  PeKaADP = AMX{68};
PeKdATP = AMX{69};  PeKdADPP = AMX{70};  PeKdADP = AMX{71};
% BeKaATD = mean([BeKaATP  BeKaADPP  BeKaADP]);
% BeKdATD = mean([BeKdATP  BeKdADPP  BeKdADP]);
% PeKaATD = mean([PeKaATP  PeKaADPP  PeKaADP]);
% PeKdATD = mean([PeKdATP  PeKdADPP  PeKdADP]);
BeKaATD = mean([BeKaATP  mean([BeKaADPP  BeKaADP])]);
BeKdATD = mean([BeKdATP  mean([BeKdADPP  BeKdADP])]);
PeKaATD = mean([PeKaATP  mean([PeKaADPP  PeKaADP])]);
PeKdATD = mean([PeKdATP  mean([PeKdADPP  PeKdADP])]);

KaATD = roundn(mean([BeKaATD PeKaATD]),-2);
KdATD = roundn(mean([BeKdATD PeKdATD]),-2);

%------
doKez = AMX{72};
KaBE = AMX{23};				% Barbed End On-Rate Scalar
KdBE = AMX{24};				% Barbed End Off-Rate Scalar
KaPE = AMX{58};				% Pointed End On-Rate Scalar
KdPE = AMX{59};				% Pointed End Off-Rate Scalar
KaTIP = (KaBE + KaPE) / 2;	% Tip Mean Ka On-Rate Scalar
KdTIP = (KdBE + KdPE) / 2;	% Tip Mean Kd Off-Rate Scalar

BeKa = KaBE * Act_uM * dT;	% Ka Barbed End ON-Rate
BeKd = KdBE * dT;			% Kd Barbed End OFF-Rate
PeKa = KaPE * Act_uM * dT;	% Ka Pointed End ON-Rate
PeKd = KdPE * dT;			% Kd Pointed End OFF-Rate
%------


%--------
if doKez
	TKa = KaTIP;
	TKd = KdTIP;
else
	TKa = KaATD;
	TKd = KdATD;
end
%---
fKa = TKa * Act_uM * dT;	% Ka Fil ON-Rate
fKd = TKd * dT;				% Kd Fil OFF-Rate







%% ----  COFILIN & GENERAL DEPOLYMERIZATION VARIABLES ----
% NOTES
%{

NA

%}
%-------
fdKd = fKd*10;				% Depoly rate when Fil is Fkd

CofR = AMX{18}* dT; % cofilin activity rate
CofS = AMX{21}; % cofilin delete Nunits
CofN = AMX{30}; % cofilin delete Nunits
delOr =AMX{14}* dT; % delete from origin rate

CofSMax = AMX{26};
TheoMaxFact = ceil(SPYhZN / nm_per_p);
ScFil = ceil(TheoMaxFact/CofSMax);



%% ----  ARP2/3 FILAMENT BRANCHING VARIABLES ----
% NOTES
%{
According to Smith & Gelles 2012:
Arp helps to nucleate new daughter filament branches 
at a rate of 2.5 to 9.7 df/(mM_Arp * s * ummf)
2.5 with no WASp and 9.7 with 300 nM WASp

At some intermediate WASp level, branching may proceed at
5 daughter filaments (df)
per mM of Arp
per second 
per micrometer of mother filament (ummf)

Example calculations with 9 uM Arp and 21600 nm of mf
2.5 df/(Arp_mM*s*ummf) * ((9/1000) * dT * (21600/1000))
9.7 df/(Arp_mM*s*ummf) * ((9/1000) * dT * (21600/1000))
5.0 df/(Arp_mM*s*ummf) * ((9/1000) * dT * (21600/1000))

5 * ((9/1000) * dT * (21600/1000))
5 * ((Arp_uM/1000) * dT * (nmmf/1000))



Friederich reports an Arp2/3 association rate of:
5.4e-4 um^-3 S^-1 (crossref from Carlsson 2004) but use a value of
1e-5 b/µM*s as a baseline parameter in their simulation software ActSimChem
ASRT:	1e-5 b/µM*s

These 3 values give branching rates that differ orders of magnitude:
5 * (9/1000) * 1	% .045
1e-5 * 9 * 1		% .00009
5.4e-4 * (9^3)		% 0.39

9000
But the Gelles and Carlsson value are fairly close when the umf is 9
5 * (9/1000) * 9	% 0.40
5.4e-4 * (9^3)		% 0.39

which would be the case if the Factin concentration was 51 uM
(9000/2.7) / 1e-16 *(1/6e23)*1e6

5.4e-4 * Arp_uM^3 * dT

%}
%-------

Arp_Sc = AMX{17} * dT;					% Arp empirical branching rate scalar
Arp_uM = AMX{33};						% Arp uM
nmmf = Actin(:,20);						% Length of filaments (nm)

FArpN = NStFils;						% #of F-Arp
GArpN = Arp_uM / 1e6 * 6e23 * SpyV;		% #of G-Arp
GArpN0 = ceil(GArpN);					% #of G-Arp (starting)
Arp_uM = GArpN / SpyV *(1/6e23)*1e6;	% Check uM_Act == Act_uM

ArpAdd = AMX{29};						% Add X units to new branches
ARPmax = AMX{22};						% Maximum Arp branches

ArpBR  = Arp_Sc * ((Arp_uM/1000) .* (nmmf/1000));	% Arp branching rate (per fil)


% MATH - Arp Branching
%{

Tnmmf = sum(Actin(:,20));				% Length of filaments (nm) Combined
ArpBRT = Arp_Sc * ((Arp_uM/1000) .* (Tnmmf/1000));	% Arp rate scalar    (total)

%ArpN = 1e3;
%ArpOn = 5;
%ArpOff = 1;
%Arp_uM = ArpN / SpyVu / mol;	% 1.6 - 16 uM
%Arp_PR = ArpOn * Arp_uM * dT;
%Arp_DR = ArpOff * dT;
%ArpR = AMX{17}* dT;	% Arp activity rate
%ArpR0 = ArpR;			% Arp activity rate
%ArpScalar = AMX{21};	% Arp filament length scalar
%}
%----------------------------------------




%% ----  THYMOSIN ACTIN-SEQUESTERING VARIABLES ----
% CONCENTRATION NOTES
%{

%---
%}
%-------

uM0_Thymosin     = AMX{73};      % Thymosin concentration (uM) total
N0_Thymosin = uM0_Thymosin / 1e6 * MOL * SpyV;    % Thymosin count in spine (N)
Thymosin_uM = N0_Thymosin / SpyV *(1/MOL)*1e6;    % Thymosin concentration in spine (uM)


ThymBoundFrac   = AMX{74}/100;  % Thymosin fraction bound to Actin
THYM_uM     = Thymosin_uM * (1-ThymBoundFrac);  % Thymosin        monomers (uM)
THYM_ACT_uM = Thymosin_uM * ThymBoundFrac;      % Thymosin-Actin  dimers (uM)
THYM_N      = THYM_uM / 1e6 * MOL * SpyV;       % Thymosin        monomers (N)
THYM_ACT_N  = THYM_ACT_uM / 1e6 * MOL * SpyV;   % Thymosin-Actin  dimers (N)
       
% -------
% THYMOSIN REACTION RATE NOTES
%{
The Law of Mass Action (*LMA) describes the rate at which chemicals collide
and interact to form different chemical combinations. When two different
chemicals can collide to form a dimer product, and the dimer can dissociate
reversably back into the individual component reactants, is described as:

        Ka>
T + A <---->  TA
       <Kd  

Where
    T : thymosin                (thymosin monomers)
    A : actin                   (Gactin monomers)
    TA: thymosin-actin          (thymosin-actin dimers)
    Ka: forward rate constant
    Kd: backward rate constant


dTA/dt = Ka[T][A] - Kd[TA]      (LMA forward reaction: TA accumulation)
dA/dt  = Kd[TA] - Ka[T][A]      (LMA reverse reaction: A accumulation)
dT/dt  = Kd[TA] - Ka[T][A]      (LMA reverse reaction: T accumulation)

%}
%-------

Ka_Thymosin     = AMX{75};      % Thymosin Ka (basal)
Kd_Thymosin     = AMX{76};      % Thymosin Kd (basal)
Ka_Thymosin_LTP = AMX{77};      % Thymosin Ka (LTP)
Kd_Thymosin_LTP = AMX{78};      % Thymosin Kd (LTP)

KaThy = (dT*Ka_Thymosin * THYM_uM * Act_uM); % Thymosin+Actin Association
KdThy = (dT*Kd_Thymosin * THYM_ACT_uM);      % Thymosin+Actin Dissociation

%% ----------------------------------------



%----------------------------------------
% ADVANCED PARAMETERS
%----------------------------------------
%{

		Abbreviations
-------------------------------------
(+)end: Barbed end of filament
(-)end: Pointed end of filament
mono: 	monomer, free protein monomer not associated with a filament
proto: 	protomer, polymerized protein or protein associated with a filament 
Pi: 	Phosphate
ATG: 	G-actin + ATP (same as ATM)
ADG: 	G-actin + ADP (same as ADM)
ATM: 	G-actin + ATP
ADM: 	G-actin + ADP 
ATF: 	F-actin proto + ATP
APF: 	F-actin proto + ADP-Pi
ADF: 	F-actin proto + ADP
FTB: 	(+)ends terminating in ATP-actin
FPB: 	(+)ends terminating in ADP-Pi-actin
FDB: 	(+)ends terminating in ADP-actin
FTP: 	(-)ends terminating in ATP-actin
FPP: 	(-)ends terminating in ADP-Pi-actin
FDP: 	(-)ends terminating in ADP-actin
CBM: 	(+)end capping mono
CBF: 	(+)end capper proto
CPM, 	(-)end capping mono
CPF: 	(-)end capper proto
FOM: 	Formin mono
FOF: 	Formin proto at (+)end
ARM: 	Arp in mono
ARF: 	Arp proto
FRP: 	Arp proto at (-)end (complex has bound actin at growth end)
FRB: 	Arp proto at (+)end (complex has no bound actin at growth end)



SNUC:	Spontanious Nucleation
FNUC:	Formin nucleation
CBNU:	Nucleation by barbed cap
CPNU:	Nucleation by pointed cap
ASTB:	Poly on rate of ATG at (+)end
ASDB:	Poly on rate of ADG at (+)end
ASTP:	Poly on rate of ATG at (-)end
ASDP:	Poly on rate of ATG at (-)end
DITB:	DPoly off rate of ATF at (+)end
DIPB:	DPoly off rate of APF at (+)end
DIDB:	DPoly off rate of ADF at (+)end
DITP:	DPoly off rate of ATF at (-)end
DIPP:	DPoly off rate of APF at (-)end
DIDP:	DPoly off rate of ADF at (-)end
TTOP:	ATP hydrolysis (ATF -> APF)
PTOD:	P release (APF -> ADF)
DTOT:	Recharge of G-actin by ATF
ASRT:	Arp association rate
DIRP:	Arp dissociation rate
FRGM:	Spontanious fragmentation
ANNL:	Annealing
ASFB:	Association of formin to (+)end
DIFB:	Dissocaition of formin
FASB:	Formin-aided filament growth
ASCB:	Capping of (+)end
DICB:	Uncapping of (-)end
ASCP:	Capping of (+)end
DICP:	Uncapping of (-)end


		Default Rates
-------------------------------------
SNUC:	1e-8	p/(uM^2*s^1)
FNUC:	7e-5	p/(uM^3*s^1)
CBNU:	1e-5	p/(uM^3*s^1)
CPNU:	1e-5	p/(uM^3*s^1)
ASTB: 	11.5	p/µM*s
ASDB: 	3.8		p/µM*s
ASTP: 	1.3		p/µM*s
ASDP: 	0.16	p/µM*s
DITB:	0.90	p/s
DIPB:	0.90	p/s
DIDB:	1.50	p/s
DITP:	0.19	p/s
DIPP:	0.19	p/s
DIDP:	0.26	p/s
FRGM:	1.8e-8	p/s
ANNL:	1e-8	p/µM*s
TTOP:	0.30	p/s
PTOD:	2.6e-3	p/s
DTOT:	1e-2	p/s
ASRT:	1.0e-5	p/µM*s
DIRP:	1.0e-3	p/s
ASCB:	3.0		p/µM*s
DICB:	4e-4	p/s
ASCP:	1.0		p/µM*s
DICP:	1e-2	p/s
ASFB:	3.0		p/µM*s
DIFB:	1e-4	p/s
FASB:	110		p/µM*s



		Rate Calculation Notes
-------------------------------------
exponents of µ^-1 == 1/µ
exponents of µ^-2 == 1/µ/µ
exponents of µ^-2 == 1/µ/µ/µ
[note these are sequential divisions, not 1/(µ/µ)]
Thus if something has a rate of [0.1 µM^-3 s^-1]
to calculate the current rate do [p=.1	uM=uM	s=dT]

	p/(uM^-3)/(s^-1)   or   p*(uM^3)*(s^1)
-------------------------------------



		Default Rates
-------------------------------------
NUCLEATIONS
SNUC:	1e-8	p/(uM^2*s^1)
FNUC:	7e-5	p/(uM^3*s^1)
CBNU:	1e-5	p/(uM^3*s^1)
CPNU:	1e-5	p/(uM^3*s^1)

ACTIN ASSOCIATIONS
ASTB: 	11.5	p/µM*s
ASDB: 	3.8		p/µM*s
ASTP: 	1.3		p/µM*s
ASDP: 	0.16	p/µM*s

ACTIN DISSOCIATIONS
DITB:	0.90	p/s
DIPB:	0.90	p/s
DIDB:	1.50	p/s
DITP:	0.19	p/s
DIPP:	0.19	p/s
DIDP:	0.26	p/s

FRAGMENTATION/ANNEALING
FRGM:	1.8e-8	p/s
ANNL:	1e-8	p/µM*s

ATP/ADP
TTOP:	0.30	p/s
PTOD:	2.6e-3	p/s
DTOT:	1e-2	p/s

ARP2/3
ASRT:	1.0e-5	p/µM*s
DIRP:	1.0e-3	p/s

CAPPING
ASCB:	3.0		p/µM*s
DICB:	4e-4	p/s
ASCP:	1.0		p/µM*s
DICP:	1e-2	p/s

FORMIN
ASFB:	3.0		p/µM*s
DIFB:	1e-4	p/s
FASB:	110		p/µM*s








Abbreviations Full Description: 
F-actin, Filamentous actin; 
nSRF model, non-structurally-resolved filament model; 
SRF model, structurally-resolved filament model; 
Pi, Phosphate; 
ATM: Globular actin (monomeric form) with incorporated ATP; 
ADM: Globular actin (monomeric form) with incorporated ADP; 
ATF: Filamentous actin protomer (F-actin) with incorporated ATP; 
APF: Filamentous actin protomer (F-actin) with incorporated ADP-Pi; 
ADF: Filamentous actin protomer (F-actin) with incorporated ADP; 
FTB: Barbed ends of filaments, terminating by ATP-actin; 
FPB: Barbed ends of filaments, terminating by ADP-Pi-actin; 
FDB: Barbed ends of filaments, terminating by ADP-actin; 
FTP: Pointed ends of filaments, terminating by ATP-actin; 
FPP: Pointed ends of filaments, terminating by ADP-Pi-actin; 
FDP: Pointed ends of filaments, terminating by ADP-actin; 
CBM: Barbed-end capping protein (capper) in monomer (free) form; 
CBF: Barbed-end capper bound to filament; CPM, Pointed-end capper in monomer form; 
CPF: Pointed-end capper bound to filament; 
FOM: Formin in monomer (free) form; 
FOF: Formin, bound to filament barbed ends; 
ARM: Arp2/3 in monomer (free) form; 
ARF: Arp2/3 associated with filament; 
FRP: Arp2/3 associated with filament pointed end (fil terminates at ARP2/3); 
FRB: Arp2/3 associated with filament with no bound actins (fil terminates at ARP2/3); 

%}
%----------------------------------------








%----------------------------------------
% LTP related variables

doActLTP = AMX{36};
ActLTPstart = AMX{37};
ActLTPend = AMX{38};
ArpRLTP = ArpBR*AMX{31};
%-------

doActT = AMX{48};
doArpT = AMX{49};


T2start = AMX{50};
T3start = AMX{51};
T4start = AMX{52};
T5start = AMX{53};

GActT2 = GActinN0 * AMX{54};
GArpT2 = GArpN0 * AMX{54};
GActT3 = GActinN0 * AMX{55};
GArpT3 = GArpN0 * AMX{55};
GActT4 = GActinN0 * AMX{56};
GArpT4 = GArpN0 * AMX{56};
GActT5 = GActinN0 * AMX{57};
GArpT5 = GArpN0 * AMX{57};

doThymT = 1;
KaThymT2 = Ka_Thymosin_LTP;
KdThymT2 = Kd_Thymosin_LTP;
KaThymT3 = Ka_Thymosin;
KdThymT3 = Kd_Thymosin;

%{
doActT = 0;
doArpT = 0;
T2start = 40000;
T3start = 80000;
T4start = 120000;
T5start = 140000;

GActT2 = GActinN0 + (GActinN0 / 4);
GArpT2 = GArpN0 + (GArpN0 / 4);

GActT3 = GActinN0;
GArpT3 = GArpN0;

GActT4 = GActinN0 - (GActinN0 / 4);
GArpT4 = GArpN0 - (GArpN0 / 8);

GActT5 = GActinN0 - (GActinN0 / 2);
GArpT5 = GArpN0 - (GArpN0 / 4);
%}
%----------------------------------------


%----------------------------------------
% PREALLOCATE COUNTERS
nT = 1:Nsteps;
NumFils = nT;	% Number of branch filaments
NumFAct = nT;	% Number FActins
NumGAct = nT;	% Number GActins
NumCOFACT = nT;	% Number of Depoly events
NumARPACT = nT;	% Number of Poly events
NumdelFi = nT;

ArpBRsum = nT;
NumAct_uM = nT;
NumfKa = nT;
NumFnowFtot = nT;
NumFArp = nT;	% Number of Arp in Filaments
NumGArp = nT;	% Number of Free Arp

NumHullV = zeros(3,round(Nsteps/LivePlotMod));

Num_nmmf = nT;
NumFded = nT;

NumGActinN0 = nT;	% Number Unbound Actins
NumTHYM = nT;       % Number Unbound Thymosin
NumTHYMACT = nT;	% Number Thymosin+GActins
%----------------------------------------



%-------------------------------------------%
% Linear Range Transformation
%-------------------------------------------%
% Equation to linear transform range
% [a b] to range [c d] and solve for [x]:
% F(x) = c*(1-(x-a)/(b-a)) + d*((x-a)/(b-a))
%-------------------------------------------%
%{
global sja;
global sjb;
global sjc;
global sjd;

SJKb = StartMonos;

sja=AMX{23};
sjb=SJKb*AMX{24};
sjc=AMX{25};
sjd=AMX{26};

linsc = @(jx) (sjc*(1-(jx-sja)/(sjb-sja)) + sjd*((jx-sja)/(sjb-sja)));

global ska;
global skb;
global skc;
global skd;
ska=AMX{48};
skb=AMX{49};
skc=AMX{50};
skd=AMX{51};

sigsc = @(jx) (1/(1+ exp(-((skc*(1-(jx-ska)/(skb-ska)) + skd*((jx-ska)/(skb-ska)))))));


%}
%----------------------------------------

%----------------------------------------
% MATH TIDBITS
% miscfun = @(xyz,tta) ([1 0 0; 0 cos(tta) sin(tta); 0 -sin(tta) cos(tta)] * xyz);
%-------------------------------------------%
%{
%-------------------------------------------%

keyboard
%%




close all
clear ArpCurve1
rbg1=.01;rbg2=.99;rbg3=.01;

TheoMaxFact = ceil(SPYhZN / nm_per_p);
ScFil = ceil(TheoMaxFact/2);
Nac = TheoMaxFact;
CofR = .002;
CofS = 30;


figure(44)
for CofS = 10:10:100

	nn = 1;
	for Flength = 0:1:Nac
	ArpCurve1(nn) = CofR * exp((Flength-CofS)/ScFil);
	nn=nn+1;
	end

ph=plot(ArpCurve1); hold on;
set(ph,'Color',[rbg1 rbg2 rbg3])
rbg1=rbg1+.08;
rbg2=rbg2-.08;
rbg3=rbg3+.08;
drawnow; 
hold on
pause(.2);
end






%%
close all
clear ArpCurve1
rbg1=.01;rbg2=.99;rbg3=.01;

ska=0;
skb=700;
skc=-7;
skd=7;

sigsc = @(jx) (1/(1+ exp(-((skc*(1-(jx-ska)/(skb-ska)) + skd*((jx-ska)/(skb-ska)))))));


figure(44)
for ArpScalr = 10:20:200

	nn = 1;
	for Flength = 0:1:skb
	ArpCurve1(nn) = ArpR * sigsc(Flength+ArpScalr);
	%ArpCurve1(nn) = ArpR * (1/(exp(1/170*(1190-14*(Flength+ArpScalr)))+1));
	nn=nn+1;
	end

ph=plot(ArpCurve1); hold on;
set(ph,'Color',[rbg1 rbg2 rbg3])
rbg1=rbg1+.08;
rbg2=rbg2-.08;
rbg3=rbg3+.08;
drawnow; 
hold on
pause(.2);
end


%%

close all
clear ArpCurve1
rbg1=.01;rbg2=.99;rbg3=.01;

ska=0;
skb=700;
skc=-7;
skd=7;

figure(46)
for ArpScalr = 10:20:200

	nn = 1;
	for Flength = 0:1:skb
	ArpCurve1(nn) = ArpR * (1/(1+ exp(-(7*(Flength+ArpScalr)/300-7))));
	nn=nn+1;
	end

ph=plot(ArpCurve1); hold on;
set(ph,'Color',[rbg1 rbg2 rbg3])
rbg1=rbg1+.08;
rbg2=rbg2-.08;
rbg3=rbg3+.08;
drawnow; pause(.2);
end







%%
% (ArpR * (1/(exp(1/170*(1190-14*(Flength+ArpScalr)))+1))) > rv(2) && ...



%-------------------------------------------%

clear ArpCurve1
close all;
rbg1=.01;rbg2=.99;rbg3=.01;
sja=0;
sjb=300*1;
sjc=0;
sjd=1;
linsc = @(jx) (sjc*(1-(jx-sja)/(sjb-sja)) + sjd*((jx-sja)/(sjb-sja)));

figure(44)
for CofSc = -30:10:30

	nn = 1;
	for Flength = 0:1:sjb
	ArpCurve1(nn) = CofR * exp(linsc(Flength-CofSc));
	nn=nn+1;
	end

ph=plot(ArpCurve1); hold on;
set(ph,'Color',[rbg1 rbg2 rbg3])
rbg1=rbg1+.08;
rbg2=rbg2-.08;
rbg3=rbg3+.08;
drawnow; pause(.5);
end

%-------------------------------------------%
clear ArpCurve1


rbg1=.1;rbg2=.9;rbg3=.1;
figure(43)
for sj = 0:.5:3
	
	linsc = @(jx) (sjc*(1-(jx-sja)/(sjb-sja)) + sjd*((jx-sja)/(sjb-sja)));

	nn = 1;
	for FLEG = 1:1:300
	ArpCurve1(nn) = ArpR * exp(linsc(FLEG-ArpScalar));
	nn=nn+1;
	end


ph=plot(ArpCurve1); hold on;
set(ph,'Color',[rbg1 rbg2 rbg3])
rbg1=rbg1+.1;
rbg2=rbg2-.1;
rbg3=rbg3+.1;

end


%-------------------------------------------%
% sjd
% as sjd increases from 0 to 4, the curve shifts from linear to exponential
% as such, a good value for sjd = 3

% sjc
% as sjc decreases from 0 to -4, the curve also becomes more exponential
% a good value for sjc = -1

% logistic sigmoid function (from -5 to 5)
% f(x) = 1 / (1 + e^-x)


clear ArpCurve1
rbg1=.1;rbg2=.9;rbg3=.1;
sja= 0;
sjb= 300;
sjc= -5;
sjd= 5;

figure(43)
for sj = 0:.5:3


linsc = @(jx) (sjc*(1-(jx-sja)/(sjb-sja)) + sjd*((jx-sja)/(sjb-sja)));

	nn = 1;
	for FLEG = 1:2:300
	ArpCurve1(nn) = (ArpR*(1/(1+ exp(-linsc(FLEG)))));
	nn=nn+1;
	end

ph=plot(ArpCurve1); hold on;
set(ph,'Color',[rbg1 rbg2 rbg3])
rbg1=rbg1+.1;
rbg2=rbg2-.1;
rbg3=rbg3+.1;
end


%}
%{

% 5e8 in cells; 1e4 synapses per neuron; 5e8/1e4 = 5e4

Polymerization Rate
(+)end: .012 N/µM*ms
** (12 N/mM*ms) (12 N/µM*s)
** thus at 1 µM free ATP-actin, .012 subunits will be added to the (+)end per ms
** at .1 mM free ATP-actin, 1.2 subunits will be added to the (+)end per ms

Depolymerization Rate
(+)end: 1.4 N/s
(-)end: 0.8 N/s
** dissociation is independent of free actin concentration

To find the critical concentration (Cc) for growth we set the two rate 
equations equal to each other:

12.0/µM*s = 1.4/s
12/1.4 = 1/µM
(+)Cc = .12 µM
(-)Cc = .6 µM





Actin monomer size:
- 5.5 nm x 5.5 nm x 3.5 nm

Factin filaments
- ?-helix composed of 2 strands of subunits
- 28 subunits (14 in each strand) in 1 full 360? turn 
- 180? turn: 36 nm
- 360? turn: 72 nm
- 28 subunits per 72 nm
--------
* For 1 filament to span a 1000 nm spine would require:
-- (mean spine length: 1.0 µm or 1000 nm)
-- each monomer spans 5.1 nm
-- every 5.1 nm requires 2 monomers (cuz double helix)
-- 13.89 turns (~14 turns)
-- 388.89 actin monomers (~400 monomers)
-- 195 monomers per strand
--------

The ATP-binding cleft of actin faces the (+) end of the filament
(+)end grows
(-)end shrinks

Polymerization Rate
(+)end: ~12.0 subunits/µM*s
(-)end: ~1.3 subunits/µM*s
** thus if there is 1 µM of free ATP-Gactin then 12 subunits will be added 
to the (+)end per second and 1.3 subunits will be added to the (-)end every second

Depolymerization Rate
(+)end: ~1.4 subunits/s
(-)end: ~0.8 subunits/s
** dissociation is independent of free Gactin concentration

To find the critical concentration (Cc) for growth we set the two rate 
equations equal to each other:

12.0/µM*s = 1.4/s
12/1.4 = 1/µM
(+)Cc = .12 µM
(-)Cc = .6 µM

Thus when the free actin concentration >.12 µM filaments will grow at the (+) end 
and when >.6 µM filaments will grow at the (-) end too.


----
There are 2 isoforms of actin
- ?-actin & ?-actin

?-actin is enriched in dendritic spines and builds filament stress fibers

G-actin = monomeric "globular" actin
F-actin = filamentous actin

Each actin molecule contains a Mg2+ ion complexed with ATP or ADP

----
Cytosolic concentration of actin in cells ranges from .1 to .5 mM

Actin makes up 1-5% of all cellular protein
(this suggests there is 10 mM total proteins in cells)


A typical cell contains around 5e8 actin molecules

The average number of synapses per neuron was 1e4
http://www.ncbi.nlm.nih.gov/pubmed/2778101

5e8 / 1e4 = 5e4
A typical spine contains 5e4 actin molecules


Total spines volume averaged 0.09 ?m^3 and ranged from 0.01 to 0.38 ?m^3 (n = 133). 
The mode peak value was 0.06 ?m^3. 
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2518053/

0.09 ?m^3 = 9e-17 L
5e4 / 9e-17 L = 5.55e20
5.55 N/L * (1 L / 6e23 N) = 0.000925 M = .9 mM
.5 - .9 mM actin / spine

%----------------------------------------
%		STATIC FINAL FIGURES SETUP
%----------------------------------------
%Fh1 = FigSetup(33);
%}
%----------------------------------------



%----------------------------------------
%	  ANIMATED REAL-TIME FIGURE SETUP
%----------------------------------------
rot0 = [5 0];
rot1 = rot0;
azel0 = [-32 12];
sz = [5e-3 6e-3 1.5 1.4];
Fh2Live = FigSetup(2,sz);
set(gcf,'Color',[1,1,1])
AxLims = [-dims(4) dims(4) -dims(5) dims(5) 0 dims(2)].*1.2;



%--------------------
%{
%--------------------
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%--------------------

nT=1;
Fh2L = FigSetup(2,sz);
AxLims = [-dims(4) dims(4) -dims(5) dims(5) 0 dims(2)].*1.2;
azel = [-32 12];
rot=[0 0];


%--------------------
figure(Fh2L)
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
	%hold on;
%--------------------
set(gca,'XTickLabel', sprintf('%.1f|',nT),'FontSize',10)
hold off;
%--------------------
figure(Fh2L)
subplot('Position',[.03 .03 .65 .95]), 

ph11a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
	axis(AxLims); axis vis3d;
	view(azel+rot)
	xlabel('X');ylabel('Y');zlabel('Z');
	set(gcf,'Color',[1,1,1])
	grid off
	hold on;
ph11b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
	axis(AxLims)
	set(ph11a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
	set(ph11b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
	hold on;
ph11c = plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]');
	axis(AxLims); axis vis3d;
	view(azel+rot)
	set(ph11c,'LineStyle','-','Color',[.7 .7 .7],'LineWidth',.1);
	hold on;


PhCx = {ph11c,ph11a,ph11b,ph12a,ph12b};

LiveAPlot(PhCx,nT,Actin,PSDTips,SPYTips,rot1,azel0,Fh2L);
rot1 = rot1 + rot0;
%}
%----------------------------------------


doProfile=0;
if doProfile
	profile on
end

doTicToc=0;
if doTicToc
    tic
end


aN = 1;
progressbar(0) % Initilize progress bar
%%
%===============================================================%
%						MAIN OUTER LOOP
for nT = 1:Nsteps
%---------------------------------------------------------------%

	%----------------------------------------------
	if doActT
	if nT==T2start; GActinN0 = GActT2; end;
	if nT==T3start; GActinN0 = GActT3; end;
	if nT==T4start; GActinN0 = GActT4; end;
	if nT==T5start; GActinN0 = GActT5; end;
	end
	if doArpT
	if nT==T2start; GArpN0 = GArpT2; end;
	if nT==T3start; GArpN0 = GArpT3; end;
	if nT==T4start; GArpN0 = GArpT4; end;
	if nT==T5start; GArpN0 = GArpT5; end;
	end
	if doThymT
	if nT==T2start; Ka_Thymosin = KaThymT2; end;
	if nT==T2start; Kd_Thymosin = KdThymT2; end;
	if nT==T3start; Ka_Thymosin = KaThymT3; end;
	if nT==T3start; Kd_Thymosin = KdThymT3; end;
	end
	
	NFact = numel(Actin(:,1));	% Current #of Filaments ("Nfi" used below loop)
	if (mod(nT,LivePlotMod)==0); disp([nT NFact]); end;
	ACTdepoly = 0;
	COFdepoly = 0;
	ARPpoly = 0;
	ACTpoly = 0;
	%----------------------------------------------
	
	
	%----------------------------------------------
	Ori_xyz = [Actin(:,3) Actin(:,6) Actin(:,9)];
	Tip_xyz = [Actin(:,4) Actin(:,7) Actin(:,10)];

	% radial distance to spine shaft membrane
	XYtipLoc = sqrt(Tip_xyz(:,1).^2 + Tip_xyz(:,2).^2);
	XYoriLoc = sqrt(Ori_xyz(:,1).^2 + Ori_xyz(:,2).^2);

	ZtipInHead = Tip_xyz(:,3) >= SPYhZS;
	ZtipInNeck = Tip_xyz(:,3) < SPYhZS;
	ZoriInHead = Ori_xyz(:,3) >= SPYhZS;
	ZoriInNeck = Ori_xyz(:,3) < SPYhZS;

	XYneckOut = XYtipLoc > SPYnXY;		% Logical array of filaments beyond SPYnXY
	XYheadOut = XYtipLoc > SPYhXY;		% Logical array of filaments beyond SPYhXY
	ZtopOut = Tip_xyz(:,3) > SPYhZN;	% Logical array of filaments above SPYhZN
	ZbotOut = Tip_xyz(:,3) < 0;			% Logical array of filaments below zero
	
	XYOneckOut = XYoriLoc > SPYnXY;		% Logical array of filaments beyond SPYnXY
	XYOheadOut = XYoriLoc > SPYhXY;		% Logical array of filaments beyond SPYhXY
	ZOtopOut = Ori_xyz(:,3) > SPYhZN;	% Logical array of filaments above SPYhZN
	ZObotOut = Ori_xyz(:,3) < 0;		% Logical array of filaments below zero

	TipOut =((XYneckOut & ZtipInNeck)+(XYheadOut & ZtipInHead)+ZtopOut+ZbotOut)>0;
	OriOut =((XYOneckOut & ZoriInNeck)+(XYOheadOut & ZoriInHead)+ZOtopOut+ZObotOut)>0;

	TipOK = ~TipOut;
	%----------------------------------------------

    Actin(:,20) = Actin(:,1) .* nm_per_p;
	

    
    ActNon0 = (Actin(:,1)>0);	% assures no negative actin values
	Actin(:,1) = Actin(:,1) .* (Actin(:,1)>0);
	Actin(:,1) = Actin(:,1) .* (Actin(:,10) > -20);
	
	
	% MCMC requires a random filament at each decision step
	rN1 = randi(NFact,1,NFact); % actin polymerization
    rN2 = randi(NFact,1,NFact); % arp branching
    rN3 = randi(NFact,1,NFact); % actin depolymerization
    rN4 = randi(NFact,1,NFact); % cofilin severing
	%=============================================================================%
	%						MAIN INNER LOOP
	for fN=1:NFact
	%-----------------------------------------------------------------------------%
	
        rv=rand(9,1); % Generate a few random vaules from uniform{0:1}


        aN = rN1(fN); % Get random filament
		%=================================================================%
		% POLYMERIZATION
        if  (fKa > rv(1)) && TipOK(aN) && (GActinN>1)
		%----------------------------------------
        
            if fKa > NFact
                Actin(aN,1) = Actin(aN,1) + 2;
                ACTpoly = ACTpoly + 2;
            else
                Actin(aN,1) = Actin(aN,1) + 1;
                ACTpoly = ACTpoly + 1;
            end
		
		%----------------------------------------
        end
		%=================================================================%
		



		aN = rN2(fN);
		%=================================================================%
		% BRANCHING
 		if  (ArpBR(aN) > rv(2)) && (GActinN>ArpAdd)
		%----------------------------------------
			fNf = NFact+1;
			%-------------------
			
			Actin(fNf,1) = ArpAdd;		% create branch: add N(ArpAdd) subunits to new branch
			Actin(fNf,11)=Actin(aN,12);	% tag branch with MomID
			Actin(fNf,12) = TagN;		% tag branch with ID
			Actin(fNf,14) = nT;			% tag branch with Born time (nT)
			TagN = TagN+1;
			
			Nmono = Actin(aN,1);		% N monomers in mother filament
			Rmm = ceil(Nmono * rand);	% Random mother monomer (== vector_length) along segment
						
			% TIP of current branch		(not actual tip, just rotational point) 
			Ct_x = (Rmm) .* nm_per_p * sin(Actin(aN,8)) * cos(Actin(aN,2)) + Actin(aN,3);
			Ct_y = (Rmm) .* nm_per_p * sin(Actin(aN,8)) * sin(Actin(aN,2)) + Actin(aN,6);
			Ct_z = (Rmm) .* nm_per_p * cos(Actin(aN,8)) + Actin(aN,9);

			% ORIGIN of current branch
			Co_x = Actin(aN,3);	% X origin (old branch)
			Co_y = Actin(aN,6);	% Y origin (current branch)
			Co_z = Actin(aN,9);	% Z origin (current branch)
			
			%=====================================================%	
			Po = [Co_x;Co_y;Co_z];
			Pt = [Ct_x;Ct_y;Ct_z];
			Pv = Pt-Po;
			%-------------------
			tL = sqrt(sum((Pv).^2));		% Length of vector PoPt (aka Pv)
			Pu = (Pv) ./ tL;				% Unit vector of Pv

			tTheta = acos(Pu(3))  +TPi;		% angle theta
			tPhi = atan2(Pu(2),Pu(1)) +PPi;	% angle phi

			x = sin(tTheta) * cos(tPhi);
			y = sin(tTheta) * sin(tPhi);
			z = cos(tTheta);
			Pr = [x;y;z]+Pt;
			%-------------------
			
			Otta = Ov(randi(26,1));		% Random rotational angle (theta)
			
			%Pn = RotateVertex(Pr(1),Pr(2),Pr(3),Pt(1),Pt(2),Pt(3),Po(1),Po(2),Po(3),...
			%			 Pv(1),Pv(2),Pv(3),Otta);

            % tictoc: 23.2298   24.2325
            dv = [Pt(1);Pt(2);Pt(3)] - [Po(1);Po(2);Po(3)];
            u=dv(1);  v=dv(2);  w=dv(3);
            Pn = RotateVec(Pr(1),Pr(2),Pr(3),Pt(1),Pt(2),Pt(3),u,v,w,Otta);
            % Pn = LineRota(Pr(1),Pr(2),Pr(3),Pt(1),Pt(2),Pt(3),u,v,w,Otta);
            %--------
			


            %-------------------
			Po2 = Pt;
			Pt2 = Pn;
			Pv2 = Pt2-Po2;
			tL2 = sqrt(sum((Pv2).^2));		% Length of vector PoPt (aka Pv)
			Pu2 = (Pv2) ./ tL2;				% Unit vector of Pv
			tTheta2 = acos(Pu2(3));			% angle theta	+TPi;
			tPhi2 = atan2(Pu2(2),Pu2(1));	% angle phi		+PPi;
			%=====================================================%
			%if ((Pv(3)<0) || (nT>20000)); keyboard; end;
            
			% New branch Angles
			Actin(fNf,2) = tPhi2;
			Actin(fNf,8) = tTheta2;
			
			% New branch Origin
			Actin(fNf,3) = Po2(1);
			Actin(fNf,6) = Po2(2);
			Actin(fNf,9) = Po2(3);
			
			% New branch Tip
			Actin(fNf,4) = (ArpAdd) .* nm_per_p * sin(tTheta2) * cos(tPhi2) + Po2(1);
			Actin(fNf,7) = (ArpAdd) .* nm_per_p * sin(tTheta2) * sin(tPhi2) + Po2(2);
			Actin(fNf,10) = (ArpAdd) .* nm_per_p * cos(tTheta2) + Po2(3);
			
            
			%---------------------------------------------%
			ARPpoly = ARPpoly + ArpAdd;
			
		%----------------------------------------
		end
		%=================================================================%

		
		
		%=================================================================%
		% PSD-BRANCHING 
		%{
		Xtip1 = Tip_xyz(:,1) >= 0;
		Ytip1 = Tip_xyz(:,2) >= 0;
		Xtip2 = Tip_xyz(:,1) <= 250;
		Ytip2 = Tip_xyz(:,2) <= 250;
		Ztip1 = Tip_xyz(:,3) >= (SPYhZS+200);
		PSDtips = (sum([Xtip1 Ytip1 Xtip2 Ytip2 Ztip1],2))==5;

		ZtopOu = Tip_xyz(:,3) > SPYhZN+50;	% Logical array of filaments above SPYhZN
		TipOu =((XYneckOut & ZtipInNeck)+(XYheadOut & ZtipInHead)+ZtopOu+ZbotOut)>0;
		TipOKk = ~TipOu .* PSDtips;

		%if nT > 10000; keyboard; end;
		%----------------------------------------------
		
		
 		if  (nT>10000) && (nT<10500) && PSDtips(aN) %&& (GActinN>ArpAdd) 
		%----------------------------------------
			fNf = NFact+1;
			%-------------------
			
			Actin(fNf,1) = ArpAdd;		% create branch: add N(ArpAdd) subunits to new branch
			Actin(fNf,11)=Actin(aN,12);	% tag branch with MomID
			Actin(fNf,12) = TagN;		% tag branch with ID
			Actin(fNf,14) = nT;			% tag branch with Born time (nT)
			TagN = TagN+1;
			
			Nmono = Actin(aN,1);		% N monomers in mother filament
			Rmm = ceil(Nmono * rand);	% Random mother monomer (== vector_length) along segment
						
			% TIP of current branch		(not actual tip, just rotational point) 
			Ct_x = (Rmm) .* nm_per_p * sin(Actin(aN,8)) * cos(Actin(aN,2)) + Actin(aN,3);
			Ct_y = (Rmm) .* nm_per_p * sin(Actin(aN,8)) * sin(Actin(aN,2)) + Actin(aN,6);
			Ct_z = (Rmm) .* nm_per_p * cos(Actin(aN,8)) + Actin(aN,9);

			% ORIGIN of current branch
			Co_x = Actin(aN,3);	% X origin (old branch)
			Co_y = Actin(aN,6);	% Y origin (current branch)
			Co_z = Actin(aN,9);	% Z origin (current branch)
			
			%=====================================================%	
			Po = [Co_x;Co_y;Co_z];
			Pt = [Ct_x;Ct_y;Ct_z];
			Pv = Pt-Po;
			%-------------------
			tL = sqrt(sum((Pv).^2));		% Length of vector PoPt (aka Pv)
			Pu = (Pv) ./ tL;				% Unit vector of Pv

			tTheta = acos(Pu(3))  +TPi;		% angle theta
			tPhi = atan2(Pu(2),Pu(1)) +PPi;	% angle phi

			x = sin(tTheta) * cos(tPhi);
			y = sin(tTheta) * sin(tPhi);
			z = cos(tTheta);
			Pr = [x;y;z]+Pt;
			%-------------------
			
			Otta = Ov(randi(26,1));		% Random rotational angle (theta)
			
			Pn = RotateVertex(Pr(1),Pr(2),Pr(3),Pt(1),Pt(2),Pt(3),Po(1),Po(2),Po(3),...
						 Pv(1),Pv(2),Pv(3),Otta);

			%-------------------
			Po2 = Pt;
			Pt2 = Pn;
			Pv2 = Pt2-Po2;
			tL2 = sqrt(sum((Pv2).^2));		% Length of vector PoPt (aka Pv)
			Pu2 = (Pv2) ./ tL2;				% Unit vector of Pv
			tTheta2 = acos(Pu2(3));			% angle theta	+TPi;
			tPhi2 = atan2(Pu2(2),Pu2(1));	% angle phi		+PPi;
			%=====================================================%
			%if ((Pv(3)<0) || (nT>20000)); keyboard; end;
            
			% New branch Angles
			Actin(fNf,2) = tPhi2;
			Actin(fNf,8) = tTheta2;
			
			% New branch Origin
			Actin(fNf,3) = Po2(1);
			Actin(fNf,6) = Po2(2);
			Actin(fNf,9) = Po2(3);
			
			% New branch Tip
			Actin(fNf,4) = (ArpAdd) .* nm_per_p * sin(tTheta2) * cos(tPhi2) + Po2(1);
			Actin(fNf,7) = (ArpAdd) .* nm_per_p * sin(tTheta2) * sin(tPhi2) + Po2(2);
			Actin(fNf,10) = (ArpAdd) .* nm_per_p * cos(tTheta2) + Po2(3);
			
            
			%---------------------------------------------%
			ARPpoly = ARPpoly + ArpAdd;
			
		%----------------------------------------
		end
		
		%=================================================================%
		%=================================================================%
		% PSD POLYMERIZATION
        if  (nT>10000) && (nT<20000) && (BeKa > rv(1)) && TipOKk(aN) %&& (GActinN>1)
		%----------------------------------------
        
            if BeKa > NFact
                Actin(aN,1) = Actin(aN,1) + 2;
                ACTpoly = ACTpoly + 2;
            else
                Actin(aN,1) = Actin(aN,1) + 1;
                ACTpoly = ACTpoly + 1;
            end
		
		%----------------------------------------
		end
		%}
		%=================================================================%
		
		
		
		
		
		aN = rN3(fN);
		%=================================================================%
		% ACTIN PASSIVE DEPOLYMERIZATION		
		if (fKd > rv(3))
		%----------------------------------------
			Actin(aN,1) = Actin(aN,1)-1;
			ACTdepoly = ACTdepoly + 1;
		%----------------------------------------
		end
		

		aN = rN4(fN);
		%=================================================================%
		% COFILIN ASSISTED DEPOLYMERIZATION
		if  ( CofR * exp((Actin(aN,1)-CofS)/ScFil) ) > rv(4)
		%----------
			
			Actin(aN,1) = Actin(aN,1)-CofN;
		
			% Depoly from origin (Tag)
			if rv(5) < delOr
			Actin(aN,13) = 1;
			end
			
			COFdepoly = COFdepoly + CofN;
		%----------
		end
		
        
        aN = fN;
		%=================================================================%
		% RAPID DEPOLYMERIZATION
		%-----------------------------------
		if (Actin(aN,13) && (fdKd > rv(6))) || OriOut(aN)
		%----------
			Actin(aN,1) = Actin(aN,1)-1;
			ACTdepoly = ACTdepoly + 1;
		%----------
		end
		%=================================================================%

		
		%=================================================================%
        % THYMOSIN+ACTIN ASSOCIATION
        if (KaThy > rv(8)) && (GActinN>1) && (THYM_N>1)

            THYM_N = THYM_N - 1;
            GActinN0 = GActinN0 - 1;
            THYM_ACT_N = THYM_ACT_N + 1;
            % disp(THYM_N); disp(GActinN0); disp(THYM_ACT_N)
        end

        % THYMOSIN+ACTIN DISSOCIATION
        if (KdThy > rv(9))  && (THYM_ACT_N>1)

            THYM_N = THYM_N + 1;
            GActinN0 = GActinN0 + 1;
            THYM_ACT_N = THYM_ACT_N - 1;
            % disp(THYM_N); disp(GActinN0); disp(THYM_ACT_N)
        end
		%=================================================================%
		
	
    %aN = fN;
	%=================================================================%
	%						ADJUST RATE VALUES
	%-----------------------------------------------------------------
	if ~mod(aN,10)
	FActinN = sum(Actin(:,1));				% Current #of FActins
	GActinN = GActinN0 - FActinN;			% Current #of GActins

	Act_uM = GActinN / SpyV *(1/MOL)*1e6;	% Actin uM
	fKa = TKa * Act_uM * dT;				% Ka Fil ON-Rate
    
    FArpN = numel(Actin(:,1));              % #of F-Arp
    GArpN = GArpN0 - FArpN;                 % #of G-Arp
                
    Arp_uM = GArpN / SpyV *(1/MOL)*1e6;    % Arp uM
    nmmf = Actin(:,20);						% Length of filaments (nm)
    ArpBR = Arp_Sc * ((Arp_uM/1000) .* (nmmf/1000));
	%ArpBR = 5.4e-4 * Arp_uM^3 * dT .* (nmmf>0);


    THYM_uM = THYM_N / SpyV *(1/MOL)*1e6;         % Thymosin concentration in spine (uM)
    THYM_ACT_uM = THYM_ACT_N / SpyV *(1/MOL)*1e6; % Thymosin+Actin concentration in spine (uM)
    KaThy = (dT*Ka_Thymosin * THYM_uM * Act_uM);  % Thymosin+Actin Association
    KdThy = (dT*Kd_Thymosin * THYM_ACT_uM);       % Thymosin+Actin Dissociation

	end
	%=================================================================%
	%%
    % keyboard
    %%
	
	
	%if nT > 5000; keyboard; end;
	%-----------------------------------------------------------------------------%
	end
	%						MAIN INNER LOOP
	%=============================================================================%
	%%
	
	
	
	
	
	
	
	%=================================================================%
	%						ADJUST RATE VALUES
	%-----------------------------------------------------------------

	FActinN = sum(Actin(:,1));				% Current #of FActins
	GActinN = GActinN0 - FActinN;			% Current #of GActins
		if GActinN<0;GActinN=0;end;

	Act_uM = GActinN / SpyV *(1/MOL)*1e6;  % Actin uM
	fKa = TKa * Act_uM * dT;				% Ka Fil ON-Rate
    
    FArpN = numel(Actin(:,1));              % #of F-Arp
    GArpN = GArpN0 - FArpN;                 % #of G-Arp
		if GArpN<0;GArpN=0;end;
                
    Arp_uM = GArpN / SpyV *(1/MOL)*1e6;    % Arp uM
    nmmf = Actin(:,20);						% Length of filaments (nm)
    ArpBR = Arp_Sc * ((Arp_uM/1000) .* (nmmf/1000));

	%=================================================================%
	
	% LTP ADJUSTMENTS
	% if doActLTP && (nT >= ActLTPstart) && (nT <= ActLTPend); ArpBR = ArpRLTP; end;
	if doActLTP && (nT >= ActLTPstart) && (nT <= ActLTPend); 
		keyboard
		ArpBR = ArpRLTP; 
	end;
	
	% DELETE ORIGINAL FILAMENTS
	if doDelOrig && (nT == doDelOrigT); Actin(1:NStFils,:) = []; end; 
	%=================================================================%
	

	%==================================================%
	%	PROCESS STORE AND REMOVE DEL BRANCHES
	%--------------------------------------------------%
	
	% Filament Lifetime
	FLif = Actin(:,16) + 1;
	Actin(:,16) = FLif;
	
	% Mean Length
	Actin(:,18) = (Actin(:,18).*(FLif-1) + Actin(:,1))./FLif;
	
	%Longest Ever Length
	MaxL = Actin(:,1) > Actin(:,17);
	Actin(MaxL,17) = Actin(MaxL,1); 
	
	
	% Delete Filaments With No Actin
	delFi = (Actin(:,1)<1);
	Actin(delFi,15) = nT;	% Death Timestamp
	ACTs = cat(1,ACTs,Actin(delFi,:));
    DiedACTs = cat(1,DiedACTs,Actin(delFi,:));
	Actin(delFi,:) = [];
	
	Nfi = numel(Actin(:,1));	% Current #of Filaments ("NFact" used above loop)
	%--------------------------------------------------%
	
	
	%==================================================%
	%	COMPUTE XYZ TIP LOCATION
	%--------------------------------------------------%
	% MATH - branch XYZ tip coordinates
	Actin(:,4) = Actin(:,1) .* nm_per_p .* sin(Actin(:,8)) .* cos(Actin(:,2)) + Actin(:,3);
	Actin(:,7) = Actin(:,1) .* nm_per_p .* sin(Actin(:,8)) .* sin(Actin(:,2)) + Actin(:,6);
	Actin(:,10) = Actin(:,1) .* nm_per_p .* cos(Actin(:,8)) + Actin(:,9);
	%--------------------------------------------------%
	
	








    %{
	if mod(nT,LivePlotMod) == 0
		LivePlot(Fh2Live,nT,Actin,inPSD,rot1,azel0,AxLims);
		TriPlot(Fh2Live,nT,Actin,inPSD,rot1,azel0,AxLims,SPYhZN,SPYhZS,SPYnXY,SPYhXY);
		rot1 = rot1 + rot0;
	end
	%}
	%==================================================%
	%				LIVE PLOT
	%--------------------------------------------------%
    if doLiveTipPlot
    if ~mod(nT,LiveTipMod)
        MainLivePlot(Fh2Live,nT,Actin,inPSD,rot1,azel0,AxLims,...
                    SPYhZN,SPYhZS,SPYnXY,SPYhXY,doLiveHullPlot);
		rot1 = rot1 + rot0;
    end
    end
	%--------------------------------------------------%
	
	%==================================================%
	%			TriHull Spine Volume
	%--------------------------------------------------%
    if ~doLiveHullPlot
	if ~mod(nT,LiveHullMod)
		%HullV = TriHull(nT,Actin,SPYhZN,SPYhZS,SPYnXY,SPYhXY);
		%NumHullV(:,nT/TriHullMod) = HullV;

        ActArpPlot(Fh2Live,nT,Nsteps,NumGArp,NumFArp,NumFAct,NumGAct,NumfKa,ArpBRsum,Num_nmmf)

	end
    end
	%--------------------------------------------------%
    if ~mod(nT,LiveTipMod)
        progressbar(nT/Nsteps)
    end
	








	%if nT == 1000; pause(5); end; if nT == 5000; pause(5); end;

	%==================================================%
	%				SAVE TipMatrix
	%--------------------------------------------------%
	if SvTpMx && (nT >SaveTipsAfter) && (mod(nT,SaveTipsRate) == 0)

		ActMx = TipMatrix(Fh2Live,doTM,nT,Actin,dims,AcMx,SPYH,rot1,azel0);
            BTs{numel(BTs)+1} = ActMx;
            AFMx{numel(AFMx)+1} = Actin;

	end
	%--------------------------------------------------%
	
	
	
	%==================================================%
	%				Counters
	%--------------------------------------------------%
	if doActCounts
		
		NumFils(nT) = Nfi;						% Number of branch filaments
		NumFAct(nT) = FActinN;					% Number FActins
		NumGAct(nT) = GActinN;					% Number GActins
		NumCOFACT(nT) = COFdepoly+ACTdepoly;	% Number of Depoly events
		NumARPACT(nT) = ACTpoly+ARPpoly;		% Number of Poly events
		NumdelFi(nT) = sum(delFi);
		
		ArpBRsum(nT) = sum(ArpBR);
		NumAct_uM(nT) = Act_uM;
		NumfKa(nT) = fKa;
		NumFnowFtot(nT) = Nfi / Actin(end,12);
		NumFArp(nT) = FArpN;
		NumGArp(nT) = GArpN;
		Num_nmmf(nT) = sum(nmmf);
		NumFded(nT) = numel(DiedACTs(:,1));

        NumGActinN0(nT) = GActinN0;     % Number Unbound Actins
		NumTHYM(nT)     = THYM_N;       % Number Unbound Thymosin
		NumTHYMACT(nT)  = THYM_ACT_N;   % Number Thymosin+GActins
	end
	%--------------------------------------------------%
    
    
    
    
    %==================================================%
	%				LIVE PLOT
	%--------------------------------------------------%
	%{
    if ~mod(nT,LivePlotMod)
        if nT > 1000
			
			
			keyboard
			NumData = {NumFils, NumFAct, NumGAct, NumCOFACT, NumCOFACT, NumARPACT...
				,NumdelFi, ArpBRsum, NumAct_uM, NumfKa, NumfKa, NumFnowFtot...
				,NumFArp, NumGArp, Num_nmmf, Num_nmmf};
			
			
    FilPlotAlone(Fh2Live,nT,Actin,inPSD,rot1,azel0,AxLims,SPYhZN,SPYhZS,SPYnXY,SPYhXY);
        rot1 = rot1 + rot0;
		%progressbar(nT/Nsteps)
		
	ActArpPlot(Fh2Live,nT,Nsteps,NumGArp,NumFArp,NumFAct,NumGAct,NumfKa,ArpBRsum,Num_nmmf);
		end
	end
	
	%}
	%--------------------------------------------------%
	
	
	
	
% if nT > 80000; keyboard; end;
%-----------------------------------------------------------------------------%
end
%						MAIN OUTER LOOP
%=============================================================================%
%%

if doTicToc
    TicToc = toc;
    disp(TicToc)
    error(strcat('Loop time was: ', num2str(TicToc)))
end

progressbar(1) % closes progress bar
if doProfile
profile viewer
end


% Lifetime of remaining filaments
Actin(:,15) = (nT + (nT-Actin(:,14)));      % Artificial Death Time
Actin(:,16) = (Actin(:,15) - Actin(:,14));  % Lifetime (nTdied - nTborn)
ACTs = cat(1,ACTs,Actin);

remove0 = (ACTs(:,16)<2);   
ACTs(remove0,:) = [];      % Remove filaments that lasted less than nT<2

DiedACTs(1:2,:) = [];



%%
% TIP DENSITY ASSESSMENT
%==================================================%
%{
clear in3db; clear xxyyzz; clear in3dn;
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];

res3d = 10;
Xdim = dims(4);
Ydim = dims(4);
Zdim = dims(2);

X3d = linspace(0,Xdim,res3d);
Y3d = linspace(0,Ydim,res3d);
Z3d = linspace(0,Zdim,res3d);


for zz = 1:(res3d-1)
	xxyyzz = 1;
	zL = Z3d(zz); zH = Z3d(zz+1);
	for xx = 1:(res3d-1)
		xL = X3d(xx); xH = X3d(xx+1); 
		for yy = 1:(res3d-1)
			yL = Y3d(yy); yH = Y3d(yy+1); 

		in3db{xxyyzz,zz} = in3Dbox(ActinTips,xL,xH,yL,yH,zL,zH);
		xxyyzz = xxyyzz+1;
		end
	end
end

%==================================================%
%%

XYZsz = size(in3db);
XYsz = XYZsz(1);
Xsz = sqrt(XYsz);
Ysz = Xsz;
Zsz = XYZsz(2);

in3dn = zeros(Xsz,Ysz,Zsz);
for SZz = 1:Zsz
	for XZz = 1:Zsz
		for YZz = 1:Zsz
	in3dn(XZz,YZz,SZz) = sum(in3db{XZz+YZz-1,SZz});
		end
	end
end
sum(sum(in3dn))

[xm,ym,zm] = meshgrid(1:9,1:9,1:9);
% scatter3(xm(:),ym(:),zm(:),5,in3dn(:))

figure
set(gcf,'OuterPosition',[300 200 400 700])
scatter3(ActinTips(:,1),ActinTips(:,2),ActinTips(:,3),50,'filled')
xlim([-250 250]);ylim([-250 250]);zlim([0 1100])
view([-27.5 6]);
hold on

%%
clear cmx szmx xmx ymx zmx
for xyzm = 1:9
cmx = in3dn(:,:,xyzm);
szmx = (in3dn(:,:,xyzm)+1) .^6;
xmx = xm(:,:,xyzm);
ymx = ym(:,:,xyzm);
zmx = zm(:,:,xyzm);
scatter3(xmx(:),ymx(:),zmx(:),szmx(:),cmx(:),'filled')
hold on;
end


%%
figure
set(gcf,'OuterPosition',[300 200 400 700])
scatter3(ActinTips(:,1),ActinTips(:,2),ActinTips(:,3),50,'filled')
xlim([-250 250]);ylim([-250 250]);zlim([0 1100])
view([-27.5 6]);
% axis vis3d

% in3dn = zeros(XYsz,Zsz)
% for szz = 1:Zsz
% 	for sxy = 1:XYsz
% 	in3dn(sxy,szz) = sum(in3db{sxy,szz});
% end
% end
% sum(sum(in3dn))

%}
%==================================================%



%%
% TRIANGULATION AND CONVEX HULL
%==================================================%
%{

%%
%==================================================%
Ori_xyz = [Actin(:,3) Actin(:,6) Actin(:,9)];
Tip_xyz = [Actin(:,4) Actin(:,7) Actin(:,10)];
ATPxyz = [Ori_xyz; Tip_xyz];

% radial distance to spine shaft membrane
XYtipLoc = sqrt(ATPxyz(:,1).^2 + ATPxyz(:,2).^2);

ZtipInHead = (ATPxyz(:,3) <= SPYhZN) & (ATPxyz(:,3) >= SPYhZS);
ZtipInNeck = ATPxyz(:,3) < SPYhZS;

XYneckOut = XYtipLoc > SPYnXY;	% Logical array of filaments beyond SPYnXY
XYheadOut = XYtipLoc > SPYhXY;  % Logical array of filaments beyond SPYhXY
ZtopOut = ATPxyz(:,3) > SPYhZN;	% Logical array of filaments above SPYhZN
ZbotOut = ATPxyz(:,3) < 0;		% Logical array of filaments below zero

TipOut = ((XYneckOut & ZtipInNeck) + (XYheadOut & ZtipInHead) + ZtopOut + ZbotOut)>0;

ATPxyz(TipOut,:) = [];

ATPxyz(end:end+4,:) = 0;

szATP = size(ATPxyz,1);

ATPxyz(szATP+1,:) = [SPYnXY/2 SPYnXY/2 0];
ATPxyz(szATP+2,:) = [-SPYnXY/2 SPYnXY/2 0];
ATPxyz(szATP+3,:) = [SPYnXY/2 -SPYnXY/2 0];
ATPxyz(szATP+4,:) = [-SPYnXY/2 -SPYnXY/2 0];

%==================================================%
fig1 = figure;
set(fig1,'OuterPosition',[500 100 800 900])

AOTxyz = ATPxyz;
minAOT = -min(AOTxyz); 
AOTxyz(:,1) = AOTxyz(:,1) + minAOT(1);
AOTxyz(:,2) = AOTxyz(:,2) + minAOT(2);
AOTxyz(:,3) = AOTxyz(:,3);

AOTxyzU = unique(AOTxyz,'rows');

Uxyz = AOTxyzU;

TinH = Uxyz(:,3) >= (SPYhZS);
HTxyz = Uxyz(TinH,:);

TinN = Uxyz(:,3) < (SPYhZS);
NTxyz = Uxyz(TinN,:);

TinNH = (Uxyz(:,3) >= (SPYhZS-80)) & (Uxyz(:,3) <= (SPYhZS+50));
NHTxyz = Uxyz(TinNH,:);

DTriTip = delaunayTriangulation(HTxyz);
[FBtri,FBpoints] = freeBoundary(DTriTip);

trimesh(FBtri,FBpoints(:,1),FBpoints(:,2),FBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
hold on

NTriTip = delaunayTriangulation(NTxyz);
[NFBtri,NFBpoints] = freeBoundary(NTriTip);

trimesh(NFBtri,NFBpoints(:,1),NFBpoints(:,2),NFBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
hold on

NHTriTip = delaunayTriangulation(NHTxyz);
[NHFBtri,NHFBpoints] = freeBoundary(NHTriTip);

trimesh(NHFBtri,NHFBpoints(:,1),NHFBpoints(:,2),NHFBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.3);
hold on


scatter3(AOTxyz(:,1),AOTxyz(:,2),AOTxyz(:,3),20,'fill','r')
axis vis3d
view([-48 8]);
%==================================================%
%%



















%==================================================%
fig1 = figure;
set(fig1,'OuterPosition',[500 100 800 900])

Ori_xyz = [Actin(:,3) Actin(:,6) Actin(:,9)];
Tip_xyz = [Actin(:,4) Actin(:,7) Actin(:,10)];
AOTxyz = [Ori_xyz; Tip_xyz];
minAOT = -min(AOTxyz); 
AOTxyz(:,1) = AOTxyz(:,1) + minAOT(1);
AOTxyz(:,2) = AOTxyz(:,2) + minAOT(2);
AOTxyz(:,3) = AOTxyz(:,3) + minAOT(3);

AOTxyzU = unique(AOTxyz,'rows');
DTriTip = delaunayTriangulation(AOTxyzU);
[FBtri,FBpoints] = freeBoundary(DTriTip);

trisurf(FBtri,FBpoints(:,1),FBpoints(:,2),FBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
hold on
scatter3(AOTxyz(:,1),AOTxyz(:,2),AOTxyz(:,3),20,'fill','r')
axis vis3d
view([-36.5 6]);
%==================================================%

%==================================================%
Ori_xyz = [Actin(:,3) Actin(:,6) Actin(:,9)];
Tip_xyz = [Actin(:,4) Actin(:,7) Actin(:,10)];
AOTxyz = [Ori_xyz; Tip_xyz];
minAOT = -min(AOTxyz); 
AOTxyz(:,1) = AOTxyz(:,1) + minAOT(1);
AOTxyz(:,2) = AOTxyz(:,2) + minAOT(2);
AOTxyz(:,3) = AOTxyz(:,3) + minAOT(3);

DTriTip = delaunayTriangulation(AOTxyz);
[DThull DThullV] = convexHull(DTriTip);

trisurf(DThull,DTriTip.Points(:,1),DTriTip.Points(:,2),DTriTip.Points(:,3),...
       'FaceColor','cyan')

trimesh(DThull,DTriTip.Points(:,1),DTriTip.Points(:,2),DTriTip.Points(:,3))
axis vis3d

%==================================================%



%%
%==================================================%
fig1 = figure;
set(fig1,'OuterPosition',[500 100 800 900])

Ori_xyz = [Actin(:,3) Actin(:,6) Actin(:,9)];
Tip_xyz = [Actin(:,4) Actin(:,7) Actin(:,10)];
AOTxyz = [Ori_xyz; Tip_xyz];
minAOT = -min(AOTxyz); 
AOTxyz(:,1) = AOTxyz(:,1) + minAOT(1);
AOTxyz(:,2) = AOTxyz(:,2) + minAOT(2);
AOTxyz(:,3) = AOTxyz(:,3) + minAOT(3);

AOTxyzU = unique(AOTxyz,'rows');

Uxyz = AOTxyzU;


TinH = Uxyz(:,3) >= (SPYhZS);
HTxyz = Uxyz(TinH,:);

TinN = Uxyz(:,3) < (SPYhZS-10);
NTxyz = Uxyz(TinN,:);

TinNH = (Uxyz(:,3) >= (SPYhZS-50)) & (Uxyz(:,3) <= (SPYhZS+30));
NHTxyz = Uxyz(TinNH,:);

DTriTip = delaunayTriangulation(HTxyz);
[FBtri,FBpoints] = freeBoundary(DTriTip);

trisurf(FBtri,FBpoints(:,1),FBpoints(:,2),FBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
hold on

NTriTip = delaunayTriangulation(NTxyz);
[NFBtri,NFBpoints] = freeBoundary(NTriTip);

trisurf(NFBtri,NFBpoints(:,1),NFBpoints(:,2),NFBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
hold on

NHTriTip = delaunayTriangulation(NHTxyz);
[NHFBtri,NHFBpoints] = freeBoundary(NHTriTip);

trisurf(NHFBtri,NHFBpoints(:,1),NHFBpoints(:,2),NHFBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
hold on

scatter3(AOTxyz(:,1),AOTxyz(:,2),AOTxyz(:,3),20,'fill','r')
axis vis3d
view([-36.5 6]);
%==================================================%
%%
%}
%==================================================%


%%
%==================================================%
%		Counter Crunch and Plot
%--------------------------------------------------%
if doActCounts


doNdel    = AMS{8};
nDelSteps = AMS{9};

if doNdel
        % Ndel = SaveTipsAfter; % Deletes the first 'Ndel' steps
        Ndel = nDelSteps;       % Deletes the first 'Ndel' steps
        nT0=nT; Ns0=Nsteps;     % Save the original Nt value
        nT = nT-Ndel;           % Reduce Nt by number of deleted steps
        Nsteps = nT;

		NumFils(1:Ndel) = [];
		NumFAct(1:Ndel) = [];
		NumGAct(1:Ndel) = [];
		NumCOFACT(1:Ndel) = [];
		NumARPACT(1:Ndel) = [];
		NumdelFi(1:Ndel) = [];
		
		ArpBRsum(1:Ndel) = [];
		NumAct_uM(1:Ndel) = [];
		NumfKa(1:Ndel) = [];
		NumFnowFtot(1:Ndel) = [];
		NumFArp(1:Ndel) = [];
		NumGArp(1:Ndel) = [];
		Num_nmmf(1:Ndel) = [];
		NumFded(1:Ndel) = [];

        NumGActinN0(1:Ndel) = [];
		NumTHYM(1:Ndel) = [];
		NumTHYMACT(1:Ndel) = [];
end

	%==================================================%
	%				REMINDER OF COUNTERS
	%--------------------------------------------------%
	%{
		NumFils(nT) = Nfi;						% Number of branch filaments
		NumFAct(nT) = FActinN;					% Number FActins
		NumGAct(nT) = GActinN;					% Number GActins
		NumCOFACT(nT) = COFdepoly+ACTdepoly;	% Number of Depoly events
		NumARPACT(nT) = ACTpoly+ARPpoly;		% Number of Poly events
		NumdelFi(nT) = sum(delFi);
		
		ArpBRsum(nT) = sum(ArpBR);
		NumAct_uM(nT) = Act_uM;
		NumfKa(nT) = BeKa;
		NumFnowFtot(nT) = Nfi / Actin(end,12);
		NumFArp(nT) = FArpN;
		NumGArp(nT) = GArpN;
		Num_nmmf(nT) = sum(nmmf);
		NumFded(nT) = numel(DiedACTs(:,1));

		NumFils
		NumFAct
		NumGAct
		NumCOFACT
		NumARPACT
		NumdelFi
		
		ArpBRsum
		NumAct_uM
		NumfKa
		NumFnowFtot
		NumFArp
		NumGArp
		Num_nmmf
		NumFded

%-----------------------------------------------------------------------------------%
								ACTs(:,1:19)
% N  Xa  Xo  Xt  Ya  Yo  Yt  Za  Zo  Zt  MomID  ID  Fkd Born Died Lif MaxL MeanL Lgth
% 1  2   3   4   5   6   7   8   9   10  11     12  13  14   15   16  17   18    19
%-----------------------------------------------------------------------------------%
	%}
	%--------------------------------------------------%





    TotN = numel(NumGAct);
    HalfN = round(TotN/2);
	
	Ga_count_mean = mean(NumGAct(end-HalfN:end));
	Fa_count_mean = mean(NumFAct(end-HalfN:end));
	Ga_uM_mean = mean(NumAct_uM(end-HalfN:end));
	Fa_uM_mean = VCP(.1,0,Fa_count_mean);
	
	Total_nm_fil_mean = mean(Num_nmmf(end-HalfN:end));
	Fil_count_mean = mean(NumFils(end-HalfN:end));
	Fil_length_um_mean = Total_nm_fil_mean / Fil_count_mean / 1000;
	
	ProtoNs = Actin(:,1); %ACTs(:,1);
	mean(ProtoNs)
	
	% ActSimChem Filament List
	%{
	ActSimC = [58 45 32 6 24 34 31 99 21 18 59 28 115 87 39 68 57 36 41 86 40 62 58 31 28 39 ...
	47 35 43 18 50 60 40 45 34 16 13 8 64 3 24 37 51 26 20 16 38 31 72 23 70 45 ...
	57 30 44 12 8 20 37 94 23 11 43 44 26 57 28 90 8 58 40 30 33 30 24 64 60 79 ...
	59 64 37 59 40 37 56 27 29 51 42 49 55 54 15 43 39 74 55 49 89 77 38 16 39 ...
	24 49 25 60 11 42 65 46 71 84 60 31 59 132 31 151 60 34 74 64 34 18 100 40 ...
	17 32 41 25 29 26 55 70 35 16 21 39 21 38 55 64 73 42 75 73 69 45 12 45 48 ...
	71 66 78 100 54 53 43 60 39 34 47 64 90 60 45 57 68 62 22 27 26 52 19 67 65 ...
	33 35 70 4 19 99 30 26 62 12 61 6 26 82 99 81 89 42 7 44 31 25 57 20 33 25 ...
	86 77 60 42 27 48 26 69 56 16 81 26 15 59 41 27 60 111 42 47 46 36 71 53 54 ...
	78 39 29 103 54 69 59 36 70 68 42 37 16 31 55 42 9 40 60 25 50 58 51 12 38 ...
	28 28 67 12 41 36 14 48 36 45 77 56 14 71 58 59 74 59 41 12 104 99 67 7 109 ...
	78 23 62 98 62 51 44 29 15 61 63 99 35 26 121 47 68 54 32 56 46 49 25 30 32 ...
	75 13 65 45 53 45 105 43 36 17 82 12 37 33 20 49 36 45 45 73 53 25 27 41 34 ...
	28 76 26 67 26 78 35 84 54 67 34 80 13 105 52 41 32 41 39 29 91 54 26 9 61 ...
	57 18 47 36 74 53 15 53 18 44 49 32 40 44 58 42 50 28 88 76 35 76 65 39 52 31 63];

	mean(ActSimC)
	
	
	ActSC = xlsread('ActSCdat.xlsx');
	[ActSCnum,ActSCtxt,ActSCraw] = xlsread('ActSCdat.xlsx');
	plot(ActSC(:,3))
	title(ActSCtxt(1,3))
	figure(1)
	plot(ActSC(:,24))
	title(ActSCtxt(1,24))
	%}
	
	


%%
keyboard
%%
clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					FINAL OUTPUT (SPLINE) FIGURE 1 OF 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scsz = get(0,'ScreenSize'); pos1=[scsz(3)/10 scsz(4)/10 scsz(3)/1.2 scsz(4)/1.2];
	
    fig88 = figure(88);
	set(fig88,'Renderer','zbuffer','Units','pixels','OuterPosition',pos1,'Color',[1 1 1]);
	
	%-------------------------------------------------------------
	c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
	c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
	pause(.1); MS1 = 5;
	




%===========================================================%
% FIG1 TOP LEFT:
%===========================================================%
	Filament_Lifetime = ACTs(:,16) .* dT;
	muFlif = mean(Filament_Lifetime);
	medFlif = median(Filament_Lifetime);
	Flast = max(Filament_Lifetime);
	[coECDF_Lif, coX_Lif] = ecdf(Filament_Lifetime);
	XLmax = 3.6e3;
	%------
	Filament_Lifetimes = sort(Filament_Lifetime);
	[CdfY,CdfX,CdfLower,CdfUpper] = ecdf(Filament_Lifetimes,...
	'Function','survivor','alpha',.0001);  % compute empirical function
	[StairsXlower,StairsYlower] = stairs(CdfX,CdfLower);
	[StairsXupper,StairsYupper] = stairs(CdfX,CdfUpper);





                Hax1 = axes('Position',[.07 .57 .4 .38]);

PhTL1 = stairs(CdfX,1-CdfY);
                grid off; hold on;

hBounds = plot([StairsXlower(:); NaN; StairsXupper(:)],...
		       [1-StairsYlower(:); NaN; 1-StairsYupper(:)]);
                hold on
                set(gca,'XLim', [0 XLmax],'YLim', [0 1.1])



	%-------------------- Tic Labels ---------------------
	if verLessThan('matlab', '8.3.1');
        ax = gca; axT = axis;
        xt = get(ax,'XTick')./60; % no need to xt*dT since data is already
        set(gca,'XTickLabel', sprintf('%.1f|',xt))
	else
        ax = gca; axT = axis;
		xt = round(ax.XTick./60);
		ax.XTickLabel = xt;
	end
    %-----------------------------------------------------
    
	text((axT(2)/1.3),.08,['Mean: ' num2str(muFlif/60,3) ' min'],...
		'FontSize',12,'BackgroundColor',[1 1 1]);
	text((axT(2)/1.3),.14,['Median: ' num2str(medFlif/60,3) ' min'],...
		'FontSize',12,'BackgroundColor',[1 1 1]);
	%set(PhTL1,'color',[1 0 1])
	rectangle('Position',[Flast,1,nT,.001],'EdgeColor',[1 0 1],'LineStyle','--')
	CDFtitle = title('Empirical CDF of Filament Lifetime');
	set(CDFtitle,'FontSize',12);
	set(PhTL1,'LineStyle','-','Color', [.3 0 .6],'LineWidth',5);
	set(hBounds,'Color',[.9 .2 .2], 'LineWidth',1.5,'LineStyle','--');
	xlabel('Seconds Existed'); ylabel('CDF (+/- 99.99% CI)');
	%----------------------

	
	hold on
	Hax1_pos = get(gca,'Position'); % store position of first axes
	Hax2 = axes('Position',Hax1_pos,'Color','none',...
				  'XAxisLocation','top','YAxisLocation','right',...
				  'XTick', [],'YTick',[]);
			  
	set(Hax2,'Box','off','TickDir','out','TickLength',[.02 .02], ...
	'XMinorTick','off','YMinorTick','off','YGrid','off');
	%----------------------

	[BarY,BarX] = ecdfhist(coECDF_Lif, coX_Lif, 20);
	BarY = BarY .* 500;
	
PhTL2 = bar(Hax2,BarX, BarY,'Parent',Hax1,'FaceColor',c3);




%===========================================================%
% FIG1 TOP RIGHT:
%===========================================================%
sbpos = [.57 .57 .38 .38];
HaxTR = axes('Position',sbpos);
%----------------------
	Filament_Lifetimes = sort(Filament_Lifetime);
	%dfittool(Filament_Lifetimes)
	LegHandles = []; LegText = {};
	%axes('Position',[.1 .15 .8 .8]);
	pause(.5);
	
[CdfY,CdfX,CdfLower,CdfUpper] = ecdf(Filament_Lifetimes,...
	'Function','survivor','alpha',1e-10);  % compute empirical function
PhTR1 = stairs(CdfX,CdfY,'Color',[.3 0 .6],'LineStyle','-', 'LineWidth',1.5);
	hold on
[StairsXlower,StairsYlower] = stairs(CdfX,CdfLower);
	hold on
[StairsXupper,StairsYupper] = stairs(CdfX,CdfUpper);
	hold on
hBounds = plot([StairsXlower(:); NaN; StairsXupper(:)], [StairsYlower(:); NaN; StairsYupper(:)],...
	'Color',[.9 .2 .2],'LineStyle',':', 'LineWidth',1);
	hold on

	xlabel('CDF Survivor (seconds)'); ylabel('Survivor function');
	LegHandles(end+1) = PhTR1; LegText{end+1} = 'Filament Lifetime';
	LegHandles(end+1) = hBounds; LegText{end+1} = '99.99% confidence bounds';
	set(gca,'XLim',[-100 3600]); box on; grid off;
	hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'Location', 'NorthEast');
	set(hLegend,'Interpreter','none');
	set(PhTR1,'LineStyle','-','Color', [.3 0 .6],'LineWidth',5);
	set(hBounds,'Color',[.9 .2 .2], 'LineWidth',1.5,'LineStyle','--');

	


%===========================================================%
% FIG1 BOT LEFT:
%===========================================================%
	% Fil_ID / (Time_of_Birth * dT)
	% (e.g. if Fid=50 and ToB=100 the branch rate .5 F/sec)
	AcTags = Actin(:,12) ./ (Actin(:,14) .* dT);
	AcToB = (Actin(:,14) .* dT);
	AcTags(1:4,:) = AcTags(1:4,:)./3; %AcToB(1:15,:) = [];
	AT2TNratio = numel(AcTags)/TagN;
%----------------------
	ummf = Num_nmmf ./ 1000;
 	CompressN = 120; 
    NcT = nT / CompressN;
 	subs = floor(linspace(1,CompressN+1,nT));
 	subs(end) = CompressN;
 	Mummf = accumarray(subs',ummf,[],@mean);
%----------------------
	sbpos = [.07 .09 .4 .35];
	HaxBL = axes('Position',sbpos);
%----------------------


[PhBL] = plot(AcToB,AcTags,'b');

	%----------------------
	xt = (get(gca,'XTick'));
	xt = linspace(0,Nsteps,numel(xt)).* dT./60;
	set(gca,'XTickLabel', sprintf('%.1f|',xt))
	%-------
	hTitle  = title('Filament Branching');
	hXLabel = xlabel('Time (min)');
	hYLabel = ylabel('Fid / (ToB * dT)');
	hold on
	HaxBL_pos = get(HaxBL,'Position'); % store position of first axes
	HaxBL2 = axes('Position',HaxBL_pos,'Color','none',...
				  'XAxisLocation','top','YAxisLocation','right');
	set(HaxBL2,'Box','off','TickDir','out','TickLength',[.02 .02], ...
	'XMinorTick','off','YMinorTick','off','YGrid','off', ...
	'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
	%----------------------

[PhBL2] = line(1:numel(ArpBRsum),ArpBRsum,'Parent',HaxBL2,'Color','k');

	%-------
	hYLabel = ylabel('Arp Branch Rate (5*mM*ummf)');
	%-------
	set(PhBL,'LineStyle','-','Color',c1,'LineWidth',1,...
	'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
	set([hXLabel, hYLabel],'FontSize',10);
	set(hTitle,'FontSize',12);
	set(HaxBL,'Box','off','TickDir','out','TickLength',[.01 .01], ...
	'XMinorTick','on','YMinorTick','on','YGrid','on','FontName','Helvetica', ...
	'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
	%-------
	set(HaxBL2,'Box','off','TickDir','out','TickLength',[.01 .01], ...
	'XMinorTick','off','YMinorTick','off','YGrid','off',...
	'XTick', []);



	HaxBL3 = axes('Position',[.07 .09 .4 .35],...
	'YAxisLocation','right','Color','none');
	hold on
	
[PhBL3] = plot(Mummf,'Parent',HaxBL3,'Color','g');

	set(HaxBL3,'XTickLabel', [],'YTickLabel', [],'XTick', [],'YTick', [])



    %---------------------- Legend -----------------------
	leg1 = legend([PhBL,PhBL2,PhBL3],{'Fid:ToB','ArpBR','ummf'});
	set(leg1, 'Location','SouthEast', 'Color', [1 1 1],'FontSize',14,'Box','off');
    set(leg1, 'Position', leg1.Position .* [1 .98 1 1.4])

	%-------------------- Tic Labels ---------------------
	if verLessThan('matlab', '8.3.1');
		xt = (get(gca,'XTick'));
        xt = linspace(0,Nsteps,numel(xt)).* dT./60;
        set(gca,'XTickLabel', sprintf('%.1f|',xt))
	else	
		xt = HaxBL.XTick;
        xt = linspace(0,Nsteps,numel(xt)).* dT./60;
		HaxBL.XTickLabel = xt;
	end
    %-----------------------------------------------------



%===========================================================%
% FIG1 BOT RIGHT:
%===========================================================%
sbpos = [.57 .09 .38 .35];
Hax1 = axes('Position',sbpos);
%----------------------

CompressN = 70; NcT = nT / CompressN;
subs = floor(linspace(1,CompressN+1,nT));
subs(end) = CompressN;
NumaBeKa = accumarray(subs',NumfKa,[],@mean);
%NumaFils = accumarray(subs',NumFils,[],@mean);
%----------------------

	%sbpos = [.07 .09 .4 .35];
	ax2 = subplot('Position',sbpos);
[ph1] = plot(NumaBeKa,'r');
	%hold on
	%[ph1] = plot([NumaBeKa NumaFils],'b','Parent',ax2);


	%-------------------- Tic Labels ---------------------
	if verLessThan('matlab', '8.3.1');
		xt = (get(gca,'XTick')).* dT./60.*NcT;
        set(gca,'XTickLabel', sprintf('%.0f|',xt))
        set(legend,'Location','NorthEast');
	else	
		xt = ax2.XTick;
        xt = linspace(0,Nsteps,numel(xt)).* dT./60;
        xt = roundn(xt,-1);
		ax2.XTickLabel = xt;
	end
    %-----------------------------------------------------

	
	%------------------------------------------%
	haxes2a=axis;
	set(ax2,'YLim',[-15 haxes2a(4)]);
	%------------------------------------------%
	MS1 = 7; MS2 = 2;
	set(ph1,'LineStyle','-','Color',c2,'LineWidth',1,...
	'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
	hTitle  = title('Actin Polymerization ');
	hXLabel = xlabel('Time (min) ');
	hYLabel = ylabel('Polymerization Rate (10*uM*dT) ');
	set(gca,'FontName','Helvetica');
	set([hYLabel],'FontName','Century Gothic');
	set([hTitle, hXLabel],'FontName','AvantGarde');
	set([hXLabel],'FontSize',10); set([hYLabel],'FontSize',11);
	set(hTitle,'FontSize',12);
	set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
	'XMinorTick','off','YMinorTick','on','YGrid','off', ...
	'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
	

	sbpos = [.61 .12 .34 .15];
	NuBeKa = NumfKa(ceil(numel(NumfKa)/10):end);
	haxes6 = axes('Position',sbpos,'YAxisLocation','right','Color','none');
	%hold on
% [ph2b] = plot(NuBeKa,'Parent',haxes6,'Color','r');
[ph2b] = semilogy(NuBeKa);
	set(gca,'Position',sbpos,'YAxisLocation','right','Color','none');
	set(haxes6,'XTickLabel', [])
	hYLabel = ylabel('Poly Rate Zoom ');
	set([hYLabel],'FontName','Century Gothic');
	set([hYLabel],'FontSize',11);
	
	axis tight
	%haxes2b=axis;
	%xlim([0 (haxes2b(2)*.25)]); ylim([0 haxes2b(4)]);
	MS1 = 3;
	set(ph2b,'LineStyle','-','Color',c2,'LineWidth',1,...
	'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					FINAL OUTPUT SECOND SET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	fig21 = figure(21);
	set(fig21,'Renderer','zbuffer')
	%set(fig21,'Renderer','OpenGL')
	set(21,'Units','pixels');scsz = get(0,'ScreenSize');
	pos1 = [scsz(3)/10  scsz(4)/10  scsz(3)/1.2  scsz(4)/1.2];
	set(fig21,'OuterPosition',pos1)
	set(gcf,'Color',[1,1,1])
	%-------------------------------------------------------------
	c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
	c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
	pause(.1);
	

%===========================================================%
% FIG1 TOP LEFT: Poly & Depoly Events
%===========================================================%
%{
CompressN = 70; NcT = nT / CompressN;
subs = floor(linspace(1,CompressN+1,nT));
subs(end) = CompressN;
AcArCOFACT = accumarray(subs',NumCOFACT,[],@mean);
AcArARPACT = accumarray(subs',NumARPACT,[],@mean);
%----------------------


	sbpos = [.07 .57 .4 .38];
	
subplot('Position',sbpos);
[ph1] = plot(AcArARPACT,'b');

	leg1 = legend(ph1,'On_{poly}');
	[LEGH,OBJH,OUTH,OUTM] = legend;
	hold on
	
subplot('Position',sbpos);
[ph2] = plot(AcArCOFACT,'r');

	legend([OUTH;ph2],OUTM{:},'Off_{depoly}');
	[LEGH,OBJH,OUTH,OUTM] = legend;
	hold on
	%------------------------------------------%
	xt = (get(gca,'XTick'));
	xt = linspace(0,Nsteps,numel(xt)).* dT./60;
	set(gca,'XTickLabel', sprintf('%.0f|',xt))
	set(legend,'Location','SouthEast');
	%------------------------------------------%
	MS1 = 7;
	set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
	'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor','none');
	set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
	'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
	hTitle  = title('Actin Polymerization & Depolymerization Events');
	hXLabel = xlabel('Time (min)');
	hYLabel = ylabel('Poly & Depoly Events');
	set(gca,'FontName','Helvetica');
	set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
	set([hXLabel, hYLabel],'FontSize',10);
	set(hTitle,'FontSize',12);
	set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
	'XMinorTick','on','YMinorTick','on','YGrid','on', ...
	'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
	haxes=axis;
	% ylim([0 haxes(4)*1.2 ]);
	% xlim([0 (haxes(2)*.9)]);
%}
%===========================================================%
%===========================================================%
% FIG1 TOP LEFT: Spine Volume (from Delaunay Triangulation Convex Hull
%===========================================================%

HulVolN = NumHullV(1,:) .* 1e-9;
HulVolH = NumHullV(2,:) .* 1e-9;
NHVol = numel(HulVolN);

if NHVol >= 60;CompressN = 30; else CompressN = 30; end;
NcT = NHVol / CompressN;
subs = floor(linspace(1,CompressN+1,NHVol));
subs(end) = CompressN;
HulVolNk = accumarray(subs',HulVolN,[],@mean);
HulVolHd = accumarray(subs',HulVolH,[],@mean);

HulVolNkV = accumarray(subs',HulVolN,[],@std);
HulVolNkSEM = HulVolNkV ./ sqrt(numel(subs(subs==1)));

HulVolHdV = accumarray(subs',HulVolH,[],@std);
HulVolHdSEM = HulVolHdV ./ sqrt(numel(subs(subs==1)));
%----------------------
	sbpos = [.06 .57 .4 .38];
	
    haxSP1 = subplot('Position',sbpos);
[ph1] = plot(HulVolNk,'b');

 	hold on
		
    haxSP2 = subplot('Position',sbpos);
[ph2] = plot(HulVolHd,'r');

 	hold on


    %---------------------- Legend -----------------------
 	leg1 = legend([ph1,ph2],{'Neck Volume','Head Volume'});
 	set(leg1, 'Location','SouthEast', 'Color', [1 1 1],'FontSize',14,'Box','on');
    set(leg1, 'Position', leg1.Position .* [1 .98 1 1.4])

	%-------------------- Tic Labels ---------------------
	if verLessThan('matlab', '8.3.1');
		xt = (get(gca,'XTick'));
        xt = linspace(0,Nsteps,numel(xt)).* dT./60;
        set(gca,'XTickLabel', sprintf('%.0f|',xt))
        set(legend,'Location','SouthEast');
	else	
		xt = haxSP2.XTick;
        xt = linspace(0,Nsteps,numel(xt)).* dT./60;
        xt = roundn(xt,-1);
		haxSP2.XTickLabel = xt;
	end
    %-----------------------------------------------------

	
	
	%------------------------------------------%
	MS1 = 7;
	set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
	'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
	set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
	'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);

	hTitle  = title('Spine Volume ');
	hXLabel = xlabel('Time (min) ');
	hYLabel = ylabel('Volume (um^3) of Convex Hull  ');
% 	hYLabel = ylabel('Volume (\mum^3) _{(Delaunay Convex Hull)}');
	set(gca,'FontName','Helvetica');
	set([hYLabel],'FontName','Century Gothic');
	set([hTitle, hXLabel],'FontName','AvantGarde');
	set([hXLabel],'FontSize',10); set([hYLabel],'FontSize',11);
	set(hTitle,'FontSize',12);
	set(gca,'Box','off','TickDir','out','TickLength',[.01 .01], ...
	'XMinorTick','off','YMinorTick','off','YGrid','off', ...
	'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
	haxes=axis;
	ylim([haxes(3) haxes(4)*1.1 ]);
	%xlim([0 (haxes(2)*.875)]);

	
	%------------------------------------------%
	%axes('Position',get(gca,'Position'),'Color','none','XTick',[],'YTick',[],...
	%			  'XAxisLocation','top','YAxisLocation','right');
	%------------------------------------------%
	
	%----------------------
	ummf = Num_nmmf ./ 1000;
 	CompressN = 120; NcT = nT / CompressN;
 	subs = floor(linspace(1,CompressN+1,nT));
 	subs(end) = CompressN;
 	Mummf = accumarray(subs',ummf,[],@mean);
	%----------------------

	HaxBL3 = axes('Position',get(gca,'Position'),...
	'YAxisLocation','right','Color','none');
	hold on
	
[PhBL3] = plot(Mummf,'Parent',HaxBL3,'Color',[.5 .5 .5]);

    hold on

	
    %---------------------- Legend -----------------------
 	leg1 = legend([ph1,ph2],{'Neck Volume','Head Volume'});
 	set(leg1, 'Location','SouthEast', 'Color', [1 1 1],'FontSize',14,'Box','off');
    set(leg1, 'Position', leg1.Position .* [1 .98 1 1.4])

	%-------------------- Tic Labels ---------------------
	if verLessThan('matlab', '8.3.1');
		LGh1 = legend([OUTH(:); PhBL3],OUTM{:},'ummf ');
        set(LGh1,'Location','SouthEast');
        hold on;
	else	
		xt = HaxBL3.XTick;
        xt = linspace(0,Nsteps,numel(xt)).* dT./60;
        xt = roundn(xt,-1);
		% HaxBL3.XTickLabel = xt;
	end
    %-----------------------------------------------------


    axis tight
	haxC=axis;
	ylim([haxC(4)/5 haxC(4)*1.1 ]);
	%xlim([0 (haxC(2)*1)]);
	%set(HaxBL3,'XTickLabel', [],'YTickLabel', [],'XTick', [],'YTick', [])
	set(HaxBL3,'Box','off','TickDir','out','TickLength',[.01 .01], ...
		'XMinorTick','off','YMinorTick','off','YGrid','off','XTick', []);
	hYLabel = ylabel('Combined Filament Length (µm) ');
	set([hYLabel],'FontName','Century Gothic','FontSize',11);
	set(gca,'Box','off','TickDir','out','TickLength',[.01 .01], ...
	'XMinorTick','off','YMinorTick','off','YGrid','off', ...
	'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);


	%------------------------------------------%
	% Error Bars
	%------------------------------------------%
	%{
[ph1e] = errorbar(HulVolNk,HulVolNkV,'rx');	
[ph2e] = errorbar(HulVolHd,HulVolHdV,'bx');

	MS1 = 7;
	set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
	'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
	set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
	'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);

	set(ph1e,'LineStyle','-','Color',c1,'LineWidth',1,...
	'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
	set(ph2e,'LineStyle','-','Color',c2,'LineWidth',1,...
	'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);

	hTitle  = title('Spine Volume ');
	hXLabel = xlabel('Time (min) ');
	hYLabel = ylabel('Volume (um^3) of Convex Hull  ');
% 	hYLabel = ylabel('Volume (\mum^3) _{(Delaunay Convex Hull)}');
	set(gca,'FontName','Helvetica');
	set([hYLabel],'FontName','Century Gothic');
	set([hTitle, hXLabel],'FontName','AvantGarde');
	set([hXLabel],'FontSize',10); set([hYLabel],'FontSize',11);
	set(hTitle,'FontSize',12);
	set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
	'XMinorTick','off','YMinorTick','on','YGrid','off', ...
	'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
	haxes=axis;
	ylim([haxes(3) haxes(4)*1.1 ]);
	xlim([0 (haxes(2)*.875)]);
%}
%===========================================================%



%===========================================================%
% FIG1 TOP RIGHT: Monomeric vs Filamentous Actin
%===========================================================%
CompressN = 200; NcT = nT / CompressN;
subs = floor(linspace(1,CompressN+1,nT));
subs(end) = CompressN;
NumGA = accumarray(subs',NumGAct,[],@mean);
NumFA = accumarray(subs',NumFAct,[],@mean);

GAct_uM = NumGA ./ SpyV .*(1./6e23).*1e6;  % GActin uM
FAct_uM = NumFA ./ SpyV .*(1./6e23).*1e6;  % FActin uM

GAct_uMa = NumGAct ./ SpyV .*(1./6e23).*1e6;  % GActin uM
FAct_uMa = NumFAct ./ SpyV .*(1./6e23).*1e6;  % FActin uM
%----------------------
	sbpos = [.55 .57 .38 .38];

	haxSP1 = subplot('Position',sbpos);
[ph1] = plot(GAct_uM,'b');

	%-------
	% leg1 = legend(ph1,'G-actin');
	% [LEGH,OBJH,OUTH,OUTM] = legend;
	hold on
	%-------
	
	haxSP2 = subplot('Position',sbpos);
[ph2] = plot(FAct_uM,'r');

	%-------
	% LGh5 = legend([OUTH;ph2],OUTM{:},'F-actin');
	% [LEGH,OBJH,OUTH,OUTM] = legend;
	% set(LGh5, 'Position', [.535 .93 .05 .04], 'Color', [1 1 1]);
	hold on;
	%-------


    %---------------------- Legend -----------------------
 	leg1 = legend([ph2,ph1],{'F-actin','G-actin'});
    % set(leg1, 'Position', leg1.Position .* [.535 .93 .05 .04])
 	set(leg1, 'Location','NorthWest', 'Color', [1 1 1],'FontSize',14,'Box','off');
    set(leg1, 'Position', leg1.Position .* [1.6 .82 1 1.4])
	%-------------------- Tic Labels ---------------------
	if verLessThan('matlab', '8.3.1');
		xt = (get(gca,'XTick')).* dT./60.*NcT;
        set(gca,'XTickLabel', sprintf('%.0f|',xt))
	else	
		xt = haxSP2.XTick;
        xt = linspace(0,Nsteps,numel(xt)).* dT./60;
        xt = roundn(xt,-1);
		haxSP2.XTickLabel = xt;
	end
    %-----------------------------------------------------


	%------------------------------------------%
	MS1 = 7;
 	set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
 	'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor','none');
 	set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
 	'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
	hTitle  = title('Monomeric vs Filamentous Actin ');
	hXLabel = xlabel('Time (min) ');
	hYLabel = ylabel('Actin Concentration (µM) ');
	set(gca,'FontName','Helvetica');
	set([hYLabel],'FontName','Century Gothic');
	set([hTitle, hXLabel],'FontName','AvantGarde');
	set([hXLabel],'FontSize',10); set([hYLabel],'FontSize',11);
	set(hTitle,'FontSize',12);
	set(gca,'Box','off','TickDir','out','TickLength',[.01 .01], ...
	'XMinorTick','off','YMinorTick','off','YGrid','off', ...
	'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
	haxesM=axis;
	ylim([-haxesM(4)/15 haxesM(4)*1.1]);


	
	
	%------------------------------------------%
	nFAct_uMa = FAct_uMa(ceil(numel(FAct_uMa)/10):end);
	mnFAct_uMa = nFAct_uMa - min(nFAct_uMa);
	nGAct_uMa = GAct_uMa(ceil(numel(GAct_uMa)/10):end);
	%-----------------------
	
	
	
	%[FAct Top Zoom]
	%------------------------------------------%
	sbtop = [.59 .80 .34 .15];
	haxes4 = axes('Position',sbtop,'XAxisLocation','top','YAxisLocation','right','Color','none');
	
[ph1b] = semilogy(mnFAct_uMa);

	set(gca,'Position',sbtop,...
		'XAxisLocation','top','YAxisLocation','right','Color','none');
	set(haxes4,'XTickLabel', [],'XTick', [])
	hYLabel = ylabel('µM F-actin ');
	set([hYLabel],'FontName','Century Gothic','FontSize',11);
	axis tight
	haxA=axis;
	%xlim([0 (haxA(2)*.5)]); 
	ylim([haxA(3) haxA(4)*1.08 ]);
	%yt = (get(gca,'YTick'))+min(nFAct_uMa);
	%set(gca,'YTickLabel', sprintf('%.1f|',yt))
	set(haxes4,'TickDir','out','TickLength',[.01 .01]...
		,'Box','off','XColor',[1 1 1],'YColor',[.3 .3 .3],'LineWidth',1)
	%------------------------------------------%

	
	
	%[GAct Bot Zoom]
	%------------------------------------------%
	sbbot = [.59 .572 .34 .15];
	haxes3 = axes('Position',sbbot,'YAxisLocation','right','Color','none');
	
[ph1b] = semilogy(nGAct_uMa,'r');

	set(gca,'Position',sbbot,'YAxisLocation','right','Color','none');
	set(haxes3,'XTickLabel', [])
	hYLabel = ylabel('µM G-Actin ');
	set([hYLabel],'FontName','Century Gothic','FontSize',11);
	axis tight
	haxB=axis;
	%xlim([0 (haxB(2))]); ylim([0 haxB(4)]);
	set(haxes3,'TickDir','out','TickLength',[.01 .01],'XTick',haxB(2)...
		,'Box','off','XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1)
	
	%------------------------------------------%
	

	%set(haxes3,'Box','off','TickDir','out','TickLength',[.01 .01], ...
	%'XMinorTick','off','YMinorTick','on','YGrid','on', ...
	%'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
	%[haxes,hline1,hline2] = plotyy(t,z1,t,z2,'semilogy','plot');
%===========================================================%





%===========================================================%
% FIG1 BOTTOM LEFT: Actin Concentration (mM)
%===========================================================%
sbpos = [.06 .09 .4 .35];
Hax1 = axes('Position',sbpos);
%----------------------
CompressN = 70; NcT = nT / CompressN;
subs = floor(linspace(1,CompressN+1,nT));
subs(end) = CompressN;
NumaFils = accumarray(subs',NumFils,[],@mean);
NumaFded = accumarray(subs',NumFded,[],@mean);
NumaFded = nthroot(NumaFded,1.2);
%-------

	Filament_Lifetimes = sort(Filament_Lifetime);
	[CdfY,CdfX,CdfLower,CdfUpper] = ecdf(Filament_Lifetimes,...
	'Function','survivor','alpha',1e-10);
	[StairsXlower,StairsYlower] = stairs(CdfX,CdfLower);
	[StairsXupper,StairsYupper] = stairs(CdfX,CdfUpper);

PhTR1 = stairs(CdfX,CdfY);
	hold on
hBounds = plot([StairsXlower(:); NaN; StairsXupper(:)],[StairsYlower(:); NaN; StairsYupper(:)]);
	hold on

	%----------------------
	LegHandles = []; LegText = {};
	LegHandles(end+1) = PhTR1; LegText{end+1} = 'Filament Lifetime';
	LegHandles(end+1) = hBounds; LegText{end+1} = '99.99% confidence bounds';
	hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'Location', 'NorthWest');
	set(hLegend,'Interpreter','none');
	%------

	title('Actin Filament Turnover','FontSize',12);
	xlabel('Time (min)'); ylabel('Filament Survivor eCDF');
	set(gca,'XLim',[-100 3600]); box on; grid off;
	set(PhTR1,'LineStyle','-','Color', [.3 0 .6],'LineWidth',5);
	set(hBounds,'Color',[.5 .5 .5], 'LineWidth',1.5,'LineStyle','--');

    
	
	%----------------------
    haxt = gca;
	%----------------------
	Hax1_pos = get(gca,'Position'); % store position of first axes
	Hax2 = axes('Position',Hax1_pos,'Color','none','XTick', [],...
				  'XAxisLocation','top','YAxisLocation','right');

	set(Hax2,'Box','off','TickDir','out','TickLength',[.01 .01], ...
	'XMinorTick','off','YMinorTick','off','YGrid','off','FontName','Helvetica');
	%----------------------

[PhBR2] = line(1:numel(NumaFded),NumaFded,'Parent',Hax2,'Color','k','LineWidth',5);
	hold on
[PhBR3] = line(1:numel(NumaFils),NumaFils,'Parent',Hax2,'Color','r','LineWidth',5);
	hold on

    ylabel('Filaments','FontName','Helvetica','Color',[.3 .3 .3],'FontSize',11);


    %---------------------- Legend -----------------------
 	leg1 = legend([PhBR2,PhBR3],{'Fil Removal Count','Fil Current Number'});
    % set(leg1, 'Position', leg1.Position .* [.535 .93 .05 .04])
 	set(leg1, 'Location','NorthEast', 'Color', [1 1 1],'FontSize',14,'Box','on');
    set(leg1, 'Position', leg1.Position .* [1 .98 1 1.4])
	%-------------------- Tic Labels ---------------------

    set(Hax2,'XTickLabel', [],'XTick', [])
	if verLessThan('matlab', '8.3.1');
		xt = (get(haxt,'XTick'))./60; % no need to xt*dT since data is already
        set(gca,'XTickLabel', sprintf('%.1f|',xt))
        set(Hax1,'Box','off','TickDir','out','TickLength',[.01 .01],...
		'FontName','Helvetica','XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
	else	
		xt = haxt.XTick;
        xt = linspace(0,Nsteps,numel(xt)).* dT./60;
        xt = roundn(xt,-1);
		haxt.XTickLabel = xt;
	end
    %-----------------------------------------------------

	%-------
	axh = axis;
	text(axh(2)/3.05,axh(4)/1.1,'Lifetimes','FontSize',12,'BackgroundColor',[1 1 1]);
	text(axh(2)/2.95,axh(4)/1.17,['Mean: ' num2str(muFlif/60,3) ' min'],...
		'FontSize',12,'BackgroundColor',[1 1 1],'Color',[.3 .3 .3]);
	text(axh(2)/2.95,axh(4)/1.24,['Median: ' num2str(medFlif/60,3) ' min'],...
		'FontSize',12,'BackgroundColor',[1 1 1],'Color',[.3 .3 .3]);

%===========================================================%





%===========================================================%
% FIG1 BOTTOM RIGHT: Branching Events (Dual Axis Plot)
%===========================================================%
CompressN = 70; NcT = nT / CompressN;
subs = floor(linspace(1,CompressN+1,nT));
subs(end) = CompressN;
NFArp = accumarray(subs',NumFArp,[],@mean);
NGArp = accumarray(subs',NumGArp,[],@mean);
%----------------------


	sbpos = [.55 .09 .38 .35];
	Hax1BR = axes('Position',sbpos);
[ph1] = plot(ArpBRsum,'r');
    hold on


	%----------------------
	hXLabel = xlabel('Time (min) ');
	hYLabela = ylabel('Arp Branch Rate (5*mM*ummf) ');
	set([hYLabela],'FontName','Century Gothic');
	set([hYLabela],'FontSize',11);
	%----------------------
	set(gca,'Box','off','TickDir','out','TickLength',[.01 .01], ...
	'XMinorTick','off','YMinorTick','off','YGrid','off', ...
	'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
	haxes1 = gca; % handle to axes
	haxes1_pos = get(haxes1,'Position'); % store position of first axes
	haxes2 = axes('Position',haxes1_pos,'Color','none',...
				  'XAxisLocation','top','YAxisLocation','right');
	set(gca,'Box','off','TickDir','out','TickLength',[.01 .01], ...
	'XMinorTick','off','YMinorTick','off','YGrid','off', ...
	'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
	%----------------------
	
[ph2] = line(1:numel(NFArp),NFArp,'Parent',haxes2,'Color','k');
	hold on
[ph3] = line(1:numel(NGArp),NGArp,'Parent',haxes2,'Color','k');
    hold on



    %---------------------- Legend -----------------------
 	leg1 = legend([ph1,ph2,ph3],{'Arp Branch Rate','Num FArp','Num GArp'});
 	set(leg1, 'Location','SouthEast', 'Color', [1 1 1],'FontSize',12,'Box','on');
    set(leg1, 'Position', leg1.Position .* [1 .98 1 1.4])
	%-------------------- Tic Labels ---------------------
    set(haxes2,'XTickLabel', [],'XTick', [])
	if verLessThan('matlab', '8.3.1');
		xt = (get(haxes1,'XTick')).* dT./60;
        set(haxes1,'XTickLabel', sprintf('%.0f|',xt))
	else	
		xt = haxes1.XTick;
        xt = linspace(0,Nsteps,numel(xt)).* dT./60;
        xt = roundn(xt,-1);
		haxes1.XTickLabel = xt;
	end
    %-----------------------------------------------------


        MS1 = 7; MS2 = 2;
    set(ph1,'LineStyle','-','Color',c1,'LineWidth',3,...
        'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c11);
    set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
        'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c2);
    set(ph3,'LineStyle','-','Color',c3,'LineWidth',1,...
        'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c3);

    hTitle  = title ('\fontsize{14} Arp2/3-Mediated Filament Branching');
    hXLabel = xlabel('\fontsize{11} Time (min)');
    hYLabel = ylabel('\fontsize{12} G-Arp & F-Arp Count');
    set([hTitle, hXLabel, hYLabel],'FontName','Helvetica Neue');
    set(gca,'Box','off','TickDir','out','TickLength',[.018 .018], ...
    'XMinorTick','off','YMinorTick','on','YGrid','on', ...
    'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%===========================================================%


%%
pause(3)
% clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					CONVEX HULL VOLUME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
	fig31 = figure(31);
	set(fig31,'Renderer','painters')
	%set(fig31,'Renderer','zbuffer')
	%set(fig21,'Renderer','OpenGL')
	set(31,'Units','pixels');scnsize = get(0,'ScreenSize');
	pos1 = [scnsize(3)/9  scnsize(4)/9  scnsize(3)/2  scnsize(4)/2];
	set(fig31,'OuterPosition',pos1)
	set(gcf,'Color',[.9,.9,.9])
	%-------------------------------------------------------------
	c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
	c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
	pause(1);
	

%===========================================================%
% FIG1 TOP LEFT: Poly & Depoly Events
%===========================================================%

HulVolN = NumHullV(1,:) .* 1e-9;
HulVolH = NumHullV(2,:) .* 1e-9;

CompressN = 20; NcT = numel(HulVolN) / CompressN;
subs = floor(linspace(1,CompressN+1,numel(HulVolN)));
subs(end) = CompressN;
HulVolNk = accumarray(subs',HulVolN,[],@mean);
HulVolHd = accumarray(subs',HulVolH,[],@mean);
%----------------------


%----------------------


	%sbpos = [.07 .57 .4 .38];
	sbpos = [.13 .13 .80 .80];
	
subplot('Position',sbpos);
[ph1] = plot(HulVolNk,'b');

	leg1 = legend(ph1,'Neck Volume');
	[LEGH,OBJH,OUTH,OUTM] = legend;
	hold on
	
subplot('Position',sbpos);
[ph2] = plot(HulVolHd,'r');

	legend([OUTH;ph2],OUTM{:},'Head Volume');
	[LEGH,OBJH,OUTH,OUTM] = legend;
	hold on
	%------------------------------------------%
	xt = (get(gca,'XTick'));
	xt = linspace(0,Nsteps,numel(xt)).* dT./60;
	set(gca,'XTickLabel', sprintf('%.0f|',xt))
	set(legend,'Location','SouthEast');
	%------------------------------------------%
	MS1 = 10;
	set(ph1,'LineStyle','-','Color',c1,'LineWidth',2,...
	'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor','none');
	set(ph2,'LineStyle','-','Color',c2,'LineWidth',2,...
	'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
	hTitle  = title('Spine Volume (from Delaunay Triangulation Convex Hull)');
	hXLabel = xlabel('Time (min)');
	hYLabel = ylabel('Volume (um3)');
	set(gca,'FontName','Helvetica');
	set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
	set([hXLabel, hYLabel],'FontSize',10);
	set(hTitle,'FontSize',12);
	set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
	'XMinorTick','on','YMinorTick','on','YGrid','on', ...
	'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
	haxes=axis;
	% ylim([0 haxes(4)*1.2 ]);
	% xlim([0 (haxes(2)*.9)]);
%}
%===========================================================%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					FINAL OUTPUT SECOND SET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	fig25 = figure(25); scsz = get(0,'ScreenSize');
	pos1 = [scsz(3)/11  scsz(4)/11  scsz(3)/1.5  scsz(4)/1.5];
	set(fig25,'Renderer','zbuffer','Units','pixels','Color',[1,1,1],'OuterPosition',pos1)
	%set(gcf,'Renderer','OpenGL')
	%-------------------------------------------------------------
	c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
	c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
	pause(.1);

%===========================================================%
% FIG1 TOP RIGHT: Monomeric vs Filamentous Actin
%===========================================================%
CompressN = 200; NcT = nT / CompressN;
subs = floor(linspace(1,CompressN+1,nT));
subs(end) = CompressN;
NumGA = accumarray(subs',NumGAct,[],@mean);
NumFA = accumarray(subs',NumFAct,[],@mean);

GAct_uM = NumGA ./ SpyV .*(1./6e23).*1e6;  % GActin uM
FAct_uM = NumFA ./ SpyV .*(1./6e23).*1e6;  % FActin uM

GAct_uMa = NumGAct ./ SpyV .*(1./6e23).*1e6;  % GActin uM
FAct_uMa = NumFAct ./ SpyV .*(1./6e23).*1e6;  % FActin uM
%----------------------
	sbpos = [.08 .09 .85 .85];

	sp1 = subplot('Position',sbpos);
[ph1] = plot(GAct_uM,'b');
    
    hold on;

	%-------
% 	leg1 = legend(ph1,'G-actin');
% 	[LEGH,OBJH,OUTH,OUTM] = legend;
	%-------
	
	sp2 = subplot('Position',sbpos);
[ph2] = plot(FAct_uM,'r');

    hold on

	%-------
% 	LGh5 = legend([OUTH;ph2],OUTM{:},'F-actin');
% 	[LEGH,OBJH,OUTH,OUTM] = legend;
% 	set(LGh5, 'Position', [.1 .9 .05 .04], 'Color', [1 1 1]);
% 	hold on;
	%-------


    %---------------------- Legend -----------------------
    leg1 = legend([ph1,ph2],{'G-actin','F-actin'});
    % set(leg1, 'Position', leg1.Position .* [.535 .93 .05 .04])
    set(leg1, 'Location','NorthWest', 'Color', [1 1 1],'FontSize',12,'Box','on');
    set(leg1, 'Position', leg1.Position .* [1 .95 1.5 1.4])

    %-------------------- Tic Labels ---------------------
    % set(haxes2,'XTickLabel', [],'XTick', [])
    if verLessThan('matlab', '8.3.1');
        xt = (get(gca,'XTick')).* dT./60.*NcT;
        set(gca,'XTickLabel', sprintf('%.0f|',xt))
    else    
        xt = sp2.XTick;
        xt = linspace(0,Nsteps,numel(xt)).* dT./60;
        xt = roundn(xt,-1);
        sp2.XTickLabel = xt;
    end

    %------------------- Axis Labels --------------------
        MS1 = 7; MS2 = 2;
    set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
        'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c1);
    set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
        'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c2);
 
    hTitle  = title ('\fontsize{14} Monomeric vs Filamentous Actin');
    hXLabel = xlabel('\fontsize{11} Time (min)');
    hYLabel = ylabel('\fontsize{12} Actin Concentration (µM)');
    set([hTitle, hXLabel, hYLabel],'FontName','Helvetica Neue');
    %------
    set(gca,'Box','off','TickDir','out','TickLength',[.01 .01], ...
	'XMinorTick','off','YMinorTick','off','YGrid','off', ...
	'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
    %------
	haxesM=axis;
	ylim([-20 haxesM(4)*1.1]);
    %-----------------------------------------------------

	
	
	%------------------------------------------%
	nFAct_uMa = FAct_uMa(ceil(numel(FAct_uMa)/2):end);
	mnFAct_uMa = nFAct_uMa - min(nFAct_uMa);
	nGAct_uMa = GAct_uMa(ceil(numel(GAct_uMa)/2):end);
	%-----------------------
	
	
	%[FAct Top Zoom]
	%------------------------------------------%
	sbtop = [.49 .53 .44 .3];
	haxes4 = axes('Position',sbtop,...
		'XAxisLocation','top','YAxisLocation','right','Color','none'...
		,'TickDir','out','TickLength',[.01 .01],'YMinorTick','off','YGrid','off'...
		,'Box','off','XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1 ...
		,'XTickLabel',[],'XTick',[],'FontName','Century Gothic','FontSize',11);
		hold on
	
[ph1b] = plot(nFAct_uMa);

	axis tight; ylabel('µM F-actin ');
	%------------------------------------------%

	
	
	%[GAct Bot Zoom]
	%------------------------------------------%
        sbbot = [.49 .09 .44 .3];
	haxes3 = axes('Position',sbbot,...
		'XAxisLocation','bottom','YAxisLocation','right','Color','none'...
		,'TickDir','out','TickLength',[.01 .01],'YMinorTick','off','YGrid','off'...
		,'Box','off','XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1 ...
		,'XTickLabel',[],'XTick',[],'FontName','Century Gothic','FontSize',11);
		hold on
	
[ph1b] = plot(nGAct_uMa,'r');

	axis tight; ylabel('µM G-actin ');
	%------------------------------------------%

	%[haxes,hline1,hline2] = plotyy(t,z1,t,z2,'semilogy','plot');
%===========================================================%








%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					FINAL OUTPUT SECOND SET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh26=figure('Position',[100 450 1400 400],'Color','w');
hax1=axes('Position',[.03 .1 .28 .8],'Color','none');
hax2=axes('Position',[.37 .1 .28 .8],'Color','none');
hax3=axes('Position',[.69 .1 .28 .8],'Color','none');
%-------------------------------------------------------------

%===========================================================%
% FIG1 TOP RIGHT: Monomeric vs Filamentous Actin
%===========================================================%
CompressN = 200; NcT = nT / CompressN;
subs = floor(linspace(1,CompressN+1,nT));
subs(end) = CompressN;
NumGAN0     = accumarray(subs',NumGActinN0,[],@mean);
NumTHY      = accumarray(subs',NumTHYM,[],@mean);
NumTHYACT   = accumarray(subs',NumTHYMACT,[],@mean);

NumGAN0_uM   = NumGAN0 ./ SpyV .*(1./6e23).*1e6;  % GActin uM
NumTHY_uM    = NumTHY ./ SpyV .*(1./6e23).*1e6;  % GActin uM
NumTHYACT_uM = NumTHYACT ./ SpyV .*(1./6e23).*1e6;  % FActin uM
%----------------------

    axes(hax1)
[ph1] = plot(NumTHY_uM,'b');

	axes(hax2)
[ph2] = plot(NumTHYACT_uM,'r');

	axes(hax3)
[ph3] = plot(NumGAN0_uM,'k');

    axes(hax1); title('\fontsize{14} NumTHYM');
    axes(hax2); title('\fontsize{14} NumTHYMACT');
    axes(hax3); title('\fontsize{14} NumGAN0');
    axes(hax1); xlabel('\fontsize{11} Time (min)');
    axes(hax2); xlabel('\fontsize{11} Time (min)');
    axes(hax3); xlabel('\fontsize{11} Time (min)');
    axes(hax1); ylabel('\fontsize{12} Thymosin (µM)');
    axes(hax2); ylabel('\fontsize{12} Thymosin + Actin (µM)');
    axes(hax3); ylabel('\fontsize{12} Actin (µM)');

%===========================================================%
%%






%%
keyboard
%------------
end
%--------------------------------------------------%
%%





















%--------------------------------------------------%
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%--------------------------------------------------%



%============================================================%
%				OUTPUT FIGURES
%------------------------------------------------------------%
%{
%==================================================%
%			MATRIX SURFACE FIGURE
%--------------------------------------------------%
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%--------------------------------------------------%
PSDXYZ = [PSDTips(:,1) PSDTips(:,2) PSDTips(:,3)];
PSDXY = round([PSDTips(:,1) PSDTips(:,2)]);
PSDactMx = zeros(SPYhXY+100,SPYhXY+100);
for mxp = 1:numel(PSDXY(:,1))
PSDactMx(PSDXY(mxp,2)+SPYhXY+10, PSDXY(mxp,1)+SPYhXY+10) = 1;
end
ActMask=[1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1];
ActMx = convn(PSDactMx,ActMask,'same');
ActMx = (ActMx>0).*1.0;
%--------------------------------------------------%
figure
subplot('Position',[.08 .05 .40 .90]), 
imagesc(ActMx)
colormap(bone)
subplot('Position',[.55 .05 .40 .90]), 
scatter(PSDXY(:,1), PSDXY(:,2))
%--------------------------------------------------%
%==================================================%
%				DEPOLY FIGURE
%--------------------------------------------------%
sz = [2.5e-3 2.5e-3 1.2 2];
Fh2 = FigSetup(2,sz);
%--------------------------------------------------%
figure(Fh2)
subplot('Position',[.03 .05 .28 .90]),
plot(ActData(BeKanT:nT,1)); title('Act uM');
subplot('Position',[.35 .05 .28 .90]),
plot(ActData(BeKanT:nT,2)); title('Act PR');
subplot('Position',[.68 .05 .28 .90]),
plot(ActData(BeKanT:nT,3)); title('Act N');

figure
plot(DePData(BeKanT:nT)); title('Depoly Sum');
%--------------------------------------------------%



%==================================================%
%					FIGURE
%--------------------------------------------------%
figure
scsz = get(0,'ScreenSize');
pos = [scsz(3)/3.5  scsz(4)/5.5  scsz(3)/2  scsz(4)/1.5];
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])
%--------------------------------------------------%
plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]')
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[1,1,1])
%--------------------------------------------------%


%==================================================%
%					FIGURE
%--------------------------------------------------%
figure
pos = [scsz(3)/3  scsz(4)/5  scsz(3)/2  scsz(4)/1.5];
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])

ph11c = plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]');
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[1,1,1])
hold on;
ph11a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph11b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[1,1,1])
set(ph11a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph11b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
set(ph11c,'LineStyle','-','Color',[.7 .7 .7],'LineWidth',.1);
%--------------------------------------------------%



%==================================================%
%					FIGURE
%--------------------------------------------------%
figure
scsz = get(0,'ScreenSize');
pos = [scsz(3)/2.8  scsz(4)/4.8  scsz(3)/2  scsz(4)/1.5];
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])
%--------------------------------------------------%
set(gca,'Color',[1,1,1])
hold on;
ph12a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph12b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[1,1,1])
set(ph12a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph12b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
%--------------------------------------------------%



%==================================================%
%				FIGURE (Fh1)
%--------------------------------------------------%
figure
scsz = get(0,'ScreenSize');
pos = [scsz(3)/2.8  scsz(4)/4.8  scsz(3)/2  scsz(4)/1.5];
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])
%--------------------------------------------------%
set(gca,'Color',[1,1,1])
hold on;
ph12a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph12b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[1,1,1])
set(ph12a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph12b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
%--------------------------------------------------%
%}
%------------------------------------------------------------%

AxLP = {Fh2Live,nT,Actin,inPSD,rot1,azel0,dims};
Ax = {0,nT,Actin,dims,AcMx,SPYH,rot1,azel0};



%% --------------           Actin                 AFMx                 --------------------
%
%                       Actin = N x 20       AFMx{nT} = Actin
%-----------------------------------------------------------------------------------------%
%                                      Actin
% 1  2   3   4   5   6   7   8   9   10  11     12  13  14   15   16  17   18    19   20
% N  Xa  Xo  Xt  Ya  Yo  Yt  Za  Zo  Zt  MomID  ID  Fkd Born Died Lif MaxL MeanL Null Lgth
%-----------------------------------------------------------------------------------------%
% 1  2   3   4   5   6   7   8   9   10  11     12  13  14   15   16  17   18    19   20
% 1  2   3   4   5   6   7   8   9   10  11     12  13  14   15   16  17   18    19   20
% 1  2   3   4   5   6   7   8   9   10  11     12  13  14   15   16  17   18    19   20
% 1  2   3   4   5   6   7   8   9   10  11     12  13  14   15   16  17   18    19   20
% 1  2   3   4   5   6   7   8   9   10  11     12  13  14   15   16  17   18    19   20
%-----------------------------------------------------------------------------------------%

% DISPLAY WHAT IS CONTAINED IN EACH CELL OF AFMx
% The first row in this preview table is the column number
% those values aren't actually part of the Nx20 matrix

Tnams={'N';'Xa';'Xo';'Xt';'Ya';'Yo';'Yt';'Za';'Zo';'Zt';...
       'MomID';'ID';'Fkd';'Born';'Died';'Lif';'MaxL';'MeanL';'Null';'Lgth'};
l=[1:20; Actin(50:60,:)]; g=1:8;
T=table(l(g,1),l(g,2),l(g,3),l(g,4),l(g,5),l(g,6),l(g,7),l(g,8),l(g,9),l(g,10),l(g,11),l(g,12),...
l(g,13),l(g,14),l(g,15),l(g,16),l(g,17),l(g,18),l(g,19),l(g,20),'VariableNames',Tnams(1:20));
disp(' '); disp(T); disp(' ');

%% ----------     GATHER AND FORMAT FILAMENT DATA     ----------

PSD_Zdim = dims(7);
SPYtips = {};
PSDtips = {};
ActXYZ = {};

for nn = 1:numel(AFMx)

    Ain = AFMx{nn};

    ActXYZ{nn} = [Ain(:,3) Ain(:,4) Ain(:,6) Ain(:,7) Ain(:,9) Ain(:,10)];
    ActinTips = [ActXYZ{nn}(:,2) ActXYZ{nn}(:,4) ActXYZ{nn}(:,6)];

    [Zrow,Zcol] = find(ActinTips(:,3) < PSD_Zdim);

    SPYtps = ActinTips(Zrow,:);
    SPYtips{nn} = SPYtps;

    [Zrow,Zcol] = find(ActinTips(:,3) >= PSD_Zdim);

    PSDtps = ActinTips(Zrow,:);
    PSDtips{nn} = PSDtps;

end


%% SAVE MAT FILE TO DISK

% ActinMx.Ax = Ax;
% ActinMx.AFMx = AFMx;
ActinMx.ActXYZ = ActXYZ;
ActinMx.SPYtips = SPYtips;
ActinMx.PSDtips = PSDtips;

cd(fileparts(which('savetipdir.m')));
save('ActinMx.mat', 'ActinMx')
cd(fileparts(which(mfilename)));

%% RETURN VALUES AS VARARGOUT

varargout = {BTs,AxLP,Ax,AFMx,ActinMx};
return

end
%####################################################################%












function [pNuM] = VCP(vol,uM,pN)

% INPUTS
% vol: volume in um^3
% uM: concentration in uM
% pN: particle count
% 
% enter zero for the unknown value
% if pN is unknown enter... pN = VCP(.1,10,0)
% if uM is unknown enter... pN = VCP(.1,0,6e5)


% [use 1e-15 to convert from um^3; use 1e-24 to convert from nm^3]
V = vol * 1e-15;	% convert um^3 to L 

if pN==0
	
	pNuM = uM / 1e6 * 6e23 * V;      % particle count

elseif uM==0

	pNuM = pN / V *(1/6e23)*1e6;     % particle concentration

end



% %----------------------------------------
% % volume of cylinder: V = pi * r^2 * h
% 
% Vn = pi * XYn^2 * Zn;
% Vh = pi * XYh^2 * Zh;
% Vnm2L = (Vn+Vh) * 1e-24;
% 
% Ka = 10 * uM * dT;		% empirical actin poly rate
% Kd = 1 * dT;				% empirical actin tip depoly rate
% %----------------------------------------


end











%####################################################################%
%					PLOTTING & HELPER FUNCTIONS
%####################################################################%




%==================================================%
%					MainLivePlot
%--------------------------------------------------%
function MainLivePlot(Fh,nT,Actin,inPSD,rot,azel,AxLims,...
                      SPYhZN,SPYhZS,SPYnXY,SPYhXY,doLiveHullPlot)
%==================================================%

figure(Fh)
%%
%--------------------
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,~] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,~] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%--------------------
subplot('Position',[.04 .03 .45 .95]), 

P3x = [Actin(:,3) Actin(:,4)]';
P3y = [Actin(:,6) Actin(:,7)]';
P3z = [Actin(:,9) Actin(:,10)]';
ph11c = plot3(P3x, P3y, P3z);
% ph11c = plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]');
axis(AxLims); axis vis3d;
%view(azel)
grid off
hold on;
ph11a = scatter3(SPYTips(:,1)', SPYTips(:,2)', SPYTips(:,3)',20,'ob');
hold on;
ph11b = scatter3(PSDTips(:,1)', PSDTips(:,2)', PSDTips(:,3)',20,'or');
view(azel+rot)
set(ph11a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph11b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
set(ph11c,'LineStyle','-','Color',[.7 .7 .7],'LineWidth',.3);
hold off;
%--------------------



%==================================================%
if doLiveHullPlot

Ori_xyz = [Actin(:,3) Actin(:,6) Actin(:,9)];
Tip_xyz = [Actin(:,4) Actin(:,7) Actin(:,10)];
ATPxyz = [Ori_xyz; Tip_xyz];



%{
Ori_xyz = [Actin(:,3) Actin(:,6) Actin(:,9)];
Tip_xyz = [Actin(:,4) Actin(:,7) Actin(:,10)];
ATPxyz = [Ori_xyz; Tip_xyz];


% radial distance to spine shaft membrane
XYtipLoc = sqrt(ATPxyz(:,1).^2 + ATPxyz(:,2).^2);

ZtipInHead = ATPxyz(:,3) >= SPYhZS;
ZtipInNeck = ATPxyz(:,3) < SPYhZS;

XYneckOut = XYtipLoc > (SPYnXY+5);	% Logical array of filaments beyond SPYnXY
XYheadOut = XYtipLoc > (SPYhXY+5);  % Logical array of filaments beyond SPYhXY
ZtopOut = ATPxyz(:,3) > (SPYhZN+5);	% Logical array of filaments above SPYhZN
ZbotOut = ATPxyz(:,3) < 0;		% Logical array of filaments below zero

TipOut = ((XYneckOut & ZtipInNeck) + (XYheadOut & ZtipInHead) + ZtopOut + ZbotOut)>0;


ATPxyz(TipOut,:) = [];
%}

ATPxyz(end:end+4,:) = 0;
szATP = size(ATPxyz,1);
ATPxyz(szATP+1,:) = [SPYnXY/2 SPYnXY/2 0];
ATPxyz(szATP+2,:) = [-SPYnXY/2 SPYnXY/2 0];
ATPxyz(szATP+3,:) = [SPYnXY/2 -SPYnXY/2 0];
ATPxyz(szATP+4,:) = [-SPYnXY/2 -SPYnXY/2 0];

%--------------------
AOTxyz = ATPxyz;
minAOT = SPYhXY+20;
AOTxyz(:,1) = AOTxyz(:,1) + minAOT;
AOTxyz(:,2) = AOTxyz(:,2) + minAOT;
AOTxyz(:,3) = AOTxyz(:,3);

AOTxyzU = unique(AOTxyz,'rows');

Uxyz = AOTxyzU;

TinHead = Uxyz(:,3) >= (SPYhZS);
TxyzHead = Uxyz(TinHead,:);

TinNeck = Uxyz(:,3) < (SPYhZS-5);
TxyzNeck = Uxyz(TinNeck,:);

TinChin = (Uxyz(:,3) >= (SPYhZS-120)) & (Uxyz(:,3) <= (SPYhZS+100));
TxyzChin = Uxyz(TinChin,:);




%---
%figure(Fh)
subplot('Position',[.55 .03 .45 .95]), 
%---

TriSZh = size(TxyzHead,1);
TriSZn = size(TxyzNeck,1);
TriSZc = size(TxyzChin,1);
TriSZall = TriSZn+TriSZc+TriSZh;

HulVolH=0;HulVolN=0;HulVolC=0;

if (TriSZall>=9)
if (TriSZn<24) || (TriSZc<24) || (TriSZh<24)
	UTriTip = delaunayTriangulation(Uxyz);
	[Utri,Upoints] = freeBoundary(UTriTip);

	trimesh(Utri,Upoints(:,1),Upoints(:,2),Upoints(:,3), ...
		   'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
	hold on
	
	
% 	[HulPtsNH,HulVolNH] = convexHull(UTriTip);
% 	HulVolN = HulVolNH / 2.1;
% 	HulVolH = HulVolNH / 2.2;
% 	%trisurf(cnvxH,UTriTip.Points(:,1),UTriTip.Points(:,2),UTriTip.Points(:,3))
		
else

	TriHead = delaunayTriangulation(TxyzHead);
	[FBtHead,FBpHead] = freeBoundary(TriHead);

	trimesh(FBtHead,FBpHead(:,1),FBpHead(:,2),FBpHead(:,3), ...
		   'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
	hold on

%---

	TriNeck = delaunayTriangulation(TxyzNeck);
	[FBtNeck,FBpNeck] = freeBoundary(TriNeck);

	trimesh(FBtNeck,FBpNeck(:,1),FBpNeck(:,2),FBpNeck(:,3), ...
		   'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
	hold on

%---

	TriChin = delaunayTriangulation(TxyzChin);
	[FBtChin,FBpChin] = freeBoundary(TriChin);

	trimesh(FBtChin,FBpChin(:,1),FBpChin(:,2),FBpChin(:,3), ...
		   'FaceColor',[.2 1 .2],'FaceAlpha', 0.3);
	hold on

%---

% 	[HulPtsH,HulVolH] = convexHull(TriHead);
% 	[HulPtsN,HulVolN] = convexHull(TriNeck);
% 	[HulPtsC,HulVolC] = convexHull(TriChin);

%---
end

% varargout = {[HulVolN;HulVolH;HulVolC]};

end
%---

%scatter3(Tip_xyz(:,1),Tip_xyz(:,2),Tip_xyz(:,3),50,'fill','r')
scatter3(AOTxyzU(:,1),AOTxyzU(:,2),AOTxyzU(:,3),25,'fill',...
			'MarkerFaceColor',[.94 .4 .4]);
hold off
%---

TriLim = [0 600 0 600 0 1200];
axis(TriLim); axis vis3d;
view(azel+rot)
grid off
hold off
%%
%drawnow;

%-------------------------
end; % doLiveHullPlot
%==================================================%



%if nT >= 15000; keyboard; end;
end
%==================================================%





%==================================================%
%					ActArpPlot
%--------------------------------------------------%
function ActArpPlot(Fh,nT,Nsteps,NumGArp,NumFArp,NumFAct,NumGAct,NumfKa,ArpBRsum,Num_nmmf)
%---
figure(Fh)
%---



%-------------------------
% Arp
	Hax1 = axes('Position',[.55 .05 .4 .4]);
[ph1] = plot(NumFArp(1:nT)','r');
	hold on
[ph2] = plot(NumGArp(1:nT)','b');
	set(Hax1,'XLim', [0 Nsteps]);
	haxesM=axis;
	hold on
	
	
	Hax1b = axes('Position',[.55 .05 .4 .4]...
			,'XAxisLocation','top','YAxisLocation','right','Color','none'...
			,'XTickLabel', [],'XTick', []...
			,'XLim', [0 Nsteps],'YLim', [0 .1]);
			hold on;
[ph3] = plot(ArpBRsum(1:nT)','g');
%-------------------------	
		


%-------------------------
% Actin
	Hax2 = axes('Position',[.55 .55 .4 .4]);
[ph4] = plot(NumFAct(1:nT)','r');
	hold on
[ph5] = plot(NumGAct(1:nT)','b');
	set(Hax2,'XLim', [0 Nsteps])
%-------------------------
	
	
if nT > 10000
	keyboard
end




% ALL AVAILABLE PLOTS
%---------------
%{
[ph1] = plot(NumfKa,'r');
[ph2b] = semilogy(NuBeKa);
[ph2] = plot(HulVolHd,'r');
[ph1] = plot(HulVolNk,'b');
[PhBL3] = plot(Mummf,'Parent',HaxBL3,'Color',[.5 .5 .5]);
%----------------------
[ph1] = plot(NumGAct,'b');
[ph2] = plot(NumFAct,'r');
GAct_uM = NumGA ./ SpyV .*(1./6e23).*1e6;  % GActin uM
FAct_uM = NumFA ./ SpyV .*(1./6e23).*1e6;  % FActin uM

GAct_uMa = NumGAct ./ SpyV .*(1./6e23).*1e6;  % GActin uM
FAct_uMa = NumFAct ./ SpyV .*(1./6e23).*1e6;  % FActin uM
%----------------------
[ph1] = plot(ArpBRsum,'r');
[ph2] = line(1:numel(NFArp),NFArp,'Parent',haxes2,'Color','k');
	hold on
[ph3] = line(1:numel(NGArp),NGArp,'Parent',haxes2,'Color','k');
%}
%---------------

end
%==================================================%






%==================================================%
%				MAKE STARTING FILAMENTS
%--------------------------------------------------%
function Actin = MakeStartFils(Actin,NStFils,StartMonos,d2r,fZa,fXYo,SPYhZN,SPYhZS,SPYhXY)



% FIL-1
%--
Actin(1,1) = StartMonos;
Actin(1,2) = 1;				% Xa
Actin(1,5) = 1;				% Ya
Actin(1,8) = 91;			% Za
Actin(1,3) = 1;				% Xo
Actin(1,6) = 1;				% Yo
Actin(1,9) = 1;				% Zo
%--

if NStFils > 1
% FIL-2
%--
Actin(2,1) = StartMonos;
Actin(2,2) = -45;			% Xa
Actin(2,5) = 45;			% Ya
Actin(2,8) = d2r*-fZa;		% Za
Actin(2,3) = fXYo;			% Xo
Actin(2,6) = fXYo;			% Yo
Actin(2,9) = 1;				% Zo
%--
end

if NStFils > 2
% FIL-3
%--
Actin(3,1) = StartMonos;
Actin(3,2) = 45;			% Xa
Actin(3,5) = -45;			% Ya
Actin(3,8) = d2r*fZa;		% Za
Actin(3,3) = fXYo;			% Xo
Actin(3,6) = -fXYo;			% Yo
Actin(3,9) = 1;				% Zo
%--
end

if NStFils > 3
% FIL-4
%--
Actin(4,1) = StartMonos;
Actin(4,2) = 45;			% Xa
Actin(4,5) = 0;				% Ya
Actin(4,8) = d2r*-fZa;		% Za
Actin(4,3) = -fXYo;			% Xo
Actin(4,6) = fXYo;			% Yo
Actin(4,9) = 1;				% Zo
%--
end

if NStFils > 4
% FIL-5
%--
Actin(5,1) = StartMonos;
Actin(5,2) = 0;				% Xa
Actin(5,5) = 45;			% Ya
Actin(5,8) = d2r*fZa;		% Za
Actin(5,3) = -fXYo;			% Xo
Actin(5,6) = -fXYo;			% Yo
Actin(5,9) = 1;				% Zo
%--
end


if NStFils > 5
% FIL-6
%--
Actin(6,1)=round(SPYhZN-SPYhZS-2);
Actin(6,2) = 225*d2r;		% Xa
Actin(6,5) = 0*d2r;			% Ya
Actin(6,8) = 45*d2r;		% Za
Actin(6,3) = (SPYhXY/2)-10;	% Xo
Actin(6,6) = (SPYhXY/2)-10;	% Yo
Actin(6,9) = SPYhZS+1;		% Zo
end

if NStFils > 6
% FIL-7
%--
Actin(7,1)=round(SPYhZN-SPYhZS-2);
Actin(7,2) = 45*d2r;		% Xa
Actin(7,5) = 0*d2r;			% Ya
Actin(7,8) = 45*d2r;		% Za
Actin(7,3) = 10-(SPYhXY/2);	% Xo
Actin(7,6) = 10-(SPYhXY/2);	% Yo
Actin(7,9) = SPYhZS+1;		% Zo
end

if NStFils > 7
% % FIL-8
% %--
Actin(8,1)=round(SPYhZN-SPYhZS-2);
Actin(8,2) = 315*d2r;		% Xa
Actin(8,5) = 0*d2r;			% Ya
Actin(8,8) = 45*d2r;		% Za
Actin(8,3) = 10-(SPYhXY/2);	% Xo
Actin(8,6) = (SPYhXY/2)-10;	% Yo
Actin(8,9) = SPYhZS+1;		% Zo
end

if NStFils > 8
% FIL-9
%--
Actin(9,1)=round(SPYhZN-SPYhZS-2);
Actin(9,2) = 135*d2r;		% Xa
Actin(9,5) = 0*d2r;			% Ya
Actin(9,8) = 45*d2r;		% Za
Actin(9,3) = (SPYhXY/2)-10;	% Xo
Actin(9,6) = 10-(SPYhXY/2);	% Yo
Actin(9,9) = SPYhZS+1;		% Zo
end


end
%==================================================%






%==================================================%
%				FigSetup FUNCTION
%--------------------------------------------------%
function [varargout] = FigSetup(varargin)

scsz = get(0,'ScreenSize');


if nargin == 1 
	Fnum=varargin{1};
	pos = scsz./[2.5e-3 2.5e-3 1.5 2];
elseif nargin == 2
	Fnum=varargin{1};
	pos=scsz./varargin{2};
else
	Fnum=1;
	pos = scsz./[2.5e-3 2.5e-3 1.5 2];
end

Fh = figure(Fnum);
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])
varargout = {Fh};

end
%==================================================%


%==================================================%
%					LIVE PLOT
%--------------------------------------------------%
function LivePlot(Fh,nT,Actin,inPSD,rot,azel,AxLims)


%--------------------
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%--------------------
figure(Fh)
subplot('Position',[.03 .03 .53 .95]), 
ph11c = plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]');
axis(AxLims); axis vis3d;
view(azel)
grid off
hold on;
ph11a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph11b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
%xlabel('X');ylabel('Y');zlabel('Z');
view(azel+rot)
set(ph11a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph11b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
set(ph11c,'LineStyle','-','Color',[.7 .7 .7],'LineWidth',.1);
hold off;
%--------------------
figure(Fh)
subplot('Position',[.65 .6 .28 .38]), 
ph12a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph12b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis(AxLims)
view([0 90])
grid off
set(ph12a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph12b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
%--------------------
set(gca,'XTickLabel', sprintf('%.1f|',nT),'FontSize',10)
hold off;
%--------------------

end
%==================================================%





%==================================================%
%					TipMatrix
%--------------------------------------------------%
function [varargout] = TipMatrix(Fh,doTM,nT,Actin,dims,AcMx,SPYH,varargin)


%{
% keyboard

scsz = get(0,'ScreenSize');
basepos = scsz./[2.5e-3 2.5e-3 1.5 2];
baseazel = [-32 12];
baserot=[0 0];

if nargin >= 2
	rot=varargin{1};
	azel=varargin{2};
elseif nargin == 1 
	rot=varargin{1};
	azel = baseazel;
else
	rot=baserot;
	azel = baseazel;
end






%--------------------
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%--------------------
figure(Fh)
subplot('Position',[.08 .15 .45 .70]), 
ph11c = plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]');
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
view(azel)
grid off
set(gca,'Color',[1,1,1])
hold on;
ph11a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph11b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis([-500 500 -500 500 0 1000])
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
subplot('Position',[.6 .55 .38 .38]), 
ph12a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph12b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis([-500 500 -500 500 0 1000])
view([0 90])
grid off
set(gca,'Color',[1,1,1])
set(ph12a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph12b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
%--------------------
set(gca,'XTickLabel', sprintf('%.1f|',nT),'FontSize',10)
hold off;
%--------------------
%}


%==================================================%
% Spine Dimensions
HX = dims(4);
HY = dims(5);
SPYz = dims(7);
SPYxy = dims(8);
%-----------------------------------%
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];


[Zin,Zic] = find(ActinTips(:,3) > SPYz);
XYsq = sqrt(ActinTips(:,1).^2 + ActinTips(:,2).^2);
[XYin,XYic] = find(XYsq < (HX-SPYxy));
XYZin = intersect(XYin, Zin);
% XYt = intersect(Xrow1, Yrow1);
%[tf, loc] = ismember(Xrow1, Yrow1)


PSDTips = ActinTips(XYZin,:);	% get just the tips, just to see how they feel
PSDXY = round([PSDTips(:,1) PSDTips(:,2)]);
PSDY = PSDXY(:,2)+HY;	% Add Y radius width
PSDX = PSDXY(:,1)+HX;	% Add X radius width
SPYHX = SPYH(1)*2;		% I'm specifying a radius in SPYH(), must double
SPYHY = SPYH(2)*2;		% I'm specifying a radius in SPYH(), must double
PSDX(PSDX>SPYHX) = SPYHX-1;	% Everything here that is somehow greater SPYHX place at edge of SPYHX
PSDY(PSDY>SPYHY) = SPYHY-1;	% Everything here that is somehow greater SPYHY place at edge of SPYHY
PSDX(PSDX<1) = 1;	% Everything here that is somehow negative place at 1
PSDY(PSDY<1) = 1;	% Everything here that is somehow negative place at 1
%==================================================%

PSDX = round(PSDX/5);
PSDY = round(PSDY/5);

for mxp = 1:numel(PSDXY(:,1))
AcMx(PSDY(mxp), PSDX(mxp)) = 1;
end

% SAVE THIS (COULD BE USEFUL)
% ActMask=ones(4);
% ActMx = convn(AcMx,ActMask,'same');
% ActMx = (ActMx>0).*1.0;

ActMx = AcMx;
%===================================%
%				FIGURE
%-----------------------------------%
if doTM
    clrmap = [1 1 1; .9 .2 .2];
    figure(Fh)
    subplot('Position',[.7 .1 .28 .38]), 
    imagesc(ActMx)
    colormap(clrmap);
    hold on
    subplot('Position',[.7 .1 .28 .38]), 
    scatter(PSDX,PSDY, 'r')
    hold off
end
%-----------------------------------%


%==================================================%
% if nT >20000; keyboard; end
varargout = {ActMx};
end
%==================================================%





%----------------------------------%
%		in3Dbox
%----------------------------------%
% Tests whether particles are in a box polygon
% and returns a logical vector
%----------------------------------%
function [in3db] = in3Dbox(ActinTips,xL,xH,yL,yH,zL,zH)

xLo = ActinTips(:,1) > xL;
xHi = ActinTips(:,1) < xH;
yLo = ActinTips(:,2) > yL;
yHi = ActinTips(:,2) < yH;
zLo = ActinTips(:,3) > zL;
zHi = ActinTips(:,3) < zH;

in3db = xLo & xHi & yLo & yHi & zLo & zHi;


% %----------------------------------%
% %		inboxfun
% %----------------------------------%
% % Tests whether particles are in a box polygon
% % and returns a logical vector
% %----------------------------------%
% % function [inbox] = inboxfun(LB,RT,xyl)
% if LB(1)>RT(1)
% 	LBt=LB;
% 	RTt=RT;
% 	LB(1)=RTt(1);
% 	RT(1)=LBt(1);
% end
% if LB(2)>RT(2)
% 	LBt=LB;
% 	RTt=RT;
% 	LB(2)=RTt(2);
% 	RT(2)=LBt(2);
% end
% 
% xylLB1 = xyl(1,:) > LB(1);
% xylRT1 = xyl(1,:) < RT(1);
% xylLB2 = xyl(2,:) > LB(2);
% xylRT2 = xyl(2,:) < RT(2);
% 
% inbox = xylLB1 & xylRT1 & xylLB2 & xylRT2;


end






%==================================================%
%			MATRIX CONVOLUTION MASK
%--------------------------------------------------%
function [hkMask dT LBR] = MaskFun(S,dT,LBR,doGaussianMask,PSAsz,scsz,AMX,MSK)

%-------------------------------%
%		Mask Setup
%-------------------------------%
hkMask=ones(AMX{6});
%hkMask=[0 1 0; 1 0 1; 0 1 0];
%----------------%
if doGaussianMask
%----------------%



% %--------------------
% A  = 2;	
% x0 = 0; 
% y0 = 0; 
% sx = .2; 
% sy = .2; 
% res= 2; 
% rx=sx; 
% ry=sy;	
% 
% t = 0;
% a = cos(t)^2/2/sx^2 + sin(t)^2/2/sy^2;
% b = -sin(2*t)/4/sx^2 + sin(2*t)/4/sy^2 ;
% c = sin(t)^2/2/sx^2 + cos(t)^2/2/sy^2;
% 
% [X, Y] = meshgrid((-sx*res):(rx):(sx*res), (-sy*res):(ry):(sy*res));
% Z = A*exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;
% %--------------------





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

%---
% figure(62)
% surf(X,Y,Z);
% view(-45,30); 
% axis([-GNspr GNspr -GNspr GNspr 0 GNpk])
% %shading interp; %axis equal; %axis vis3d; 
% xlabel('x-axis');ylabel('y-axis');zlabel('z-axis')
%--------------------


%--------------------
hkMask=Z;
hk = convn(S,hkMask,'same');
hkor = hk(PSAsz+1,PSAsz+1);
%-----
%LBR(1) = hkor-sqrt(GNpk); LBR(2) = hkor+sqrt(GNpk);
%--------------------


%----------------%
% FIGURE: 3D Gaussian Distribution
fh5 = figure(5); %set(fh5,'OuterPosition',(scsz./[2e-3 2e-3 2 2]))
fh5op = get(fh5,'OuterPosition');
set(fh5,'OuterPosition',[fh5op(1) fh5op(2) fh5op(3)*1.5 fh5op(4)*1.3]);
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
%==================================================%

%####################################################################%


%==================================================%
%					TriHull
%--------------------------------------------------%
function [varargout] = TriHull(nT,Actin,SPYhZN,SPYhZS,SPYnXY,SPYhXY)
%==================================================%

Ori_xyz = [Actin(:,3) Actin(:,6) Actin(:,9)];
Tip_xyz = [Actin(:,4) Actin(:,7) Actin(:,10)];
ATPxyz = [Ori_xyz; Tip_xyz];

ATPxyz(end:end+4,:) = 0;
szATP = size(ATPxyz,1);
ATPxyz(szATP+1,:) = [SPYnXY/2 SPYnXY/2 0];
ATPxyz(szATP+2,:) = [-SPYnXY/2 SPYnXY/2 0];
ATPxyz(szATP+3,:) = [SPYnXY/2 -SPYnXY/2 0];
ATPxyz(szATP+4,:) = [-SPYnXY/2 -SPYnXY/2 0];

%--------------------
AOTxyz = ATPxyz;
minAOT = SPYhXY+20;
AOTxyz(:,1) = AOTxyz(:,1) + minAOT;
AOTxyz(:,2) = AOTxyz(:,2) + minAOT;
AOTxyz(:,3) = AOTxyz(:,3);

AOTxyzU = unique(AOTxyz,'rows');

Uxyz = AOTxyzU;

TinHead = Uxyz(:,3) >= (SPYhZS);
TxyzHead = Uxyz(TinHead,:);

TinNeck = Uxyz(:,3) < (SPYhZS-5);
TxyzNeck = Uxyz(TinNeck,:);

TinChin = (Uxyz(:,3) >= (SPYhZS-120)) & (Uxyz(:,3) <= (SPYhZS+100));
TxyzChin = Uxyz(TinChin,:);


%---

TriSZh = size(TxyzHead,1);
TriSZn = size(TxyzNeck,1);
TriSZc = size(TxyzChin,1);
TriSZall = TriSZn+TriSZc+TriSZh;

HulVolH=0;HulVolN=0;HulVolC=0;

if (TriSZall>=9)
if (TriSZn<24) || (TriSZc<24) || (TriSZh<24)
	
	UTriTip = delaunayTriangulation(Uxyz);
	%[Utri,Upoints] = freeBoundary(UTriTip);

	[HulPtsNH,HulVolNH] = convexHull(UTriTip);
	HulVolN = HulVolNH / 3.0;
	HulVolH = HulVolNH / 4.0;
		
else

	TriHead = delaunayTriangulation(TxyzHead);
	%[FBtHead,FBpHead] = freeBoundary(TriHead);

%---

	TriNeck = delaunayTriangulation(TxyzNeck);
	%[FBtNeck,FBpNeck] = freeBoundary(TriNeck);

%---

	TriChin = delaunayTriangulation(TxyzChin);
	%[FBtChin,FBpChin] = freeBoundary(TriChin);

%---

	[HulPtsH,HulVolH] = convexHull(TriHead);
	[HulPtsN,HulVolN] = convexHull(TriNeck);
	[HulPtsC,HulVolC] = convexHull(TriChin);

%---
end

end
%---

varargout = {[HulVolN;HulVolH;HulVolC]};

%if nT >= 15000; keyboard; end;
end
%==================================================%



%==================================================%
%					FilPlotAlone
%--------------------------------------------------%
function FilPlotAlone(Fh,nT,Actin,inPSD,rot,azel,AxLims,SPYhZN,SPYhZS,SPYnXY,SPYhXY)
%==================================================%

figure(Fh)
%%
%--------------------
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,~] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,~] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%--------------------
subplot('Position',[.04 .03 .45 .95]), 

P3x = [Actin(:,3) Actin(:,4)]';
P3y = [Actin(:,6) Actin(:,7)]';
P3z = [Actin(:,9) Actin(:,10)]';
ph11c = plot3(P3x, P3y, P3z);
% ph11c = plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]');
axis(AxLims); axis vis3d;
%view(azel)
grid off
hold on;
ph11a = scatter3(SPYTips(:,1)', SPYTips(:,2)', SPYTips(:,3)',20,'ob');
hold on;
ph11b = scatter3(PSDTips(:,1)', PSDTips(:,2)', PSDTips(:,3)',20,'or');
view(azel+rot)
set(ph11a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph11b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
set(ph11c,'LineStyle','-','Color',[.7 .7 .7],'LineWidth',.3);
hold off;
%--------------------



%if nT >= 15000; keyboard; end;
end
%==================================================%




%==================================================%
%					TriPlot
%--------------------------------------------------%
function TriPlot(Fh,nT,Actin,inPSD,rot,azel,AxLims,SPYhZN,SPYhZS,SPYnXY,SPYhXY)
%==================================================%
Ori_xyz = [Actin(:,3) Actin(:,6) Actin(:,9)];
Tip_xyz = [Actin(:,4) Actin(:,7) Actin(:,10)];
ATPxyz = [Ori_xyz; Tip_xyz];

% radial distance to spine shaft membrane
XYtipLoc = sqrt(ATPxyz(:,1).^2 + ATPxyz(:,2).^2);

ZtipInHead = (ATPxyz(:,3) <= SPYhZN) & (ATPxyz(:,3) >= SPYhZS);
ZtipInNeck = ATPxyz(:,3) < SPYhZS;

XYneckOut = XYtipLoc > SPYnXY;	% Logical array of filaments beyond SPYnXY
XYheadOut = XYtipLoc > SPYhXY;  % Logical array of filaments beyond SPYhXY
ZtopOut = ATPxyz(:,3) > SPYhZN;	% Logical array of filaments above SPYhZN
ZbotOut = ATPxyz(:,3) < 0;		% Logical array of filaments below zero

TipOut = ((XYneckOut & ZtipInNeck) + (XYheadOut & ZtipInHead) + ZtopOut + ZbotOut)>0;

ATPxyz(TipOut,:) = [];

ATPxyz(end:end+4,:) = 0;

szATP = size(ATPxyz,1);

ATPxyz(szATP+1,:) = [SPYnXY/2 SPYnXY/2 0];
ATPxyz(szATP+2,:) = [-SPYnXY/2 SPYnXY/2 0];
ATPxyz(szATP+3,:) = [SPYnXY/2 -SPYnXY/2 0];
ATPxyz(szATP+4,:) = [-SPYnXY/2 -SPYnXY/2 0];

%==================================================%
AOTxyz = ATPxyz;
%minAOT = -min(AOTxyz);
minAOT = SPYhXY+20;
AOTxyz(:,1) = AOTxyz(:,1) + minAOT;
AOTxyz(:,2) = AOTxyz(:,2) + minAOT;
AOTxyz(:,3) = AOTxyz(:,3);

AOTxyzU = unique(AOTxyz,'rows');

Uxyz = AOTxyzU;

TinH = Uxyz(:,3) >= (SPYhZS);
HTxyz = Uxyz(TinH,:);

TinN = Uxyz(:,3) < (SPYhZS);
NTxyz = Uxyz(TinN,:);

TinNH = (Uxyz(:,3) >= (SPYhZS-80)) & (Uxyz(:,3) <= (SPYhZS+50));
NHTxyz = Uxyz(TinNH,:);


figure(Fh)
subplot('Position',[.56 .02 .47 .55]),


if numel(HTxyz) > 12
DTriTip = delaunayTriangulation(HTxyz);
[FBtri,FBpoints] = freeBoundary(DTriTip);

trimesh(FBtri,FBpoints(:,1),FBpoints(:,2),FBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
hold on
end

if numel(NTxyz) > 12
NTriTip = delaunayTriangulation(NTxyz);
[NFBtri,NFBpoints] = freeBoundary(NTriTip);

trimesh(NFBtri,NFBpoints(:,1),NFBpoints(:,2),NFBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
hold on
end

if numel(NHTxyz) > 12
NHTriTip = delaunayTriangulation(NHTxyz);
[NHFBtri,NHFBpoints] = freeBoundary(NHTriTip);

trimesh(NHFBtri,NHFBpoints(:,1),NHFBpoints(:,2),NHFBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.3);
hold on
end

scatter3(AOTxyz(:,1),AOTxyz(:,2),AOTxyz(:,3),20,'fill','r')
TriLim = [0 600 0 600 0 1200];
axis(TriLim); axis vis3d;
view(azel+rot)
grid off
hold off

end
%==================================================%
%{
%==================================================%
%					LIVE PLOT
%--------------------------------------------------%
function TriPlot(Fh,nT,Actin,inPSD,rot,azel,AxLims,SPYhZN,SPYhZS,SPYnXY,SPYhXY)
%==================================================%
Ori_xyz = [Actin(:,3) Actin(:,6) Actin(:,9)];
Tip_xyz = [Actin(:,4) Actin(:,7) Actin(:,10)];
ATPxyz = [Ori_xyz; Tip_xyz];

% radial distance to spine shaft membrane
XYtipLoc = sqrt(ATPxyz(:,1).^2 + ATPxyz(:,2).^2);

ZtipInHead = (ATPxyz(:,3) <= SPYhZN) & (ATPxyz(:,3) >= SPYhZS);
ZtipInNeck = ATPxyz(:,3) < SPYhZS;

XYneckOut = XYtipLoc > SPYnXY;	% Logical array of filaments beyond SPYnXY
XYheadOut = XYtipLoc > SPYhXY;  % Logical array of filaments beyond SPYhXY
ZtopOut = ATPxyz(:,3) > SPYhZN;	% Logical array of filaments above SPYhZN
ZbotOut = ATPxyz(:,3) < 0;		% Logical array of filaments below zero

TipOut = ((XYneckOut & ZtipInNeck) + (XYheadOut & ZtipInHead) + ZtopOut + ZbotOut)>0;

ATPxyz(TipOut,:) = [];

ATPxyz(end:end+4,:) = 0;

szATP = size(ATPxyz,1);

ATPxyz(szATP+1,:) = [SPYnXY/2 SPYnXY/2 0];
ATPxyz(szATP+2,:) = [-SPYnXY/2 SPYnXY/2 0];
ATPxyz(szATP+3,:) = [SPYnXY/2 -SPYnXY/2 0];
ATPxyz(szATP+4,:) = [-SPYnXY/2 -SPYnXY/2 0];

%==================================================%
AOTxyz = ATPxyz;
minAOT = -min(AOTxyz); 
AOTxyz(:,1) = AOTxyz(:,1) + minAOT(1);
AOTxyz(:,2) = AOTxyz(:,2) + minAOT(2);
AOTxyz(:,3) = AOTxyz(:,3);

AOTxyzU = unique(AOTxyz,'rows');

Uxyz = AOTxyzU;

TinH = Uxyz(:,3) >= (SPYhZS);
HTxyz = Uxyz(TinH,:);

TinN = Uxyz(:,3) < (SPYhZS);
NTxyz = Uxyz(TinN,:);

TinNH = (Uxyz(:,3) >= (SPYhZS-80)) & (Uxyz(:,3) <= (SPYhZS+50));
NHTxyz = Uxyz(TinNH,:);


figure(Fh)
subplot('Position',[.7 .05 .28 .48]),


if numel(HTxyz) > 9
DTriTip = delaunayTriangulation(HTxyz);
[FBtri,FBpoints] = freeBoundary(DTriTip);

trimesh(FBtri,FBpoints(:,1),FBpoints(:,2),FBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
hold on
end

if numel(NTxyz) > 9
NTriTip = delaunayTriangulation(NTxyz);
[NFBtri,NFBpoints] = freeBoundary(NTriTip);

trimesh(NFBtri,NFBpoints(:,1),NFBpoints(:,2),NFBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.35);
hold on
end

if numel(NHTxyz) > 9
NHTriTip = delaunayTriangulation(NHTxyz);
[NHFBtri,NHFBpoints] = freeBoundary(NHTriTip);

trimesh(NHFBtri,NHFBpoints(:,1),NHFBpoints(:,2),NHFBpoints(:,3), ...
       'FaceColor',[.2 1 .2],'FaceAlpha', 0.3);
hold on
end

scatter3(AOTxyz(:,1),AOTxyz(:,2),AOTxyz(:,3),20,'fill','r')
TriLim = [0 600 0 600 0 1200];
axis(TriLim); %axis vis3d;
view(azel+rot)
grid off
hold off

end
%==================================================%
%}




%==================================================%
%					LIVE PLOT
%--------------------------------------------------%
function LiveAPlot(Ph,nT,Ac,PT,ST,rot,azel,Fh)

%Ac34 = [Ac(:,3) Ac(:,4)]';
%Ac67 = [Ac(:,6) Ac(:,7)]';
%Ac91 = [Ac(:,9) Ac(:,10)]';

%set(get(gca))



%for nPh = 1:numel(Ph{1})
%set(Ph{1}(nPh),'XData',Ac34(:,nPh),'YData',Ac67(:,nPh),'ZData',Ac91(:,nPh));
%end

set(Ph{4},'XData',ST(:,1)','YData',ST(:,2)','ZData',ST(:,3)');
set(Ph{5},'XData',PT(:,1)','YData',PT(:,2)','ZData',PT(:,3)');
view([0 90])
%set(gca,'XTickLabel', sprintf('%.1f|',nT),'FontSize',10)

set(Ph{2},'XData',ST(:,1)','YData',ST(:,2)','ZData',ST(:,3)');
set(Ph{3},'XData',PT(:,1)','YData',PT(:,2)','ZData',PT(:,3)');
view(azel+rot)
drawnow
ph11c = plot3([Ac(:,3) Ac(:,4)]', [Ac(:,6) Ac(:,7)]', [Ac(:,9) Ac(:,10)]');
set(ph11c,'LineStyle','-','Color',[.7 .7 .7],'LineWidth',.1);



end
%==================================================%



%####################################################################%
%						UNUSED FUNCTIONS
%####################################################################%

%==================================================%
%		SAVE DATA (PlotCsave)
%==================================================%
%		ADD STANDARD IJK ORTHONORMAL AXIS
%==================================================%
%		3D VECTOR LENGTH FUNCTION
%==================================================%
%		3D VECTOR LENGTH FUNCTION
%--------------------------------------------------%
%{


%---------------------------------------------%
%	ADD STANDARD IJK ORTHONORMAL AXIS
%---------------------------------------------%
function [varargout] = addaxis(varargin)

if nargin == 1
	axvals = varargin{1};
XAxO = [axvals{1}(1) 0 0]; XAxT = [axvals{1}(2) 0 0];
YAxO = [0 axvals{1}(3) 0]; YAxT = [0 axvals{1}(4) 0];
ZAxO = [0 0 axvals{1}(5)]; ZAxT = [0 0 axvals{1}(6)];




else
XAxO = [-15 0 0]; XAxT = [15 0 0];
YAxO = [0 -15 0]; YAxT = [0 15 0];
ZAxO = [0 0 -5]; ZAxT = [0 0 10];
end
[XAx YAx ZAx] = plot3prep({XAxT YAxT ZAxT},{XAxO YAxO ZAxO});
%-----------------
hAax = plot3(XAx,YAx,ZAx,'k','LineWidth',1,'MarkerSize',10);
set(hAax,{'LineWidth'},{1,1,1}'); hold on
varargout={hAax};
end
%---------------------------------------------%

%==================================================%
%	3D VECTOR LENGTH FUNCTION
%--------------------------------------------------%
function [L] = LV(V)

L = sqrt(sum(V.^2));

end
%==================================================%


%==================================================%
%	3D VECTOR LENGTH FUNCTION
%--------------------------------------------------%
function [varargout] = FigSetup2(varargin)

scsz = get(0,'ScreenSize');


if nargin == 1 
	Fnum=varargin{1};
	pos = scsz./[2.5e-3 2.5e-3 1.5 2];
elseif nargin >= 2
	Fnum=varargin{1};
	pos=scsz./varargin{2};
else
	Fnum=1;
	pos = scsz./[2.5e-3 2.5e-3 1.5 2];
end

Fh = figure(Fnum);
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])



if nargin >= 3
	fpos=varargin{3};
	set(gca,'Position',fpos)
end

varargout = {Fh};

end
%--------------------------------------------------%



%==================================================%
%			REPORT DATA (Nenex)
%--------------------------------------------------%
function [] = Nenex(NSteps,stepN,S0,S,CSen,CSex,NoffN,Soffs,Sons,doNoff)



OnPerStep = sum(CSex)/(stepN);
OffPerStep = sum(CSen)/(stepN);

OffStepRate = 1 / OffPerStep;
Offmins = OffStepRate/60;

ClusMins = stepN/60;

ClusMins / Offmins;

Sstart = sum(sum(S0));
Send = sum(sum(S));
Spct = (Send / Sstart)*100;

PperMin = 1/Offmins;
PperHr = PperMin*60;

PctClusPHr = PperHr / Sstart * 100;
NumClusPHr = PctClusPHr/100*Sstart;





disp([' The Cluster to Particle lifetime ratio (in minutes): '])
disp([' (cluster) ' num2str(ClusMins) ' : ' num2str(Offmins) ' (particle)'  ])

disp([' ' ])
disp(['with...' ])
disp([' ' num2str(OnPerStep) ' on-events per step (on average)' ])
disp([' ' num2str(OffPerStep) ' off-events per step (on average)' ])

disp(['if step = 1 second...' ])
disp([' ' num2str(OffStepRate) ' seconds between off events (on average)' ])
disp([' ' num2str(OffStepRate) ' seconds = ' num2str(Offmins) ' minutes'])

disp([' ' ])
disp(['The starting cluster size was ' num2str(Sstart) ' particles,' ])
disp(['with ' num2str(NumClusPHr) ' particles dissociating per hour.' ])
% disp([' equivalent to ' num2str(PctClusPHr) '% of the starting cluster.'])
disp(['The ending cluster size was ' num2str(Send) ' particles, '])
disp([' which is ' num2str(Spct) '%. of the starting cluster size.'])

%-------------------------------------------%
OffPerHr = OffPerStep * 3600;
OffPerHr15 = OffPerHr / 15;
OffPerHr5 = OffPerHr / 5;

FdecPctA = (OffPerHr15/Sstart*100);
FdecPctB = Send/Sstart;

Smin = 0;
Smax = 2;
Fmin = 0;
Fmax = 100-FdecPctA;
Xval = FdecPctB;
PctX = Fmin*(1-(Xval-Smin)/(Smax-Smin)) + Fmax*((Xval-Smin)/(Smax-Smin));
Pval = Fmax-PctX;
FdecPctC = FdecPctA + Pval/5;



disp([' ' ])
disp(['In 1 hour ' num2str(OffPerHr) ' saps dissociated from the lattice. Of those...' ])
disp(['  adjusting for rebinding rate (5:1 to 15:1), '...
	num2str(OffPerHr5) ' to ' num2str(OffPerHr15) ' saps exited the PSD' ])


if (100-FdecPctC) > 0
disp([' ' ])
disp(['linear scaled fluorescence change (approximation):' ])
disp([' flourescence decrease: ' num2str(FdecPctC) '%' ])
disp([' flourescence remaining: ' num2str(100-FdecPctC) '%' ])
else
disp([' estimated flourescence remaining: <1%'])
end
%-------------------------------------------%

% Equation to linear transform range [a b] to range [c d] and solve for [x]:
%				c*(1-(x-a)/(b-a)) + d*((x-a)/(b-a))
% example:
% a=-200
% b=200
% c=0
% d=1
% x=50
% c*(1-(x-a)/(b-a)) + d*((x-a)/(b-a))

%-------------------------------------------%
% Noffs,Soffs,Sons

Soffon = (Sons+Soffs)./2;

%----------------------------%
% figure(3);
% subplot(2,1,1),imagesc(Soffs);
% colormap('bone')
% title('OFF Activity Map');
% colorbar
% %----------------------------%
% figure(3);
% subplot(2,1,2),imagesc(Sons);
% colormap('bone')
% title('ON Activity Map');
% colorbar
%----------------------------%

figure(4);
imagesc(Soffon);
colormap('bone')
axis off
hT = title('  Activity Map (On/Off Events) ');
set(hT,'FontName','Arial','FontSize',16);
colorbar


%----------------------------------------------------------------------%
% FONTS
%{
fig21 = figure(21);
set(21,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig21,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%------

hTXT = text(.1,.9,0,'Century Gothic Font Name');
set(hTXT,'FontName','Century Gothic','FontSize',18);

hTXT = text(.1,.7,0,'Arial Font Name');
set(hTXT,'FontName','Arial','FontSize',18);

hTXT = text(.1,.6,0,'Microsoft Sans Serif Font Name');
set(hTXT,'FontName','Microsoft Sans Serif','FontSize',18);

hTXT = text(.1,.4,0,'Times New Roman Font Name');
set(hTXT,'FontName','Times New Roman','FontSize',18);
%}
%----------------------------------------------------------------------%



% xlabel('Average On/Off Rate')
% labels = {'High Activity','Low Activity'};
% lcolorbar(labels,'fontweight','bold');
%----------------------------%
if doNoff
Noffs = sum(NoffN);
No123 = sum(Noffs(1:3));
No4 = sum(Noffs(4));
N1234pct = (No123)/(No123+No4)*100;
disp([' Of the off events, ' num2str(N1234pct) '% was along an edge'])
end
%-------------------------------------------%


end



%==================================================%
%			FLUORESCENT PLOT (FluorPlot)
%--------------------------------------------------%
function [] = FluorPlot(NSoff,NNoff,SN0,stepN,ph10,ph11,Rev)
%=================================%
%			LIVE FRAP
%=================================%

% FRAP = (SN0-(sum(NSoff)))/SN0;
FRAPL = (SN0-NNoff)./SN0;
SoCSN = (SN0-NNoff./Rev)./SN0;

set(ph10,'XData',[1:stepN],'YData',[FRAPL]);
drawnow; hold on;
set(ph11,'XData',[1:stepN],'YData',[SoCSN]);
drawnow; hold on;

% if stepN==3900;keyboard;end;
%=================================%
end



%==================================================%
%				FIGURE SETUP FUNCTION
%--------------------------------------------------%
function [varargout] = FigSetup(varargin)

scsz = get(0,'ScreenSize');


if nargin == 1 
	Fnum=varargin{1};
	pos = scsz./[2.5e-3 2.5e-3 1.5 2];
elseif nargin == 2
	Fnum=varargin{1};
	pos=scsz./varargin{2};
else
	Fnum=1;
	pos = scsz./[2.5e-3 2.5e-3 1.5 2];
end

Fh = figure(Fnum);
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])
varargout = {Fh};

end


%==================================================%
%					LIVE PLOT
%--------------------------------------------------%
function [varargout] = LivePlot(Fh,nT,Actin,inPSD,varargin)

% keyboard

scsz = get(0,'ScreenSize');
basepos = scsz./[2.5e-3 2.5e-3 1.5 2];
baseazel = [-32 12];
baserot=[0 0];

if nargin >= 2
	rot=varargin{1};
	azel=varargin{2};
elseif nargin == 1 
	rot=varargin{1};
	azel = baseazel;
else
	rot=baserot;
	azel = baseazel;
end



%--------------------
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%--------------------
figure(Fh)
subplot('Position',[.08 .15 .45 .70]), 
ph11c = plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]');
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
view(azel)
grid off
set(gca,'Color',[1,1,1])
hold on;
ph11a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph11b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis([-500 500 -500 500 0 1000])
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
subplot('Position',[.6 .55 .38 .38]), 
ph12a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph12b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis([-500 500 -500 500 0 1000])
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
SPYnXY = dims(1);
SPYhZN = dims(2);
SPYhZS = dims(3);
SPYhXY = dims(4);
SPYhXY = dims(5);
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
PSDactMx = zeros(SPYhXY+100,SPYhXY+100);
for mxp = 1:numel(PSDXY(:,1))
PSDactMx(PSDXY(mxp,2)+SPYhXY+10, PSDXY(mxp,1)+SPYhXY+10) = 1;
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


%==================================================%
%					LIVE PLOT
%--------------------------------------------------%
function [varargout] = TipMatrix(Fh,nT,Actin,dims,AcMx,SPYH,varargin)


%{
% keyboard

scsz = get(0,'ScreenSize');
basepos = scsz./[2.5e-3 2.5e-3 1.5 2];
baseazel = [-32 12];
baserot=[0 0];

if nargin >= 2
	rot=varargin{1};
	azel=varargin{2};
elseif nargin == 1 
	rot=varargin{1};
	azel = baseazel;
else
	rot=baserot;
	azel = baseazel;
end






%--------------------
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%--------------------
figure(Fh)
subplot('Position',[.08 .15 .45 .70]), 
ph11c = plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]');
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
view(azel)
grid off
set(gca,'Color',[1,1,1])
hold on;
ph11a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph11b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis([-500 500 -500 500 0 1000])
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
subplot('Position',[.6 .55 .38 .38]), 
ph12a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph12b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis([-500 500 -500 500 0 1000])
view([0 90])
grid off
set(gca,'Color',[1,1,1])
set(ph12a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph12b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
%--------------------
set(gca,'XTickLabel', sprintf('%.1f|',nT),'FontSize',10)
hold off;
%--------------------
%}


%==================================================%
% Spine Dimensions
% SPYnXY = dims(1);
% SPYhZN = dims(2);
% SPYhZS = dims(3);
% PSDproxy = dims(6);
HX = dims(4);
HY = dims(5);

inPSD = dims(7);
%-----------------------------------%
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
PSDXY = round([PSDTips(:,1) PSDTips(:,2)]);
PSDY = PSDXY(:,2)+HY;
PSDX = PSDXY(:,1)+HX;
SPYHX = SPYH(1)*2;SPYHY = SPYH(2)*2;
% PSDXok = (PSDX<SPYHX);
% PSDYok = (PSDY<SPYHY);
% PSDXk = (PSDX .* PSDXok)+1;
% PSDYk = (PSDY .* PSDYok)+1;

PSDX(PSDX>SPYHX) = SPYHX-1;
PSDY(PSDY>SPYHY) = SPYHY-1;
PSDX(PSDX<1) = 1;
PSDY(PSDY<1) = 1;


% mxpv = 10;
for mxp = 1:numel(PSDXY(:,1))
% PSDactMx(PSDXY(mxp,2)+SPYhXY+mxpv, PSDXY(mxp,1)+SPYhXY+mxpv) = 1;
AcMx(PSDY(mxp), PSDX(mxp)) = 1;
end
ActMask=[1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1];
ActMx = convn(AcMx,ActMask,'same');
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
scatter(PSDX,PSDY, 'r')
% scatter(PSDXY(:,1),PSDXY(:,2), 'r')
hold off
%-----------------------------------%


%==================================================%
% if nT >20000; keyboard; end


varargout = {ActMx};
end
%=============================================================================%




%==================================================%
%			SAVE DATA (PlotCsave)
%--------------------------------------------------%
function [] = PlotCsave(Csave,datarate,dT,stepN,NSteps)
%-------------------------------------------%
	figure(2); hold on;
	plot(Csave)
	set(get(gca,'XLabel'),'String','Steps')
	set(get(gca,'YLabel'),'String','SAP')
	set(gca,'YLim',[0 (max(Csave)+10)])
	xt = (get(gca,'XTick'))*datarate;
	set(gca,'XTickLabel', sprintf('%.0f|',xt));
	SAPtitle = title(['    Cluster Size Over Time '...
	'    Total Steps: ' num2str(NSteps)...
	'  \bullet  Time-Step: ' num2str(dT)]);
	set(SAPtitle, 'FontSize', 12);
	
	
disp(['This cluster lasted for ' num2str(stepN) ' steps'])
disp(['of the ' num2str(NSteps) ' requested steps'])
disp(' ')

end
%=============================================================================%




%}
%####################################################################%







