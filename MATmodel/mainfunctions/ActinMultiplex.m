function [varargout] = ActinMultiplex(LBR,TIME,SIZE,DOES,REVA,AMX,MSK,AMS)
clc, close all; scsz = get(0,'ScreenSize');


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
[ActinTips AxLP] = ActinMainStage2(ActSteps,AMX,AMS);

assignin('base', 'ATs', ActinTips)
ATs = evalin('base', 'ATs');
assignin('base', 'Ax', AxLP)
Ax = evalin('base', 'Ax');

assignin('base', 'AMx', AMx)
AMx = evalin('base', 'AMx');

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

BTs = [];
Nsteps = ActSteps;
dT = AMX{47};

LivePlotMod = AMS{1};
SaveTipsAfter = AMX{19};
SaveTipsRate = AMX{20};

doActCounts = AMX{35};

doDelOrig = AMX{39};
doDelOrigT = AMX{40};







%---------------------------------------------%
%				ANGLE SETUP
%---------------------------------------------%
% unitsratio('rad', 'deg')
d2r = 1/(180/pi);
Ov = [2.0 29.7 57.4 85.1 112.8 140.5 168.2 195.8 223.5 251.2 278.9 306.6 334.3];
Ov = Ov * d2r;

%----------------------------------------
% Spine Dimensions
%----------------------------------------
%GaSz = 1;						% Actin size (5/2)
%Ga2RR = 5.1/2/GaSz;			% Ratio of sim to real

SPYnXY = AMX{10};	% Spy neck XY
SPYhZN = AMX{11};	% Spy head north
SPYhZS = AMX{12};	% Spy head south
SPYhX = AMX{13};	% Spy head X
SPYhY = AMX{14};	% Spy head Y

PSDproxy = AMX{15};
PSDproxyXY = AMX{27};
inPSD = SPYhZN - PSDproxy;

SPYH = [SPYhX SPYhY];
AcMx = zeros(SPYhY*2,SPYhX*2);

dims = [SPYnXY SPYhZN SPYhZS SPYhX SPYhY PSDproxy inPSD PSDproxyXY];
%----------------------------------------


%		INITIALIZE AND TAG STARTING FILAMENTS
%-----------------------------------------------------------%
NStFils = AMX{42};      Actin = zeros(NStFils,18);          
%   N  Xa  Xo  Xt  Ya  Yo  Yt  Za  Zo  Zt  MomID  ID  Fkd Born Died Lif MaxL MeanL sdL
%   1  2   3   4   5   6   7   8   9   10  11     12  13  14   15   16  17   18    19
%-----------------------------------------------------------%
TagN = numel(Actin(:,1));
Actin(:,12) = 1:TagN;
TagN = TagN+1;
%--------------------------------------------------------%

% Starting Length Loc & Angles
StartMonos = AMX{41};
fXYo = AMX{45};
fZo = AMX{43};
fXYa = AMX{46};
fZa = AMX{44};

% Branching Angles
TPi = d2r*AMX{9};
PPi = 0;


% FIL-1
%--
Actin(1,1) = StartMonos;
Actin(1,2) = 0;				% Xa
Actin(1,5) = 0;				% Ya
Actin(1,8) = 0;				% Za
Actin(1,3) = 0;				% Xo
Actin(1,6) = 0;				% Yo
Actin(1,9) = 0;				% Zo
%--

% FIL-2
%--
Actin(2,1) = StartMonos;
Actin(2,2) = -45;				% Xa
Actin(2,5) = 45;				% Ya
Actin(2,8) = d2r*-fZa;		% Za
Actin(2,3) = fXYo;			% Xo
Actin(2,6) = fXYo;			% Yo
Actin(2,9) = 0;				% Zo
%--

% FIL-3
%--
Actin(3,1) = StartMonos;
Actin(3,2) = 45;			% Xa
Actin(3,5) = -45;				% Ya
Actin(3,8) = d2r*fZa;		% Za
Actin(3,3) = fXYo;			% Xo
Actin(3,6) = -fXYo;			% Yo
Actin(3,9) = 0;				% Zo
%--

% FIL-4
%--
Actin(4,1) = StartMonos;
Actin(4,2) = 45;				% Xa
Actin(4,5) = 0;				% Ya
Actin(4,8) = d2r*-fZa;		% Za
Actin(4,3) = -fXYo;			% Xo
Actin(4,6) = fXYo;			% Yo
Actin(4,9) = 0;				% Zo
%--

% FIL-5
%--
Actin(5,1) = StartMonos;
Actin(5,2) = 0;				% Xa
Actin(5,5) = 45;				% Ya
Actin(5,8) = d2r*fZa;		% Za
Actin(5,3) = -fXYo;			% Xo
Actin(5,6) = -fXYo;			% Yo
Actin(5,9) = 0;				% Zo
%--




% 6th-9th Starting Filaments
Actin(6:9,1)=round(SPYhZN-SPYhZS-2);

% FIL-6
%--
Actin(6,2) = 225*d2r;		% Xa
Actin(6,5) = 0*d2r;			% Ya
Actin(6,8) = 45*d2r;		% Za
Actin(6,3) = (SPYhX/2)-10;	% Xo
Actin(6,6) = (SPYhY/2)-10;	% Yo
Actin(6,9) = SPYhZS+1;		% Zo

% FIL-7
%--
Actin(7,2) = 45*d2r;		% Xa
Actin(7,5) = 0*d2r;			% Ya
Actin(7,8) = 45*d2r;		% Za
Actin(7,3) = 10-(SPYhX/2);	% Xo
Actin(7,6) = 10-(SPYhY/2);	% Yo
Actin(7,9) = SPYhZS+1;		% Zo


% % FIL-8
% %--
% Actin(8,2) = 315*d2r;		% Xa
% Actin(8,5) = 0*d2r;			% Ya
% Actin(8,8) = 45*d2r;		% Za
% Actin(8,3) = 10-(SPYhX/2);	% Xo
% Actin(8,6) = (SPYhY/2)-10;	% Yo
% Actin(8,9) = SPYhZS+1;		% Zo

% FIL-9
%--
% Actin(7,2) = 135*d2r;		% Xa
% Actin(7,5) = 0*d2r;			% Ya
% Actin(7,8) = 45*d2r;		% Za
% Actin(7,3) = (SPYhX/2)-10;	% Xo
% Actin(7,6) = 10-(SPYhY/2);	% Yo
% Actin(7,9) = SPYhZS+1;		% Zo



%----------------------------------------
% TRIG: branch XYZ tip coordinates
Actin(:,4) = Actin(:,1) .* sin(Actin(:,8)) .* cos(Actin(:,2)) + Actin(:,3);
Actin(:,7) = Actin(:,1) .* sin(Actin(:,8)) .* sin(Actin(:,2)) + Actin(:,6);
Actin(:,10) = Actin(:,1) .* cos(Actin(:,8)) + Actin(:,9);
%----------------------------------------

ACTs = Actin;
oActin = Actin;



%----------------------------------------
% SPINE VOLUME
% A typical spine volume is <= 0.1 um^3 or 1e-16 L

Vneck = pi * SPYnXY^2 * SPYhZS;
Vhead = pi * SPYhX^2 * (SPYhZN-SPYhZS);
SpyV = (Vneck+Vhead) * 1e-22;	%dim inputs are in .1um^3 units, conver to L


% MATH - Spine Volume
%{

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

% Actin Polymerization (12 N/然*s)
% Act_Na = 1e5;
% Act_PRnT = 1000;

% Actin Depolymerization (2 N/s)	
Act_DRnT = 200;
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

Act_mM = Act_N/SpyV * (1/6e23) * (1000/1);	% 1.6e-10 

% Act_mM = Act_N / SpyVu / mol;	% 1.6e-10 
% Act_N = Act_N / SpyVu / mol;
Act_PR = 12 * Act_mM * dT;
Act_DR = 2 * dT;

volume of cylinder
V = pi * r^2 * h

volume of a sphere
V = (4*pi*r^3)/3
%}
%----------------------------------------
% Actin Variables (GActin & FActin)
%----------------------------------------
mM_Act = AMX{16};							% Actin mM
GActinN0 = mM_Act * (1/1000) * SpyV * 6e23;	% t0 N Gactin monomer units (6e4)
Act_mM = GActinN0/SpyV *(1/6e23)*(1000/1);	% Check mM_Act == Act_mM (1.6e-10)

Nfi = numel(Actin(:,1));		% current number of branches
FActinN = sum(Actin(:,1));		% current number FActins
GActinN = GActinN0 - FActinN;	% current number GActins

Act_PR = 12 * Act_mM * dT;		% empirical actin poly rate
Act_DR = 2 * dT;				% empirical actin tip depoly rate
Act_DRo = 8 * dT;				% empirical actin origin depoly rate

%FActinN0 = FActinN;				% t0 number FActins
%Act_PR0 = Act_PR;				% t0 poly rate
%----------------------------------------



%----------------------------------------
% Cofilin Variables
%----------------------------------------
CofR = AMX{18}* dT; % cofilin activity rate
CofN = AMX{30}; % cofilin delete Nunits
CofS = AMX{34}; % cofilin delete Nunits
delOr =AMX{52}* dT; % delete from origin rate
%----------------------------------------



%----------------------------------------
% Arp Variables
%----------------------------------------
ArpR = AMX{17}* dT;			% Arp activity rate
ArpR0 = ArpR;			% Arp activity rate
ArpScalar = AMX{21};	% Arp filament length scalar
ArpAdd = AMX{29};		% Add X units to new branches
ARPmax = AMX{22};		% Maximum Arp branches


GArpN = AMX{33};	% G-Arp free monomer num
FArpN = NStFils;	% F-Arp in filaments num


% MATH - Arp Branching
%{
%ArpN = 1e3;
%ArpOn = 5;
%ArpOff = 1;
%Arp_uM = ArpN / SpyVu / mol;	% 1.6 - 16 uM
%Arp_PR = ArpOn * Arp_uM * dT;
%Arp_DR = ArpOff * dT;
%}
%----------------------------------------



%----------------------------------------
% LTP related variables

doActLTP = AMX{36};
ActLTPstart = AMX{37};
ActLTPend = AMX{38};
ArpRLTP = ArpR*AMX{31};

%----------------------------------------




%-------------------------------------------%
% Linear Range Transformation
%-------------------------------------------%
% Equation to linear transform range
% [a b] to range [c d] and solve for [x]:
% F(x) = c*(1-(x-a)/(b-a)) + d*((x-a)/(b-a))
%-------------------------------------------%
global sja;
global sjb;
global sjc;
global sjd;

sja=AMX{23};
sjb=StartMonos*AMX{24};
sjc=AMX{25};
sjd=AMX{26};

linsc = @(jx) (sjc*(1-(jx-sja)/(sjb-sja)) + sjd*((jx-sja)/(sjb-sja)));

global ska;
global skb;
global skc;
global skd;
ska=AMX{48};
skb=StartMonos*AMX{49};
skc=AMX{50};
skd=AMX{51};

sigsc = @(jx) (1/(1+ exp(-((skc*(1-(jx-ska)/(skb-ska)) + skd*((jx-ska)/(skb-ska)))))));



% miscfun = @(xyz,tta) ([1 0 0; 0 cos(tta) sin(tta); 0 -sin(tta) cos(tta)] * xyz);
%-------------------------------------------%
%{
%-------------------------------------------%

keyboard


close all
clear ArpCurve1
rbg1=.01;rbg2=.99;rbg3=.01;

ska=0;
skb=170;
skc=-7;
skd=7;

sigsc = @(jx) (1/(1+ exp(-((skc*(1-(jx-ska)/(skb-ska)) + skd*((jx-ska)/(skb-ska)))))));


figure(44)
for ArpScalr = 30:10:80

	nn = 1;
	for Flength = 0:1:StartMonos
	%ArpCurve1(nn) = ArpR * sigsc(Flength+ArpScalr);
	ArpCurve1(nn) = ArpR * (1/(exp(1/170*(1190-14*(Flength+ArpScalr)))+1));
	nn=nn+1;
	end

ph=plot(ArpCurve1); hold on;
set(ph,'Color',[rbg1 rbg2 rbg3])
rbg1=rbg1+.08;
rbg2=rbg2-.08;
rbg3=rbg3+.08;
drawnow; pause(.2);
end

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
	for Flength = 0:1:StartMonos
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



%----------------------------------------
% MATH TIDBITS
%{

% 5e8 in cells; 1e4 synapses per neuron; 5e8/1e4 = 5e4

Polymerization Rate
(+)end: .012 N/然*ms
** (12 N/mM*ms) (12 N/然*s)
** thus at 1 然 free ATP-actin, .012 subunits will be added to the (+)end per ms
** at .1 mM free ATP-actin, 1.2 subunits will be added to the (+)end per ms

Depolymerization Rate
(+)end: 1.4 N/s
(-)end: 0.8 N/s
** dissociation is independent of free actin concentration

To find the critical concentration (Cc) for growth we set the two rate 
equations equal to each other:

12.0/然*s = 1.4/s
12/1.4 = 1/然
(+)Cc = .12 然
(-)Cc = .6 然





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
-- (mean spine length: 1.0 痠 or 1000 nm)
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
(+)end: ~12.0 subunits/然*s
(-)end: ~1.3 subunits/然*s
** thus if there is 1 然 of free ATP-Gactin then 12 subunits will be added 
to the (+)end per second and 1.3 subunits will be added to the (-)end every second

Depolymerization Rate
(+)end: ~1.4 subunits/s
(-)end: ~0.8 subunits/s
** dissociation is independent of free Gactin concentration

To find the critical concentration (Cc) for growth we set the two rate 
equations equal to each other:

12.0/然*s = 1.4/s
12/1.4 = 1/然
(+)Cc = .12 然
(-)Cc = .6 然

Thus when the free actin concentration >.12 然 filaments will grow at the (+) end 
and when >.6 然 filaments will grow at the (-) end too.


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


%}
%----------------------------------------


%----------------------------------------
%		STATIC FINAL FIGURES SETUP
%----------------------------------------
Fh1 = FigSetup(33);


%----------------------------------------
%	  ANIMATED REAL-TIME FIGURE SETUP
%----------------------------------------
rot0 = [5 0];
rot1 = rot0;
azel0 = [-32 12];
sz = [2.5e-3 2.5e-3 1.5 1.4];
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

%===============================================================%
%						MAIN OUTER LOOP
for nT = 1:Nsteps
%---------------------------------------------------------------%



	Nfi = numel(Actin(:,1));	% Current number of Filaments (mature >= 1000 filaments)
	if (mod(nT,LivePlotMod)==0); disp(Nfi); end;
	ACTdepoly = 0;
	COFdepoly = 0;
	ARPpoly = 0;
	ACTpoly = 0;
	
	% radial distance to spine shaft membrane
	ActinXYh = sqrt(Actin(:,4).^2 + Actin(:,7).^2);

	Xspi = ActinXYh > SPYnXY;			% Logical array of filaments beyond SPYnXY
	Yspi = ActinXYh > SPYnXY;			% Logical array of filaments beyond SPYnXY
	Zspi = (abs(Actin(:,10)) > SPYhZN);	% Logical array of filaments beyond SPYhZN

	Xhed = ActinXYh >= SPYhX;
	Yhed = ActinXYh >= SPYhY;
	Zhed = (abs(Actin(:,10)) <= SPYhZN) & (abs(Actin(:,10)) >= SPYhZS);
		
	ActMem = [Xspi Yspi Zspi Xhed Yhed Zhed];
	
 	Act0 = (Actin(:,1)>0);	% assures no negative actin values
	Actin(:,1) = Actin(:,1) .* (Actin(:,1)>0);
	Actin(:,1) = Actin(:,1) .* (Actin(:,10) > -20);


	NFact = numel(Actin(:,1));	% current number of actin filaments
	

	%=============================================================================%
	%						MAIN INNER LOOP
	for aN=1:NFact
	%-----------------------------------------------------------------------------%
	
	rv=rand(6,1); % Generate a few random vaules from uniform{0:1}

	
	
		%=================================================================%
		% POLYMERIZATION
		if  (Act_PR > rv(1)) && ...
			(((~sum(ActMem(aN,:)) && Act0(aN)) || Zhed(aN))) && ...
			(~Xhed(aN) && ~Yhed(aN))
		%----------------------------------------
		
			Actin(aN,1) = Actin(aN,1) + 1;
			
			ACTpoly = ACTpoly + ceil(Act_PR);
			
		%----------------------------------------
		end
		%=================================================================%
		
		
		
		
		
		%=================================================================%
		% BRANCHING
		
		%(ArpR * sigsc(Actin(aN,1)+ArpScalar)) > rv(2) && ...
		
		if  (Nfi < ARPmax) && ...
			(ArpR * (1/(exp(1/170*(1190-14*(Actin(aN,1)+ArpScalar)))+1))) > rv(2) && ...
			((~sum(ActMem(aN,:)) && Act0(aN)) || Zhed(aN)) && ...
			(~Xhed(aN) && ~Yhed(aN))
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
			Ct_x = (Rmm) * sin(Actin(aN,8)) * cos(Actin(aN,2)) + Actin(aN,3);
			Ct_y = (Rmm) * sin(Actin(aN,8)) * sin(Actin(aN,2)) + Actin(aN,6);
			Ct_z = (Rmm) * cos(Actin(aN,8)) + Actin(aN,9);

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
			
			Otta = Ov(randi(13,1));		% Random rotational angle (theta)
			
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
			
			% New branch Angles
			Actin(fNf,2) = tPhi2;
			Actin(fNf,8) = tTheta2;
			
			% New branch Origin
			Actin(fNf,3) = Po2(1);
			Actin(fNf,6) = Po2(2);
			Actin(fNf,9) = Po2(3);
			
			% New branch Tip
			Actin(fNf,4) = (ArpAdd) * sin(tTheta2) * cos(tPhi2) + Po2(1);
			Actin(fNf,7) = (ArpAdd) * sin(tTheta2) * sin(tPhi2) + Po2(2);
			Actin(fNf,10) = (ArpAdd) * cos(tTheta2) + Po2(3);
			
			
			%---------------------------------------------%
			ARPpoly = ARPpoly + ArpAdd;
			
		%----------------------------------------
		end
		%=================================================================%


		
		
		
		%=================================================================%
		% ACTIN PASSIVE DEPOLYMERIZATION
		
		if (Act_DR > rv(3))
		%----------------------------------------
			Actin(aN,1) = Actin(aN,1)-1;
			ACTdepoly = ACTdepoly + 1;
		%----------------------------------------
		end
		
		
		
		
		%=================================================================%
		% COFILIN ASSISTED DEPOLYMERIZATION
		
		%if  ( CofR * exp(linsc(Actin(aN,1)-CofS)) ) > rv(4)
		if  ( CofR * exp((Actin(aN,1)-CofS)/200) ) > rv(4)
		%----------
			
			Actin(aN,1) = Actin(aN,1)-CofN;
		
			% Depoly from origin (Tag)
			if rv(5) < delOr
			Actin(aN,13) = 1;
			end
			
			COFdepoly = COFdepoly + CofN;
		%----------
		end
		
		%-----------------------------------
		if Actin(aN,13) && (Act_DRo > rv(6))
		%----------
			Actin(aN,1) = Actin(aN,1)-1;
			ACTdepoly = ACTdepoly + 1;
		%----------
		end
		%=================================================================%

		
		
		
	
	%-----------------------------------------------------------------------------%
	end
	%						MAIN INNER LOOP
	%=============================================================================%
	
	
	%=================================================================%
	% ADJUST RATE VALUES
	%----------------------
	FArpN = Nfi;					% Current number of branch filaments
	FActinN = sum(Actin(:,1));		% Current number FActins
	GActinN = GActinN0 - FActinN;	% Current number GActins

	ArpR = ((GArpN-FArpN)/GArpN)*ArpR0;

	Act_mM = GActinN/SpyV * (1/6e23) * (1000/1);
	Act_PR = 12 * (Act_mM * dT);

	%----------------------
	% LTP ADJUSTMENTS
	%----------------------
	if doActLTP
	if (nT >= ActLTPstart) && (nT <= ActLTPend)
		ArpR = ArpRLTP;
	end
	end
	
	%----------------------
	% DELETE ORIGINAL FILAMENTS
	%----------------------
	if doDelOrig
	if (nT == doDelOrigT)
		Actin(1:NStFils,:) = [];
	end
	end
	%=================================================================%
	

	
	%==================================================%
	%	PROCESS STORE AND REMOVE DEL BRANCHES
	%--------------------------------------------------%
	
	Actin(:,16) = Actin(:,16) + 1; %Lifetime
	
	MaxL = Actin(:,1) > Actin(:,17);
	Actin(MaxL,17) = Actin(MaxL,1); %Longest Ever Length
	
	
	% Mean Length
	Actin(:,18) = (Actin(:,18).*(Actin(:,16)-1) + Actin(:,1))./Actin(:,16);
	
	delFi = (Actin(:,1)<2);
	Actin(delFi,15) = nT;
	ACTs = cat(1,ACTs,Actin(delFi,:));
	Actin(delFi,:) = [];
	%--------------------------------------------------%
	
	
	%==================================================%
	%	USE ANGLES TO COMPUTE TIP LOCATION
	%--------------------------------------------------%
	% MATH - branch XYZ tip coordinates
	Actin(:,4) = Actin(:,1) .* sin(Actin(:,8)) .* cos(Actin(:,2)) + Actin(:,3);
	Actin(:,7) = Actin(:,1) .* sin(Actin(:,8)) .* sin(Actin(:,2)) + Actin(:,6);
	Actin(:,10) = Actin(:,1) .* cos(Actin(:,8)) + Actin(:,9);
	%--------------------------------------------------%
	
	

	
	%==================================================%
	%				LIVE PLOT
	%--------------------------------------------------%
	if mod(nT,LivePlotMod) == 0
				
		LivePlot(Fh2Live,nT,Actin,inPSD,rot1,azel0,AxLims);
		rot1 = rot1 + rot0;
	end
	%--------------------------------------------------%
	
	
	%==================================================%
	%				SAVE TipMatrix
	%--------------------------------------------------%
	if nT >SaveTipsAfter;
	if mod(nT,SaveTipsRate) == 0
		ActMx = TipMatrix(Fh2Live,nT,Actin,dims,AcMx,SPYH,rot1,azel0);
		BTs{numel(BTs)+1} = ActMx;
	end
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
		
		NumArpR(nT) = ArpR;
		NumAct_mM(nT) = Act_mM;
		NumAct_PR(nT) = Act_PR;
		
	end
	%--------------------------------------------------%
	
%if nT == 30000; keyboard; end;
%{
% ActData(nT,:) = [Act_mM, Act_PR, Act_N];
% DePData(nT,:) = DePSum;
% DePSum = 0;
% if nT == 125; keyboard; end;
%}
%-----------------------------------------------------------------------------%
end
%						MAIN OUTER LOOP
%=============================================================================%

if doProfile
profile viewer
end

% Lifetime of remaining filaments
Actin(:,15) = (nT - Actin(:,14))*2;
Actin(:,16) = (Actin(:,15) - Actin(:,14));
ACTs = cat(1,ACTs,Actin);

remove0 = (ACTs(:,16)<2);
ACTs(remove0,:) = [];



%==================================================%
%		Counter Crunch and Plot
%--------------------------------------------------%
if doActCounts
%%
%----------------------------------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6]; c5= [.01 .9 .01];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7]; c55=[.01 .9 .01];
%----------------------------------------------------------------------%
%The plot shows the probability for each nonnegative integer when ? = 5.
%------------------------------------------%
% Filament Lifetime CDF
%------------------------------------------%
Filament_Lifetime = ACTs(:,16) .* dT;
muFlif = mean(Filament_Lifetime);
Flast = max(Filament_Lifetime);
[coECDF_Lif, coX_Lif] = ecdf(Filament_Lifetime);
XLmax = 3.6e3;
%------
	cdfexit = figure(69);
	set(cdfexit,'Color',[1,1,1])
	set(cdfexit,'OuterPosition',[400,100,600,700])
axes('Position',[.08 .48 .86 .48]);
PmeanFlif = cdfplot(Filament_Lifetime);
	set(gca,'XLim', [0 XLmax],'YLim', [0 1.1])
	axT = axis;
	xt = (get(gca,'XTick'))./60;
	set(gca,'XTickLabel', sprintf('%.0f|',xt))
	text((axT(2)/2.5),.05,'Minutes Existed','FontSize',12)
	set(PmeanFlif,'color',[1 0 1])
	rectangle('Position',[Flast,1,nT,.001],'EdgeColor',[1 0 1],'LineStyle','--')
	CDFtitle = title([...
	'    Empirical CDF of Filament lifetime (mean lifetime was:' int2str(muFlif/60) ' min)']);
	set(CDFtitle,'FontSize',12);
	set(PmeanFlif,'LineStyle','-','Color', [.9 .2 .2],'LineWidth',2);
axes('Position',[.12 .06 .78 .32]);
ecdfhist(coECDF_Lif, coX_Lif, 50)
	axis([0 XLmax/2 0 100]); axis 'auto y'; title('Lifetime ECDF Histogram');
	xt = (get(gca,'XTick'))./60; set(gca,'XTickLabel', sprintf('%.0f|',xt))
	xlabel('minutes'); ylabel('CDF hist  (range 0-1) ');
%------------------------------------------%
Filament_Lifetimes = sort(Filament_Lifetime);
%dfittool(Filament_Lifetimes)

% Prepare figure
F70 = figure(70);
set(F70,'OuterPosition',[400,100,600,700])
set(gcf,'Color',[1,1,1])
LegHandles = []; LegText = {};
% --- Plot data originally in dataset "meanlife data"
axes('Position',[.08 .48 .86 .48]);
pause(.5)
[CdfY,CdfX,CdfLower,CdfUpper] = ecdf(Filament_Lifetimes,...
	'Function','survivor','alpha',1e-10);  % compute empirical function
hLine = stairs(CdfX,CdfY,'Color',[.3 0 .6],'LineStyle','-', 'LineWidth',1.5);
hold on
[StairsXlower,StairsYlower] = stairs(CdfX,CdfLower);
hold on
[StairsXupper,StairsYupper] = stairs(CdfX,CdfUpper);
hold on
hBounds = plot([StairsXlower(:); NaN; StairsXupper(:)], [StairsYlower(:); NaN; StairsYupper(:)],...
	'Color',[.9 .2 .2],'LineStyle','-.', 'LineWidth',1);
hold on

xlabel('CDF Survivor (seconds)');
ylabel('Survivor function')
LegHandles(end+1) = hLine;
LegText{end+1} = 'Filament Lifetime';
LegHandles(end+1) = hBounds;
LegText{end+1} = '99.99% confidence bounds';
set(gca,'XLim',[-100 3600]);
box on;
grid on;
hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'Location', 'NorthEast');
set(hLegend,'Interpreter','none');

% Ratio of all Tags generated to current number of Tags
AcTags = Actin(:,12);
AT2TNratio = numel(AcTags)/TagN;
%---
axes('Position',[.12 .08 .78 .30]);
[ph1] = plot(AcTags,'b');
leg1 = legend(ph1,'Arp Rate');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
set(legend,'Location','SouthEast');
%------------------------------------------%
MS1 = 5;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
hTitle  = title('Arp Branching Rate');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Branching Events');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set(hTitle,'FontSize',12);
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%------------------------------------------%
% NOTES
%{
meanFlif = ACTs(:,16) .* dT;
muFlif = mean(meanFlif);
Flast = max(meanFlif);
[coECDF_Lif, coX_Lif] = ecdf(meanFlif);
XLmax = 2e3;
%------
	cdfexit = figure(69);
	set(cdfexit,'OuterPosition',[400,100,700,800])
axes('Position',[.06 .55 .90 .40]);
PmeanFlif = cdfplot(meanFlif);
	set(gca,'XLim', [0 XLmax],'YLim', [0 1.1])
	axT = axis;
	xt = (get(gca,'XTick'))./60;
	set(gca,'XTickLabel', sprintf('%.0f|',xt))
	text((axT(2)/2.1),.05,'Minutes Existed','FontSize',12)
	set(PmeanFlif,'color',[1 0 1])
	rectangle('Position',[Flast,1,nT,.001],'EdgeColor',[1 0 1],'LineStyle','--')
	CDFtitle = title([...
	'    Empirical CDF of Filament lifetime (mean lifetime was:' int2str(muFlif/60) ' min)']);
	set( CDFtitle,'FontSize',12,'FontWeight','bold');
axes('Position',[.06 .30 .40 .20]);
plot(sort(meanFlif),1:numel(meanFlif));
	axis([0 XLmax 0 100]); axis 'auto y'; title('Filament Lifetimes');
	xt = (get(gca,'XTick'))./60; set(gca,'XTickLabel', sprintf('%.0f|',xt))
axes('Position',[.55 .30 .40 .20]);
ecdfhist(coECDF_Lif, coX_Lif)
	axis([0 XLmax 0 100]); axis 'auto y'; title('Lifetime ECDF Histogram');
	xt = (get(gca,'XTick'))./60; set(gca,'XTickLabel', sprintf('%.0f|',xt))
axes('Position',[.06 .04 .40 .20]);
plot(sort(meanFlif),1:numel(meanFlif));
	axis([0 XLmax 0 100]); axis 'auto y'; title('Filament Lifetimes');
	xt = (get(gca,'XTick'))./60; set(gca,'XTickLabel', sprintf('%.0f|',xt))
axes('Position',[.55 .04 .40 .20]);
hist(meanFlif);%ecdfhist(coX_G1, coX_G1)
	axis([0 XLmax 0 100]); axis 'auto y'; title('Lifetime CDF Histogram');
	xt = (get(gca,'XTick'))./60; set(gca,'XTickLabel', sprintf('%.0f|',xt))
%----------------------------

% FilTurnover is the ratio between the average number of
% filaments throughout the simulation, and the number of
% filaments that were completely eliminated in the
% course of the simulation
% if this ratio is less than 1, there is a possibility that
% all original filaments were eliminated.
% The smaller this ratio, the greater the likelihood
%}

SNumFils = sum(NumFils) / nT;
SNumdelFi = sum(NumdelFi);
FilTurnover = SNumFils / SNumdelFi;

CompressN = 100;
subsCOFACT = floor(linspace(1,CompressN+1,numel(NumCOFACT)));
subsCOFACT(numel(NumCOFACT)) = CompressN;
Ncomp = numel(NumCOFACT) / CompressN;
AcArCOFACT = accumarray(subsCOFACT',NumCOFACT) ./ Ncomp;
CumSumCOFACT = cumsum(AcArCOFACT);

subsARPACT = floor(linspace(1,CompressN+1,numel(NumARPACT)));
subsARPACT(numel(NumARPACT)) = CompressN;
AcArARPACT = accumarray(subsARPACT',NumARPACT) ./ Ncomp;
CumSumARPACT = cumsum(AcArARPACT);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					FINAL OUTPUT (SPLINE) FIGURE 1 OF 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------------------%
fig21 = figure(21);
set(21,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig21,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6]; c5= [.01 .9 .01];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7]; c55=[.01 .9 .01];
pause(1);
%===========================================================%
% FIG1 TOP LEFT: Poly & Depoly Events
%===========================================================%
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
MS1 = 5; MS2 = 2;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
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
%===========================================================%


%===========================================================%
% FIG1 TOP RIGHT: Monomeric vs Filamentous Actin
%===========================================================%
sbpos = [.55 .57 .38 .38];
subplot('Position',sbpos);
[ph1] = plot(NumGAct,'b');
leg1 = legend(ph1,'G-actin');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
subplot('Position',sbpos);
[ph2] = plot(NumFAct,'r');
legend([OUTH;ph2],OUTM{:},'F-actin');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%------------------------------------------%
xt = (get(gca,'XTick')).* dT./60;
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','SouthEast');
%------------------------------------------%
MS1 = 5; MS2 = 2;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
hTitle  = title('Monomeric vs Filamentous Actin');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Actin Molecules');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set(hTitle,'FontSize',12);
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%===========================================================%



%===========================================================%
% FIG1 BOTTOM LEFT: Actin Concentration (mM)
%===========================================================%
sbpos = [.07 .09 .4 .35];
subplot('Position',sbpos);
[ph1] = plot(NumAct_mM,'b');
leg1 = legend(ph1,'Conc. mM');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
subplot('Position',sbpos);
[ph2] = plot(NumAct_PR,'r');
legend([OUTH;ph2],OUTM{:},'Poly Rate (12*mM*dT)');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%------------------------------------------%
xt = (get(gca,'XTick')).* dT./60;
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthEast');
%------------------------------------------%
MS1 = 5; MS2 = 2;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
hTitle  = title('Actin Polymerization');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Actin Concentration (mM)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set(hTitle,'FontSize',12);
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%===========================================================%




%===========================================================%
% FIG1 BOTTOM RIGHT: Branching Events
%===========================================================%
%----------------------------
% Dual Axis Plot
%----------------------------
sbpos = [.55 .09 .38 .35];
subplot('Position',sbpos);
[ph1] = plot(NumArpR,'r');
leg1 = legend(ph1,'Arp Rate');
[LEGH,OBJH,OUTH,OUTM] = legend;
set(legend,'Location','SouthEast');
haxes1 = gca; % handle to axes
%set(haxes1,'XColor','r','YColor','r')
hold on
haxes1_pos = get(haxes1,'Position'); % store position of first axes
haxes2 = axes('Position',haxes1_pos,...
              'XAxisLocation','top',...
              'YAxisLocation','right',...
              'Color','none');
%subplot('Position',sbpos);
[ph2] = line(1:numel(AcTags),AcTags,'Parent',haxes2,'Color','k');
LGh1 = legend([OUTH;ph2],OUTM{:},'Actin Tags');
[LEGH,OBJH,OUTH,OUTM] = legend;
set(LGh1,'Location','SouthEast');
hold on
%------------------------------------------%
xt = (get(haxes1,'XTick')).* dT./60;
set(haxes1,'XTickLabel', sprintf('%.0f|',xt))
%------------------------------------------%
MS1 = 5; MS2 = 2;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
%hTitle  = title('Arp Branching Rate');
%hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Branching Events');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set(hTitle,'FontSize',12);
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
%===========================================================%




%{
% %===========================================================%
% % FIG3 TOP RIGHT: Synapse Particle Counts
% %===========================================================%
% sbpos = [.55 .57 .4 .38]; ptype = 4;
% cOLOR = [c1; c2; c3; c4];
% cOLOR = repmat(cOLOR,15,1);
% %------------------------------------------%
% itemN = 16; 
% [ph1 hax1 mu1] = CIplot(reDATAGluRdata,applered,5,sbpos,itemN,1);
% legend(ph1,'G2-PSA1');
% [LEGH,OBJH,OUTH,OUTM] = legend;
% hold on
% %ph1 = plot(mu1);
% %hold on
% %---
% itemN = 17; 
% [ph2 hax2 mu2] = CIplot(reDATAGluRdata,oceanblue,5,sbpos,itemN,1);
% legend([OUTH;ph2],OUTM{:},'G2-PSA2');
% [LEGH,OBJH,OUTH,OUTM] = legend;
% hold on
% %ph2 = plot(mu2);
% %hold on
% %---
% itemN = 18; 
% [ph3 hax3 mu3] = CIplot(reDATAGluRdata,neongreen,5,sbpos,itemN,1);
% legend([OUTH;ph3],OUTM{:},'G2-PSD1');
% [LEGH,OBJH,OUTH,OUTM] = legend;
% hold on
% %ph3 = plot(mu3);
% %hold on
% %---
% itemN = 19; 
% [ph4 hax4 mu4] = CIplot(reDATAGluRdata,hotpink,5,sbpos,itemN,1);
% legend([OUTH;ph4],OUTM{:},'G2-PSD2');
% [LEGH,OBJH,OUTH,OUTM] = legend;
% hold on
% %ph4 = plot(mu4);
% %hold on
% %------------------------------------------%
% xt = (get(gca,'XTick'))*100;
% set(gca,'XTickLabel', sprintf('%.0f|',xt))
% %------------------------------------------%
% MS1 = 5; LnW=1;
% set(ph1,'LineStyle','-','Color',applered,'LineWidth',LnW,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',applered,'MarkerFaceColor',applered);
% set(ph2,'LineStyle','-','Color',oceanblue,'LineWidth',LnW,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',oceanblue,'MarkerFaceColor',oceanblue);
% set(ph3,'LineStyle',':','Color',neongreen,'LineWidth',LnW,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',neongreen,'MarkerFaceColor',neongreen);
% set(ph4,'LineStyle',':','Color',hotpink,'LineWidth',LnW,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',hotpink,'MarkerFaceColor',hotpink);
% 
% hTitle  = title('Distribution of Particles with Brownian Motion');
% hXLabel = xlabel('Time');
% hYLabel = ylabel('Particles (+/- SEM)');
% set(gca,'FontName','Helvetica');
% set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
% set([hXLabel, hYLabel],'FontSize',10);
% set( hTitle,'FontSize',12);
% set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
% 'XMinorTick','on','YMinorTick','on','YGrid','on', ...
% 'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
% haxes=axis;
% %ylim([0 haxes(4)*1.2 ]);
% % xlim([0 (haxes(2)*.9)]);
% %======================================================================%
%%
%}

%%
%------------
end
%--------------------------------------------------%
keyboard


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
PSDactMx = zeros(SPYhY+100,SPYhX+100);
for mxp = 1:numel(PSDXY(:,1))
PSDactMx(PSDXY(mxp,2)+SPYhY+10, PSDXY(mxp,1)+SPYhX+10) = 1;
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
plot(ActData(Act_PRnT:nT,1)); title('Act uM');
subplot('Position',[.35 .05 .28 .90]),
plot(ActData(Act_PRnT:nT,2)); title('Act PR');
subplot('Position',[.68 .05 .28 .90]),
plot(ActData(Act_PRnT:nT,3)); title('Act N');

figure
plot(DePData(Act_PRnT:nT)); title('Depoly Sum');
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
Ax = {Fh2Live,nT,Actin,dims,AcMx,SPYH,rot1,azel0};




varargout = {BTs,AxLP,Ax};

end
%####################################################################%







%####################################################################%
%					PLOTTING & HELPER FUNCTIONS
%####################################################################%

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
subplot('Position',[.03 .03 .65 .95]), 
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
subplot('Position',[.7 .55 .28 .38]), 
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



%==================================================%
%					TipMatrix
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
clrmap = [1 1 1; .9 .2 .2];
figure(Fh)
subplot('Position',[.7 .1 .28 .38]), 
imagesc(ActMx)
colormap(clrmap);
hold on
subplot('Position',[.7 .1 .28 .38]), 
scatter(PSDX,PSDY, 'r')
hold off
%-----------------------------------%


%==================================================%
% if nT >20000; keyboard; end
varargout = {ActMx};
end
%==================================================%


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
SPYhX = dims(4);
SPYhY = dims(5);
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
PSDactMx = zeros(SPYhY+100,SPYhX+100);
for mxp = 1:numel(PSDXY(:,1))
PSDactMx(PSDXY(mxp,2)+SPYhY+10, PSDXY(mxp,1)+SPYhX+10) = 1;
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
% PSDactMx(PSDXY(mxp,2)+SPYhY+mxpv, PSDXY(mxp,1)+SPYhX+mxpv) = 1;
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







