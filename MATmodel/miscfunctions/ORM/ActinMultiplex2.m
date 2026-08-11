function [varargout] = ActinMultiplex2(LBR,TIME,SIZE,MODS,DOES,REVA,GLU,GT,GTab,doTs)
clc, close all; scsz = get(0,'ScreenSize');

%-------------------------------%
% Steps-Related Parameters
%-------------------------------%
if doTs(1)
NSteps = doTs(6);
elseif doTs(2)
NSteps = doTs(7);
elseif doTs(3)
NSteps = doTs(8);
elseif doTs(4)
NSteps = doTs(9);
elseif doTs(5)
NSteps = doTs(10);
else
NSteps = TIME(1);
end

%-------------------------------%
% Timing-Related Parameters
%-------------------------------%
dT = TIME(2);
datarate = TIME(3);
viewTime = TIME(4);
FluorTime=DOES(4);


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
doGaussianMask = GT(4);
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


%-------------------------------%
% Data Collection
%-------------------------------%
Soffs = zeros(SYNsz);
Sons = zeros(SYNsz);
Csave = zeros(1,(NSteps/datarate));
NoffN = [];



%-------------------------------%
% Mask Setup
%-------------------------------%
[hkMask dT LBR] = MaskFun(S,dT,LBR,doGaussianMask,PSAsz,scsz);




%===============================================%
%				AMPAR STUFF
%===============================================%

doAMPARs = GLU(1);
amparate = GLU(3);
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




loadActinTips = 0;
generateActinTips = 1;

%---
%{.
if loadActinTips
ATdat = which('ATdata1.mat');
load(ATdat);
ActinTips = ATs;
TipCellN = numel(ActinTips);
TipCellA = TipCellN - 200;
end
%}
%---


%---
%{.
if generateActinTips
ActSteps = 60000;
ActinTips = ActinMainStage(ActSteps);
assignin('base', 'ATs', ActinTips)

end
%}

if ~loadActinTips
	ATs = evalin('base', 'ATs');
end
%---



ActinTips = ATs;
TipCellN = numel(ActinTips);
TipCellA = TipCellN - 200;


ACTINp = ActinTips{TipCellA};
ACTINp(:,251:end) = [];
ACTINp(251:end,:) = [];

% ACTINp(:,601:end) = [];
% ACTINp(601:end,:) = [];
% ACTINp(:,1:200) = [];
% ACTINp(1:200,:) = [];




% imagesc(ActinTips{140})
%-------------------------------%
% AMask=[0 1 0; 1 0 1; 0 1 0];
% AMask=[1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1];
% hkMask=[1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1];
AMask=ones(5);
hkMask=ones(5);
ActinBranches = 50;
ActinSlotsV = 2;
ACTINpos = zeros(PSDsz);
ACTINrand=randi([1 (PSDsz*PSDsz)],1,ActinBranches);
ACTINpos(ACTINrand)=ActinSlotsV;
% ACTINp=padarray(ACTINpos,[PSAsz PSAsz], 0);
%-------------------------------%
Ahk = convn(ACTINp,AMask,'same');
S = (Ahk>0).*1.0;


%================================================%
%				FIGURE SETUP
%------------------------------------------------%
if doPlot
%--------

%----------------------------------------
%	  FIGURE SETUP
%----------------------------------------
sz = [2.5e-3 2.5e-3 1.6 1.2];
Fh1 = FigSetup(11,sz);
Ph1 = imagesc(S);
colormap('bone')
end
%---------------------------
if doFluorPlot
%---------------
fh10=figure(14); colormap('bone'); 
set(fh10,'OuterPosition',(scsz./[5e-3 5e-3 3 2.5]))
%----
ph10 = plot([1 2],[1 .9]);
axis([0 NSteps 0 1.1]); axis manual; hold on;
ph11 = plot([1 2],[1 .9],'r');
axis([0 NSteps 0 1.1]); axis manual; hold on;
%---------------
end
%-------------------------------------------------%





%=========================================================================%
%=========================================================================%
for stepN = 1:NSteps
%-----------------------------------------------%


%------------
if mod(stepN,20) == 0
TipCellA = TipCellA+1;
ACTINp = ActinTips{TipCellA};
ACTINp(:,801:end) = [];
ACTINp(801:end,:) = [];
Ahk = convn(ACTINp,AMask,'same');
% S = (Ahk>0).*1.0;
end
%------------




%------------
Pmx = rand(size(S));
Soc = (S>0);
Sno = ~Soc;
%---
% Aoc = (ACTINp>0);
% Ah = convn(Aoc,hkMask,'same');
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

%======================================%
		S = (Soc-Soff) + Son;
%======================================%

%{
%====================================%
%		Counters & Plots
%====================================%

Soffs = Soffs+Soff;
Sons = Sons+Son;

Noff = numel(find(Soff));
Non = numel(find(Son));

CSen(stepN) = Noff;
CSex(stepN) = Non;


%-------
if doNoff
if Noff >= 1;
	Ntic = Ntic+1;
	Foff = hk .* Soff;
	[~,~,Fv] = find(Foff);
	NoffN(Ntic,1) = sum(Fv==1);
	NoffN(Ntic,2) = sum(Fv==2);
	NoffN(Ntic,3) = sum(Fv==3);
	NoffN(Ntic,4) = sum(Fv==4);
end;
end
%-------
%}

%-------
if doPlot
if mod(stepN, viewTime) == 0
set(Ph1,'CData',S);
drawnow
end
end
%-------

%{
%-------
if doFluorPlot
	NSoff(stepN) = Noff;
	NNoff(stepN) = sum(NSoff);
if mod(stepN, FluorTime) == 0
FluorPlot(NSoff,NNoff,SN0,stepN,ph10,ph11,Rev)
end
end
%-------


%-------
if mod(stepN, datarate) == 0
Csave(stepN/datarate) = sum(S(:));
end
%-------


%-------
if sum(S) < 1
PlotCsave(Csave,datarate,dT,stepN,NSteps)
varargout = {NSteps,stepN,S0,S,CSen,CSex,NoffN,Soffs,Sons,doNoff};
return;
end
%-------
%}


%-----------------------------------------------%
end % end main loop
%=========================================================================%
%=========================================================================%

PlotCsave(Csave,datarate,dT,stepN,NSteps)

varargout = {NSteps,stepN,S0,S};

%---------------------------------------%
end % end main function
%=========================================================================%





%=============================================================================%
%==================================================%
%			MATRIX CONVOLUTION MASK
%--------------------------------------------------%
function [hkMask dT LBR] = MaskFun(S,dT,LBR,doGaussianMask,PSAsz,scsz)

%-------------------------------%
%		Mask Setup
%-------------------------------%
hkMask=[0 1 0; 1 0 1; 0 1 0];
%----------------%
if doGaussianMask
%----------------%
A = 2;	x0=0; y0=0;	sx = .2; sy = .2; rx=sx; ry=sy;	res=2;

t = 0;
a = cos(t)^2/2/sx^2 + sin(t)^2/2/sy^2;
b = -sin(2*t)/4/sx^2 + sin(2*t)/4/sy^2 ;
c = sin(t)^2/2/sx^2 + cos(t)^2/2/sy^2;

[X, Y] = meshgrid((-sx*res):(rx):(sx*res), (-sy*res):(ry):(sy*res));
Z = A*exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;

hkMask=Z;
hk = convn(S,hkMask,'same');
hkor = hk(PSAsz+1,PSAsz+1);
LBR(1) = hkor-sqrt(A); LBR(2) = hkor+sqrt(A);
% dT=dT*2;

%----------------%
% 3D Gaussian Distribution
fh5 = figure(5); set(fh5,'OuterPosition',(scsz./[2e-3 2e-3 2 2]))
%----------------%
figure(fh5)
subplot('Position',[.05 .05 .40 .90]); ph5 = imagesc(hkMask); 
ph5 = surf(X,Y,Z);
view(0,90); axis equal; drawnow;
xlabel('x-axis');ylabel('y-axis');zlabel('z-axis')
subplot('Position',[.55 .05 .40 .90]); ph6 = surf(X,Y,Z);
view(-13,22); shading interp; 
xlabel('x-axis');ylabel('y-axis');zlabel('z-axis')
%----------------%
end
%-------------------------------%

end
%--------------------------------------------------%

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
%=============================================================================%







%=============================================================================%
function [varargout] = ActinMainStage(ActSteps,varargin)
clc, close all; clear all; scsz = get(0,'ScreenSize');
BTs = [];
Nsteps = 50000;
% Nsteps = ActSteps;
dT = 1/1000;


%---------------------------------------------%
%				ANGLE SETUP
%---------------------------------------------%
% unitsratio('rad', 'deg')
rad2dg = 180/pi;
dg2rad = 1/(180/pi);
Ovec = [1.0 28.7 56.4 84.1 111.8 139.5 167.2 194.8 222.5 250.2 277.9 305.6 333.3];
Ovec = Ovec * dg2rad;



%----------------------------------------------------------------------------%
							Actin = zeros(5,12);
%[Nact	Xang	Xorg	Xtip	Yang	Yorg	Ytip	Zang	Zorg	Ztip	Lact	OrO ]%
%[1		2		3		4		5		6		7		8		9		10		11		12	]%
%----------------------------------------------------------------------------%
GaSize = 2.5;		% 5.1 / 2;	% Actin Size
Actin(:,1) = 300;	% N monomers in 5 Starting Filaments


%{
% (90 deg = pi/2) (70: 7*pi/18) (60: pi/3) (45: pi/4) (30: pi/6)
TPi = (pi/2);
% TPi = (7*pi/18);
PPi = 0;
ArpX = 70;
ArpY = 13;
ArpZ = 0;

% Act = struct([]);
Act.Nm = 0;


%=============================================================================%
Po = [1.5;.5;.5];
Pt = [2;-1.5;1.5];
Dv = Pt-Po;
% Pr = Pt+(Veq.*sign(Pt+Si));
%----------------------------------------------------------------------------%
Nfils = 1;
for fN = 1:Nfils; 
	Act(fN).Nm = 270;
	Act(fN).Nm = Act(fN).Nm .* GaSize;
	Act(fN).Tta = TPi;
	Act(fN).Phi = PPi;
	Act(fN).ArpX = ArpX; 
	Act(fN).ArpY = ArpY; 
	Act(fN).ArpZ = ArpZ;
	Act(fN).Po = [1.5;.5;.5];
	Act(fN).Pt = [2;-1.5;1.5];
	Act(fN).Dv = Act(fN).Pt - Act(fN).Po;
	Act(fN).Pr = [0;0;0];
	Act(fN).Pe = [0;0;0];	
	Act(fN).Dr = zeros(3,14);
end;
%----------------------------------------------------------------------------%
TL0 = sqrt(sum((Pt - Po).^2));
Txyz = (Pt - Po) ./ TL0;
TL = sqrt(sum((Txyz).^2));
Ttta = acos(Txyz(3)/TL) + TPi;
Tphi = atan(Txyz(2)/Txyz(1)) + PPi;
LTx = (TL) * sin(Ttta) * cos(Tphi);
LTy = (TL) * sin(Ttta) * sin(Tphi);
LTz = (TL) * cos(Ttta);

Pe = [LTx;LTy;LTz]+Pt;
Pr = Pe;

for fN = 1:Nfils; 
Act(fN).Pe = Pe;
Act(fN).Pr = Pr;
end


fPo = Act(fN).Po;
fPt = Act(fN).Pt;
fPr = Act(fN).Pr;
fDv = Act(fN).Dv;
fDr = [];

%----------------------------------------------------------------------------%
% TPi = (pi/2);
TPi = (7*pi/18);
PPi = 0;
%----------------------------------------------------------------------------%
for nLp = 1:3
Nfils = 50;
%=============================================================================%
for fN = 1:Nfils;
%=============================================================================%


Fdots = 14;
for Fbranch = 1:Fdots
	%Rmono = ceil(200 * rand);		% Random monomer along segment
	%Rmono13 = mod(Rmono,13)+1;		% Get monomer repeat among the 13 rotational axis angles
	Rmono13 = mod(Fbranch,13)+1;			% Get monomer repeat among the 13 rotational axis angles
	tta = Ovec(Rmono13);			% Rotational angle of new branch
	fDr(:,Fbranch) = RotateVertex(fPr(1),fPr(2),fPr(3),fPt(1),fPt(2),fPt(3),fPo(1),fPo(2),fPo(3),...
								  fDv(1),fDv(2),fDv(3),tta);
	Act(fN).Dr(:,Fbranch) = fDr(:,Fbranch);
end
%----------------------------------------------------------------------------%
	Act(fN).Po = fPo;
	Act(fN).Pt = fPt;
	Act(fN).Pr = fPr;
	Act(fN).Dv = fDv;
	Act(fN).Dr = fDr;

	%---------------------------------------------%
	Rfil = ceil(200 * rand);
	Pfil = mod(Rfil,Fdots)+1;
	%-----
	fPo = fPt;
	fPt = fDr(:,Pfil); % fPt = fPr;
	fDv = fPt-fPo;
	%---------------------------------------------%
	TL0 = sqrt(sum((fPt - fPo).^2));
	Txyz = (fPt - fPo) ./ TL0;
	TL = sqrt(sum((Txyz).^2));
	Ttta = acos(Txyz(3)/TL) + TPi;
	Tphi = atan2(Txyz(2),Txyz(1)) + PPi;
	TL=2;
	LTx = (TL) * sin(Ttta) * cos(Tphi);
	LTy = (TL) * sin(Ttta) * sin(Tphi);
	LTz = (TL) * cos(Ttta);
	Pe = [LTx;LTy;LTz]+fPt;
	%-----
	fPr = Pe;
	%---------------------------------------------%

%=============================================================================%
end
%=============================================================================%




%---------------------------------------------%
%				FIGURE SETUP
%---------------------------------------------%
sz = [2.5e-3 2.5e-3 1.6 1.2];
fpos = [.05 .05 .95 .95];
Fh25 = FigSetup2(25,sz,fpos);
%---------------------------------------------%

[QZQ FilN] = size(Act);
for PLoop = 1:FilN-1
	
Po = Act(PLoop).Po;
Pt = Act(PLoop).Pt;
Pr = Act(PLoop).Pr;
Dv = Act(PLoop).Dv;
Dr = Act(PLoop).Dr;






%---------------------------------------------%
%		PLOT ORIGIN->TIP LINE
%---------------------------------------------%
figure(Fh25);
[XMx YMx ZMx] = plot3prep({Pt},{Po});
hA1 = plot3(XMx,YMx,ZMx,'LineWidth',2,'MarkerSize',10);
%----------------
hold on
set(hA1,{'Color'},{'r'}')
grid on;
xlabel('X');ylabel('Y');zlabel('Z');
axis equal
hold on;
%---------------------------------------------%
% axvals = {[-10 10 -10 10 -10 10]./2};
axvals = {[-10 10 -10 10 -10 10].*2};
addaxis(axvals);
axis(axvals{1})
%---------------------------------------------%





%---------------------------------------------%
%		PLOT ORIGIN->TIP LINE
%---------------------------------------------%
figure(Fh25);
[XMx YMx ZMx] = plot3prep({Pt Pr},{Po Pt});
hA1 = plot3(XMx,YMx,ZMx,'LineWidth',2,'MarkerSize',10);
%----------------
hold on
set(hA1,{'Color'},{'r','b'}')
grid on;
xlabel('X');ylabel('Y');zlabel('Z');
axis equal
hold on;
%---------------------------------------------%
% axvals = {[-10 10 -10 10 -10 10]./2};
axvals = {[-10 10 -10 10 -10 10]};
addaxis(axvals);
axis(axvals{1})
%---------------------------------------------%

%---------------------------------------------%
%	LINE: Pr-Pt Circle to Point
%---------------------------------------------%
PtMx = repmat(Pt,1,Fdots);
PtPrMxX = [PtMx(1,:); Dr(1,:)];
PtPrMxY = [PtMx(2,:); Dr(2,:)];
PtPrMxZ = [PtMx(3,:); Dr(3,:)];

hA2 = plot3(PtPrMxX,PtPrMxY,PtPrMxZ,'c');
% FMx.Fdots = {.5};
% set(hA2,{'LineWidth'},repmat(FMx.Fdots,1,Fdots)');
% hold on;
%---------------------------------------------%


end

end

keyboard
%}

Actin(:,11) = Actin(:,1); % Length of 5 Starting Filaments

% Branching Angles
Pi90 = (pi/2);
Pi70 = (7*pi/18);
Pi30 = pi/6;
Pi45 = pi/4;
Pi20 = 20*dg2rad;

TPi = Pi20;
PPi = 0;

% Angle of 5 Starting Filaments
Actin(:,2) = 0;	% X angle
Actin(:,5) = 0;	% Y angle
Actin(:,8) = 0;	% Z angle

		
% Origin of 5 Starting Filaments
Actin(2,3) = 20;		% X origin
Actin(2,6) = 20;		% Y origin
Actin(2,9) = 0;			% z origin
Actin(3,3) = 20;		% X origin
Actin(3,6) = -20;		% Y origin
Actin(3,9) = 0;			% z origin
Actin(4,3) = -20;		% X origin
Actin(4,6) = 20;		% Y origin
Actin(4,9) = 0;			% z origin
Actin(5,3) = -20;		% X origin
Actin(5,6) = -20;		% Y origin
Actin(5,9) = 0;			% z origin

% TRIG: Filament Tip Locations
% x = Actin(:,1) .* sin(TPi) .* cos(PPi) + Actin(:,3);
% y = Actin(:,1) .* sin(TPi) .* sin(PPi) + Actin(:,6);
% z = Actin(:,1) .* cos(TPi) + Actin(:,9);
% Actin(:,4) = x;
% Actin(:,7) = y;
% Actin(:,10) = z;

% MATH - branch XYZ tip coordinates
Actin(:,4) = Actin(:,1) .* sin(Actin(:,8)) .* cos(Actin(:,2)) + Actin(:,3);
Actin(:,7) = Actin(:,1) .* sin(Actin(:,8)) .* sin(Actin(:,2)) + Actin(:,6);
Actin(:,10) = Actin(:,1) .* cos(Actin(:,8)) + Actin(:,9);



%----------------------------------------
% Spine Dimensions
SPYneckXY = round(150 / GaSize);
SPYheadZN = round(1000 / GaSize); % North Z dim
SPYheadZS = round(700 / GaSize); % South Z dim
SPYheadX = round(300 / GaSize); % X dim
SPYheadY = round(300 / GaSize); % Y dim

PSDproxy = round(100 / GaSize);
inPSD = SPYheadZN - PSDproxy;

SPYH = [SPYheadX SPYheadY];
AcMx = zeros(SPYheadY*2+10,SPYheadX*2+10);

dims = [SPYneckXY SPYheadZN SPYheadZS SPYheadX SPYheadY PSDproxy...
	inPSD];
%----------------------------------------


%----------------------------------------
% Quantitative Values (mol, volume, rate parameters)

% Units
mM = 1e-3;
uM = 1e-6;
upM = 1e3;
dnM = 1e-3;

%{
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
%}

mol = 6e23;		% in N

% Spine volume (0.1 um^3) in L and uL
SpyV = 1e-16;	% in L
SpyVu = .1e-9;	% in uL

% Actin Polymerization (12 N/然*s)
Act_Na = 1e5;
% Act_Nb = 1e6;
% Act_N = Act_Na;
Act_PRnT = 1000;

% Actin Depolymerization (2 N/s)	
Act_DRnT = 200;
DePSum = 0;

% Cofilin Depoly
CofR = .0005;
CofN = 40;
CofS = 200;

% Actin Poly Math
% Cytosolic concentration of actin in cells ranges from .1 to .5 mM
% Given a spine volume of 'SpyV' we can find how many actin monomers
% are in an average spine:
%
% .1 mM (mmol/L) * (1 mol / 1000 mmol) * SpyV = 1e-17 mol
% 1e-17 mol * (6e23 units / 1 mol) = 6e3 monomer units

Act_N = .1 * (1/1000) * SpyV * 6e23; % 6e3 monomer units

% we can check our math starting with a set number of actin monomers
% and calculate the spine molarity (6e3 monomer units as an example):
% 
% 6e3 units/SpyV * (1 mol / 6e23 units) * (1000 mmol / 1 mol)
% 6e3/SpyV*(1/6e23) 

Act_mM = Act_N/SpyV * (1/6e23) * (1000/1);	% 1.6e-10 

% Act_mM = Act_N / SpyVu / mol;	% 1.6e-10 
% Act_N = Act_N / SpyVu / mol;
Act_N = Act_N;
Act_PR = 12 * Act_mM * dT;
Act_DR = 2 * dT;



% Arp Branching Rate
ArpRa = .0005;
ArpRb = .0005;
ArpR = ArpRa;
ArpST = 3000;
ArpScalar = 30;
ArpDec = .001;

% Arp Branching Math
ArpN = 1e3;
ArpOn = 5;
ArpOff = 1;
Arp_uM = ArpN / SpyVu / mol;	% 1.6 - 16 uM
Arp_PR = ArpOn * Arp_uM * dT;
Arp_DR = ArpOff * dT;

%-------------------------------------------%

% Equation to linear transform range [a b] to range [c d] and solve for [x]:
%				c*(1-(x-a)/(b-a)) + d*((x-a)/(b-a))
% example:
global sja;
global sjb;
global sjc;
global sjd;

sja=0;
sjb=350;
sjc=-.3;
sjd=3;

linsc = @(jx) (sjc*(1-(jx-sja)/(sjb-sja)) + sjd*((jx-sja)/(sjb-sja)));
% miscfun = @(xyz,tta) ([1 0 0; 0 cos(tta) sin(tta); 0 -sin(tta) cos(tta)] * xyz);
%-------------------------------------------%


%----------------------------------------


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
%----------------------------------------





%===============================================================%
%						MAIN OUTER LOOP
for nT = 1:Nsteps
%---------------------------------------------------------------%

	% if nT == Act_PRnT; Act_N = Act_Nb; end;
	if nT == ArpST; ArpR = ArpRb; end;
	
	% radial distance to spine shaft membrane
	ActinXYh = sqrt(Actin(:,4).^2 + Actin(:,7).^2);

	Xspi = ActinXYh > SPYneckXY;			% Logical array of filaments beyond SPYneckXY
	Yspi = ActinXYh > SPYneckXY;			% Logical array of filaments beyond SPYneckXY
	Zspi = (abs(Actin(:,10)) > SPYheadZN);	% Logical array of filaments beyond SPYheadZN

	Xhed = ActinXYh >= SPYheadX;
	Yhed = ActinXYh >= SPYheadY;
	Zhed = (abs(Actin(:,10)) <= SPYheadZN) & (abs(Actin(:,10)) >= SPYheadZS);
		
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
		if Act_PR > rv(1)
		if ((~sum(ActMem(aN,:)) && Act0(aN)) || Zhed(aN))
		if (~Xhed(aN) && ~Yhed(aN))

			Actin(aN,1) = Actin(aN,1) + ceil(Act_PR);
			
			Act_N = Act_N - ceil(Act_PR);

		end
		end
		end
		%=================================================================%
		
		
		
		
		
		%=================================================================%
		% BRANCHING
 		% if ( ArpR * (Actin(aN,1) / ArpScalar ) ) > rv(2)
		if ( ArpR * exp(linsc(Actin(aN,1)-ArpScalar)) ) > rv(2)
		if ((~sum(ActMem(aN,:)) && Act0(aN)) || Zhed(aN))
		if (~Xhed(aN) && ~Yhed(aN))
			fNf = NFact+1;
			
			Actin(fNf,1) = 20;		% create branch: add 20 subunits to new branch
			
			Nmono = Actin(aN,1);		% N monomers in mother filament
			Rmono = ceil(Nmono * rand); % Random monomer (== vector_length) along segment
			Rmono13 = mod(Rmono,13)+1;	% Get monomer repeat among the 13 rotational axis angles
			Rang = Ovec(Rmono13);		% Rotational angle of new branch
			Actin(fNf,12) = Rang;		% Store this rotational angle

			% New branch XYZ origin coordinates
			Ox = (Rmono) * sin(Actin(aN,8)) * cos(Actin(aN,2)) + Actin(aN,3);
			Oy = (Rmono) * sin(Actin(aN,8)) * sin(Actin(aN,2)) + Actin(aN,6);
			Oz = (Rmono) * cos(Actin(aN,8)) + Actin(aN,9);
			Actin(fNf,3) = Ox;	% X origin
			Actin(fNf,6) = Oy;	% Y origin
			Actin(fNf,9) = Oz;	% Z origin
			
			% New branch XYZ temporary tip coordinates
			TOx = (Rmono+1) * sin(Actin(aN,8)) * cos(Actin(aN,2)) + Actin(aN,3);
			TOy = (Rmono+1) * sin(Actin(aN,8)) * sin(Actin(aN,2)) + Actin(aN,6);
			TOz = (Rmono+1) * cos(Actin(aN,8)) + Actin(aN,9);
			
			%=====================================================%	
 			% fPo = [Actin(NFact,3);Actin(NFact,6);Actin(NFact,9)];
			fPo = [Ox;Oy;Oz];
			fPt = [TOx;TOy;TOz];
			fDv = fPt-fPo;
			%---------------------------------------------%
			TL0 = sqrt(sum((fPt - fPo).^2));
			Txyz = (fPt - fPo) ./ TL0;
			TL = sqrt(sum((Txyz).^2));
			Ttta = acos(Txyz(3)/TL) + TPi;
			Tphi = atan2(Txyz(2),Txyz(1)) + PPi;
			TL=1;
			LTx = (TL) * sin(Ttta) * cos(Tphi);
			LTy = (TL) * sin(Ttta) * sin(Tphi);
			LTz = (TL) * cos(Ttta);
			fPr = [LTx;LTy;LTz]+fPt;
			%---------------------------------------------%
			% keyboard

			Rmono = ceil(200 * rand);		% Random monomer along segment
			Rmono13 = mod(Rmono,13)+1;		% Get monomer repeat among the 13 rotational axis angles
			tta = Ovec(Rmono13);			% Rotational angle of new branch

			fDr = RotateVertex( fPr(1),fPr(2),fPr(3),...
								fPt(1),fPt(2),fPt(3),...
								fPo(1),fPo(2),fPo(3),...
								fDv(1),fDv(2),fDv(3),...
								tta );

			fPo = fPt;
			fPt = fDr;
			fDv = fPt-fPo;
			%---------------------------------------------%
			TL0 = sqrt(sum((fPt - fPo).^2));
			Txyz = (fPt - fPo) ./ TL0;
			TL = sqrt(sum((Txyz).^2));
			Ttta = acos(Txyz(3)/TL) + TPi;
			Tphi = atan2(Txyz(2),Txyz(1)) + PPi;
			TL=10;
			LTx = (TL) * sin(Ttta) * cos(Tphi);
			LTy = (TL) * sin(Ttta) * sin(Tphi);
			LTz = (TL) * cos(Ttta);
			fPr = [LTx;LTy;LTz]+fPt;
			%=======================================================%
			
			% Store angle
			Actin(fNf,2) = Tphi;
			Actin(fNf,5) = 0;
			Actin(fNf,8) = Ttta;
			
			% New branch XYZ tip coordinates
			Actin(fNf,4) = LTx  + Ox;	% X tip
			Actin(fNf,7) = LTy  + Oy;	% Y tip
			Actin(fNf,10) = LTz + Oz;	% Z tip
			
			ArpR = ArpR - (ArpR*ArpDec);
			% ArpN = ArpN - 1;
			% ArpRate = ArpRate*(ArpMol/ArpMol0);
			
		end
		end
		end
		%=================================================================%


		
		
		
		%=================================================================%
		% DEPOLYMERIZATION
		if nT > Act_DRnT
		if Act_DR > rv(3);
			
			Actin(aN,1) = Actin(aN,1)-ceil(Act_DR) .* (Actin(aN,1)>0);
			
			DePSum = DePSum +1;
			
			Act_N = Act_N + ceil(Act_DR);
			
		end
		end
		%---------------
		if nT > CofS
		if CofR > rv(3)
			Actin(aN,1) = Actin(aN,1)-CofN;
			
			Actin(aN,1) = Actin(aN,1) .* (Actin(aN,1)>0);
			
			Act_N = Act_N + CofN;
			% Act_N = Act_N + (CofN*(Actin(aN,1)>0));
		end
		end
		%=================================================================%
		
		
		
		
		
		%=================================================================%
		% ADJUST RATE VALUES
		
		Act_mM = Act_N/SpyV * (1/6e23) * (1000/1);
		Act_PR = 12e3 * (Act_mM * dT);
		
		% Act_mM = Act_N / SpyVu / mol;	% .17 uM
		% Act_N = Act_mM * (1/1000) * SpyV * 6e23;
		% Act_PR = 12 * Act_uM * dT;
		
		if Act_N < 0; Act_N = 0; Act_PR = 0; end;
		%=================================================================%
		


	%-----------------------------------------------------------------------------%
	end
	%						MAIN INNER LOOP
	%=============================================================================%
	
	
	NoAct = find(Actin(:,1)<1);
	Actin(NoAct,:) = [];
	if numel(NoAct); 
		ArpR = ArpR + (ArpR*ArpDec*numel(NoAct)); 
	end;
	
	Actin(:,1) = Actin(:,1);	% Length of Factin segments
	
	%{
% 	%=================================================%
% 	NFact = numel(Actin(:,1));
% 	for NFa = 1:NFact
% 	%---------------------------------------------%
% 	fPo = [Actin(NFact,3);Actin(NFact,6);Actin(NFact,9)];
% 	fPt = [Actin(NFact,4);Actin(NFact,7);Actin(NFact,10)];
% 	%---------------------------------------------%
% 	TL0 = sqrt(sum((fPt - fPo).^2));
% 	Txyz = (fPt - fPo) ./ TL0;
% 	TL = sqrt(sum((Txyz).^2));
% 	Ttta = acos(Txyz(3)/TL);
% 	Tphi = atan2(Txyz(2),Txyz(1));
% 	TL1 = Actin(NFa,1).*GaSize;
% 	LTx = (TL1) * sin(Ttta) * cos(Tphi);
% 	LTy = (TL1) * sin(Ttta) * sin(Tphi);
% 	LTz = (TL1) * cos(Ttta);
% 	Pe = [LTx;LTy;LTz]+fPo;
% 	%-----
% 	fPr = Pe;
% 	%---------------------------------------------%
% 	end
% 	%=================================================%
	%}
	
	
	% MATH - branch XYZ tip coordinates
	Actin(:,4) = Actin(:,1) .* sin(Actin(:,8)) .* cos(Actin(:,2)) + Actin(:,3);
	Actin(:,7) = Actin(:,1) .* sin(Actin(:,8)) .* sin(Actin(:,2)) + Actin(:,6);
	Actin(:,10) = Actin(:,1) .* cos(Actin(:,8)) + Actin(:,9);
	
	ActAngs = [Actin(:,2) Actin(:,5) Actin(:,8)];

	
	
	%==================================================%
	%				LIVE PLOT
	%--------------------------------------------------%
	if mod(nT,500) == 0
		LivePlot(Fh2Live,nT,Actin,inPSD,rot1,azel0,dims);
		rot1 = rot1 + rot0;
		% disp(nT); disp(Act_N);disp(Act_PR);disp(Act_mM);
		
	end
	%--------------------------------------------------%
	% if mod(nT,11500) == 0; keyboard; end;
	% if mod(nT,2) == 0; keyboard; end;
	
	%==================================================%
	%				SAVE TipMatrix
	%--------------------------------------------------%
	if nT >4000;
	if mod(nT,100) == 0
		ActMx = TipMatrix(Fh2Live,nT,Actin,dims,AcMx,SPYH,rot1,azel0);
		BTs{numel(BTs)+1} = ActMx;
	end
	end
	%--------------------------------------------------%


	
	
	
	
% ActData(nT,:) = [Act_mM, Act_PR, Act_N];
% DePData(nT,:) = DePSum;
% DePSum = 0;

% if nT == 125; keyboard; end;

%-----------------------------------------------------------------------------%
end
%						MAIN OUTER LOOP
%=============================================================================%




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
PSDactMx = zeros(SPYheadY+100,SPYheadX+100);
for mxp = 1:numel(PSDXY(:,1))
PSDactMx(PSDXY(mxp,2)+SPYheadY+10, PSDXY(mxp,1)+SPYheadX+10) = 1;
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

varargout = {BTs};

end
%=============================================================================%

%{.
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
% SPYneckXY = dims(1);
% SPYheadZN = dims(2);
% SPYheadZS = dims(3);
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
% PSDactMx(PSDXY(mxp,2)+SPYheadY+mxpv, PSDXY(mxp,1)+SPYheadX+mxpv) = 1;
AcMx(PSDY(mxp), PSDX(mxp)) = 1;
end
ActMask=[1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1];
ActMx = convn(AcMx,ActMask,'same');
ActMx = (ActMx>0).*1.0;
%===================================%
%				FIGURE
%-----------------------------------%
figure(Fh)
subplot('Position',[.7 .1 .28 .38]), 
imagesc(ActMx)
colormap(bone)
hold on
% subplot('Position',[.55 .05 .40 .40]), 
subplot('Position',[.7 .1 .28 .38]), 
scatter(PSDX,PSDY, 'r')
% scatter(PSDXY(:,1),PSDXY(:,2), 'r')
hold off
%-----------------------------------%


%==================================================%
% if nT >20000; keyboard; end


varargout = {ActMx};
end
%=============================================================================%
%}

















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



%{
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
% SPYneckXY = dims(1);
% SPYheadZN = dims(2);
% SPYheadZS = dims(3);
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
% PSDactMx(PSDXY(mxp,2)+SPYheadY+mxpv, PSDXY(mxp,1)+SPYheadX+mxpv) = 1;
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
%}









