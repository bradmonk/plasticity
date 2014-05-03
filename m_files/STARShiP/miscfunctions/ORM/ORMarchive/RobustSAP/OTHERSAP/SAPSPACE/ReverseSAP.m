function [] = ReverseSAP(Lon,Qon,Ron,Loff,Qoff,Roff)
clc, close all, clear all

%-------------------------------%
% Parameters
%-------------------------------%
NSteps = 5000;
dT = .01;
datarate = 10;
viewTime = 10;
pauseTime = .0;

PSDsz = 8;
PSAsz = 4;


%-------------------------------%
% Presets
%-------------------------------%

hkMask=[0 1 0; 1 0 1; 0 1 0];
SYNsz = PSDsz+(PSAsz*2);
S=padarray(ones(PSDsz),[PSAsz PSAsz], 0);
S0=S;

Soffs = zeros(SYNsz);
Sons = zeros(SYNsz);
Csave = zeros(1,(NSteps/datarate));
Ntic = 0; doNoff = 0; NoffN = [];

%-------------------------------%
% Lon = On Energy (lower = more on events)
% Qon = On Neighbor-Independant Rate (lower = more on) (new growth @<10) 
% Ron = On Neighbor-Dependant Rate (higher = more on) (cluster fill-in @>10) 
%--
% Loff = Off Energy (higher = more off events)
% Qoff = Off Neighbor-Inependant Rate (uniform off)  (lower = more off)
% Roff = Off Neighbor-Dependant Rate (edge off) (higher = more off)
%--------
% PROTOTYPE
% Lon = 1.8;	% (lower = faster ON rate)
% Qon = 20;		% (lower = faster ON lattice) (CONV INDY)
% Ron = 20;		% (higher = faster ON cluster) (CONV DEP)
% %--
% Loff = 1.8;	% (higher = faster OFF rate)
% Qoff = 1.2;	% (lower = faster OFF lattice) (CONV INDY)
% Roff = 3.0;	% (higher = faster OFF cluster) (CONV DEP)
% %--------
% Lon = 1.8;
% Loff = 1.8;
% %-- INDY (low=fast)
% Qon = 20; (CONV INDY)
% Qoff = 1; (CONV INDY)
% %-- DEP (high=fast)
% Ron = 4; (CONV DEP)
% Roff = 1; (CONV DEP)
% %--------
%--------
Lon = 2;	%(lo=fst)
Loff = 2;	%(hi=fst)
%-- INDY (low=fast)
Qon = 10;
Qoff = 1;
%-- DEP (high=fast)
Ron = 4;
Roff = 1;
%--------

% Qon: growth starts at <13, can go upwards to infinite with no impact
% Qoff: can be 0 and cluster is ok, >4 and no dots ever leave center
% If Roff is 3x Ron cluster breaks down quickly



%-------------------------------%
% ALTERNATIVE SCENARIOS
%{
%-------------------------------%
% UNSTABLE CLUSTER
% Lon = 1.8;	% 1.8 On Energy (lower = more on events)
% Qon = 32;		% 32 On Neighbor-Independant Rate (new growth) (lower = more on)
% Ron = 3;		% 20 On Neighbor-Dependant Rate (cluster fill-in) (higher = more on)
% 
% Loff = 1.8;	% 1.8 Off Energy (higher = more off events)
% Qoff = 1.2;	% 1.2 Off Neighbor-Inependant Rate (uniform off)  (lower = more off)
% Roff = 3;		% 3.0 Off Neighbor-Dependant Rate (edge off) (higher = more off)
%-------------------------------%

%-------------------------------%
% ULTRASTABLE CLUSTER
Lon = 1.8;	% On Energy (lower = more on events)
Qon = 32;	% On Neighbor-Independant Rate (new growth) (lower = more on)
Ron = 20;	% On Neighbor-Dependant Rate (cluster fill-in) (higher = more on)

Loff = 1.8;	% Off Energy (higher = more off events)
Qoff = 1.2;	% Off Neighbor-Inependant Rate (uniform off)  (lower = more off)
Roff = 3.0;	% Off Neighbor-Dependant Rate (edge off) (higher = more off)
%-------------------------------%
%}
%-------------------------------%

%===============================================%
for stepN = 1:NSteps
%-----------------------------------------------%

Pmx = rand(size(S));
Soc = (S>0);
Sno = ~Soc;
hk = convn(Soc,hkMask,'same');

Lhon = hk-Lon;
Pon = 1 ./ (1+exp((-1*Qon)*Lhon));
Pkon = Sno .* ( Ron * dT * Pon );
Son = (Pkon>Pmx);

Lhoff = (-hk)+Loff;
Poff = 1 ./ (1+exp((-1*Qoff)*+Lhoff));
Pkoff = Soc .* ( Roff * dT * Poff );
Soff = (Pkoff>Pmx);

S = (Soc-Soff) + Son;



%-------------------------------%
% Parameter Notes
%{
% Pen = (S) * mu * dT
% Pex = 1/(1+exp(-B*hk))
% Pkex = (1-S) ( rho * r * dT * Pex)
% 
% hk:	cluster force (2D field mask:Mx convelution)
% Le:	repulsion constant (lattice repulsion force)	[hi = shrink]
% Lhk:	(cluster force) - (repulsion constant)
% B:	slope of Pex(Lhk) function						[hi = shrink]
% rho:	probability an endo pool sap is available		[hi = grow]
% r:	transition rate from endo to empty site			[hi = grow]
% dT:	time-step interval
% mu:	internalization rate							[hi = shrink]
% S:	Surface SAP matrix (Soc:occupied | Snoc:not occupied)



% ORIGINAL SHOUVAL SETUP
% CONSTANT OFF, NEIGHBOR-DEPENDANT ON SETUP
hkMask=[0 1 0; 1 0 1; 0 1 0];
S=padarray(ones(7),[4 4], 0);
B = 60;			% hi = shrink
Le = 1.9;		% hi = shrink
rho = 0.9;		% hi = grow
r = 10;			% hi = grow
mu = .95;		% hi = shrink
Pmx = rand(size(S));
Soc = (S>0);
Sno = ~Soc;
hk = convn(Soc,hkMask,'same');
Lhk = hk-Le;
Pen = Soc .* mu * dT;
Pex = 1 ./ (1+exp((-1*B)*Lhk));
Pkex = Sno .* ( rho * r * dT * Pex );
Sen = (Pen>Pmx);
Sex = (Pkex>Pmx);
S = (Soc-Sen) + Sex;


% CONSTANT ON, NEIGHBOR-DEPENDANT OFF SETUP
Lhk = (-hk)+Le;
Pex = Sno .* mu * dT;
Pen = 1 ./ (1+exp((-1*B)*+Lhk));
Pken = Soc .* ( rho * r * dT * Pen );
% Pex = 1 ./ (1+exp((-1*B)*Lhk));
% Pkex = Sno .* ( rho * r * dT * Pex );
Sen = (Pken>Pmx);
Sex = (Pex>Pmx);
S = (Soc-Sen) + Sex;


% NEIGHBOR-DEPENDANT OFF & ON SETUP
hkMask=[0 1 0; 1 0 1; 0 1 0];
S=padarray(ones(7),[4 4], 0);
B = 10;
Le = 1.9;		
rho = 0.9;		
r = 10;			
% mu = .95;		
Lon = 1.8;
Loff = 1.9;
% for...
Pmx = rand(size(S));
Soc = (S>0);
Sno = ~Soc;
hk = convn(Soc,hkMask,'same');
Lhon = hk-Lon;
Pon = 1 ./ (1+exp((-1*B)*Lhon));
Pkon = Sno .* ( rho * r * dT * Pon );
Son = (Pkon>Pmx);
Lhoff = (-hk)+Loff;
Poff = 1 ./ (1+exp((-1*B)*+Lhoff));
Pkoff = Soc .* ( rho * r * dT * Poff );
Soff = (Pkoff>Pmx);
S = (Soc-Soff) + Son;
%...end

%}
%-------------------------------%
%------------------------%
% Counters & Plots
%------------------------%

Soffs = Soffs+Soff;
Sons = Sons+Son;

Noff = numel(find(Soff));
Non = numel(find(Son));

CSen(stepN) = Noff;
CSex(stepN) = Non;


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


if mod(stepN, viewTime) == 0
PlotClusFun1(S);
pause(pauseTime);
end

if mod(stepN, datarate) == 0
Csave(stepN/datarate) = sum(S(:));
end

if sum(S) < 1
PlotCsave(Csave,datarate,dT,stepN,NSteps)
Nenex(NSteps,stepN,S0,S,CSen,CSex,NoffN,Soffs,Sons,doNoff)
return;
end


%-----------------------------------------------%
end % end main loop
%===============================================%

PlotCsave(Csave,datarate,dT,stepN,NSteps)
% 	set(gcf, 'PaperPositionMode', 'auto');
% 	print -depsc2 -tiff SClusterPlot


Nenex(NSteps,stepN,S0,S,CSen,CSex,NoffN,Soffs,Sons,doNoff)


% profile viewer
%---------------------------------------%
end % end main function
%---------------------------------------%



function [] = PlotClusFun1(S)
%-------------------------------------------%
	figure(1);
	imagesc(S);
    colormap('bone')
    title('2D Particle Map');
    title('PSD1');
end


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
title('Activity Map : (ON + OFF) / 2 ');
colorbar
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


