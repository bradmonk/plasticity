function [varargout] = OFFDEPSAP(QRvals,params)
clc, close all

%-------------------------------%
% Parameters
%-------------------------------%
NSteps = params(1);
dT = params(2);
datarate = 10;
viewTime = 10;
pauseTime = .0;



PSDsz = 8;
PSAsz = 4;

doAMPARs = 0;
AMPARN = 3;
amparate = 50;
G1RT = .1;

doFluorPlot=1;
FluorTime=100;
SN0 = PSDsz^2;

%================================================%
%				FIGURE SETUP
%------------------------------------------------%
Fh1 = figure(1);
Flh = figure(10);

set(Flh,'Units','pixels');	scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/7  scnsize(4)/2  scnsize(3)/5  scnsize(4)/3];
set(Flh,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])

% P1h = bar(1,.4);
% axis([0 2 0 1]) 
% axis manual

P1h = plot([1 2],[1 .9]);
axis([0 NSteps 0 1.1]) 
axis manual
hold on
%----------------------------------------------------------------------%


% IndOffRate = 200;
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

%===============================================%
%{
% Lon = On Energy (lower = more on events)
% Bon = On Neighbor-Independant Rate (lower = more on) (new growth @<10) 
% Ron = On Neighbor-Dependant Rate (higher = more on) (cluster fill-in @>10) 
%--
% Loff = Off Energy (higher = more off events)
% Boff = Off Neighbor-Inependant Rate (uniform off)  (lower = more off)
% Roff = Off Neighbor-Dependant Rate (edge off) (higher = more off)
%--------
% PROTOTYPE
% Lon = 1.8;	% (lower = faster ON rate)
% Bon = 20;		% (lower = faster ON lattice) (CONV INDY)
% Ron = 20;		% (higher = faster ON cluster) (CONV DEP)
% %--
% Loff = 1.8;	% (higher = faster OFF rate)
% Boff = 1.2;	% (lower = faster OFF lattice) (CONV INDY)
% Roff = 3.0;	% (higher = faster OFF cluster) (CONV DEP)
% Bon: growth starts at <13, can go upwards to infinite with no impact
% Boff: can be 0 and cluster is ok, >4 and no dots ever leave center
% If Roff is 3x Ron cluster breaks down quickly
% %--------
% Lon = 1.8;
% Loff = 1.8;
% %-- INDY (low=fast)
% Bon = 20; (CONV INDY)
% Boff = 1; (CONV INDY)
% %-- DEP (high=fast)
% Ron = 4; (CONV DEP)
% Roff = 1; (CONV DEP)
% %--------
%-- B: INDY (low=fast) R: DEP (high=fast)
% SHOUVAL RATE CONSTANTS
dT = .01
Lon = 1.5;
Bon = 60;
Ron = 9.5;
Loff = 2;
Boff = 0;
Roff = 2;
%}
%--------
Lon = QRvals(5);	%(lo=fst)
Loff = QRvals(6);	%(hi=fst)
%-- INDY (low=fast)
Bon = QRvals(1);
Boff = QRvals(3);
%-- DEP (high=fast)
Ron = QRvals(2);
Roff = QRvals(4);
%--------

%===============================================%
%===============================================%
for stepN = 1:NSteps
%-----------------------------------------------%

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

% if mod(stepN, IndOffRate)==0
% RIoff = Soc .* (rand(size(S)) .* .1);
% SRIoff = (RIoff>Pmx);
% end


%====================================%
if doAMPARs
if mod(stepN, amparate) == 0
%-------------------------------%
SG1oc = S .* 0.0;
%-------------------------------%
RPOS=randi([4 12],1,AMPARN);CPOS=randi([4 12],1,AMPARN);
SG1oc(RPOS,CPOS)=1;
SG1oc(RPOS,CPOS+1)=1;
SG1oc(RPOS+1,CPOS)=1;
SG1oc(RPOS+1,CPOS+1)=1;
%-------------------------------%
Gex = Sno .* SG1oc .* G1RT + Pkon;
Gp = 1 ./ (1+exp((-1*1)*hk)) + .001;
Gpex = Gex .* Gp;
%---

Son = (Gpex>Pmx);
%-------------------------------%
end; end;
%====================================%

%-------------------------------%
S = (Soc-Soff) + Son;
%-------------------------------%

% if mod(stepN, IndOffRate)==0
% S = (Soc-Soff-SRIoff) + Son;
% end


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



if doFluorPlot
	NSoff(stepN) = Noff;
	NNoff(stepN) = sum(NSoff);
if mod(stepN, FluorTime) == 0
FluorPlot(NSoff,NNoff,SN0,stepN,Flh,P1h)
end
end

if mod(stepN, datarate) == 0
Csave(stepN/datarate) = sum(S(:));
end

if sum(S) < 1
PlotCsave(Csave,datarate,dT,stepN,NSteps)
Nenex(NSteps,stepN,S0,S,CSen,CSex,NoffN,Soffs,Sons,doNoff)
varargout = {NSteps,stepN,S0,S,CSen,CSex,NoffN,Soffs,Sons,doNoff};
return;
end


%-----------------------------------------------%
end % end main loop
%===============================================%

PlotCsave(Csave,datarate,dT,stepN,NSteps)
% 	set(gcf, 'PaperPositionMode', 'auto');
% 	print -depsc2 -tiff SClusterPlot


Nenex(NSteps,stepN,S0,S,CSen,CSex,NoffN,Soffs,Sons,doNoff)

varargout = {NSteps,stepN,S0,S,CSen,CSex,NoffN,Soffs,Sons,doNoff};
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

function [] = FluorPlot(NSoff,NNoff,SN0,stepN,Flh,P1h)
%=================================%
%			LIVE FRAP
%=================================%
% pause(.5)

% if stepN==900;keyboard;end;



FRAP = (SN0-(sum(NSoff)/15))/SN0;
FRAPL = (SN0-NNoff./15)./SN0;

set(P1h,'XData',[1:stepN],'YData',[FRAPL]);
drawnow

% set(P1h,'XData',[(stepN-100) stepN],'YData',[1 FRAP]);  %# Update the y data of the line
% drawnow 

% hold off
% set(0,'CurrentFigure',Flh)
% subplot(3,3,[7 9]),
% bar(numel(NSoff),FRAP)
% % axis manual
% % axis([0 100])
% hold on
%=================================%
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


