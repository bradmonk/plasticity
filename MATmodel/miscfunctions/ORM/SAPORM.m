function [varargout] = SAPORM(LBR,TIME,SIZE,MODS,DOES,REVA,GLU,GT,GTab,doTs)
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

%================================================%
%				FIGURE SETUP
%------------------------------------------------%
if doPlot
%--------
Fh1 = figure(1);
Ph1 = imagesc(S);
colormap('bone')
end
%---------------------------
if doFluorPlot
%---------------
fh10=figure(10); colormap('bone'); 
set(fh10,'OuterPosition',(scsz./[5e-3 5e-3 3 2.5]))
%----
ph10 = plot([1 2],[1 .9]);
axis([0 NSteps 0 1.1]); axis manual; hold on;
ph11 = plot([1 2],[1 .9],'r');
axis([0 NSteps 0 1.1]); axis manual; hold on;
%---------------
end
%-------------------------------------------------%

%===============================================%
for stepN = 1:NSteps
%-----------------------------------------------%

Pmx = rand(size(S));
Soc = (S>0);
Sno = ~Soc;
hk = convn(Soc,hkMask,'same');

Pon = 1 ./ (1+exp((hk-Lon).*(-Bon)));
Pkon = Sno .* ( Ron * dT * Pon );
Son = (Pkon>Pmx);

Poff = 1 ./ (1+exp(((-hk)+Loff).*(-Boff)));
Pkoff = Soc .* ( Roff * dT * Poff );
Soff = (Pkoff>Pmx);


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


%-------
if doPlot
if mod(stepN, viewTime) == 0
set(Ph1,'CData',S);
drawnow
end
end
%-------


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
% Nenex(NSteps,stepN,S0,S,CSen,CSex,NoffN,Soffs,Sons,doNoff)
varargout = {NSteps,stepN,S0,S,CSen,CSex,NoffN,Soffs,Sons,doNoff};
return;
end
%-------



%-----------------------------------------------%
end % end main loop
%===============================================%

PlotCsave(Csave,datarate,dT,stepN,NSteps)

% Nenex(NSteps,stepN,S0,S,CSen,CSex,NoffN,Soffs,Sons,doNoff)

varargout = {NSteps,stepN,S0,S,CSen,CSex,NoffN,Soffs,Sons,doNoff};

%---------------------------------------%
end % end main function
%---------------------------------------%


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
subplot('Position',[.05 .05 .40 .90]); ph5 = imagesc(hkMask); % ph5 = surf(X,Y,Z);
view(0,90); axis equal; drawnow;
xlabel('x-axis');ylabel('y-axis');zlabel('z-axis')
subplot('Position',[.55 .05 .40 .90]); ph6 = surf(X,Y,Z);
view(-13,22); shading interp; 
xlabel('x-axis');ylabel('y-axis');zlabel('z-axis')
%----------------%
end
%-------------------------------%

end



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


