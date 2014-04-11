function [varargout] = clustersim(LBR,TIME,SIZE,DOES,GLU,GT,GTab)
clc, close all; scsz = get(0,'ScreenSize');
% If function receives no inputs (run as stand-alone)
if nargin < 1 

LBR = [2     2    10     1    15     4];
TIME = [3600 .0014 10 10 0 20];
SIZE = [10 4];
DOES = [1 1 0 100];
GLU = [0 10 1 10 -5];
GT = [10 -5 0];
GTab = [0 1 0;1 1 1;0 1 0];
% SVALS=padarray(ones(PSDsz),[PSAsz PSAsz], 0);

end

%-------------------------------%
% Parameters
%-------------------------------%
NSteps = TIME(1);
dT = TIME(2);
datarate = TIME(3);
viewTime = TIME(4);

PSDsz = SIZE(1);
PSAsz = SIZE(2);
SN0 = PSDsz^2;
SYNsz = PSDsz+(PSAsz*2);

doPlot = DOES(1);
doGaussianMask = DOES(2);
doFluorPlot=DOES(3);
FluorTime=DOES(4);




%-------------------------------------------------%
%				Mask Setup
%-------------------------------------------------%
hkMask=[0 1 0; 1 0 1; 0 1 0];
%----------------------------------------%
if doGaussianMask
%----------------------------------------%
A = .76;	x0=0; y0=0;	sx = .2; sy = .2; rx=sx; ry=sy;	sxy=3;

t = 0;
a = cos(t)^2/2/sx^2 + sin(t)^2/2/sy^2;
b = -sin(2*t)/4/sx^2 + sin(2*t)/4/sy^2 ;
c = sin(t)^2/2/sx^2 + cos(t)^2/2/sy^2;

[X, Y] = meshgrid((-sx*sxy):(rx):(sx*sxy), (-sy*sxy):(ry):(sy*sxy));
Z = A*exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;

hkMask=Z;

%----------------------------------------%
% 3D Gaussian Distribution
fh5 = figure(5); set(fh5,'OuterPosition',(scsz./[2e-3 2e-3 2 2]))
%----------------------------------------%
figure(fh5)
subplot('Position',[.05 .05 .40 .90]); ph5 = imagesc(hkMask); % ph5 = surf(X,Y,Z);
view(0,90); axis equal; drawnow;
xlabel('x-axis');ylabel('y-axis');zlabel('z-axis')
subplot('Position',[.55 .05 .40 .90]); ph6 = surf(X,Y,Z);
view(-13,22); shading interp; 
xlabel('x-axis');ylabel('y-axis');zlabel('z-axis')
%----------------------
for t = 0:pi/100:pi
	a = cos(t)^2/2/sx^2 + sin(t)^2/2/sy^2;
	b = -sin(2*t)/4/sx^2 + sin(2*t)/4/sy^2 ;
	c = sin(t)^2/2/sx^2 + cos(t)^2/2/sy^2;

	[X, Y] = meshgrid((-sx*3):(rx):(sx*3), (-sy*3):(ry):(sy*3));
	Z = A*exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;

	set(ph5,'CData',hkMask);
	drawnow;
	set(ph6,'XData',X,'YData',Y,'ZData',Z);
	drawnow;
end
%-------------------------------------------------%
end % doGaussianMask
%-------------------------------------------------%

%-------------------------------%
% Presets
%-------------------------------%
S=padarray(ones(PSDsz),[PSAsz PSAsz], 0);
S0=S;

Soffs = zeros(SYNsz);
Sons = zeros(SYNsz);
Csave = zeros(1,(NSteps/datarate));
doNoff = 0; NoffN = [];



%===============================================%
%				AMPAR STUFF
%===============================================%
% This section will eventually be populated by code to parse the 
% actual surface particle locations durring diffusion
doAMPARs = GLU(1);
AMPARN = GLU(4);

GhkMask=GTab;
GTon = GT(1);
GToff = GT(2);
LTPv = GT(3);

if (doAMPARs && doPlot)
Fh101 = figure(101);
set(Fh101,'OuterPosition',(scsz./[2e-1 .2 4 4]))
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
%---------------------------


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
if doAMPARs;
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
end;
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
FluorPlot(NNoff,SN0,stepN,ph10,ph11,Rev)
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

%-----------------------------------------------%
end % end main loop
%===============================================%

PlotCsave(Csave,datarate,dT,stepN,NSteps)
varargout = {NSteps,stepN,S0,S,CSen,CSex,NoffN,Soffs,Sons,doNoff};

%---------------------------------------%
end % end main function
%---------------------------------------%



%===============================================%
%				Subfunctions
%===============================================%
function [] = FluorPlot(NNoff,SN0,stepN,ph10,ph11,Rev)

FRAPL = (SN0-NNoff)./SN0;
SoCSN = (SN0-NNoff./Rev)./SN0;

set(ph10,'XData',[1:stepN],'YData',[FRAPL]);
drawnow; hold on;
set(ph11,'XData',[1:stepN],'YData',[SoCSN]);
drawnow; hold on;

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




