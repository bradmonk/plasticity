function [] = SimpleSAP()
clc, close all, clear all

NSteps = 50000;
dT = .01;
datarate = 10;
viewTime = 100;
pauseTime = .0;


% Interesting Parameters:
%{
% Brads Longest Lasting (dT=1.0)
B = 80;
Le = 1.21;
rho = 0.9;		
r = 10;
mu = .15;

% Shouval (dT=0.01)
B = 60;
Le = 1.2;
rho = 0.95;		
r = 10;
mu = .9;
%}

Csave = zeros(1,(NSteps/datarate));

hkMask=[0 1 0; 1 0 1; 0 1 0];
S=padarray(ones(7),[4 4], 0);
B = 50;			% hi = shrink
Le = 1.5;		% hi = shrink
rho = 0.95;		% hi = grow
r = 10;			% hi = grow
mu = 1;		% hi = shrink
%-------------------------------%
for stepN = 1:NSteps
%-------------------------------%
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
%-------------------------------%

Pmx = rand(size(S));
Soc = (S>0);

hk = convn(Soc,hkMask,'same');
Lhk = hk-Le;

Pen = (Soc) .* mu * dT;
Pex = 1 ./ (1+exp((-1*B)*Lhk));
Pkex = (1-Soc) .* ( rho * r * dT * Pex );

Sen = (Pen>Pmx);
Sex = (Pkex>Pmx);

S = (Soc-Sen) + Sex;

%------------------------%
% if mod(stepN, 4000) == 0; keyboard; end;
CSen(stepN) = numel(find(Sen));
CSex(stepN) = numel(find(Sex));
% CSex(stepN) = numel(find(Soc))-numel(find(Sex));
%------------------------%

if mod(stepN, viewTime) == 0
PlotClusFun1(S);
pause(pauseTime);
end

if mod(stepN, datarate) == 0
Csave(stepN/datarate) = sum(S(:));
end

if sum(S) < 1
PlotCsave(Csave,datarate,dT,stepN,NSteps)
Nenex(NSteps,stepN,dT,CSen,CSex)
return;
end


%-------------------------------%
end % end main loop
%-------------------------------%

PlotCsave(Csave,datarate,dT,stepN,NSteps)
% 	set(gcf, 'PaperPositionMode', 'auto');
% 	print -depsc2 -tiff SClusterPlot


Nenex(NSteps,stepN,dT,CSen,CSex)



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
	set(get(gca,'XLabel'),'String','Time (min)')
	set(get(gca,'YLabel'),'String','SAP')
	set(gca,'YLim',[0 (max(Csave)+10)])
	xt = (get(gca,'XTick'))*datarate*(dT)/(60);
	set(gca,'XTickLabel', sprintf('%.1f|',xt));
	SAPtitle = title(['    Cluster Size Over Time '...
	'    Total Steps: ' num2str(NSteps)...
	'  \bullet  Step Duration: ' num2str(dT) ' s']);
	set(SAPtitle, 'FontSize', 12);
	
	
disp('This cluster lasted...')
disp(['seconds: ' num2str(stepN*dT)])
disp(['minutes: ' num2str(stepN*dT/60)])
disp(['hours: ' num2str(stepN*dT/60/60)])
disp(['days: ' num2str(stepN*dT/60/60/24)])
end


function [] = Nenex(NSteps,stepN,dT,CSen,CSex)

disp([ 'time step was: ' num2str(dT) ])

disp([ num2str(sum(CSen)/(stepN)) ' S.en per step' ])
disp([ num2str(sum(CSex)/(stepN)) ' S.ex per step' ])

enPERstep = sum(CSen)/(stepN);
LifeRatio = stepN / enPERstep;

disp([ num2str(LifeRatio) ' Cluster Lifetime : SAP Lifetime' ])


% disp([ num2str((stepN)/sum(CSen)) ' steps between S.en' ])
% disp([ num2str((stepN)/sum(CSex)) ' steps between S.ex' ])

% disp([ num2str(sum(CSen)/(stepN)*60*60) ' S.en per step' ])
% disp([ num2str(sum(CSex)/(stepN)*60*60) ' S.ex per step' ])

% disp([ num2str(sum(CSen)/(stepN*dT)) ' S.en per second' ])
% disp([ num2str(sum(CSex)/(stepN*dT)) ' S.ex per second' ])

% disp([ num2str((stepN*dT)/sum(CSen)) ' seconds between S.en' ])
% disp([ num2str((stepN*dT)/sum(CSex)) ' seconds between S.ex' ])

end


% Fold and Save These
%{
Pmx = rand(size(S));
Soc = (S>0);

hk = convn(Soc,hkMask,'same');
Lhk = Le - hk;
% Lhk = hk-Le;

Pen = Soc .* (1-exp(-mu*dT));
Pex = 1 ./ (1+exp(B*Lhk));
Pkex = (1-Soc) .* ( rho * r * dT * Pex );

Sex = (Pkex>Pmx);
Sen = (Pen>Pmx);

S = (Soc-Sen) .* S + Sex;
%}
%{
Pmx = rand(size(S));
Soc = (S>0);

hk = convn(Soc,hkMask,'same');
% Lhk = Le - hk;
Lhk = hk-Le;

Pen = (Soc) .* mu * dT;
Pex = 1 ./ (1+exp((-1*B)*Lhk));
Pkex = (1-Soc) .* ( rho * r * dT * Pex );

Sen = (Pen>Pmx);
Sex = (Pkex>Pmx);

S = (Soc-Sen) .* S + Sex;
%}
%{
%% See what Beta and Le does to the slope:
clc, close all, clear all
Le = 1.5;
Lhk = (0:4) - Le;
B = 1:10:50;

for i = 1:5
Pex(i,:) = 1 ./ (1 + exp(-B .* Lhk(i) ) );
end

figure(1)
plot(Lhk,Pex)
xt = (get(gca,'XTick'))+Le;
set(gca,'XTickLabel', sprintf('%.1f|',xt))
%%
%}

