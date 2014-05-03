function [] = ReverseSAP2()
clc, close all, clear all

NSteps = 5000;
dT = .001;
datarate = 10;
viewTime = 10;
pauseTime = .0;


Csave = zeros(1,(NSteps/datarate));
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
hkMask=[0 1 0; 1 0 1; 0 1 0];
S=padarray(ones(7),[4 4], 0);
S0=S;

% Lon = 1.5;	% On Energy (lower = more on events)
% Qon = 14;	% On Neighbor-Independant Rate (new growth) (lower = more on)
% Ron = 22;	% On Neighbor-Dependant Rate (cluster fill-in) (higher = more on)
% 
% Loff = 2.0;	% Off Energy (higher = more off events)
% Qoff = 2;	% Off Neighbor-Inependant Rate (uniform off)  (lower = more off)
% Roff = 12;	% Off Neighbor-Dependant Rate (edge off) (higher = more off)


Lon = 1.8;	% On Energy (lower = more on events)
Qon = 32;	% On Neighbor-Independant Rate (new growth) (lower = more on)
Ron = 20;	% On Neighbor-Dependant Rate (cluster fill-in) (higher = more on)

Loff = 1.8;	% Off Energy (higher = more off events)
Qoff = 1.8;	% Off Neighbor-Inependant Rate (uniform off)  (lower = more off)
Roff = 3;	% Off Neighbor-Dependant Rate (edge off) (higher = more off)

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


% if mod(stepN, 501) == 0; keyboard; end;
%------------------------%
% Counters & Plots
%------------------------%

CSen(stepN) = numel(find(Soff));
CSex(stepN) = numel(find(Son));

if mod(stepN, viewTime) == 0
PlotClusFun1(S);
pause(pauseTime);
end

if mod(stepN, datarate) == 0
Csave(stepN/datarate) = sum(S(:));
end

if sum(S) < 1
PlotCsave(Csave,datarate,dT,stepN,NSteps)
Nenex(NSteps,stepN,S0,CSen,CSex)
return;
end


%-----------------------------------------------%
end % end main loop
%===============================================%

PlotCsave(Csave,datarate,dT,stepN,NSteps)
% 	set(gcf, 'PaperPositionMode', 'auto');
% 	print -depsc2 -tiff SClusterPlot


Nenex(NSteps,stepN,S0,CSen,CSex)



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
	
	
disp('This cluster lasted for...')
disp([ num2str(stepN) ' steps'])
disp(['of the ' num2str(NSteps) ' requested steps'])
disp(' ')

end


function [] = Nenex(NSteps,stepN,S0,CSen,CSex)


enPERstep = sum(CSen)/(stepN);
LifeRatio = stepN / enPERstep;
disp(['Cluster:Particle lifetime ratio: ' num2str(LifeRatio)  ])

OnPerStep = sum(CSex)/(stepN);
OffPerStep = sum(CSen)/(stepN);
OffStepRate = 1 / OffPerStep;
Offmins = OffStepRate/60;

Sstart = sum(sum(S0));
PperMin = 1/Offmins;
PperHr = PperMin*60;

PctClusPHr = PperHr / Sstart * 100;
NumClusPHr = PctClusPHr/100*Sstart;

disp([' ' ])
disp(['with...' ])
disp([' ' num2str(OnPerStep) ' on-events per step (on average)' ])
disp([' ' num2str(OffPerStep) ' off-events per step (on average)' ])

disp([' ' ])
disp(['if step = 1 second...' ])
disp([' ' num2str(OffStepRate) ' seconds between off events (on average)' ])
disp([' ' num2str(OffStepRate) ' seconds = ' num2str(Offmins) ' minutes'])

disp([' ' ])
disp(['The starting cluster size was ' num2str(Sstart) ' particles' ])
disp(['A total of ' num2str(NumClusPHr) ' particles dissociated per hour' ])
disp([' equivalent to ' num2str(PctClusPHr) ' percent of the starting cluster'])



end


