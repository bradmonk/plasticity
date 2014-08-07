function [varargout] = ORMATHFUN(LBR,TIME,SIZE,MODS,DOES,REVA,GLU,GT,GTab,doTs)
%============================================%
clc, close all; scsz = get(0,'ScreenSize');
%============================================%
%{
%-- This rate can be calculated using two similar equations
%- 1. Exponential Decay Rate (solution to differential equation)
%- 2. Exponential Mean Lifetime (if decay elements are discrete units)

% Nt = No * exp(-Lam*dt);	% [1]

% Nt = No * exp(-dt/Tau);	% [2]
% 
% 
% Lam (?) : exponential decay constant
% Tau (?) : exponential mean element lifetime
% Nt      : quantity remaining at time t
% No      : quantity at start
% dt      : elapsed time from start
% 
% - ? and ? are proportional:
% 
% ? = 1/?
% ? = 1/?

%- ? is the time at which the population of elements is reduced to:
%- 1/e = 0.368 * initial value.
%- If No=100, at time dt=?, Nt=36.8
%}
%============================================%


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
S0=S; So=sum(sum(S0));

%-------------------------------%
% doItems
%-------------------------------%
doGaussianMask = GT(4);
doRev = REVA(1);
doPlot = DOES(1);
doFluorPlot=DOES(3);
doNoff = 0; 

%-------------------------------%
% Run Models
%-------------------------------%
OnRM = MODS(1);
OffRM = MODS(2);
ORM = MODS(3);

%-------------------------------%
% Revolving Particle Entry
%-------------------------------%
Rev=REVA(2);
if doRev
MuT = dT*Rev;
else
MuT = dT;
end
%--------


%-------------------------------%
% Mask Setup
%-------------------------------%
hkMask=[0 1 0; 1 0 1; 0 1 0];

% [hkMask dT LBR] = MaskFun(S,dT,LBR,doGaussianMask,PSAsz,scsz);
% MuT=dT;



%-------------------------------%
% Data Collection
%-------------------------------%
Soffs = zeros(SYNsz);
Sons = zeros(SYNsz);
Csave = zeros(1,(NSteps/datarate));
NoffN = [];


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


%======================================================================%
%======================================================================%
for re=1:TIME(6)
%----------------
S = S0;

	%----------------
	for t = 1:NSteps
	%----------------
		Pmx = rand(size(S));
		Soc = (S>0);
		Sno = ~Soc;
		hk = convn(Soc,hkMask,'same');

		%====================================%
		if ORM
			Pon = 1 ./ (1+exp((hk-Lon).*(-Bon)));
			Pkon = Sno .* ( Ron * MuT * Pon );
			Son = (Pkon>Pmx);

			Poff = 1 ./ (1+exp(((-hk)+Loff).*(-Boff)));
			Pkoff = Soc .* ( Roff * MuT * Poff );
			Soff = (Pkoff>Pmx);
			
			%-------------------------------%
			if doAMPARs
			%----------
			SG1oc = zeros(SYNsz);
			GRPOS=randi([1 (SYNsz*SYNsz)],1,AMPARN);
			SG1oc(GRPOS)=1;
			GRhk = convn(SG1oc,GhkMask,'same');
			GRk=(GRhk.*GTon); GSk=(GRhk.*GToff);
			%----------
			Gon = Pkon+(GRk.*(Pkon+LTPv));
			Goff = Pkoff+(GSk.*Pkoff);

			Son = (Gon>Pmx);
			Soff = (Goff>Pmx);
			%----------
			end
			%-------------------------------%
					S = (Soc-Soff) + Son;
			%====================================%
			
		elseif OffRM
			Poff = 1 ./ (1+exp(((-hk)+Loff).*(-Boff)));
 			% Poff = 1 ./ (1+exp((hk-Loff).*Boff)); % this is equivalent to Poff above
			Pkoff = Soc .* ( Roff * MuT * Poff );
			Soff = (Pkoff>Pmx);
			
			S = (Soc-Soff);
			
		elseif OnRM
			Poff = 1 ./ (1+exp(((-hk)+Loff).*(-Boff)));
			Pkoff = Soc .* (Roff * MuT);
			Soff = (Pkoff>Pmx);
			
			S = (Soc-Soff);
		end
		
		
		% INSIDE COUNTERS
		SPoff = Soc .* Poff;
		SnSum(t) = sum(sum(S));
		Noff(t) = sum(Soff(:));	
		
		
	% LIVE PLOT
	if doPlot
	set(Ph1,'CData',S);
	drawnow
	end
	%----------------
	end %%%%%%%%%%%%%
	%----------------


% OUTSIDE COUNTERS
if ~OnRM
Pms(re) = mean(Poff(:));
SPms(re) = mean(nonzeros(SPoff));
end


RNoff(re,:) = Noff;
SNoff(re) = sum(Noff(:));

SnSums(re,:) = SnSum;
Sns(re) = sum(sum(S));

%-------------------
end %%%%%%%%%%%%%%%%
%======================================================================%
%======================================================================%







%==========================================================%
%		OnRM (Fixed Off) Calculations
%==========================================================%
if OnRM  %%%%%%%%%%%%
%--------------------
Sn = mean(Sns);				% Particles Remaining
Sd = So-Sn;					% Particles Diminished

Lam = (MuT*Roff);			% Lamda Calculation
Tau = 1/(MuT*Roff);			% Tau Calculation
LSn = So * exp(-Lam*t);		% Particles Remaining Lamda Calc
TSn = So * exp(-t/Tau);		% Particles Remaining Tau Calc (adjusted)
%-------------------


%-------------------
FSo = 1.0;
FSn = Sn/So;		% Fluor Remaining (Simulated)
FSd = 1-FSn;		% Fluor Diminished (Simulated)
FCn = TSn/So;		% Fluor Remaining (Calculated)
FCd = 1-FCn;		% Fluor Diminished (Calculated)

%-------------------
disp('    Sn    LSn   TSn')
disp(round([Sn  LSn  TSn]))
disp('   MuT    Lam     Tau')
fprintf('%0.5f ',MuT),fprintf(' %0.5f ',Lam),fprintf(' %0.1f\n',Tau)
LBRo = 1.00001*[Loff,Boff,Roff];
disp('    Loff      Boff      Roff')
disp(LBRo)
%------------------------------------------%
figure(9); set(gcf,'Position',[600 400 800 400])
%-----
axes('Position',[.05 .05 .55 .9])
plot(mean(SnSums))
axis([0 t 0 So])
%------------------------------------------%
figure(9);
axes('Position',[.67 .05 .30 .9])
fhb = bar([FSd FCd; FSn FCn]);
axis([.5 2.5 0 1]);
set(gca,'YLim',[0 1]);
set(get(gca,'YLabel'),'String','Fluorescence')
xt = {'Loss','Remain'};
set(gca,'XTickLabel', sprintf('%s|',xt{1:2}))
legend(fhb,{'sim','calc'});
%------------------------------------------%




%======================================================================%
SnTempM=[];
for nn = 0:59
SnTemp1 = circshift(SnSums,[0 (-60*nn)]);
SnTemp2 = SnTemp1(:,1:60);
SnTempM(:,nn+1) = mean(SnTemp2,2);
end
CellSnSums = {SnTempM'};
%======================================================================%
%					FINAL OUTPUT (SPLINE) FIGURE
%======================================================================%
fig21 = figure(21);
pos1 = [scsz(3)/3  scsz(4)/5  scsz(3)/2.5  scsz(4)/2];
set(fig21,'OuterPosition',pos1)
%----------------------------------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6]; c11=[.9 .3 .3];
%===========================================================%
% FIG1 TOP LEFT: GluR Subtypes in Synapses
%===========================================================%
sbpos = [.09 .13 .85 .80]; ptype = 4;
cOLOR = [c1; c2; c3; c4; c1; c2; c3; c4];
itemN = 1; 
[ph1 hax1] = CIplot(CellSnSums,sbpos,itemN,cOLOR,ptype);
hold on
%------------------------------------------%
MS1 = 5;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
hTitle  = title('Fluorescence Loss Over 60 Minutes');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Fluorescence Level (+/- SEM)');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel, hTitle],'FontSize',12);
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
axis([0 t/60 0 So])
%======================================================================%


%--------------------
end; %%%  OnRM   %%%%
%==========================================================%





%==========================================================%
%		OffRM (Dependent Off) Calculations
%==========================================================%
if OffRM  %%%%%%%%%%%%
%--------------------
Pm = mean(Pms);				% mean(Poff(:));
SPm = mean(SPms);			% mean(nonzeros(SPoff));
SNo = mean(SNoff);			% sum(Noff(:));

Sn = mean(Sns);				% Particles Remaining
Sd = So-Sn;					% Particles Diminished

Lam = (MuT*Roff);			% Lamda Calculation
Tau = 1/(MuT*Roff * SPm);	% Tau Calculation (adjusted)
LSn = So * exp(-Lam*t);		% Particles Remaining Lamda Calc
TSn = So * exp(-t/Tau);		% Particles Remaining Tau Calc (adjusted)

Lam0 = (MuT*Roff);			% Lamda Calculation
Tau0 = 1/(MuT*Roff);		% Tau Calculation
LSn0 = So * exp(-Lam0*t);	% Particles Remaining Lamda Calc
TSn0 = So * exp(-t/Tau0);	% Particles Remaining Tau Calc (adjusted)
%------------------------------------------%

FSo = 1.0;
FSn = Sn/So;		% Fluor Remaining (Simulated)
FSd = 1-FSn;		% Fluor Diminished (Simulated)
FCn = TSn/So;		% Fluor Remaining (reCalculated)
FCd = 1-FCn;		% Fluor Diminished (reCalculated)
FCn0 = TSn0/So;		% Fluor Remaining (Calculated)
FCd0 = 1-FCn0;		% Fluor Diminished (Calculated)

%------------------------------------------%

disp('    Sn    LSn   TSn    LSn0   TSn0')
disp(round([Sn  LSn  TSn  LSn0  TSn0]))
disp('   MuT    Lam     Tau    Lam0     Tau0')
fprintf('%0.5f ',MuT),fprintf(' %0.5f ',Lam),fprintf(' %0.1f ',Tau)
fprintf(' %0.5f ',Lam0),fprintf(' %0.1f\n',Tau0)
LBRo = 1.00001*[Loff,Boff,Roff];
disp('    Loff      Boff      Roff')
disp(LBRo)

%------------------------------------------%
figure(9); set(gcf,'Position',[600 400 800 400])
%-----
axes('Position',[.05 .05 .55 .9])
plot(mean(SnSums))
axis([0 t 0 So])
%------------------------------------------%
figure(9);
axes('Position',[.67 .05 .30 .9])
fhb = bar([FSd FCd FCd0; FSn FCn FCn0]);
axis([.5 2.5 0 1]);
set(gca,'YLim',[0 1]);
set(get(gca,'YLabel'),'String','Fluorescence')
xt = {'Loss','Remain'};
set(gca,'XTickLabel', sprintf('%s|',xt{1:2}))
legend(fhb,{'sim','calc','calc0'});
%------------------------------------------%


%======================================================================%
SnTempM=[];
for nn = 0:59
SnTemp1 = circshift(SnSums,[0 (-60*nn)]);
SnTemp2 = SnTemp1(:,1:60);
SnTempM(:,nn+1) = mean(SnTemp2,2);
end
CellSnSums = {SnTempM'};
%======================================================================%
%					FINAL OUTPUT (SPLINE) FIGURE
%======================================================================%
fig21 = figure(21);
pos1 = [scsz(3)/3  scsz(4)/5  scsz(3)/2.5  scsz(4)/2];
set(fig21,'OuterPosition',pos1)
%----------------------------------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6]; c11=[.9 .3 .3];
%===========================================================%
% FIG1 TOP LEFT: GluR Subtypes in Synapses
%===========================================================%
sbpos = [.09 .13 .85 .80]; ptype = 4;
cOLOR = [c1; c2; c3; c4; c1; c2; c3; c4];
itemN = 1; 
[ph1 hax1] = CIplot(CellSnSums,sbpos,itemN,cOLOR,ptype);
hold on
%------------------------------------------%
MS1 = 5;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
hTitle  = title('Fluorescence Loss Over 60 Minutes');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Fluorescence Level (+/- SEM)');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel, hTitle],'FontSize',12);
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
axis([0 t/60 0 So])
%======================================================================%


%--------------------
end; %%%  OffRM   %%%%
%==========================================================%





%==========================================================%
%		ORM Calculations
%==========================================================%
if ORM   %%%%%%%%%%%%
%--------------------

Sn = So - mean(SNoff);

MRN = mean(RNoff);
subs1 = 1:36;
subs2 = repmat(subs1,[100 1]);
subs3 = reshape(subs2,1,3600);
MRNoff = (accumarray(subs3', MRN))';
CSNoff = cumsum(MRNoff); % last value (36) should match mean(SNoff)



%------------------------------------------%
Lam = (MuT*Roff);			% Lamda Calculation
LSn = So * exp(-Lam*t);		% Particles Remaining Lamda Calc

fXoff = (mean(MRNoff) * (t/100))/100;

%------------------------------------------%

FSd = (mean(SNoff))/100;
FSn = (So - mean(SNoff))/100;

FCn = (So * (1-fXoff))/100;
FCd = (So * (fXoff))/100;

%------------------------------------------%
figure(9); set(gcf,'Position',[600 400 800 400])
%-----
ah9a = axes('Position',[.05 .05 .55 .9]);
plot(So-CSNoff); hold on;
axis([0 t/100 0 So])
xt = (get(gca,'XTick'))*100;
set(gca,'XTickLabel', sprintf('%.0f|',xt))
%------------------------------------------%
figure(9);
axes('Position',[.67 .05 .30 .9])
fhb = bar([FSd FCd; FSn FCn]);
axis([.5 2.5 0 1]);
set(gca,'YLim',[0 1]);
set(get(gca,'YLabel'),'String','Fluorescence')
xt = {'Loss','Remain'};
set(gca,'XTickLabel', sprintf('%s|',xt{1:2}))
legend(fhb,{'sim','calc'});
%------------------------------------------%



%=====================================================%
if doRev
%---------------------------------------------%
SoCSN = So-cumsum(log10(MRNoff));
CSrev = So-(CSNoff ./ Rev);
%---------------------------------------------%
axes(ah9a)
% could also do plot(ah9a,CSrev);
%---------------------------------------------% 
ph1 = plot(SoCSN);
leg1 = legend(ph1,'log10');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%---
ph2 = plot(CSrev);
legend([OUTH;ph2],OUTM{:},'cumsum / Rev');
[LEGH,OBJH,OUTH,OUTM] = legend;
axis([0 t/100 0 So])
xt = (get(gca,'XTick'))*100;
set(gca,'XTickLabel', sprintf('%.0f|',xt))
hold on
%---
set(ph1,'LineStyle','-','Color',[.9 .2 .2],'LineWidth',3,...
'Marker','none','MarkerSize',6,'MarkerFaceColor',[.9 .3 .3]);
set(ph2,'LineStyle',':','Color',[.2 .4 .6],'LineWidth',3,...
'Marker','none','MarkerSize',6,'MarkerFaceColor',[.3 .5 .7]);
%---------------------------------------------%
end
%=====================================================%


% Mean ending cluster size
disp(Sns); disp(mean(Sns));
%======================================================================%
cumsumRNoff = cumsum(RNoff,2);
csRNoff = So-cumsumRNoff;

SnTempM=[];
for nn = 0:59
SnTemp1 = circshift(csRNoff,[0 (-60*nn)]);
SnTemp2 = SnTemp1(:,1:60);
SnTempM(:,nn+1) = mean(SnTemp2,2);
end

CellSnSums = {SnTempM'};

%------
cumsumOFFS = cumsum(RNoff,2);
fluorLeft = So-cumsumOFFS;

SnTempM=[];
for nn = 0:59
SnT1 = circshift(RNoff,[0 (-60*nn)]);
SnT2 = SnT1(:,1:60);
SnT3(:,nn+1) = sum(SnT2,2);
end
SnT4 = So-(cumsum(SnT3,2)./Rev);

CellSnSumsLog10 = {SnT4'};
%------


%======================================================================%
%					FINAL OUTPUT (SPLINE) FIGURE
%======================================================================%
fig21 = figure(21);
pos1 = [scsz(3)/3  scsz(4)/5  scsz(3)/2.5  scsz(4)/2];
set(fig21,'OuterPosition',pos1)
%----------------------------------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6]; c11=[.9 .3 .3];
%===========================================================%
% FIG1 TOP LEFT: GluR Subtypes in Synapses
%===========================================================%
sbpos = [.09 .13 .85 .80]; ptype = 4;
cOLOR = [c1; c2; c3; c4; c1; c2; c3; c4];
itemN = 1; 
[ph1 hax1] = CIplot(CellSnSums,sbpos,itemN,cOLOR,ptype);
hold on
itemN = 1; cOLOR = [c2; c3; c4; c1; c2; c3; c4; c1];
[ph2 hax2] = CIplot(CellSnSumsLog10,sbpos,itemN,cOLOR,ptype);
% legend([OUTH;ph2],OUTM{:},'SnSums');
% [LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%------------------------------------------%
MS1 = 5;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c2);
hTitle  = title('Fluorescence Loss Over 60 Minutes');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Fluorescence Level (+/- SEM)');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel, hTitle],'FontSize',12);
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
axis([0 t/60 0 So])
%======================================================================%



%----------------------------------------------------------%
end; %%%  ORM   %%%%
%==========================================================%






% keyboard

%============================================%
varargout = {[FSd FSn]};
%============================================%
end










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



function varargout = CIplot(varargin)



dataINo = varargin{1};
sbpos   = varargin{2};
itemN   = varargin{3};
cOLOR   = varargin{4};
ptype	= varargin{5};



dataIN = dataINo;

%==============================================%
Mu = [mean(dataIN{itemN},2)]';
Sd = [std(dataIN{itemN},0,2)]';
Se = [Sd./sqrt(numel(dataIN{itemN}(1,:)))];
y_Mu = Mu;
x_Mu = 1:(size(y_Mu,2));
e_Mu = Se;
xx_Mu = 1:0.1:max(x_Mu);
% yy_Mu = spline(x_Mu,y_Mu,xx_Mu);
% ee_Mu = spline(x_Mu,e_Mu,xx_Mu);
yy_Mu = interp1(x_Mu,y_Mu,xx_Mu,'pchip');
ee_Mu = interp1(x_Mu,e_Mu,xx_Mu,'pchip');
p_Mu = polyfit(x_Mu,Mu,3);
x2_Mu = 1:0.1:max(x_Mu);
y2_Mu = polyval(p_Mu,x2_Mu);


%===========================================================%
plottype = ptype;

switch plottype
    case 1 
		subplot('Position',sbpos),...
        confplot(x_Mu,y_Mu,e_Mu);
		hold on;
		hax = gca;
		ph = plot(x_Mu,y_Mu,'o',xx_Mu,yy_Mu);
    case 2
		subplot('Position',sbpos),...
        ph = errorbar(Mu,Se,'x');
		set(ph,'LineStyle','none','LineWidth',1,'Marker','o',...
					'MarkerSize',5,'MarkerEdgeColor',[.2 .2 .2]);
		hax = gca;
    case 3 
		subplot('Position',sbpos),...
        ph = plot(x_Mu,y_Mu,'o',xx_Mu,yy_Mu);
		hax = gca;
	case 4 
		XT_Mu = xx_Mu';
		YT_Mu = yy_Mu';
		ET_Mu = ee_Mu';
		subplot('Position',sbpos),...
        [hl, hp] = boundedline(XT_Mu,YT_Mu, ET_Mu,...
		'cmap',cOLOR(itemN,:),'alpha','transparency', 0.4);
		ph = hl;
		po = hp;
		hax = gca;
	case 5 
		subplot('Position',sbpos),...
        ciplot((yy_Mu-ee_Mu),(yy_Mu+ee_Mu),xx_Mu,'b');
		% lighting gouraud; alpha(.4); 
		lighting phong; alpha(.4);
		hold on;
		hax = gca;
		ph = plot(x_Mu,y_Mu,'o',xx_Mu,yy_Mu);
    otherwise
        warning('Unexpected plot type.');
end
%===========================================================%



nargchk(0, 2, nargout);

if nargout >= 1
    varargout{1} = ph;
end

if nargout >= 2
    varargout{2} = hax;
end

if nargout >= 3
    varargout{3} = po;
end
% varargout={dataIN;hl};
% varargout={hl;hp};
return
%===========================================================%
%===========================================================%




%===========================================================%
%			ERROR BAR & CONFIDENCE ENVELOPE PLOTS
%===========================================================%
FIGgcf = gcf;
set(FIGgcf,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.2  scnsize(4)/1.7];
set(FIGgcf,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])

% use: subplot('Position',sbpos),...

%==============================================%
% confplot
% confplot(X,Y,SEM)
%------------------
% plots jagged line and jagged CI envelope
% SEM is distance from Y at each X
% X,Y,SEM must all have same numel
% Y = mean values; X = time; SEM of Y values
%------------------

%{
figure(FIGgcf);
subplot('Position',sbpos),...
confplot(x_Mu,y_Mu,e_Mu);
hold on
%}

figure(1);
confplot(x_Mu,y_Mu,e_Mu);
hold on
hTitle  = title ('confplot');

%==============================================%




%==============================================%
% errorbar
% errorbar(Y,SEM,'ok')
% errorbar(X,Y,SEM,'ok')
%-----------------------
% plots error bars 
% SEM is distance from Y at each X
% X,Y,SEM must all have same numel
% Y = mean values; X = time; SEM of Y values
%-----------------------

%{
subplot('Position',sbpos),...
FIGerb = errorbar(Mu,Se,'x');
set(FIGerb,'LineStyle','none','LineWidth',1,'Marker','o',...
'MarkerSize',5,'MarkerEdgeColor',[.2 .2 .2]);
hold on
%}

figure(1);
FIGerb = errorbar(Mu,Se,'x');
set(FIGerb,'LineStyle','none','LineWidth',1,'Marker','o',...
'MarkerSize',5,'MarkerEdgeColor',[.2 .2 .2]);
hTitle  = title ('errorbar');



%==============================================%




%==============================================%
% spline
%-----------------------
% Smoothed Plot Line (SPLine)
%
% 1) create fine grain x-axis timesteps
%	 >> X=timesteps; Y=Mu;
%    >> XX = 1:0.1:max(X);
%
% 2) use spline(X,Y,XX)
%    >> YY = spline(X,Y,XX);
%
% 3) plot smooth line
%	 >> plot(X,Y,'o',XX,YY);
% 
%-----------------------

%{
subplot('Position',sbpos),...
plot(x_Mu,y_Mu,'o',xx_Mu,yy_Mu);
hold on
%}

figure(1);
plot(x_Mu,y_Mu,'o',xx_Mu,yy_Mu);
hTitle  = title ('spline');

%==============================================%



%==============================================%
% boundedline
% boundedline(xx, yy, ee)
%-----------------------
% plots boundedline
%-----------------------

%{
subplot('Position',sbpos),...
[hl, hp] = boundedline(xx_Mu, yy_Mu, ee_Mu);
hold on
lighting gouraud; 
alpha(.4)


% XT_Mu = xx_Mu';
% YT_Mu = yy_Mu';
% ET_Mu = ee_Mu';
% %XT_Mu(2,:) = sqrt(xx_Mu)
% YT_Mu(2,:) = YT_Mu-(YT_Mu./2)
% ET_Mu(2,:) = ET_Mu-(ET_Mu./2)
% XT_Mu = XT_Mu';
% YT_Mu = YT_Mu';
% ET_Mu = ET_Mu';
% figure(2);
% [hl, hp] = boundedline(XT_Mu,YT_Mu, ET_Mu,...
% 'cmap', [.4 .1 .3; .1 .8 .7],'alpha','transparency', 0.07);
% 'cmap',lines(4),'alpha','transparency', 0.07);
% lighting gouraud; alpha(.4);
%}

figure(1);
[hl, hp] = boundedline(xx_Mu, yy_Mu, ee_Mu);
hTitle  = title ('boundedline');
lighting gouraud; 
alpha(.4)

%==============================================%


%==============================================%
% ciplot
% ciplot((YY-EE),(YY+EE),XX,'b')
%-----------------------
% confidence interval plot (CIPlot)
%
% 1) create fine grain x-axis timesteps
%	 >> X=timesteps; Y=Mu;
%    >> XX = 1:0.1:max(X);
%
% 2) use spline function to get YY and EE
%    >> YY = spline(X,Y,XX);
%    >> EE = spline(X,Y,XX);
%
% 3) plot
%	 >> ciplot((YY-EE),(YY+EE),XX,'b')
% 
%-----------------------

%{
subplot('Position',sbpos),...
ciplot((yy_Mu-ee_Mu),(yy_Mu+ee_Mu),xx_Mu,'b')
%ciplot((y_Mu-e_Mu),(y_Mu+e_Mu),x_Mu,'b')
% lighting gouraud; alpha(.4);
%}

figure(1);
ciplot((yy_Mu-ee_Mu),(yy_Mu+ee_Mu),xx_Mu,'b')
hTitle  = title ('ciplot');
lighting gouraud; 
alpha(.4)



%==============================================%






varargout={dataIN};


end





%======================================================================%
%					FINAL OUTPUT (SPLINE) FIGURE 1 OF 2
%======================================================================%
%{

SnTempM=[];
for nn = 0:59
SnTemp1 = circshift(SnSums,[0 (-60*nn)]);
SnTemp2 = SnTemp1(:,1:60);
SnTempM(:,nn+1) = mean(SnTemp2,2);
end

CellSnSums = {SnTempM'};

CellSnMeans = {mean(SnTempM)'};

%======================================================================%
%					FINAL OUTPUT (SPLINE) FIGURE
%======================================================================%
fig21 = figure(21);
set(21,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/2  scnsize(4)/2];
set(fig21,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
%===========================================================%
% FIG1 TOP LEFT: GluR Subtypes in Synapses
%===========================================================%
sbpos = [.09 .09 .85 .85]; ptype = 4;
cOLOR = [c1; c2; c3; c4; c1; c2; c3; c4];
itemN = 1; 
[ph1 hax1] = CIplot(CellSnSums,sbpos,itemN,cOLOR,ptype);
% leg1 = legend(ph1,'SnSums');
% [LEGH,OBJH,OUTH,OUTM] = legend;
hold on
% itemN = 1;
% [ph2 hax2] = CIplot(CellSnSums,sbpos,itemN,cOLOR,ptype);
% legend([OUTH;ph2],OUTM{:},'SnSums');
% [LEGH,OBJH,OUTH,OUTM] = legend;
% hold on
%------------------------------------------%
MS1 = 5; MS2 = 2;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
% set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
hTitle  = title('Fluorescence Loss Over 60 Minutes');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('Fluorescence Level (+/- SEM)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
axis([0 t/60 0 So])
%======================================================================%








%======================================================================%
%					FINAL OUTPUT (SPLINE) FIGURE
%======================================================================%
fig21 = figure(21);
set(21,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig21,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%

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



%}



%======================================================================%
%					FINAL OUTPUT (SPLINE) FIGURE 1 OF 2
%======================================================================%
%{
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



%}



%============================================%
%			MATH
%--------------------------------------------%
%{
%--------------------------------------------%
% In Malenka (2008), the PSD95 FRAP data shows that after 60 min,
% fluorescence goes down to 64 percent. The other data points are:

10 min: 82%
30 min: 78%
45 min: 66%
60 min: 64%

% This can be used to find the exponential decay rate (lamda) and the 
% exponential mean element lifetime (tau):

	exp(-Lam*t) = .64;
	exp(-t/Tau) = .64;

%-----
% Solving for Tau at 60 min...

	exp(-60/Tau) = .64

% Take the reciprocal of both sides:

	exp(60/Tau) = 1/.64

% Remove the exponent by taking the natural log of both sides:

	60/Tau = log(1/.64)
	60/Tau = 0.446

% Bring the variable (Tau) to the numerator by taking reciprocal:

	Tau/60 = 1/0.446

% Multiply both sides by 60

	Tau = 60/0.446

% Solve for Tau

	Tau = 134.5

%-----

Thus the exponential mean element lifetime (Tau) is 134.5 min.
If there were 100 PSD95 molecules at baseline (No) we can solve for how
many there were at 60 min (Nt) using:

	Nt = No * exp(-t/Tau)

% Solving for Nt...

	Nt = 100 * exp(-60/134.5)
	Nt = 64

As we expected [64/100 = .64], there was 64 percent of fluor remaining at 60 min. 

Tau from all time-point measurements 
10 min (82%): 50.4
30 min (77%): 114.8
45 min (66%): 108.3
60 min (64%): 134.5

Conversion to seconds
600  s (10m 82%): 3023
1800 s (30m 77%): 6886
2700 s (45m 66%): 6498
3600 s (60m 64%): 8066
(theoretically, these should all be the same)

% Using a tau of 7000 seems reasonable:
100 * exp(-3600/7000) % 60% remaining @ 1hour

% For the simulation, and in real life, we must assume that fluor molecules removed from
the cluster can reattach. This happens at a rate of 15:1 (fluor:new). Given this
adjustment factor, we can set a constant to achieve the desired off rate. 

Mu = 1/Tau
Tau = 1/Mu
1/(Mu*x) = 7000

% If Mu = 1
1/x = 7000
x = .000143

x*15 = 0.0021

Thus, the scalar should be around 0.002
Which is what I've been using all along...


%--------------------------------------------%
%}
%============================================%


%=============================%
%{
%============================================%
hkN = 3;
%--------------------------------------------%
Pk = dT * Roff / (1+ exp(((-hkN)+Loff) * (-Boff))) ;
LFTime = dT / Pk;
%--------------------------------------------%
PkoffSho = ((2 * Roff * exp(-Boff * (hkN-Loff))) / (1 + exp(-Boff * (hkN-Loff)))) * dT;
LFTimeSho = dT / PkoffSho;
%============================================%
disp([LFTime LFTimeSho LFTimeSho*2])
% Shouval's last 2* longer because of (2 .* Roff) - Why does he do that ???
% Also, I'm not exactly sure why we are doing 1+exp

% find the average number of neighbors in the ORM using:
NBRS(stepN,1) = find(hk,1);
NBRS(stepN,2) = find(hk,2);
NBRS(stepN,3) = find(hk,3);
NBRS(stepN,4) = find(hk,4);

%}
%=============================%

%=============================%
%{
clc, clear all, close all


Boff = 0;

for n=1:10
	
Nsteps = 80000;
dT = .01;
Lon = 1.5;
Loff = 1.5;
Bon = 0;
% Boff = 1;
Ron = .01;
Roff = 1;

Boff = Boff+1;

QRvals = [Bon Ron Boff Roff Lon Loff];
params = [Nsteps dT];


	
for m=1:10
[NSteps,stepN,S0,S,CSen,CSex,NoffN,Soffs,Sons,doNoff] = SAPORM(QRvals,params);

% NStepsOUT(m) = Nsteps;
stepNOUT(m) = stepN;
% S0OUT{m} = S0;
% SOUT{m} = S;
% CSenOUT{m} = CSen;
% CSexOUT{m} = CSex;
% NoffNOUT(m) = NoffN;
% SoffsOUT{m} = Soffs;
% SonsOUT{m} = Sons;
% doNoffOUT(m) = doNoff;

end

SLife(n,:) = stepNOUT;

end

SSLife = SLife.*dT;
SLF = {log10(SSLife)};

Boff = 0;
for m = 1:10
Boff=Boff+1;
NBRS = 1;
Pkoff = (2 .* Roff .* exp(-Boff * (NBRS-Loff))) ./ (1 + exp(-Boff * (NBRS-Loff))) .* dT;
LF1Time(m) = dT / Pkoff;
end

Boff = 0;
for m = 1:10
Boff=Boff+1;
NBRS = 2;
Pkoff = (2 .* Roff .* exp(-Boff * (NBRS-Loff))) ./ (1 + exp(-Boff * (NBRS-Loff))) .* dT;
LF2Time(m) = dT / Pkoff;
end

Boff = 0;
for m = 1:10
Boff=Boff+1;
NBRS = 3;
Pkoff = (2 .* Roff .* exp(-Boff * (NBRS-Loff))) ./ (1 + exp(-Boff * (NBRS-Loff))) .* dT;
LF3Time(m) = dT / Pkoff;
end

Boff = 0;
for m = 1:10
Boff=Boff+1;
NBRS = 4;
Pkoff = (2 .* Roff .* exp(-Boff * (NBRS-Loff))) ./ (1 + exp(-Boff * (NBRS-Loff))) .* dT;
LF4Time(m) = dT / Pkoff;
end

LogLF1Time = [];
LogLF2Time = [];
LogLF3Time = [];
LogLF4Time = [];

LogLF1Time = abs(log10(LF1Time));
LogLF2Time = log10(LF2Time);
LogLF3Time = log10(LF3Time);
LogLF4Time = log10(LF4Time);

LoogLF1Time = [LogLF1Time;(LogLF1Time.*.99);(LogLF1Time.*1.01);LogLF1Time;(LogLF1Time.*.99)];
LoogLF2Time = [LogLF2Time;(LogLF2Time.*.99);(LogLF2Time.*1.01);LogLF2Time;(LogLF2Time.*.99)];
LoogLF3Time = [LogLF3Time;(LogLF3Time.*.99);(LogLF3Time.*1.01);LogLF3Time;(LogLF3Time.*.99)];
LoogLF4Time = [LogLF4Time;(LogLF4Time.*.99);(LogLF4Time.*1.01);LogLF4Time;LogLF4Time];

LgLF1Time = {LoogLF1Time'};
LgLF2Time = {LoogLF2Time'};
LgLF3Time = {LoogLF3Time'};
LgLF4Time = {LoogLF4Time'};

%======================================================================%
%							FIGURE SETUP
%----------------------------------------------------------------------%
fig21 = figure(21);
set(21,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/3.1  scnsize(4)/1.5];
set(fig21,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
%===========================================================%
% FIG1 TOP LEFT: GluR Subtypes in Synapses
%===========================================================%
cOLOR = [c1; c2; c3; c4; c11; c22; c33; c44];
sbpos = [.1 .1 .8 .8]; ptype = 4;
%---
itemN = 1;  
[ph1 hax1] = CIplot(SLF,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'Cluster Sim');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%---
itemN = 1;  
[ph2 hax2] = CIplot(LgLF1Time,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'1-Neighbor');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%---
itemN = 1;
[ph3 hax3] = CIplot(LgLF2Time,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph3],OUTM{:},'2-Neighbors');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%---
itemN = 1;
[ph4 hax4] = CIplot(LgLF3Time,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph4],OUTM{:},'3-Neighbors');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%---
itemN = 1;
[ph5 hax5] = CIplot(LgLF4Time,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph5],OUTM{:},'4-Neighbors');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%---
%------------------------------------------%
set(legend,'Location','NorthWest');
set(legend,'FontSize',16);
legend('off')
%------------------------------------------%
% %------------------------------------------%
% set(get(gca,'XLabel'),'String','Time (min)')
% set(get(gca,'YLabel'),'String','Occupied Slots')
% xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
% set(gca,'XTickLabel', sprintf('%.0f|',xt))
% set(legend,'Location','NorthWest');
% %------------------------------------------%
MS1 = 7; MS2 = 2; clr1=[.1 .1 .1];clr2=[.1 .1 .1];
set(ph1,'LineStyle','-','Color',clr1,'LineWidth',1,...
'Marker','s','MarkerSize',7,'MarkerEdgeColor',clr1,'MarkerFaceColor',clr1);
set(ph2,'LineStyle',':','Color',c2,'LineWidth',4,...
'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(ph3,'LineStyle','-.','Color',c3,'LineWidth',4,...
'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
set(ph4,'LineStyle','-','Color',c4,'LineWidth',4,...
'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
set(ph5,'LineStyle','--','Color',c1,'LineWidth',3,...
'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
hTitle  = title('Cluster vs. Particle Lifetimes');
hXLabel = xlabel('Beta');
hYLabel = ylabel('Lifetime (log steps +/- SEM) ');
set(gca,'FontName','Century Gothic');
set([hTitle, hXLabel, hYLabel],'FontName','Arial');
set([hXLabel, hYLabel],'FontSize',16);
% set( hTitle,'FontSize',18,'FontWeight','bold');
set(hTitle,'FontSize',18);
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
% haxes=axis;
% ylim([0 haxes(4)*1.2 ]);
% xlim([0 (haxes(2))]);
%======================================================================%




%======================================================================%
% FONT NAME TESTERS
%{
%----------------------------------------------------------------------%
fig21 = figure(21);
set(21,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig21,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%------


hTXT = text(.1,.9,0,'Century Gothic Font Name');
set(hTXT,'FontName','Century Gothic','FontSize',18);
% set(hTXT,'FontWeight','bold');

hTXT = text(.1,.8,0,'Book Antiqua Font Name');
set(hTXT,'FontName','Book Antiqua','FontSize',18);
% set(hTXT,'FontWeight','bold');

hTXT = text(.1,.7,0,'Arial Font Name');
set(hTXT,'FontName','Arial','FontSize',18);
% set(hTXT,'FontWeight','bold');

hTXT = text(.1,.6,0,'Microsoft Sans Serif Font Name');
set(hTXT,'FontName','Microsoft Sans Serif','FontSize',18);
% set(hTXT,'FontWeight','bold');

hTXT = text(.1,.5,0,'Calibri Font Name');
set(hTXT,'FontName','Calibri','FontSize',18);
% set(hTXT,'FontWeight','bold');

hTXT = text(.1,.4,0,'Times New Roman Font Name');
set(hTXT,'FontName','Times New Roman','FontSize',18);
% set(hTXT,'FontWeight','bold');

hTXT = text(.1,.3,0,'Verdana Font Name');
set(hTXT,'FontName','Verdana','FontSize',18);
% set(hTXT,'FontWeight','bold');

hTXT = text(.1,.2,0,'xxxxxx Font Name');
set(hTXT,'FontName','xxxxxx','FontSize',18);
% set(hTXT,'FontWeight','bold');

hTXT = text(.1,.1,0,'xxxxxx TT Font Name');
set(hTXT,'FontName','xxxxxx','FontSize',18);
% set(hTXT,'FontWeight','bold');
%}
%----------------------------------------------------------------------%


%}
%=============================%




