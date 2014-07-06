function [output] = EXITSPINE(pbox,um)
format compact;format short;close all;
doLivePlot = 0;
%=========================================================%
% INPUT PARAMETERS FROM GUI

%{.
if 1-exist('pbox','var')
pbox=[30 200 4 0 0.1 1000 50 50 3 2.4 0.1 3 2.4 0.1 0.2 0.2 1 1 2 1000 2 1000 7 7 1 1 1 1];
um=[3 6 .4 .4 .8 .8];
Krate=[2 1000 2 1000];
end
%}
Krate=[2 1000 2 1000];


GraphTime = pbox(1);
AllowedTime = pbox(2);
loops = pbox(3);
doSteadyState = pbox(4);
Scale = pbox(5);			% scale of model
t = pbox(6)/1000;			% time step

GR1Ndots = pbox(7);
GR2Ndots = pbox(8);
GluR1Ndots = GR1Ndots;
GluR2Ndots = GR2Ndots;

DGR1spy = pbox(9);
kGR1 = pbox(10);
DGR1psd = pbox(11);

DGR2spy = pbox(12);
kGR2 = pbox(13);
DGR2psd = pbox(14);

LsGR1psd = pbox(15);
LsGR2psd = pbox(16);


addTime1 = 4/DGR1psd;
addTime2 = 4/DGR2psd;
AllowedTime = AllowedTime+addTime1+addTime2;


%===========================%
% Slots
%---------------------------%  
useGluR1slots = pbox(17);
useGluR2slots = pbox(18);

KonGR1= Krate(1);
KoffGR1= Krate(2);
KonGR2= Krate(3);
KoffGR2= Krate(4);

S1sum= pbox(23)*pbox(23);
S2sum= pbox(24)*pbox(24);
SAP5= [pbox(25) pbox(26) pbox(27) pbox(28)];

GR1slotN = round(S1sum/5*SAP5(1));
GR2slotN = round(S2sum/5*SAP5(3));
GR1slotNo = GR1slotN;
GR2slotNo = GR2slotN;

GluR1_TdwellPSD = zeros(2,GluR1Ndots);
GluR2_TdwellPSD = zeros(2,GluR2Ndots);


%--Decrease Open Slots Subroutine--%
% G2DOSS = 10*round(1/DGR2psd);
% G1DOSS = 10*round(1/DGR2psd);
G1DOSS = 100;
G2DOSS = 100;

stepN = 1;
%===========================%

%%
%=========================================================%
% AMPAR PARTICLE VECTORS
%---------------------------------------------------------%
MaxTime = (AllowedTime*60);     % 200 min = 12,000 s
GraphTime = (GraphTime*60);     % 30 min = 1,800 s
Nsteps = MaxTime;

GR2Sdots = GR2Ndots;			% saves GluR2 Starting particle count
GR1Sdots = GR1Ndots;			% saves GluR1 Starting particle count
%=========================================================%



%===========================================%
% TIME AND PARTICLE COUNT VARIABLES
%-------------------------------------------%
max2graph = GraphTime*2;
exittime = MaxTime * ones(1,GR2Ndots);
exittime2 = MaxTime * ones(1,GR1Ndots);


PSD1n = 0; PSD2n = 0;
PSDTn = PSD1n+PSD2n;
PSD1CaT = 0; PSD2CaT = 0;
%===========================================%


%%
%=========================================================%
% PSD & FIELD MATRIX SETUP
%---------------------------------------------------------%
DENXum = um(1); % DENXum
DENYum = um(2); % DENYum
PSD1um = um(3); % PSD1um
PSD2um = um(4); % PSD2um
PERI1um = um(5); % PERI1um
PERI2um = um(6); % PERI2um


PSD1um = .4; % PSD1um
PSD2um = .4; % PSD2um
PERI1um = .6; % PERI1um
PERI2um = .6; % PERI2um
DENXum = PSD1um+(PERI1um*2);
DENYum = PSD1um+(PERI1um*2);

DENXumS = round(DENXum/Scale);
DENYumS = round(DENYum/Scale);
PSD1umS = round(PSD1um/Scale);
PSD2umS = round(PSD2um/Scale);
PERI1umS = round(PERI1um/Scale);
PERI2umS = round(PERI2um/Scale);


[Yrows Xcols YrPr2bot XcPr2rit YrPr2top XcPr2lft...
YrPr1bot XcPr1rit YrPr1top XcPr1lft PSDfield...
PSDSZE Zfield XWIDE YHIGH...
PSD1WH PSD2WH PERI1WH PERI2WH SPYN1WH SPYN2WH...
XYLTpr1 XYRTpr1 XYLBpr1 XYRBpr1 XYLTpr2 XYRTpr2 XYLBpr2 XYRBpr2...
XYLTp1 XYRTp1 XYLBp1 XYRBp1 XYLTp2 XYRTp2 XYLBp2 XYRBp2...
XYBOXpr1 XYBOXpr2 XYBOXp1 XYBOXp2...
SPYN1xv SPYN1yv SPYN2xv SPYN2yv...
PSD1xv PSD1yv PSD2xv PSD2yv...
PERI1xv PERI1yv PERI2xv PERI2yv...
fPSD1 fPSD2...
] = FieldFun2(DENXumS,DENYumS,PSD1umS,PSD2umS,PERI1umS,PERI2umS);

LFTPSA = XYLBpr1(1); BOTPSA = XYLBpr1(2);
RITPSA = XYRTpr1(1); TOPPSA = XYRTpr1(2);
LFTPSD = XYLBp1(1); BOTPSD = XYLBp1(2);
RITPSD = XYRTp1(1); TOPPSD = XYRTp1(2);

%=========================================================%
% PUT DOTS IN BOX IN RANDOM LOCATIONS
%---------------------------------------------------------%

GluR2x = randi([(LFTPSA+2),(RITPSA-2)],[1,GR2Ndots]);
GluR2y = randi([(BOTPSA+2),(TOPPSA-2)],[1,GR2Ndots]);
GluR2xy = [GluR2x;GluR2y];

GluR1x = randi([(LFTPSA+2),(RITPSA-2)],[1,GR1Ndots]);
GluR1y = randi([(BOTPSA+2),(TOPPSA-2)],[1,GR1Ndots]);
GluR1xy = [GluR1x;GluR1y];

GluR2xyl = GluR2xy;
GluR1xyl = GluR1xy;

GluR2xyds = preXYDSGR2(GR2Ndots, kGR2);
GluR1xyds = preXYDSGR1(GR1Ndots, kGR1);

[GluR2xyl] = preMOVEGLUR2(stepN,GluR2xyds,GluR2xyl,XYLBpr1,XYRTpr1);
[GluR1xyl] = preMOVEGLUR1(stepN,GluR1xyds,GluR1xyl,XYLBpr1,XYRTpr1);

preMAINPLOT(GluR1xyl,GluR2xyl,...
	XYLBpr1,XYRTpr1,XYLBp1,XYRTp1,...
	PSD1WH,PSD2WH,PERI1WH,PERI2WH,SPYN1WH,SPYN2WH,XWIDE,YHIGH);


%=========================================================%
% DO STEADY STATE - RUN SIM (NO EXIT) FOR 100 STEPS
%---------------------------------------------------------%
if doSteadyState
stepN = 1;
for Nt = 1:100
	  
GluR2xyds = preXYDSGR2(GR2Ndots, kGR2);
GluR1xyds = preXYDSGR1(GR1Ndots, kGR1);

[GluR2xyl] = preMOVEGLUR2(stepN,GluR2xyds,GluR2xyl,XYLBpr1,XYRTpr1);
[GluR1xyl] = preMOVEGLUR1(stepN,GluR1xyds,GluR1xyl,XYLBpr1,XYRTpr1);

preMAINPLOT(GluR1xyl,GluR2xyl,...
	XYLBpr1,XYRTpr1,XYLBp1,XYRTp1,...
	PSD1WH,PSD2WH,PERI1WH,PERI2WH,SPYN1WH,SPYN2WH,XWIDE,YHIGH);

stepN = stepN+1;
end
end %if doSteadyState
%=========================================================%

PreGluR2xyl = GluR2xyl;
PreGluR1xyl = GluR1xyl;

TexitGR1 = [];
TexitGR2 = [];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					MAIN LOOPS START
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%============================================================%
%						GLUR2 LOOP
%============================================================%


TexitGR2all= zeros(loops,GR2Ndots);

%====================%
for re = 1:loops
%====================%
	GluR2xyl = PreGluR2xyl;
	TexitGR2 = [];
	GR2slotN = GR2slotNo;
	GR2RSOT = [];
	GluR2_TdwellPSD = zeros(2,GluR2Ndots);
%===============%
stepN = 1;
for Nt = 1:Nsteps
%===============%
    
  
    GluR2xyds = XYDSGR2(GR2Ndots, kGR2);
	
	
	if useGluR2slots
	[GluR2xyds GluR2_TdwellPSD G2FSLOTS] = SLOTSG2(GluR2xyds,GluR2_TdwellPSD,...
    KonGR2,KoffGR2,GR2slotN);
	end
	
	
	
	%===============================%
	GR2dotsLeft = size(GluR2xyl,2);
	if GR2dotsLeft >=1
		
	G2INPSD = INBOXFun(XYLBp1,XYRTp1,GluR2xyl);
	G2INPSA = INBOXFun(XYLBpr1,XYRTpr1,GluR2xyl);
	
	[GluR2xyl TexitGR2 GluR2_TdwellPSD] =...
	MOVEGLUR2(stepN,GluR2xyds,GluR2xyl,GluR2_TdwellPSD,LsGR2psd,...
	TexitGR2,t,XYLBpr1,XYRTpr1,G2INPSD,G2INPSA,GR2dotsLeft);
	
	end %if GR2dotsLeft >=1
	%===============================%
	% if stepN == 20; keyboard; end
	
	
	
	%--Decrease Open Slots Subroutine--%
	if useGluR2slots
	if mod(stepN, G2DOSS) == 0
	%-----------------------------
		% Proportion of new dots
		GR2PoND = numel(TexitGR2)/GR2Ndots;
	
		% Number of open slots
		GR2NoOS = GR2slotN - G2FSLOTS; 
	
		% Most likely number of open slots filled by new dots
		NOS=0:GR2NoOS;
		y = binopdf(NOS,GR2NoOS,GR2PoND);
		[x,i]=max(y);
		NOSFND = NOS(i);
	
		GR2slotN=GR2slotN-NOSFND; % current number of slots
			% if GR2slotN<2;GR2slotN=1;end
			
		GR2RSOT(:,stepN/G2DOSS) = [stepN;GR2slotN]; % remaining slots over time
	%-----------------------------	
	% if mod(stepN, 2002) == 0; keyboard; end
	
	if stepN<=6000
	SOTylim = [0 GR2slotNo+5];
	SOTxlim = [0 Nsteps/2];

	SLOTPLOT2(GluR2xyl,...
	GR2RSOT,SOTxlim,SOTylim,stepN,G2DOSS,...
	XYLBpr1,XYRTpr1,XYLBp1,XYRTp1,...
	PSD1WH,PSD2WH,PERI1WH,PERI2WH,SPYN1WH,SPYN2WH,XWIDE,YHIGH)

	end
	
	end
	end %if useGluR2slots
	%--Decrease Open Slots Subroutine--%
	
	
	
  %===============%	
  stepN = stepN+1;
  if mod(stepN, 1000) == 0
  stepN
  end
%===============%
end %for Nt = 1:Nsteps
%===============%


TexitGR2all(re,:) = TexitGR2;
if useGluR2slots
GR2RSOTall(re,:) = GR2RSOT(2,:);
end %if useGluR2slots
%====================%
end %for re = 1:loops
%====================%
TexitGR2 = mean(TexitGR2all,1);
if useGluR2slots
GR2RSOTMEAN = mean(GR2RSOTall,1);
GR2RSOTMEANT = padarray(GR2RSOTMEAN,[1 0],0,'pre');
GR2RSOTMEANT(1,:) = GR2RSOT(1,:);
GLURX = 2;
SLOTPLOTEND(GR2RSOTMEANT,SOTxlim,SOTylim,GLURX);
end %if useGluR2slots
%============================================================%
%						GLUR1 LOOP
%============================================================%

TexitGR1all= zeros(loops,GR1Ndots);

for re = 1:loops
	GluR1xyl = PreGluR1xyl;
	TexitGR1 = [];
	GR1slotN = GR1slotNo;
	GR1RSOT = [];
	GluR1_TdwellPSD = zeros(2,GluR1Ndots);
%===============%
stepN = 1;
for Nt = 1:Nsteps
%===============%

  
    GluR1xyds = XYDSGR1(GR1Ndots, kGR1);
	
	if useGluR1slots
	[GluR1xyds GluR1_TdwellPSD G1FSLOTS] = SLOTSG1(GluR1xyds,GluR1_TdwellPSD,...
	KonGR1,KoffGR1,GR1slotN);
	end
    
    %===============================%
	GR1dotsLeft = size(GluR1xyl,2);
	if GR1dotsLeft >=1
		
	G1INPSD = INBOXFun(XYLBp1,XYRTp1,GluR1xyl);
	G1INPSA = INBOXFun(XYLBpr1,XYRTpr1,GluR1xyl);
	
	[GluR1xyl TexitGR1 GluR1_TdwellPSD] =...
    MOVEGLUR1(stepN,GluR1xyds,GluR1xyl,GluR1_TdwellPSD,LsGR1psd,...
    TexitGR1,t,XYLBpr1,XYRTpr1,G1INPSD,G1INPSA,GR1dotsLeft);
	
	end %if GR2dotsLeft >=1
	%===============================%
	% if stepN == 20; keyboard; end
	
	
	
	%--Decrease Open Slots Subroutine--%
	if useGluR1slots
	if mod(stepN, G1DOSS) == 0
	%-----------------------------
		% Proportion of new dots
		GR1PoND = numel(TexitGR1)/GR1Ndots;
	
		% Number of open slots
		GR1NoOS = GR1slotN - G1FSLOTS; 
	
		% Most likely number of open slots filled by new dots
		NOS=0:GR1NoOS;
		y = binopdf(NOS,GR1NoOS,GR1PoND);
		[x,i]=max(y);
		NOSFND = NOS(i);
	
		GR1slotN=GR1slotN-NOSFND; % current number of slots
			% if GR1slotN<1;GR1slotN=1;end
		
		GR1RSOT(:,stepN/G1DOSS) = [stepN;GR1slotN];	% remaining slots over time
	%-----------------------------
	
	
	if stepN<=6000
	SOTylim = [0 GR1slotNo+5];
	SOTxlim = [0 Nsteps/2];

	SLOTPLOT2(GluR1xyl,...
	GR1RSOT,SOTxlim,SOTylim,stepN,G1DOSS,...
	XYLBpr1,XYRTpr1,XYLBp1,XYRTp1,...
	PSD1WH,PSD2WH,PERI1WH,PERI2WH,SPYN1WH,SPYN2WH,XWIDE,YHIGH);
	end
	
	end
	end %if useGluR1slots
	%--Decrease Open Slots Subroutine--%
	
	
	
	
	
  %===============%	
  stepN = stepN+1;
  if mod(stepN, 1000) == 0
  stepN
  end
%===============%
end %for Nt = 1:Nsteps
%===============%

TexitGR1all(re,:) = TexitGR1;
if useGluR1slots
GR1RSOTall(re,:) = GR1RSOT(2,:);
end %if useGluR1slots
%====================%
end %for re = 1:loops
%====================%
TexitGR1 = mean(TexitGR1all,1);
if useGluR1slots
GR1RSOTMEAN = mean(GR1RSOTall,1);
GR1RSOTMEANT = padarray(GR1RSOTMEAN,[1 0],0,'pre');
GR1RSOTMEANT(1,:) = GR1RSOT(1,:);
GLURX = 1;
SLOTPLOTEND(GR1RSOTMEANT,SOTxlim,SOTylim,GLURX);
end %if useGluR1slots
%---------------------------------------------------------%
%					MAIN LOOPS END
%=========================================================%

exittime = TexitGR2;
exittime2 = TexitGR1;

GluR2exT30 = [exittime];
GluR2exT30 = sort(GluR2exT30);

GluR1exT30 = [exittime2];
GluR1exT30 = sort(GluR1exT30);

GluR2exT30(numel(GluR2exT30)) = []
GluR1exT30(numel(GluR1exT30)) = []


GluR1last = GluR1exT30(find(GluR1exT30,1,'last'));
GluR2last = GluR2exT30(find(GluR2exT30,1,'last'));




%=========================================================%
%					FINAL GRAPHIC OUTPUTS
%=========================================================%
scD = DGR1spy*Scale;
scPSD1 = DGR1psd*Scale;
scPSD2 = DGR2psd*Scale;

GraphT = GraphTime/t;
xlim = [0 GraphT];
ylim = [0 1.1];


[coECDF_G2, coX_G2] = ecdf(exittime);
[coECDF_G1, coX_G1] = ecdf(exittime2);




figure(69)
set(gcf,'OuterPosition',[400,100,800,800])
cdfexit = gcf;
figure(cdfexit)
subplot(5,4,[1 12]),GluR2plot = cdfplot(GluR2exT30);
axis([xlim,ylim]);
xt = (get(gca,'XTick'))*(t)/(60);
set(GluR2plot,'color',[1 0 1])
rectangle('Position',[GluR2last,1,GraphT,.001],'EdgeColor',[1 0 1],'LineStyle','--')
hold on
subplot(5,4,[1 12]),cdfplot(GluR1exT30);
set(gca,'XTickLabel', sprintf('%.1f|',xt))
text((GraphT/2),.05,'Minute Exited','FontSize',12)
rectangle('Position',[GluR1last,1,GraphT,.001],'EdgeColor','b','LineStyle','--')
CDFtitle = title(['CDF Exited'...
	'    D: SPI = ' num2str(scD)...
	'    GR2 = ' num2str(scPSD2)...
	'    GR1 = ' num2str(scPSD1)...
	'  \bullet  N: GR2 =' int2str(GR2Sdots)...
	'    GR1 =' int2str(GR1Sdots)]);
set(CDFtitle, 'FontSize', 12);
leg1=legend('  GluR2', 'GluR1');
set(leg1,'Location','SouthEast');
figure(cdfexit)
subplot(5,4,[13 14]),plot(sort(exittime),1:numel(exittime));
title('GluR2 Exit Times');
subplot(5,4,[15 16]),ecdfhist(coECDF_G2, coX_G2)
title('GluR2 CDF Histogram');
subplot(5,4,[17 18]),plot(sort(exittime2),1:numel(exittime2));
set(get(gca,'XLabel'),'String','Time-Step')
title('GluR1 Exit Times');
subplot(5,4,[19 20]),hist(exittime2);%ecdfhist(coX_G1, coX_G1)
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Exit Minute')
title('GluR1 CDF Histogram');








for k=1:2, varargout(k) = {eval(['GluR' int2str(k) 'exT30'])}; end
output = varargout;
%-------------############################------------------%
end %		  ##    END MAIN FUNCTION   ##
%-------------############################------------------% 
%%



%=========================================================%
%					FUNCTION SUBROUTINES
%=========================================================%
%-------------------------------------------%
% GLUR1 BROWNIAN STEP GENERATOR
%-------------------------------------------%
function GluR1xyds = preXYDSGR1(GR1Ndots, k)
    GluR1xyds = (k * randn(2,GR1Ndots));
end
 
%-------------------------------------------%
% GLUR2 BROWNIAN STEP GENERATOR
%-------------------------------------------%
function GluR2xyds = preXYDSGR2(GR2Ndots, k)
    GluR2xyds = (k * randn(2,GR2Ndots));
end

%-------------------------------------------%
% GLUR1 BROWNIAN STEP GENERATOR
%-------------------------------------------%
function GluR1xyds = XYDSGR1(GR1Ndots, k)
    GluR1xyds = (k * randn(2,GR1Ndots));
end
 
%-------------------------------------------%
% GLUR2 BROWNIAN STEP GENERATOR
%-------------------------------------------%
function GluR2xyds = XYDSGR2(GR2Ndots, k)
    GluR2xyds = (k * randn(2,GR2Ndots));
end




%===================================%
% MOVE PARTICLES MAIN FUNCTION
%===================================%
function [GluR2xyl TexitGR2 GluR2_TdwellPSD] =...
	MOVEGLUR2(stepN,GluR2xyds,GluR2xyl,GluR2_TdwellPSD,LsGR2psd,...
	TexitGR2,t,XYLBpr1,XYRTpr1,G2INPSD,G2INPSA,GR2dotsLeft)


GluR2xyds(:,G2INPSD) = GluR2xyds(:,G2INPSD)*(LsGR2psd);
GluR2_TdwellPSD(1,G2INPSD) = GluR2_TdwellPSD(1,G2INPSD)+t;
GluR2_TdwellPSD(1,~G2INPSD) = 0;

for i = 1:2
for j = 1:GR2dotsLeft
	GluR2xyl(i,j) = GluR2xyl(i,j)+GluR2xyds(i,j);
end
end

%{
%===================================%
GluR1xyl = GluR1xyl+GluR1xyds;
 
for j = 1:GluR1Ndots 
    if GluR1xyl(1,j)>(XWIDE) || GluR1xyl(1,j)<(0)
            GluR1xyl(1,j) = uint8(sign(GluR1xyl(1,j)))*(XWIDE);
    elseif GluR1xyl(2,j)>(YHIGH) || GluR1xyl(2,j)<(-YHIGH)
            GluR1xyl(2,j) = sign(GluR1xyl(2,j))*(YHIGH);
    end    
        
end
%===================================%
GluR2xyl = GluR2xyl+GluR2xyds;

for j = 1:GluR2Ndots 
	if GluR2xyl(1,j)>(XWIDE) || GluR2xyl(1,j)<(0)
            GluR2xyl(1,j) = uint8(sign(GluR2xyl(1,j)))*(XWIDE);
	elseif GluR2xyl(2,j)>(YHIGH) || GluR2xyl(2,j)<(-YHIGH)
            GluR2xyl(2,j) = sign(GluR2xyl(2,j))*(YHIGH);
	end	   
        
end
%}

LFT = XYLBpr1(1); BOT = XYLBpr1(2);
RIT = XYRTpr1(1); TOP = XYRTpr1(2);
%===========================%
    for j = 1:GR2dotsLeft  
		if GluR2xyl(1,j)<=(LFT)
			GluR2xyl(1,j) = LFT;
		end
		
		if GluR2xyl(1,j)>=(RIT)
			GluR2xyl(1,j) = RIT;
		end
		
		if GluR2xyl(1,j)>=(RIT) && GluR2xyl(2,j)<=(TOP+4)
			GluR2xyl(1,j) = 111;
			TexitGR2(numel(TexitGR2)+1) = stepN;
		end
		
	end
	
	for j = 1:GR2dotsLeft
		if GluR2xyl(2,j)<=(BOT)
			GluR2xyl(2,j) = BOT;
		end
		
		if GluR2xyl(2,j)>=(TOP)
			GluR2xyl(2,j) = TOP;
		end
	end
%===========================%

delthis = find(GluR2xyl(1,:)>100);
GluR2xyl(:,delthis) = [];
% countdots = size(GluR2xyl,2); 

GluR2_TdwellPSD(:,delthis) = [];



%===================================%
% if stepN == 500; keyboard; end
end
%===================================%
%===================================%
function [GluR1xyl TexitGR1 GluR1_TdwellPSD] =...
    MOVEGLUR1(stepN,GluR1xyds,GluR1xyl,GluR1_TdwellPSD,LsGR1psd,...
    TexitGR1,t,XYLBpr1,XYRTpr1,G1INPSD,G1INPSA,GR1dotsLeft)
 
 
GluR1xyds(:,G1INPSD) = GluR1xyds(:,G1INPSD)*(LsGR1psd);
GluR1_TdwellPSD(1,G1INPSD) = GluR1_TdwellPSD(1,G1INPSD)+t;
GluR1_TdwellPSD(1,~G1INPSD) = 0;
 
for i = 1:2
for j = 1:GR1dotsLeft
    GluR1xyl(i,j) = GluR1xyl(i,j)+GluR1xyds(i,j);
end
end
 
%{
%===================================%
GluR1xyl = GluR1xyl+GluR1xyds;
 
for j = 1:GluR1Ndots 
    if GluR1xyl(1,j)>(XWIDE) || GluR1xyl(1,j)<(0)
            GluR1xyl(1,j) = uint8(sign(GluR1xyl(1,j)))*(XWIDE);
    elseif GluR1xyl(2,j)>(YHIGH) || GluR1xyl(2,j)<(-YHIGH)
            GluR1xyl(2,j) = sign(GluR1xyl(2,j))*(YHIGH);
    end    
        
end
%===================================%
GluR1xyl = GluR1xyl+GluR1xyds;
 
for j = 1:GluR1Ndots 
    if GluR1xyl(1,j)>(XWIDE) || GluR1xyl(1,j)<(0)
            GluR1xyl(1,j) = uint8(sign(GluR1xyl(1,j)))*(XWIDE);
    elseif GluR1xyl(2,j)>(YHIGH) || GluR1xyl(2,j)<(-YHIGH)
            GluR1xyl(2,j) = sign(GluR1xyl(2,j))*(YHIGH);
    end    
        
end
%}
 
LFT = XYLBpr1(1); BOT = XYLBpr1(2);
RIT = XYRTpr1(1); TOP = XYRTpr1(2);
%===========================%
    for j = 1:GR1dotsLeft  
        if GluR1xyl(1,j)<=(LFT)
            GluR1xyl(1,j) = LFT;
        end
        
        if GluR1xyl(1,j)>=(RIT)
            GluR1xyl(1,j) = RIT;
        end
        
        if GluR1xyl(1,j)>=(RIT) && GluR1xyl(2,j)<=(TOP+4)
            GluR1xyl(1,j) = 111;
            TexitGR1(numel(TexitGR1)+1) = stepN;
        end
        
    end
    
    for j = 1:GR1dotsLeft
        if GluR1xyl(2,j)<=(BOT)
            GluR1xyl(2,j) = BOT;
        end
        
        if GluR1xyl(2,j)>=(TOP)
            GluR1xyl(2,j) = TOP;
        end
    end
%===========================%
 
delthis = find(GluR1xyl(1,:)>100);
GluR1xyl(:,delthis) = [];
% countdots = size(GluR1xyl,2); 
 
GluR1_TdwellPSD(:,delthis) = [];
 
 
 
%===================================%
% if stepN == 500; keyboard; end
end
%===================================%



%-------------------------------------------%
% INBOXFun
%-------------------------------------------%
function [inbox] = INBOXFun(LB,RT,xyl)

%{
inboxfun(LB,RT,xyl)

The inboxfun() function takes the right top (RT) and left bottom (LB)
X Y coordinates of a square in the form of a 1-row 2-col vector.
It also takes a 2xN matrix of point locations 
and determines if those points are inside the square

xyl = [x x x ...Nx;
	   y y y ...Ny];

EXAMPLE:

RT = [2 2];
LB = [0 0];
xyl = [0:4;0:4];

inbox = inboxfun(LB,RT,varargin)

inbox =
     0     1     0     0     0




if ~nargin
	RT = [2 2];
	LB = [0 0];
	xyl = [0:4;0:4];
end

%}

if LB(1)>RT(1)
	LBt=LB;
	RTt=RT;
	LB(1)=RTt(1);
	RT(1)=LBt(1);
end
if LB(2)>RT(2)
	LBt=LB;
	RTt=RT;
	LB(2)=RTt(2);
	RT(2)=LBt(2);
end


sz = numel(xyl(1,:));
inbox = zeros(1,sz);

for s = 1:sz
	
	if xyl(1,s) > LB(1) && xyl(1,s) < RT(1) &&...
	   xyl(2,s) > LB(2) && xyl(2,s) < RT(2)
		
		inbox(s) = 1;
	end
   
end

inbox = inbox>0;

end



%-------------------------------------------%
% PLOT Particle Motion
%-------------------------------------------%
function [] = SLOTPLOT2(GLURXYL,...
	SOT,SOTxlim,SOTylim,stepN,DOSS,...
	XYLBpr1,XYRTpr1,XYLBp1,XYRTp1,...
	PSD1WH,PSD2WH,PERI1WH,PERI2WH,SPYN1WH,SPYN2WH,XWIDE,YHIGH)
%-------------------------------------------%


%=================================%
%       SLOT LINE PLOT
%=================================%
hold off
figure(1)
subplot(3,3,[7 9]),
plot(SOT(1,:),SOT(2,:))
axis manual
axis([SOTxlim SOTylim])
%=================================%

%=================================%
%       MAIN 2D PLOT
%---------------------------------%
LFTPSA = XYLBpr1(1); BOTPSA = XYLBpr1(2);
RITPSA = XYRTpr1(1); TOPPSA = XYRTpr1(2);
LFTPSD = XYLBp1(1); BOTPSD = XYLBp1(2);
RITPSD = XYRTp1(1); TOPPSD = XYRTp1(2);

%---
xlim = [-2 (XWIDE+2)];
ylim = [-2 (XWIDE+2)];
%---
%---------------------------------%
%     MAIN DOTS
%----------------------%
figure(1)
hold off
subplot(3,3,[1 6]),
AMPARPlot = gscatter(GLURXYL(1,:),GLURXYL(2,:));
axis([xlim, ylim]);
set(AMPARPlot,'marker','.','markersize',[6],'color',[1 0 0])
rectangle('Position',[LFTPSA,BOTPSA,SPYN1WH,SPYN1WH])
rectangle('Position',[LFTPSD,BOTPSD,PSD1WH,PSD1WH])
rectangle('Position',[RITPSA,BOTPSA,1,4],'Curvature',[.2,.1],...
	'EdgeColor','w','FaceColor','w','LineWidth',1,'LineStyle','--')
%======================%
hold off



end




%-------------------------------------------%
% SLOTPLOTEND
%-------------------------------------------%
function [] = SLOTPLOTEND(SOT,SOTxlim,SOTylim,GLURX)
%-------------------------------------------%

if GLURX == 1
figure(43)
subplot(2,1,1),
plot(SOT(1,:),SOT(2,:))
axis manual
axis([SOTxlim SOTylim])
title('BLEACHED GluR1 IN SLOTS');
end

if GLURX == 2
figure(43)
subplot(2,1,2),
plot(SOT(1,:),SOT(2,:))
axis manual
axis([SOTxlim SOTylim])
title('BLEACHED GluR2 IN SLOTS');
end

end


%-------------------------------------------%
% PLOT Particle Motion
%-------------------------------------------%
function [] = MAINPLOT(GluR1xyl, GluR2xyl, Fcol, PSDLOC, PSDSZE,...
	row1F,row1L,col1F,col1L)
%-------------------------------------------%
xlim = [0 Fcol];
ylim = [0 Fcol];
%---

%---
PSDLOC1 = PSDLOC;
PSDLOC1 = [PSDLOC1(2,:); PSDLOC1(1,:)];

SYN1SZ = (PSDSZE(1,1) + (PSDSZE(2,1)*2))-1;
PSD1SZ = (PSDSZE(1,1))-1;

P1INs = PSDSZE(2,1);
P1INx = (col1L);
P1INy = (row1F);

P1OUTx=PSDLOC(2,1);
P1OUTy=PSDLOC(1,1);

%---

%=================================%
%       MAIN 2D PLOT
%---------------------------------%
%     MAIN DOTS
%----------------------%
figure(1) 
AMPARPlot = gscatter(GluR2xyl(1,:),GluR2xyl(2,:));
axis([xlim, ylim]);
set(AMPARPlot,'marker','.','markersize',[6],'color',[1 0 0])
hold on
AMPARPlot = gscatter(GluR1xyl(1,:),GluR1xyl(2,:));
axis([xlim, ylim]);
set(AMPARPlot,'marker','.','markersize',[6],'color',[0 1 1])
hold off


rectangle('Position',[P1OUTx,P1OUTy,SYN1SZ,SYN1SZ])
rectangle('Position',[P1INx,P1INy,1,4],'Curvature',[.2,.1],...
	'EdgeColor','w','FaceColor','w','LineWidth',1,'LineStyle','--')
%======================%

end



%%			  DENDRITIC FIELD STEP TOOLS
%-------------##########################------------------%
%		FIELD MAP FOR DENDRITE AREA AND PSD AREAS
%-------------##########################------------------%
function [Frow Fcol row1L col1L row1F col1F PSDfield...
PSDLOC PSDSZE] = FieldFun(fsizeX, fsizeY, PSD1size, periPSD1size)

Syn1Size = PSD1size+(periPSD1size*2);

PSD1padX = round((fsizeX - Syn1Size)/2);
PSD1padY = round(((fsizeY - Syn1Size)/2)/2);

fPSD1 = ones(PSD1size);						% PSD1 SIZE
fPSD1 = padarray(fPSD1,[periPSD1size periPSD1size], 2);	% PSD1 SIZE
pfPSD1 = padarray(fPSD1,[PSD1padY PSD1padX], 0);	% PAD PSD1 [Y-rows X-cols]

[pfPSD1R pfPSD1C] = size(pfPSD1);


PSDfield = pfPSD1;					% CONCAT PSD FIELDS
% figure(5); imagesc(pfPSD1)

[Frow Fcol] = size(PSDfield);           % Fcol: #cols | Frow: #cols x2
[row1F col1F] = find(pfPSD1,1,'first'); % PSD1 1st row & col
[row1L col1L] = find(pfPSD1,1,'last');  % PSD1 Lst row & col
P1TL = [row1F col1F]';					% PSD1 top R location
P1BR = [row1L col1L]';					% PSD1 btm L location
PSDLOC = [P1TL P1BR];					% ALL PSD TL BR locations
PSDSZE = [PSD1size PSD1size; periPSD1size periPSD1size];


end



%%			4. DENDRITIC FIELD MAPPING TOOLS
%-------------##########################------------------%
%		FIELD MAP FOR DENDRITE AREA AND PSD AREAS
%-------------##########################------------------%
% FieldFun
function [Yrows Xcols YrPr2bot XcPr2rit YrPr2top XcPr2lft...
YrPr1bot XcPr1rit YrPr1top XcPr1lft PSDfield...
PSDSZE Zfield XWIDE YHIGH...
PSD1WH PSD2WH PERI1WH PERI2WH SPYN1WH SPYN2WH...
XYLTpr1 XYRTpr1 XYLBpr1 XYRBpr1 XYLTpr2 XYRTpr2 XYLBpr2 XYRBpr2...
XYLTp1 XYRTp1 XYLBp1 XYRBp1 XYLTp2 XYRTp2 XYLBp2 XYRBp2...
XYBOXpr1 XYBOXpr2 XYBOXp1 XYBOXp2...
SPYN1xv SPYN1yv SPYN2xv SPYN2yv...
PSD1xv PSD1yv PSD2xv PSD2yv...
PERI1xv PERI1yv PERI2xv PERI2yv...
fPSD1 fPSD2...
] = FieldFun2(fsizeX,fsizeY,PSD1size,PSD2size,periPSD1size,periPSD2size)


Syn1Size = PSD1size+(periPSD1size*2);
Syn2Size = PSD2size+(periPSD2size*2);

PSD1padX = round((fsizeX - Syn1Size)/2);
PSD1padY = round(((fsizeY - Syn1Size)/2)/2);
PSD2padX = round((fsizeX - Syn2Size)/2);
PSD2padY = round(((fsizeY - Syn2Size)/2)/2);

fPSD1 = ones(PSD1size+1);						% PSD1 SIZE
fPSD2 = ones(PSD2size+1);						% PSD2 SIZE
fPSD1 = padarray(fPSD1,[periPSD1size periPSD1size], 2);	% PSD1 SIZE
fPSD2 = padarray(fPSD2,[periPSD2size periPSD2size], 2);	% PSD2 SIZE
pfPSD1 = padarray(fPSD1,[PSD1padY PSD1padX], 0);	% PAD PSD1 [Y-rows X-cols]
pfPSD2 = padarray(fPSD2,[PSD2padY PSD2padX], 0);	% PAD PSD2 [Y-rows X-cols]

[pfPSD1R pfPSD1C] = size(pfPSD1);
[pfPSD2R pfPSD2C] = size(pfPSD2);


if pfPSD1C > pfPSD2C
		pfPSD2 = padarray(pfPSD2,[0 (pfPSD1C - pfPSD2C)], 0,'post');
elseif pfPSD1C < pfPSD2C
		pfPSD1 = padarray(pfPSD1,[0 (pfPSD2C - pfPSD1C)], 0,'post');
end


%=========================================================%
PSDfield = cat(1, pfPSD1, pfPSD2);				% CONCAT PSD FIELDS
[Yrows Xcols] = size(PSDfield);					% Even value, total field
[YrPr1top XcPr1lft] = find(pfPSD1,1,'first');	% PSD1 1st row & col
[YrPr1bot XcPr1rit] = find(pfPSD1,1,'last');	% PSD1 Lst row & col
[YrPr2top XcPr2lft] = find(pfPSD2,1,'first');	% PSD2 1st row & col
[YrPr2bot XcPr2rit] = find(pfPSD2,1,'last');	% PSD2 Lst row & col
PSDSZE = [PSD1size PSD2size; periPSD1size periPSD2size];
YrowsHIGH = Yrows/2;							% Rows high from Y=0
XcolsWIDE = Xcols;								% Cols wide from X=0
YHIGH = YrowsHIGH;								% Rows high from Y=0
XWIDE = XcolsWIDE;								% Cols wide from X=0
%=========================================================%
Y1RowsUpTop = (YrowsHIGH-YrPr1top);
Y1RowsUpBot = (YrowsHIGH-YrPr1bot);
X1ColsOvrL = XcPr1lft;
X1ColsOvrR = XcPr1rit;

Y2RowsDnTop = (YrowsHIGH-YrPr2bot)*-1;
Y2RowsDnBot = (YrowsHIGH-YrPr2top)*-1;
X2ColsOvrL = XcPr2lft;
X2ColsOvrR = XcPr2rit;
%=========================================================%
PSD1WH = PSD1size;
PSD2WH = PSD2size;
PERI1WH = periPSD1size;
PERI2WH = periPSD2size;
SPYN1WH = PSD1WH+(PERI1WH*2);
SPYN2WH = PSD2WH+(PERI2WH*2);
%=========================================================%
XYLTpr1 = [X1ColsOvrL Y1RowsUpTop];
XYRTpr1 = [X1ColsOvrR Y1RowsUpTop];
XYLBpr1 = [X1ColsOvrL Y1RowsUpBot];
XYRBpr1 = [X1ColsOvrR Y1RowsUpBot];
XYLTpr2 = [X2ColsOvrL Y2RowsDnTop];
XYRTpr2 = [X2ColsOvrR Y2RowsDnTop];
XYLBpr2 = [X2ColsOvrL Y2RowsDnBot];
XYRBpr2 = [X2ColsOvrR Y2RowsDnBot];
XYLTp1 = [(XYLTpr1(1,1)+PERI1WH) (XYLTpr1(1,2)-PERI1WH)];
XYRTp1 = [(XYRTpr1(1,1)-PERI1WH) (XYRTpr1(1,2)-PERI1WH)];
XYLBp1 = [(XYLBpr1(1,1)+PERI1WH) (XYLBpr1(1,2)+PERI1WH)];
XYRBp1 = [(XYRBpr1(1,1)-PERI1WH) (XYRBpr1(1,2)+PERI1WH)];
XYLTp2 = [(XYLTpr2(1,1)+PERI2WH) (XYLTpr2(1,2)-PERI2WH)];
XYRTp2 = [(XYRTpr2(1,1)-PERI2WH) (XYRTpr2(1,2)-PERI2WH)];
XYLBp2 = [(XYLBpr2(1,1)+PERI2WH) (XYLBpr2(1,2)+PERI2WH)];
XYRBp2 = [(XYRBpr2(1,1)-PERI2WH) (XYRBpr2(1,2)+PERI2WH)];
%=========================================================%
% rectangle('Position',[x,y,w,h])
XYBOXpr1 = [XYLBpr1 SPYN1WH SPYN1WH];
XYBOXpr2 = [XYLBpr2 SPYN2WH SPYN2WH];
XYBOXp1 = [XYLBp1 PSD1WH PSD1WH];
XYBOXp2 = [XYLBp2 PSD2WH PSD2WH];
%=========================================================%
Zfield = changem(PSDfield,[1 2 3],[0 2 1]);
%=========================================================%
xlim = [0 XcolsWIDE];
ylim = [-YrowsHIGH YrowsHIGH];
%---
SNSZ = get(0,'ScreenSize');
fig99 = figure(99);

figure(fig99)
subplot(5,5,[3 20]), subplot('Position',[.45 .40 .5 .59]),...
scatter(1:4,1:4,5,[1 0 0]);
axis([xlim, ylim]);
grid on
set(get(gca,'XLabel'),'String','GluR2:Blue \bullet \bullet GluR1:Red')
rectangle('Position',[XYBOXpr1(1),XYBOXpr1(2),XYBOXpr1(3),XYBOXpr1(4)])
rectangle('Position',[XYBOXpr2(1),XYBOXpr2(2),XYBOXpr2(3),XYBOXpr2(4)])
rectangle('Position',[XYBOXp1(1),XYBOXp1(2),XYBOXp1(3),XYBOXp1(4)])
rectangle('Position',[XYBOXp2(1),XYBOXp2(2),XYBOXp2(3),XYBOXp2(4)])
%---
text(XYLTpr1(1),XYLTpr1(2),...
strcat(num2str(XYLTpr1(1)), '\bullet',num2str(XYLTpr1(2)),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(XYRTpr1(1),XYRTpr1(2),...
strcat('\leftarrow ', num2str(XYRTpr1(1)), '\bullet',num2str(XYRTpr1(2))),...
'FontSize',12,'HorizontalAlignment','left');
text(XYLBpr1(1),XYLBpr1(2),...
strcat(num2str(XYLBpr1(1)), '\bullet',num2str(XYLBpr1(2)),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(XYRBpr1(1),XYRBpr1(2),...
strcat('\leftarrow ', num2str(XYRBpr1(1)), '\bullet',num2str(XYRBpr1(2))),...
'FontSize',12,'HorizontalAlignment','left');
%---
text(XYLTpr2(1),XYLTpr2(2),...
strcat(num2str(XYLTpr2(1)), '\bullet',num2str(XYLTpr2(2)),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(XYRTpr2(1),XYRTpr2(2),...
strcat('\leftarrow ', num2str(XYRTpr2(1)), '\bullet',num2str(XYRTpr2(2))),...
'FontSize',12,'HorizontalAlignment','left');
text(XYLBpr2(1),XYLBpr2(2),...
strcat(num2str(XYLBpr2(1)), '\bullet',num2str(XYLBpr2(2)),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(XYRBpr2(1),XYRBpr2(2),...
strcat('\leftarrow ', num2str(XYRBpr2(1)), '\bullet',num2str(XYRBpr2(2))),...
'FontSize',12,'HorizontalAlignment','left');
%---
text(XYLTp1(1),XYLTp1(2),...
strcat(num2str(XYLTp1(1)), '\bullet',num2str(XYLTp1(2)),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(XYRTp1(1),XYRTp1(2),...
strcat('\leftarrow ', num2str(XYRTp1(1)), '\bullet',num2str(XYRTp1(2))),...
'FontSize',12,'HorizontalAlignment','left');
text(XYLBp1(1),XYLBp1(2),...
strcat(num2str(XYLBp1(1)), '\bullet',num2str(XYLBp1(2)),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(XYRBp1(1),XYRBp1(2),...
strcat('\leftarrow ', num2str(XYRBp1(1)), '\bullet',num2str(XYRBp1(2))),...
'FontSize',12,'HorizontalAlignment','left');
%---
text(XYLTp2(1),XYLTp2(2),...
strcat(num2str(XYLTp2(1)), '\bullet',num2str(XYLTp2(2)),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(XYRTp2(1),XYRTp2(2),...
strcat('\leftarrow ', num2str(XYRTp2(1)), '\bullet',num2str(XYRTp2(2))),...
'FontSize',12,'HorizontalAlignment','left');
text(XYLBp2(1),XYLBp2(2),...
strcat(num2str(XYLBp2(1)), '\bullet',num2str(XYLBp2(2)),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(XYRBp2(1),XYRBp2(2),...
strcat('\leftarrow ', num2str(XYRBp2(1)), '\bullet',num2str(XYRBp2(2))),...
'FontSize',12,'HorizontalAlignment','left');
%---
text(XcolsWIDE,0,...
strcat(num2str(XcolsWIDE), '\bullet',num2str(0),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(0,0,...
strcat('\leftarrow ', num2str(0), '\bullet',num2str(0)),...
'FontSize',12,'HorizontalAlignment','left');
text(XcolsWIDE,YrowsHIGH,...
strcat(num2str(XcolsWIDE), '\bullet',num2str(YrowsHIGH),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');
text(0,-YrowsHIGH,...
strcat('\leftarrow ', num2str(0), '\bullet',num2str(-YrowsHIGH)),...
'FontSize',12,'HorizontalAlignment','left');
%---

figure(fig99)
subplot(5,5,[1 6]), imagesc(Zfield);     % <--FIG--##
colormap('bone')
title('PSD1');
figure(fig99)
subplot(5,5,[2 7]), imagesc(PSDfield); % <--FIG--##
title('PSD2');

set(gcf,'Position',...
	[SNSZ(3)/1.5	SNSZ(4)/5		SNSZ(3)/3	SNSZ(4)/1.5]);


%=========================================================%

SPYN1xv = [XYLTpr1(1) XYRTpr1(1) XYRBpr1(1) XYLBpr1(1) XYLTpr1(1)]';
SPYN1yv = [XYLTpr1(2) XYRTpr1(2) XYRBpr1(2) XYLBpr1(2) XYLTpr1(2)]';
SPYN2xv = [XYLTpr2(1) XYRTpr2(1) XYRBpr2(1) XYLBpr2(1) XYLTpr2(1)]';
SPYN2yv = [XYLTpr2(2) XYRTpr2(2) XYRBpr2(2) XYLBpr2(2) XYLTpr2(2)]';

PSD1xv = [XYLTp1(1) XYRTp1(1) XYRBp1(1) XYLBp1(1) XYLTp1(1)]';
PSD1yv = [XYLTp1(2) XYRTp1(2) XYRBp1(2) XYLBp1(2) XYLTp1(2)]';
PSD2xv = [XYLTp2(1) XYRTp2(1) XYRBp2(1) XYLBp2(1) XYLTp2(1)]';
PSD2yv = [XYLTp2(2) XYRTp2(2) XYRBp2(2) XYLBp2(2) XYLTp2(2)]';
 
[PERI1xv, PERI1yv] = polybool('xor', PSD1xv, PSD1yv, SPYN1xv, SPYN1yv);
[PERI2xv, PERI2yv] = polybool('xor', PSD2xv, PSD2yv, SPYN2xv, SPYN2yv);



%=================================%
%{
plot(SPYN1xv,SPYN1yv,GluR2xyl(1,G2INSPYN1),GluR2xyl(2,G2INSPYN1),...
    'r+',GluR2xyl(1,~G2INSPYN1),GluR2xyl(2,~G2INSPYN1),'bo');
hold on
plot(SPYN2xv,SPYN2yv,GluR2xyl(1,G2INSPYN2),GluR2xyl(2,G2INSPYN2),...
    'r+',GluR2xyl(1,~G2INSPYN2),GluR2xyl(2,~G2INSPYN2),'bo');
hold off

plot(PSD1xv,PSD1yv,GluR2xyl(1,G2INPSD1),GluR2xyl(2,G2INPSD1),...
	'r+',GluR2xyl(1,~G2INPSD1),GluR2xyl(2,~G2INPSD1),'bo');
hold on
plot(PSD2xv,PSD2yv,GluR2xyl(1,G2INPSD2),GluR2xyl(2,G2INPSD2),...
	'r+',GluR2xyl(1,~G2INPSD2),GluR2xyl(2,~G2INPSD2),'bo');
hold off

plot(PERI1xv,PERI1yv,GluR2xyl(1,G2INPERI1),GluR2xyl(2,G2INPERI1),...
    'r+',GluR2xyl(1,~G2INPERI1),GluR2xyl(2,~G2INPERI1),'bo');
hold on
plot(PERI2xv,PERI2yv,GluR2xyl(1,G2INPERI2),GluR2xyl(2,G2INPERI2),...
    'r+',GluR2xyl(1,~G2INPERI2),GluR2xyl(2,~G2INPERI2),'bo');
hold off
%}
%{
plot(SPYN1xv,SPYN1yv,GluR1xyl(1,G1INSPYN1),GluR1xyl(2,G1INSPYN1),...
    'r+',GluR1xyl(1,~G1INSPYN1),GluR1xyl(2,~G1INSPYN1),'bo');
hold on
plot(SPYN2xv,SPYN2yv,GluR1xyl(1,G1INSPYN2),GluR1xyl(2,G1INSPYN2),...
    'r+',GluR1xyl(1,~G1INSPYN2),GluR1xyl(2,~G1INSPYN2),'bo');
hold off
 
plot(PSD1xv,PSD1yv,GluR1xyl(1,G1INPSD1),GluR1xyl(2,G1INPSD1),...
    'r+',GluR1xyl(1,~G1INPSD1),GluR1xyl(2,~G1INPSD1),'bo');
hold on
plot(PSD2xv,PSD2yv,GluR1xyl(1,G1INPSD2),GluR1xyl(2,G1INPSD2),...
    'r+',GluR1xyl(1,~G1INPSD2),GluR1xyl(2,~G1INPSD2),'bo');
hold off
 
plot(PERI1xv,PERI1yv,GluR1xyl(1,G1INPERI1),GluR1xyl(2,G1INPERI1),...
    'r+',GluR1xyl(1,~G1INPERI1),GluR1xyl(2,~G1INPERI1),'bo');
hold on
plot(PERI2xv,PERI2yv,GluR1xyl(1,G1INPERI2),GluR1xyl(2,G1INPERI2),...
    'r+',GluR1xyl(1,~G1INPERI2),GluR1xyl(2,~G1INPERI2),'bo');
hold off
%}
%=================================%
end






%-------------------------------------------%
% preGR2STEP
%-------------------------------------------%
function [GluR2xyl GR2Ndots] = preGR2STEP(stepN, GR2Ndots,GluR2xyds,GluR2xyl,LsGR2psd,...
	Frow,Fcol,row1F,col1F,row1L,col1L,PSDLOC,PSDSZE,GR2Sdots)



GR2Ndots = GR2Ndots;


%===========================%
% Add Step To Location

if GR2Ndots >=1
for j = 1:GR2Ndots
if GluR2xyl(1,j)>=(col1F+3) && GluR2xyl(1,j)<=(col1L-3) &&...
          GluR2xyl(2,j)>=(row1F+11) && GluR2xyl(2,j)<=(row1L-0)
	  
            GluR2xyds(:,j) = GluR2xyds(:,j)*(LsGR2psd);
end
end
end

for i = 1:2
for j = 1:GR2Ndots
	GluR2xyl(i,j) = GluR2xyl(i,j)+GluR2xyds(i,j);
end
end


%===========================%
% Keep In Box

    for j = 1:GR2Ndots  
		if GluR2xyl(1,j)<=(col1F)
			GluR2xyl(1,j) = col1F;
		end
		
		if GluR2xyl(1,j)>=(col1L)
			GluR2xyl(1,j) = col1L;
		end
		
	end
	
	for j = 1:GR2Ndots
		if GluR2xyl(2,j)<=(row1F)
			GluR2xyl(2,j) = row1F;
		end
		
		if GluR2xyl(2,j)>=(row1L)
			GluR2xyl(2,j) = row1L;
		end
	end
%===========================%

	
end
%-------------------------------------------%
% preGR1STEP
%-------------------------------------------%
function [GluR1xyl GR1Ndots] = preGR1STEP(stepN, GR1Ndots,GluR1xyds,GluR1xyl,LsGR1psd,...
    Frow,Fcol,row1F,col1F,row1L,col1L,PSDLOC,PSDSZE,GR2Sdots)
 

%===========================%
% Add Step To Location

for j = 1:GR1Ndots
if GluR1xyl(1,j)>=(col1F+3) && GluR1xyl(1,j)<=(col1L-3) &&...
          GluR1xyl(2,j)>=(row1F+11) && GluR1xyl(2,j)<=(row1L-0)
	  
            GluR1xyds(:,j) = GluR1xyds(:,j)*(LsGR1psd);
end
end

for i = 1:2
for j = 1:GR1Ndots
	GluR1xyl(i,j) = GluR1xyl(i,j)+GluR1xyds(i,j);
end
end
 
%===========================%
% Keep In Box

    for x = 1:GR1Ndots  
        if GluR1xyl(1,x)<=(col1F)
            GluR1xyl(1,x) = col1F;
        end
        
        if GluR1xyl(1,x)>=(col1L)
            GluR1xyl(1,x) = col1L;
		end
        
    end
    
    for y = 1:GR1Ndots
        if GluR1xyl(2,y)<=(row1F)
            GluR1xyl(2,y) = row1F;
        end
        
        if GluR1xyl(2,y)>=(row1L)
            GluR1xyl(2,y) = row1L;
        end
	end
%===========================%
    
end
%-------------------------------------------%
% PreMAINPLOT
%-------------------------------------------%
function [] = PreMAINPLOT(GluR1xyl, GluR2xyl, Fcol, PSDLOC, PSDSZE,...
	row1F,row1L,col1F,col1L)
%-------------------------------------------%
xlim = [0 Fcol];
ylim = [0 Fcol];
%---

%---
PSDLOC1 = PSDLOC;
PSDLOC1 = [PSDLOC1(2,:); PSDLOC1(1,:)];

% SYN1SZ = (PSDSZE(1,1) + (PSDSZE(2,1)*2))-1;
% PSD1SZ = (PSDSZE(1,1))-1;

SYN1SZ = (PSDSZE(1,1) + (PSDSZE(2,1)*2))-1;
PSD1SZ = (PSDSZE(1,1))-1;

P1INs = PSDSZE(2,1);
% P1INx = PSDLOC1(1,1) + P1INs;
% P1INy = PSDLOC1(2,1) + P1INs;
P1INx = (col1L);
P1INy = (row1F);

P1OUTx=PSDLOC(2,1);
P1OUTy=PSDLOC(1,1);

%---

%=================================%
%       MAIN 2D PLOT
%---------------------------------%
%     MAIN DOTS
%----------------------%
figure(70) 
AMPARPlot = gscatter(GluR2xyl(1,:),GluR2xyl(2,:));
axis([xlim, ylim]);
set(AMPARPlot,'marker','.','markersize',[6],'color',[1 0 0])
hold on
AMPARPlot = gscatter(GluR1xyl(1,:),GluR1xyl(2,:));
axis([xlim, ylim]);
set(AMPARPlot,'marker','.','markersize',[6],'color',[0 1 1])
hold off


rectangle('Position',[P1OUTx,P1OUTy,SYN1SZ,SYN1SZ])
rectangle('Position',[P1INx,P1INy,1,4],'Curvature',[.2,.1],...
	'EdgeColor','w','FaceColor','w','LineWidth',1,'LineStyle','--')
%======================%





end




%===================================%
% preMOVEGLUR
%===================================%
function [GluR2xyl] = preMOVEGLUR2(stepN,GluR2xyds,GluR2xyl,XYLBpr1,XYRTpr1)


% GluR2xyds(:,G2INPSD) = GluR2xyds(:,G2INPSD)*(LsGR2psd);

GluR2xyl = GluR2xyl+GluR2xyds;

GR2dots=size(GluR2xyl,2);

LFT = XYLBpr1(1); BOT = XYLBpr1(2);
RIT = XYRTpr1(1); TOP = XYRTpr1(2);
%===========================%
    for j = 1:GR2dots 
		if GluR2xyl(1,j)<=(LFT)
			GluR2xyl(1,j) = LFT;
		end
		
		if GluR2xyl(1,j)>=(RIT)
			GluR2xyl(1,j) = RIT;
		end
	end
	
	for j = 1:GR2dots
		if GluR2xyl(2,j)<=(BOT)
			GluR2xyl(2,j) = BOT;
		end
		
		if GluR2xyl(2,j)>=(TOP)
			GluR2xyl(2,j) = TOP;
		end
	end
%===========================%

%===================================%
% if stepN == 500; keyboard; end
end
%===================================%
function [GluR1xyl] = preMOVEGLUR1(stepN,GluR1xyds,GluR1xyl,XYLBpr1,XYRTpr1)
 
 
GluR1xyl = GluR1xyl+GluR1xyds;

GR1dots=size(GluR1xyl,2);

LFT = XYLBpr1(1); BOT = XYLBpr1(2);
RIT = XYRTpr1(1); TOP = XYRTpr1(2);
%===========================%
    for j = 1:GR1dots
        if GluR1xyl(1,j)<=(LFT)
            GluR1xyl(1,j) = LFT;
        end
        
        if GluR1xyl(1,j)>=(RIT)
            GluR1xyl(1,j) = RIT;
		end        
    end
    
    for j = 1:GR1dots
        if GluR1xyl(2,j)<=(BOT)
            GluR1xyl(2,j) = BOT;
        end
        
        if GluR1xyl(2,j)>=(TOP)
            GluR1xyl(2,j) = TOP;
        end
    end
%===========================%

 
%===================================%
% if stepN == 500; keyboard; end
end
%===================================%
% preMAINPLOT
%===================================%
function [] = preMAINPLOT(GluR1xyl,GluR2xyl,...
	XYLBpr1,XYRTpr1,XYLBp1,XYRTp1,...
	PSD1WH,PSD2WH,PERI1WH,PERI2WH,SPYN1WH,SPYN2WH,XWIDE,YHIGH)
%-------------------------------------------%


%=================================%
%       MAIN 2D PLOT
%---------------------------------%
LFTPSA = XYLBpr1(1); BOTPSA = XYLBpr1(2);
RITPSA = XYRTpr1(1); TOPPSA = XYRTpr1(2);
LFTPSD = XYLBp1(1); BOTPSD = XYLBp1(2);
RITPSD = XYRTp1(1); TOPPSD = XYRTp1(2);

%---
xlim = [0 XWIDE];
ylim = [0 XWIDE];
%---
figure(70) 
AMPARPlot = gscatter(GluR2xyl(1,:),GluR2xyl(2,:));
axis([xlim, ylim]);
set(AMPARPlot,'marker','.','markersize',[6],'color',[1 0 0])
hold on
AMPARPlot = gscatter(GluR1xyl(1,:),GluR1xyl(2,:));
axis([xlim, ylim]);
set(AMPARPlot,'marker','.','markersize',[6],'color',[0 1 1])
rectangle('Position',[LFTPSA,BOTPSA,SPYN1WH,SPYN1WH])
rectangle('Position',[LFTPSD,BOTPSD,PSD1WH,PSD1WH])
rectangle('Position',[RITPSA,BOTPSA,1,4],'Curvature',[.2,.1],...
	'EdgeColor','w','FaceColor','w','LineWidth',1,'LineStyle','--')
%======================%
hold off

end
%===================================%



%-------------------------------------------%
%		SLOTS SLOTS SLOTS SLOTS SLOTS
%-------------------------------------------%
% ----SLOTS GluR1----
function [GluR1xyds GluR1_TdwellPSD G1FSLOTS] = SLOTSG1(GluR1xyds,GluR1_TdwellPSD,...
	KonGR1,KoffGR1,GR1slotN)

PSD1slotN = GR1slotN;

[r,c,v] = find(GluR1_TdwellPSD);
rcv = [v r c]; % VALUE ROW COLUMN
[r1] = find(rcv(:,2)<2);
p1r = rcv(r1,:);
p1rs = sortrows(p1r,-1);
if size(p1rs,1) < PSD1slotN
	p1rs = vertcat(p1rs,zeros(PSD1slotN,3));
end

SLPSD1 = p1rs(1:PSD1slotN,:); %SLOT LOCATION PSD1
SLPSD1LOCkon = SLPSD1((SLPSD1(:,1)>KonGR1),3);
SLPSD1LOCkoff = SLPSD1((SLPSD1(:,1)>(KoffGR1+KonGR1)),3);
GluR1xyds(:,SLPSD1LOCkon) = GluR1xyds(:,SLPSD1LOCkon)*(.001);
GluR1xyds(:,SLPSD1LOCkoff) = GluR1xyds(:,SLPSD1LOCkoff)/(.001);
G1FSLOTS = numel(SLPSD1LOCkon)-numel(SLPSD1LOCkoff);
%================%
end
%===================================%
% ----SLOTS GluR2----
function [GluR2xyds GluR2_TdwellPSD G2FSLOTS] = SLOTSG2(GluR2xyds,GluR2_TdwellPSD,...
    KonGR2,KoffGR2,GR2slotN)
 
PSD1slotN = GR2slotN;
 
[r,c,v] = find(GluR2_TdwellPSD);
rcv = [v r c]; % VALUE ROW COLUMN
[r1] = find(rcv(:,2)<2);
p1r = rcv(r1,:);
p1rs = sortrows(p1r,-1);
if size(p1rs,1) < PSD1slotN
    p1rs = vertcat(p1rs,zeros(PSD1slotN,3));
end 

SLPSD1 = p1rs(1:PSD1slotN,:); %SLOT LOCATION PSD1
SLPSD1LOCkon = SLPSD1((SLPSD1(:,1)>KonGR2),3);
SLPSD1LOCkoff = SLPSD1((SLPSD1(:,1)>(KoffGR2+KonGR2)),3);
GluR2xyds(:,SLPSD1LOCkon) = GluR2xyds(:,SLPSD1LOCkon)*(.001);
GluR2xyds(:,SLPSD1LOCkoff) = GluR2xyds(:,SLPSD1LOCkoff)/(.001);

G2FSLOTS = numel(SLPSD1LOCkon)-numel(SLPSD1LOCkoff);
%================%
end
%===================================%




