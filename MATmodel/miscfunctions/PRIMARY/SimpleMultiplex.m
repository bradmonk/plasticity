function [] = SimpleMultiplex()
clc; close all; clear all;


% USER ENTERED VALUES
NSteps = 500;			% number of steps (loop time)
Ndots = 100;			% number of particles
NLoops= 100;
OLoopN = 10;

Dsurf = 1.0;		% D coef ES
DboxA = 0.1;		% D coef SYN1
DboxB = 0.5;		% D coef SYN2

Dbox = fliplr(linspace(DboxA,DboxB,OLoopN));

%===========================================%
for OLoop = 1:OLoopN
%===========================================%

Dsurf = Dsurf;
Dbox1 = Dbox(OLoop);
Dbox2 = DboxB;

% Choquet found that diffusion rates for...
% ES: .1	SYN_Base: .05	SYN_LTP: .01

Scale = 1;				% scale of model (life:model) 1/10
TimeStep = 1000;		% time step (ms)
	
% BASE DIFFUSION RATES EQUATIONS
Sc = Scale;					% scale of model (life:model)
t = TimeStep/1000;			% time step (ms)
dm = 2;                     % dimensions
Da = Dbox1*t/Sc;		% Diffusion Rate A (D = L² / 2d*t)
Db = Dbox2*t/Sc;		% Diffusion Rate B
Dr = Da/Db;					% Ratio of Da:Ds (1/Ls)^2;
Dn = Da/Dr;					% new D after scaling L
k = sqrt(dm*Da);			% stdev of D's step size distribution
L = sqrt(2*dm*Da);			% average diagonal (2D) step size
Lx = L/sqrt(2);             % average linear (1D) step size
Ls = 1/sqrt(Dr);			% scales Lx values for Dn
MSD = 2*dm*Da;				% mean squared displacement

Ds = Dsurf*t/Sc;		% Diffusion Rate of Surf (D = L² / 2d*t)
Db1 = Dbox1*t/Sc;		% Diffusion Rate Box1
Db2 = Dbox2*t/Sc;		% Diffusion Rate Box2

k = sqrt(dm*Ds);		% stdev of Dsurf step size distribution

Lsb1 = 1/sqrt(Ds/Db1);	% scales Lx values for Dn
Lsb2 = 1/sqrt(Ds/Db2);	% scales Lx values for Dn





XYL = zeros(2,Ndots);		% XY particle locations
XYS = zeros(2,Ndots);		% XY step sizes
XYLp = zeros(2,NSteps);		% preallocate matrix for trace dot
XYL(1,:) = XYL(1,:);
XYL(2,:) = XYL(2,:);


POLYSz = [100 100];			% size of polygon dendrite enclosure (XY in µm)
POLYSz = POLYSz./Scale;		% scale enclosures 
XWIDE = POLYSz(1)/2;		% half X enclosure size (will double below)
YHIGH = POLYSz(2)/2;		% half X enclosure size (will double below)

% [POLYSz(1)/2 POLYSz(2)/2 POLYSz(1) POLYSz(2)]
BOARDER = [-POLYSz(1)/2 -POLYSz(2)/2 POLYSz(1) POLYSz(2)];

% Box1 location [X Y W H]
BOXxy=20; BOXwh=20;
BOXLOC1 = [-10 15 BOXwh BOXwh]; 	
BOXLOC2 = [-10 -35 BOXwh BOXwh];

Box1LB = BOXLOC1(1:2);
Box1RT = Box1LB + BOXLOC1(3:4);
Box2LB = BOXLOC2(1:2);
Box2RT = Box2LB + BOXLOC2(3:4);


doTrace=0;
doLiveScatter = 0;
doLiveBar = 0;

if doLiveScatter
%================================================%
%               FIGURE SETUP
%------------------------------------------------%
Flh = figure(1);
set(Flh,'Units','pixels');  scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/2  scnsize(4)/3  scnsize(3)/2  scnsize(4)/1.8];
set(Flh,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%--------
xlim = [-XWIDE*1.01 XWIDE*1.01]; ylim = [-YHIGH*1.01 YHIGH*1.01];
%---
subplot('Position',[.05 .05 .65 .9])
Ph1 = scatter(XYL(1,:),XYL(2,:),5,[0 0 1]);
set(Ph1,'Marker','o','SizeData',60,'LineWidth',.5,...
	'MarkerFaceColor',[.95 .1 .1],'MarkerEdgeColor','none')
axis([xlim, ylim]);
% axis off
set(gca,'XTick',[],'YTick',[],'Color',[.8 .8 .8])
hold on

% pause(1)
% SPh2 = subplot(1,1,1);
% Ph2 = plot(XYLp(1,:),XYLp(2,:));
% axis([xlim, ylim]);
% axis off
% hold(SPh2,'on');
rectangle('Position',BOXLOC1); hold on;
rectangle('Position',BOXLOC2); hold on;
rectangle('Position',BOARDER); hold on;

pause(1)

Ah3 = subplot('Position',[.75 .05 .22 .9]);
Ph3 = bar([1 1]);
set(Ah3,'Ylim',[0 100])

pause(1)
%-------------------------------------------------%
end %if doLivePlot


SMxE = zeros(50);


%===========================================%
for Ln = 1:NLoops
mm=1;
%===========================================%

	%===========================================%
	for Nt = 1:NSteps 
	%===========================================%


	XYS = (k * randn(2,Ndots));	% generates step sizes
	
	% Sticky Functions
	%{.
	% [XYS] = STICKYFUN(XYL,XYS,Box1LB,Box2LB,Ls);
	Box1 = inboxfun(Box1LB,Box1RT,XYL);
	Box2 = inboxfun(Box2LB,Box2RT,XYL);
	XYS(:,Box1) = XYS(:,Box1)*(Lsb1);
	XYS(:,Box2) = XYS(:,Box2)*(Lsb2);
	%}
	
	XYL = XYL+XYS;				% adds step to location
	
	% Keep everything inside enclosure  %
	[XYL] = ENCLOSE(Nt,XYL,XWIDE,YHIGH,Ndots);

	%[G1xy] = dopartmx(XYL,SMxE);
	%figure(1)
	%imagesc(G1xy)
	
	if doTrace
	XYLp(:,Nt) = XYL(:,1);		% save step of first dot (for trace)
	DUALPLOT(Nt,XYL,Ph1,XYLp,Ph2)
	XL(Nt,:) = XYL(1,:);
	YL(Nt,:) = XYL(2,:);
	end
	
	if doLiveScatter
	%figure(SPh2)
	set(Ph1,'XData',XYL(1,:),'YData',XYL(2,:));
	drawnow
	end
	

	
	if mod(Nt,100)==0	
	Box1N(mm) = sum(Box1>0);
	Box2N(mm) = sum(Box2>0);
		if doLiveBar
		set(Ph3,'YData',[Box1N(mm) Box2N(mm)]);
		drawnow
		end
	mm=mm+1;
	end

	
	%===========================================%
	end % for Nt = 1:Nsteps 
	%===========================================%




Boxy1{Ln} = Box1N';
Boxy2{Ln} = Box2N';
Boxy3{Ln} = ((Box1N+1) ./ (Box2N+1))';
Boxy0{Ln} = (Ndots - (Box1N+Box2N))';

MeanBox0(Ln) = mean(Boxy0{Ln});
MeanBox1(Ln) = mean(Boxy1{Ln});
MeanBox2(Ln) = mean(Boxy2{Ln});
MeanRatio(Ln) = mean(Boxy3{Ln});


if mod(Ln,10)==0;disp(Ln); end;
%===========================================%
end
%===========================================%

Mbox0(OLoop) = mean(MeanBox0);
Mbox1(OLoop) = mean(MeanBox1);
Mbox2(OLoop) = mean(MeanBox2);
MRb1b2(OLoop) = mean(MeanRatio);

DboxD(:,OLoop) = [Dsurf; Dbox1; Dbox2];
DboxR(OLoop) = Dbox2/Dbox1;

%===========================================%
end % OLoop
%===========================================%


DPLOT(DboxD,DboxR,MRb1b2,Mbox0,Mbox1,Mbox2,NSteps)

FINALPLOT(Boxy1,Boxy2,Boxy3,Ds,Db1,Db2)



%==============================================================%
end %varargout={Boxy1,Boxy2};
%==============================================================%







function [G1xy] = dopartmx(XYL,SMxE)


GR1c = round(XYL(1,:))+55;
GR1r = round(XYL(2,:))+55;
% GR1xy(1,:) = GR1c;
% GR1xy(2,:) = GR1r;

G1xy = SMxE;
for xy = 1:numel(GR1c)
G1xy(GR1r(xy),GR1c(xy)) = 1;
end

end



%-------------------------------------------%
% STICKYFUN: keep particles inside polygon
%-------------------------------------------%
function [XYS] = STICKYFUN(XYL,XYS,Box1LB,Box2LB,Ls)

INPSD1 = inboxfun(Box1LB,Box1RT,XYL);
INPSD2 = inboxfun(Box2LB,Box2RT,XYL);

XYS(:,INPSD1) = XYL(:,INPSD1)*(Ls);
XYS(:,INPSD2) = XYL(:,INPSD2)*(Ls);

end



%-------------------------------------------%
% ENCLOSE: keep particles inside polygon
%-------------------------------------------%
function [XYL] = ENCLOSE(Nt,XYL,XWIDE,YHIGH,Ndots)

	for j = 1:Ndots 
		if XYL(1,j)>(XWIDE) || XYL(1,j)<(-XWIDE)
			XYL(1,j) = sign(XYL(1,j))*(XWIDE);
		elseif XYL(2,j)>(YHIGH) || XYL(2,j)<(-YHIGH)
			XYL(2,j) = sign(XYL(2,j))*(YHIGH);
		end
	end

end



%===================================%
%		inboxfun
%===================================%
% Tests whether particles are in a box polygon
% and returns a logical vector
%===================================%
function [inbox] = inboxfun(LB,RT,xyl)

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

xylLB1 = xyl(1,:) > LB(1);
xylRT1 = xyl(1,:) < RT(1);
xylLB2 = xyl(2,:) > LB(2);
xylRT2 = xyl(2,:) < RT(2);

inbox = xylLB1 & xylRT1 & xylLB2 & xylRT2;

end
%===================================%



%-------------------------------------------%
% DUALPLOT: plots diffusion
%-------------------------------------------%
function [] = DUALPLOT(Nt,XYL,Ph1,XYLp,SPh2)
%-------------------------------------------%

LLg = 3;
if Nt > LLg
xp = (XYLp(1,(Nt-LLg):Nt));
yp = (XYLp(2,(Nt-LLg):Nt));
else
xp = (XYLp(1,Nt));
yp = (XYLp(2,Nt));
end

set(Ph1,'XData',XYL(1,:),'YData',XYL(2,:));
drawnow

plot(xp,yp,'Parent',SPh2);
drawnow

 
end


%-------------------------------------------%
% FINALPLOT: plots steady-state averages
%-------------------------------------------%
function [] = DPLOT(DboxD,DboxR,MRb1b2,Mbox0,Mbox1,Mbox2,NSteps)


%%
%------------------------------------------%
Fh = figure;
set(Fh,'OuterPosition',[200 200 1400 500])
% Color
neongreen = [.1 .9 .1];
liteblue = [.2 .9 .9];
hotpink=[.9 .1 .9];
%------------------------------------------%
axes('Position',[.07 .14 .4 .78]);
%------------------------------------------%
% [ph1] = line(1:numel(DboxR),DboxR);
[ph1] = plot(DboxR);
legend(ph1, 'Set Ratio');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
[ph2] = plot(MRb1b2);
legend([OUTH;ph2],OUTM{:},'Sim Ratio');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%------------------------------------------%
% xt = (get(gca,'XTick'));
% xt = linspace(0,NSteps,numel(xt)); %.* dT./60
% set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 9; c1=neongreen; c2=hotpink; c3=liteblue;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',2,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c1);
set(ph2,'LineStyle',':','Color',c2,'LineWidth',2,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c2);
hTitle  = title('Set D-coef Ratio (s2:s1) vs Resulting SSn Ratio');
hXLabel = xlabel('Trial');
hYLabel = ylabel('Particles (+/- SEM)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12);
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
ylim([0 haxes(4)*1.2 ]);
% xlim([0 (haxes(2)*.9)]);
%======================================================================%

%------------------------------------------%
axes('Position',[.55 .14 .4 .78]);
%------------------------------------------%
% [ph1] = line(1:numel(DboxR),DboxR);
[ph1] = plot(Mbox1);
legend(ph1, 'SYN1 N');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
[ph2] = plot(Mbox2);
legend([OUTH;ph2],OUTM{:},'SYN2 N');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
[ph3] = plot(Mbox0);
legend([OUTH;ph3],OUTM{:},'ES N');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%------------------------------------------%
% xt = (get(gca,'XTick'));
% xt = linspace(0,NSteps,numel(xt)); %.* dT./60
set(gca,'XTickLabel', sprintf('%.2f|',DboxR))
set(legend,'Location','NorthEast');
%------------------------------------------%
MS1 = 9; c1=neongreen; c2=hotpink; c3=liteblue;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',2,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c1);
set(ph2,'LineStyle',':','Color',c2,'LineWidth',2,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c2);
set(ph3,'LineStyle','-.','Color',c3,'LineWidth',2,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c3);
hTitle  = title('Particle SSn in ES, S1 & S2 at Different D-coef Ratios');
hXLabel = xlabel('D-coef Ratio s2:s1');
hYLabel = ylabel('Particles (+/- SEM)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12);
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
ylim([0 haxes(4)*1.2 ]);
% xlim([0 (haxes(2)*.9)]);
%======================================================================%
%%



end



%-------------------------------------------%
% FINALPLOT: plots steady-state averages
%-------------------------------------------%
function [] = FINALPLOT(data1,data2,data3,Ds,Db1,Db2)



doRatFig = 0;
doRawFig = 1;

%%

if doRawFig
Fh = figure;
set(Fh,'OuterPosition',[200 200 650 700])

% Color
neongreen = [.1 .9 .1];
liteblue = [.2 .9 .9];
hotpink=[.9 .1 .9];
%------------------------------------------%
[ph1 hax1 mu1] = CIplot(data3,neongreen,5);
legend(ph1, horzcat('Ratio _{ES: ', num2str(Ds), ' um²/s}'));
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
ph11 = plot(mu1);
hold on
[ph2 hax2 mu2] = CIplot(data1,hotpink,5);
legend([OUTH;ph2],OUTM{:},strcat('Synapse-1 _{', num2str(Db1), ' um²/s}'));
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
ph22 = plot(mu2);
hold on
pause(.5)
[ph3 hax3 mu3] = CIplot(data2,liteblue,5);
legend([OUTH;ph3],OUTM{:},strcat('Synapse-2 _{', num2str(Db2), ' um²/s}'));
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
ph33 = plot(mu3);
hold on
%------------------------------------------%
xt = (get(gca,'XTick'))*100;
set(gca,'XTickLabel', sprintf('%.0f|',xt))
%------------------------------------------%
MS1 = 9; c1=neongreen; c2=hotpink; c3=liteblue;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',6,'MarkerEdgeColor',c1,'MarkerFaceColor','none');
set(ph11,'LineStyle','none','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c1);
set(ph22,'LineStyle','none','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c2);
set(ph33,'LineStyle','none','Color',c3,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c3);
hTitle  = title('Distribution of Particles with Brownian Motion');
hXLabel = xlabel('Time');
hYLabel = ylabel('Particles (+/- SEM)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12);
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
ylim([0 haxes(4)*1.2 ]);
% xlim([0 (haxes(2)*.9)]);
%======================================================================%
end
%%


end



%-------------------------------------------%
% FINALPLOT: plots steady-state averages
%-------------------------------------------%
function [] = FINALPLOT2(data1,data2,data3)



doRatFig = 1;
doRawFig = 1;

%%

if doRawFig
Fh = figure;
set(Fh,'OuterPosition',[200 200 900 600])

% Color
neongreen = [.1 .9 .1];
liteblue = [.2 .9 .9];
hotpink=[.9 .1 .9];
%------------------------------------------%
[ph1 hax1 mu1] = CIplot(data1,neongreen,5);
legend(ph1,'Box1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
ph1 = plot(mu1);
hold on
[ph2 hax2 mu2] = CIplot(data2,hotpink,5);
legend([OUTH;ph2],OUTM{:},'Box2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
ph2 = plot(mu2);
hold on
%------------------------------------------%
xt = (get(gca,'XTick'))*100;
set(gca,'XTickLabel', sprintf('%.0f|',xt))
%------------------------------------------%
MS1 = 8; c1=neongreen; c2=hotpink;
set(ph1,'LineStyle','none','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c1);
set(ph2,'LineStyle','none','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c2);
hTitle  = title('Distribution of Particles with Brownian Motion');
hXLabel = xlabel('Time');
hYLabel = ylabel('Particles (+/- SEM)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12);
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
ylim([0 haxes(4)*1.2 ]);
% xlim([0 (haxes(2)*.9)]);
%======================================================================%
end


if doRatFig
Fh = figure;
set(Fh,'OuterPosition',[200 200 900 600])

% Color
neongreen = [.1 .9 .1];
liteblue = [.2 .9 .9];
hotpink=[.9 .1 .9];
%------------------------------------------%
[ph1 hax1 mu1] = CIplot(data3,neongreen,5);
legend(ph1,'Box1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
ph1 = plot(mu1);
% hold on
% [ph2 hax2 mu2] = CIplot(data2,hotpink,5);
% legend([OUTH;ph2],OUTM{:},'Box2');
% [LEGH,OBJH,OUTH,OUTM] = legend;
% hold on
% ph2 = plot(mu2);
% hold on
%------------------------------------------%
xt = (get(gca,'XTick'))*100;
set(gca,'XTickLabel', sprintf('%.0f|',xt))
%------------------------------------------%
MS1 = 8; c1=neongreen; c2=hotpink;
set(ph1,'LineStyle','none','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c1);
% set(ph2,'LineStyle','none','Color',c2,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c2);
hTitle  = title('Distribution of Particles with Brownian Motion');
hXLabel = xlabel('Time');
hYLabel = ylabel('Ratio (+/- SEM)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12);
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
ylim([1 haxes(4)*1.2 ]);
% xlim([0 (haxes(2)*.9)]);
%======================================================================%
end
%%


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
%%


end
