%============================================%
clc, close all; scsz = get(0,'ScreenSize');
%============================================%

%-------------------------------------------%
% COLOR
%-------------------------------------------%
RGB = [255 20 147] ./ 255;
neongreen = [.1 .9 .1];
liteblue = [.2 .9 .9];
hotpink=[.9 .1 .9];

%-------------------------------------------%
% ULTIMATE PLOT SETUP
%-------------------------------------------%
clc, close all;
XDATA = 1:10; YDATA=11:20; ZDATA=rand(1,10);
CDATA = padarray(randi([1,10],[5,5]),[2 2],0);
nn=10;

scsz = get(0,'ScreenSize'); scsx=scsz(3); scsy=scsz(4);
Fh1 = figure(1);
set(gcf,'Position',[scsx/4 scsy/4 scsx/2 scsy/2]);

axes('Position',[.05 .55 .40 .40]);
Ph1 = plot(XDATA,YDATA,'LineWidth',3);

axes('Position',[.05 .05 .40 .40]);
Ph2 = scatter(XDATA,YDATA);

axes('Position',[.55 .55 .40 .40]);
Ph3 = plot3(XDATA,YDATA,ZDATA,'r');

axes('Position',[.55 .05 .40 .40]);
Ph4 = imagesc(CDATA);

%------------------------------%
% UPDATING FIGURES IN FOR-LOOPS
%------------------------------%
figure(Fh1)

set(Ph1,'XData',(1:nn),'YData',YDATA);
drawnow; hold on;

set(Ph2,'XData',XDATA,'YData',YDATA);
drawnow; hold on;

set(Ph3,'XData',XDATA,'YData',YDATA,'ZData',ZDATA);
drawnow;

set(Ph4,'CData',CDATA);
drawnow
%------------------------------%


%-----
set(gca, 'Box', 'on'),
set(gca,'xticklabel',[]),
set(gca,'yticklabel',[]),
set(gca,'zticklabel',[]),
%-----

%-----
set(gca,'YLim',[0 1]);
set(get(gca,'YLabel'),'String','MyYLabel')
xt = {'Loss','Remain'};
set(gca,'XTickLabel', sprintf('%s|',xt{1:2}))
legend(fhb,{'sim','calc'});
%-----


%---------------------------------------------%
Ph10 = plot3(PoRMx2DvX,PoRMx2DvY,PoRMx2DvZ,'r');
FMx.Fdots = {2};
set(Ph10,{'LineWidth'},repmat(FMx.Fdots,1,Fdots)')
hold on
%---------------------------------------------%
[XMx YMx ZMx] = scatter3prep({Po00 Po Pt Dv Pr v0 v1});
Ph99 = scatter3(XMx,YMx,ZMx,'Tag','TMxH');
%----------------
child_handles = get(Ph99,'Children');
mrkr = {'d';'d';'o';'x';'s';'v';'v'};
set(child_handles,{'Marker'},mrkr)
hold on
colr = {[0 0 0]; [1 0 0]; [0 0 0]; [0 1 0]; [0 .8 .8]; [1 .5 .1]; [.7 .5 .7]};
set(child_handles,{'MarkerFaceColor'},colr)
hold on
set(child_handles,'LineWidth',1,'MarkerSize',11)
%---------------------------------------------%
% anotext({Po00, 'Po00'},{Po, 'Po'},{Pt, 'Pt'},{Pr, 'Pr'},{v1, 'v1'},{v0, 'v0'})





set(gcf,'WindowStyle','docked') 
PLh = light('Position',[.5 .5 .5],'Style','infinite');
set(hp,'FaceLighting','gouraud','AmbientStrength',0.9)

S2Ph2 = scatter(PSDX,PSDY);
set(S2Ph2,'Marker','d','SizeData',70,'LineWidth',1.5,...
	'MarkerFaceColor','none','MarkerEdgeColor',[.8 .08 .14])

%---------------------------------------------%
% child_handles
%---------------------------------------------%
PoRMx2DvX = [PoMX(1,:); RMX2(1,:)];
PoRMx2DvY = [PoMX(2,:); RMX2(2,:)];
PoRMx2DvZ = [PoMX(3,:); RMX2(3,:)];
%---------------------------------------------%
Ph10 = plot3(PoRMx2DvX,PoRMx2DvY,PoRMx2DvZ,'r');
FMx.Fdots = {2};
set(Ph10,{'LineWidth'},repmat(FMx.Fdots,1,Fdots)')
hold on
%---------------------------------------------%
[XMx YMx ZMx] = scatter3prep({Po00 Po Pt Dv Pr v0 v1});
Ph99 = scatter3(XMx,YMx,ZMx,'Tag','TMxH');
%----------------
child_handles = get(Ph99,'Children');
mrkr = {'d';'d';'o';'x';'s';'v';'v'};
set(child_handles,{'Marker'},mrkr)
hold on
colr = {[0 0 0]; [1 0 0]; [0 0 0]; [0 1 0]; [0 .8 .8]; [1 .5 .1]; [.7 .5 .7]};
set(child_handles,{'MarkerFaceColor'},colr)
hold on
set(child_handles,'LineWidth',1,'MarkerSize',11)
%---------------------------------------------%
% anotext({Po00, 'Po00'},{Po, 'Po'},{Pt, 'Pt'},{Pr, 'Pr'},{v1, 'v1'},{v0, 'v0'})



%================================================%
%				Camera Movement
%------------------------------------------------%

%==================================================%
%				PLOT FIGURE AUTOROTATION
%--------------------------------------------------%
function [varargout] = RotaFig(varargin)


if nargin == 3
	do2D=varargin{1};
	pauseT1=varargin{2};
	pauseT2=varargin{3};
elseif nargin == 2
	do2D=varargin{1};
	pauseT1=varargin{2};
	pauseT2=.01;
elseif nargin == 1 
	do2D=varargin{1};
	pauseT1=1;
	pauseT2=.01;
else
	do2D = 1;
	pauseT1=1;
	pauseT2=.01;
end


%----------------------
% figure(Fh1)
camP0 = campos;
camT0 = camtarget;

nfram = 50;
cxp = linspace(-10,25,nfram);
for cf = 1:nfram
    campos([cxp(cf),-25,25])
    drawnow
	pause(pauseT2)
end
% figure(Fh1)
camP1 = campos;
camT1 = camtarget;
%----------------------



%----------------------
% camT = camT1+5;
camT = camT1;
cTx = 5;
cTy = -5;
cTz = 10;
nfram = 20;
cxt = fliplr(linspace(camT(1)-cTx,camT(1),nfram));
cyt = fliplr(linspace(camT(2)-cTy,camT(2),nfram));
czt = fliplr(linspace(camT(3)-cTz,camT(3),nfram));
cxt = [cxt fliplr(cxt)];
cyt = [cyt fliplr(cyt)];
czt = [czt fliplr(czt)];

for cf = 1:nfram*2
	camtarget([cxt(cf),cyt(cf),czt(cf)])
    drawnow
	pause(pauseT2)
end
for cf = 1:nfram*2
	camtarget([-cxt(cf),-cyt(cf),-czt(cf)])
    drawnow
	pause(pauseT2)
end
%----------------------


%----------------------
if do2D
view(0,0)
pause(1)
view(0,90)
pause(1)
view(90,0)
pause(1)
view(18,28)
end
%----------------------
varargout = {campos, camtarget};
end





%================================================%
%			imagesc LOOP FIGURE SETUP
%------------------------------------------------%
% Outside Loop
Fh1 = figure(1);
set(Fh1,'OuterPosition',(scsz./[2e-1 .2 4 4]))
Ph1 = imagesc(S);
% Inside Loop
set(Ph1,'CData',S);
drawnow



% imagesc extras
colormap('bone')
axis off
hT = title('  Activity Map (On/Off Events) ');
set(hT,'FontName','Arial','FontSize',16);
colorbar
% xlabel('Average On/Off Rate')
% labels = {'High Activity','Low Activity'};
% lcolorbar(labels,'fontweight','bold');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					Camera Movement
%=============================================================%
sz = [2.5e-3 2.5e-3 1.6 1.2];
Fh1 = FigSetup(1,sz);
%------------
figure(Fh1)
plot3(xp',yp',zp')
axis([-15 15 -15 15 0 30])
xlabel('X');ylabel('Y');zlabel('Z');
grid on


% figure(Fh1)
camP0 = campos;
camT0 = camtarget;

nfram = 50;
cxp = linspace(-10,25,nfram);
for cf = 1:nfram
    campos([cxp(cf),-25,25])
    drawnow
	pause(.05)
end


% figure(Fh1)
camP1 = campos;
camT1 = camtarget;


camT = camT1+5;
nfram = 50;
cxt = fliplr(linspace(-camT(1),camT(1),nfram));
cyt = fliplr(linspace(-camT(2),camT(2),nfram));
czt = fliplr(linspace(-camT(3),camT(3),nfram));
cxt = [cxt fliplr(cxt)];
cyt = [cyt fliplr(cyt)];
czt = [czt fliplr(czt)];

for cf = 1:nfram*2
	camtarget([cxt(cf),cyt(cf),czt(cf)])
    drawnow
	pause(.05)
end

% figure(Fh1)
campos(camP1)
camtarget(camT1)

% Set or get the value of the camera view angle
% camva
% to half or double the zoom use: camzoom(.5) or camzoom(2)
% camzoom(.5)

% vis3d
% campos([-10,-25,25])
% camzoom(.12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









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


%-------------------------------------------%
if nargin > 0 
spikes=varargin{1};
dospikes=1;
SStep=1;
SStepMod=50;
SpikeT = numel(spikes(:,1));
FStep=1;

else
dospikes=0;
spikes=0;
end
%-------------------------------------------%

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



%================================================%
%				FIGURE SETUP
%------------------------------------------------%
% FRAP = (SN0-(sum(NSoff)))/SN0;
FRAPL = (SN0-NNoff)./SN0;
SoCSN = (SN0-NNoff./Rev)./SN0;

set(ph10,'XData',[1:stepN],'YData',[FRAPL]);
drawnow; hold on;
set(ph11,'XData',[1:stepN],'YData',[SoCSN]);
drawnow; hold on;

% if stepN==3900;keyboard;end;
%=================================%

%-------------------------------------------------%



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
disp([' which is ' num2str(Spct) '%. size.'])
%-------------------------------------------%





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




%===============================%
% Homeostatic Bar Graph
%-------------------------------%
if doPLOT3D == 0 && doHOMEO ==1
if mod(stepN, 100) == 0
PSDHSn = [SPY1N (SPY1N+SPY2N) SPY2N];
% if stepN == 10
HSline = [HSlow HSlow HSlow HShi HShi HShi];
HSshape = [0 4 0 0 4 4]; 
figure(1);
subplot(5,5,[1 7]);
plot(HSshape,HSline,'LineWidth',3)
hold on
subplot(5,5,[1 7]);
bar(1:3, [PSDHSn']);
set(gca,'Ylim',[0 (HShi+20)])
hold off
end
end
%===============================%



%===============================%
% Live Particle Totals Bar Graph
%-------------------------------%
if doPLOT3D == 0 && doHOMEO == 0
	if mod(stepN, 100) == 0
	figure(1)
	subplot(5,5,[1 7]);
	bar(1:4, [EST'], .8, 'EdgeColor',[1 0.5 0.5]);
	set(gca,'Ylim',[0 GluR2Ndots],'XLim',[.5 4.5])
	title('AMPAR#:   EC    -   SPINE1   -     SPINE2  -    SPINES');
	end
end

if doPLOT3D == 1
	if mod(stepN, 100) == 0
	figure(1)
	subplot(5,5,[21 22]);
	bar(1:4, [EST'], .8, 'EdgeColor',[1 0.5 0.5]);
	set(gca,'Ylim',[0 GluR2Ndots],'XLim',[.5 4.5])
	title('AMPAR#:   EC    -   SPINE1   -     SPINE2  -    SPINES');
	end
end



if doGluR1	
	subplot(5,5,[16 17]), bar(1:4, [G1FSLOTS'], .8, 'EdgeColor',[1 0.5 0.5]);
		set(gca,'Ylim',[0 ymax],'XLim',[.5 4.5])
        title('GLUR1 IN SLOTS');
		xt = [G1PSD1slotN G1PSD2slotN G1PERI1slotN G1PERI2slotN];
		set(gca,'XTickLabel', sprintf('%i|',xt))
		%set(get(gca,'XLabel'),'String','\bullet \bullet  PSD  \bullet  \bullet PERI')
end
if doGluR2
	subplot(5,5,[21 22]), bar(1:4, [G2FSLOTS'], .8, 'EdgeColor',[1 0.5 0.5]);
		set(gca,'Ylim',[0 ymax],'XLim',[.5 4.5])
        title('GLUR2 IN SLOTS');
		xt = [G2PSD1slotN G2PSD2slotN G2PERI1slotN G2PERI2slotN];
		set(gca,'XTickLabel', sprintf('%i|',xt))
		set(get(gca,'XLabel'),'String',...
			strcat('PSD \bullet \bullet \bullet \bullet \bullet \bullet',...
			'\bullet \bullet \bullet \bullet \bullet \bullet PERI'))
end


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
%===============================%



%-------------------------------------------%
% PLOT Particle Motion
%-------------------------------------------%
% function [] = MAINPLOT(G1Ph1, G2Ph1, xyl)
%-------------------------------------------%

set(G2Ph1,'XData',xyl(1,:),'YData',xyl(2,:));
drawnow

set(G1Ph1,'XData',xyl(1,:),'YData',xyl(2,:));
drawnow

% end
%-------------------------------------------%



%===========================================%
% PREP GRAPHICS
%-------------------------------------------%
%{
% fig1 = figure(1);
% fig1box(fig1);
% datafig = figure(88);
% figure('Position',[1 SCREENSIZE(4)/8 SCREENSIZE(3)/8 SCREENSIZE(4)/8])
% set(fig1,'Renderer','painters')
% set(fig1,'DoubleBuffer','off')
% set(fig1,'WindowStyle','modal')
%}

figure(1)
fig1 = figure(1);
set(gcf,'OuterPosition',[300,100,900,800])


xlim = [0 XWIDE];
ylim = [-YHIGH YHIGH];

subplot(5,5,[3 25]), 
G2Ph1 = scatter(GluR2xyl(1,:),GluR2xyl(2,:),5,[0 0 1]);
hold on;
subplot(5,5,[3 25]), 
G1Ph1 = scatter(GluR1xyl(1,:),GluR1xyl(2,:),5,[1 0 0]);
axis([xlim, ylim]);
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(get(gca,'XLabel'),'String','GluR2:Blue \bullet \bullet GluR1:Red')
rectangle('Position',[XYBOXpr1(1),XYBOXpr1(2),XYBOXpr1(3),XYBOXpr1(4)])
rectangle('Position',[XYBOXpr2(1),XYBOXpr2(2),XYBOXpr2(3),XYBOXpr2(4)])
rectangle('Position',[XYBOXp1(1),XYBOXp1(2),XYBOXp1(3),XYBOXp1(4)])
rectangle('Position',[XYBOXp2(1),XYBOXp2(2),XYBOXp2(3),XYBOXp2(4)])
hold off;





%-- S1 Cluster Plots --
subplot(5,5,1); S1Ph1 = imagesc(S1);
colormap('hot');
title('SYNAPSE-1 CLUSTERS');

subplot(5,5,6),
S1Ph2 = imagesc(S1);
title('hk');

subplot(5,5,11),
S1Ph3 = imagesc(S1);
title('Pkon');

subplot(5,5,16),
S1Ph4 = imagesc(S1);
title('Sno');

subplot(5,5,21),
S1Ph5 = imagesc(S1);
title('TETHERED G1');

%-- S2 Cluster Plots --
subplot(5,5,2); 
S2Ph1 = imagesc(S2);
colormap('hot');
title('SYNAPSE-2 CLUSTERS');


subplot(5,5,7),
S2Ph2 = imagesc(S2);
title('hk');


subplot(5,5,12),
S2Ph3 = imagesc(S2);
title('Pkon');


subplot(5,5,17),
S2Ph4 = imagesc(S2);
title('Sno');


subplot(5,5,22),
S2Ph5 = imagesc(S2);
title('TETHERED G1');


%===========================================%



%%
%======================================================================%
%					SAVE FIGURE 1 PLOTS
%======================================================================%
set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,'STARShiP2','png');
saveas(gcf, ['outputfigs/FIGURE1.png']);
%======================================================================%


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



%===========================================================%
% FIG1 TOP RIGHT: GluR Subtypes in PSD vs PSA
%===========================================================%
sbpos = [.55 .77 .4 .18]; ptype = 4;
cOLOR = [c1; c2; c3; c4];
itemN = 1; 
[ph1 hax1] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'PSD');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 2;
[ph2 hax2] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'PSA');
[LEGH,OBJH,OUTH,OUTM] = legend;
%------------------------------------------%
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','Occupied Slots')
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 5;
set(ph1,'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
hTitle  = title ('GluR Subtypes in PSD vs PSA');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('GluR1/2');
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
sbpos = [.55 .57 .4 .18]; ptype = 4;
cOLOR = [c1; c2; c3; c4];
itemN = 3; 
[ph1 hax1] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'PSD');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 4;
[ph2 hax2] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'PSA');
[LEGH,OBJH,OUTH,OUTM] = legend;
%------------------------------------------%
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','Occupied Slots')
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 5;
set(ph1,'LineStyle','-','Color', c3,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
set(ph2,'LineStyle','-','Color',c4,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('-');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('GluR2/3');
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



%===========================================================%
% FIG1 BOTTOM LEFT: GluR Subtypes in Synapse 1v2
%===========================================================%
sbpos = [.05 .28 .4 .18]; ptype = 4;
cOLOR = [c1; c2; c3; c4; c1; c2; c3; c4];
itemN = 5; 
[ph1 hax1] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'Synapse 1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 6;
[ph2 hax2] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'Synapse 2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 2;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
hTitle  = title ('GluR Subtypes in Synapse 1v2');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('GluR1/2');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
% ylim([0 haxes(4)]);
% xlim([0 (haxes(2)*.9)]);
%===========================================================%
sbpos = [.05 .09 .4 .18]; ptype = 4;
cOLOR = [c1; c2; c3; c4; c1; c2; c1; c2];
itemN = 7; 
[ph1 hax1] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'Synapse 1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 8;
[ph2 hax2] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'Synapse 2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 2;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
hTitle  = title ('-');
hXLabel = xlabel('Time (min)');
hYLabel = ylabel('GluR2/3');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
haxes=axis;
% ylim([0 haxes(4)]);
% xlim([0 (haxes(2)*.9)]);
%======================================================================%





%===========================================================%
% FIG1 BOTTOM RIGHT: Synapse Particle Counts
%===========================================================%
sbpos = [.55 .09 .4 .38]; ptype = 4;
cOLOR = [c1; c2; c1; c2; c3; c2; c3; c4];
itemN = 3; 
[ph1 hax1] = CIenvFun(reDATAdataset,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'Synapse 1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 4;
[ph2 hax2] = CIenvFun(reDATAdataset,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'Synapse 2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 5; 
[ph3 hax3] = CIenvFun(reDATAdataset,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph3],OUTM{:},'Synapse Total');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
%------------------------------------------%
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','Occupied Slots')
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 5;
set(ph1,'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(ph3,'LineStyle','-','Color',c3,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
hTitle  = title ('Synapse Particle Counts');
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


%======================================================================%
set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,'STARShiP1','png');
saveas(gcf, ['outputfigs/STARShiP1.png']);
% printpreview
%======================================================================%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%======================================================================%
%					FINAL OUTPUT (SPLINE) FIGURE 2 OF 2
%======================================================================%
%							FIGURE SETUP
%----------------------------------------------------------------------%
fig12 = figure(12);
set(12,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig12,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
cOLOR = [c1; c2; c3; c4; c1; c2; c3; c4];



%===========================================================%
% FIG2 TOP LEFT: Synaptic AMPARs
%===========================================================%
sbpos = [.055 .57 .4 .38]; ptype = 4;
cOLOR = [c1; c4; c1; c2; c3; c3];
itemN = 3; 
[ph1 hax1] = CIenvFun(reDATAdataset,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'Synapse 1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 4;
[ph2 hax2] = CIenvFun(reDATAdataset,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'Synapse 2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 5; 
%------------------------------------------%
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','Occupied Slots')
% xt = (get(gca,'XTick'));
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 2;
set(ph1,'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
% set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(ph2,'LineStyle','-.','Color',c2,'LineWidth',3,...
'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c2);
hTitle  = title ('Synaptic AMPARs');
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
% ylim([haxes(3)/1.1 haxes(4)*1.1]);
% xlim([0 (haxes(2)*.9)]);
%======================================================================%




%===========================================================%
% WAS FIG1 BOTTOM LEFT: GluR Subtypes in Synapse 1v2
%===========================================================%
%===========================================================%
% FIG2 TOP RIGHT: GluR Subtypes in Synapse 1v2
%===========================================================%
% sbpos = [.05 .09 .4 .38]; ptype = 4;
sbpos = [.55 .57 .4 .38]; ptype = 4;
cOLOR = [c1; c2; c3; c4; c1; c2; c3; c4];
itemN = 5; 
[ph1 hax1] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'G1/2 Synapse 1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 6;
[ph2 hax2] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'G1/2 Synapse 2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 7; 
[ph3 hax3] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph3],OUTM{:},'G2/3 Synapse 1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 8;
[ph4 hax4] = CIenvFun(reDATAGluRdata,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph4],OUTM{:},'G2/3 Synapse 2');
hold on
%------------------------------------------%
set(get(gca,'XLabel'),'String','Time (min)')
set(get(gca,'YLabel'),'String','Occupied Slots')
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 5;
set(ph1,'LineStyle','-','Color', c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(ph3,'LineStyle','-','Color',c3,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
set(ph4,'LineStyle','-','Color',c4,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('GluR Subtypes in Synapse 1v2');
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
% ylim([0 haxes(4)]);
% xlim([0 (haxes(2)*.9)]);
%======================================================================%



%===========================================================%
% FIG2 BOTTOM RIGHT: Occupied Slots
%===========================================================%
sbpos = [.55 .09 .4 .38]; ptype = 4;
cOLOR = [c1; c2; c3; c4; c1; c2; c3; c4];
itemN = 1; 
[ph1 hax1] = CIenvFun(reDATAG1SLOTSDATA,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'G1/2 PSD 1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 2;
[ph2 hax2] = CIenvFun(reDATAG1SLOTSDATA,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'G1/2 PSD 2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
cOLOR = [c3; c4; c3; c4; c1; c2; c3; c4];
itemN = 1; 
[ph3 hax3] = CIenvFun(reDATAG2SLOTSDATA,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph3],OUTM{:},'G2/3 PSD 1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 2;
[ph4 hax4] = CIenvFun(reDATAG2SLOTSDATA,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph4],OUTM{:},'G2/3 PSD 2');
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 2;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(ph3,'LineStyle','-','Color',c3,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c3,'MarkerFaceColor',c33);
set(ph4,'LineStyle','-','Color',c4,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c4,'MarkerFaceColor',c44);
hTitle  = title ('Occupied Slots');
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
% ylim([0 haxes(4)*1.1 ]);
% xlim([0 (haxes(2)*.9)]);
%======================================================================%



%===========================================================%
% FIG2 BOTTOM LEFT: Synaptic SAP Cluster Size
%===========================================================%
sbpos = [.055 .09 .4 .38]; ptype = 4;
cOLOR = [c1; c2; c3; c4; c1; c2; c3; c4];
itemN = 1; 
[ph1 hax1] = CIenvFun(reDATASAPdata,sbpos,itemN,cOLOR,ptype);
leg1 = legend(ph1,'Synapse 1');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
itemN = 2;
[ph2 hax2] = CIenvFun(reDATASAPdata,sbpos,itemN,cOLOR,ptype);
legend([OUTH;ph2],OUTM{:},'Synapse 2');
[LEGH,OBJH,OUTH,OUTM] = legend;
hold on
xt = (get(gca,'XTick'))*AveOver*DATARATE*(t)/(60);
set(gca,'XTickLabel', sprintf('%.0f|',xt))
set(legend,'Location','NorthWest');
%------------------------------------------%
MS1 = 5; MS2 = 2;
set(ph1,'LineStyle','-','Color',c1,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c1,'MarkerFaceColor',c11);
% set(ph2,'LineStyle',':','Color',c2,'LineWidth',1,...
% 'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
set(ph2,'LineStyle','-.','Color',c2,'LineWidth',3,...
'Marker','none','MarkerSize',MS1,'MarkerEdgeColor',c2);
hTitle  = title ('Synaptic SAP Cluster Size');
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
% ylim([haxes(3)/1.1 haxes(4)*1.1 ]);
% xlim([0 (haxes(2)-3)]);
% ylim([min(min(reDATASAPdata{2}))/2 max(max(reDATASAPdata{2}))*1.2]);
% xlim([0 (haxes(2)*.9)]);
%======================================================================%


%======================================================================%
set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,'STARShiP2','png');
saveas(gcf, ['outputfigs/STARShiP2.png']);




%-------------------------------------------%
% GAUSSIAN TETHERING COLORMAP
%-------------------------------------------%
function [] = SLOTMAPFUN(G1STBASE,G1RTBASE,G1STLTP,G1RTLTP,G1BSMu,G1LSMu,...
G2STBASE,G2RTBASE,G2STLTP,G2RTLTP,G2BSMu,G2LSMu,GT,GTab)

saps = 0:4;
tails = 0:.5:4;
offsap = G1BSMu;

%param = [saps+tails(1) saps+tails(2) saps+tails(3) saps+tails(4) saps+tails(5)]-6

param = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G1BSMu;

cdfMx = 1-cdf('norm',0,param,1);
cdfMx2 = [linspace(0,1,numel(tails)); cdfMx];

cdfMxpad = padarray(cdfMx,[1 1],0,'pre');
cdfMxpad(:,1)=[0 saps]';
cdfMxpad(1,:)=[0 tails];

ProbCDF = cdfMx2(1,:);

	
%=========================================================%
% FIGURE SETUP
%---------------------------------------------------------%
fig99 = figure(99);
figure(fig99)
set(99,'Units','pixels') 
scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/1.8  scnsize(4)/7  scnsize(3)/2.3  scnsize(4)/1.2];
% OuterPosition(left, bottom, width, height)
set(fig99,'OuterPosition',pos1)
fig99 = figure(99);
figure(fig99)
set(gcf,'Color',[.9,.9,.9])
%=========================================================%
yt = [-1 saps];
xt = tails; 

subplot(5,5,[16 22]), subplot('Position',[.075 .052 .7 .3]),...
imagesc(cdfMx2)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
set(gca,'YTickLabel', sprintf('%.1f|',yt))
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(.8,1.8,strcat('\color[rgb]{.7 .7 .5}Gaussian \phi = ',...
	num2str(offsap)),'FontSize',18,'HorizontalAlignment','left');
text(6,1,strcat('REFERENCE ROW'),'FontSize',14.5,'HorizontalAlignment','center');
% text(0,1,strcat('P(on)= ', num2str(ProbCDF(:)')),'FontSize',14.5,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
colorbar('location','North')
text(.5,1,strcat('-------------------------------------------------------'),...
	'FontSize',28,'HorizontalAlignment','left',...
	'BackgroundColor',[1 1 1]);
hold on

G1STBL = plot((G1STBASE*2+1),2:6,...
	'MarkerFaceColor',[.9 .3 .5],'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',8,'Marker','d');
hold on
G1STLL = plot((G1STLTP*2+1),2:6,'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',24,'Marker','+');
hold on
G2STBL = plot((G2STBASE*2+1),2:6,...
	'MarkerFaceColor',[.2 .7 .8],'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',8,'Marker','o');
hold on
G2STLL = plot((G2STLTP*2+1),2:6,'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',18,'Marker','x');
hold on
STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');







%=========================================================%
% CDF SETUP		% P(G|S) = P(x??):{x?N(?+?,1)}
%---------------------------------------------------------%
saps = 0:4; tails = 0:.5:4;

G1BSMuTSmx = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G1BSMu;

G1LSMuTSmx = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G1LSMu;
	 
G2BSMuTSmx = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G2BSMu;

G2LSMuTSmx = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G2LSMu;
	 
%--------------------------------------%
% G1BSMu - G1 BASE 
%--------------------------------------%
G1BSMuCDFmx = 1-cdf('norm',0,G1BSMuTSmx,1);
G1BSMuCDF = [linspace(0,1,numel(tails)); G1BSMuCDFmx];

G1BSMuCDFmxpad = padarray(G1BSMuCDFmx,[1 1],0,'pre');
G1BSMuCDFmxpad(:,1)=[0 saps]';
G1BSMuCDFmxpad(1,:)=[0 tails];
G1BSMuCDFp = G1BSMuCDF(1,:);
%--------------------------------------%
%--------------------------------------%
% G1LSMu - G1 LTP 
%--------------------------------------%
G1LSMuCDFmx = 1-cdf('norm',0,G1LSMuTSmx,1);
G1LSMuCDF = [linspace(0,1,numel(tails)); G1LSMuCDFmx];

G1LSMuCDFmxpad = padarray(G1LSMuCDFmx,[1 1],0,'pre');
G1LSMuCDFmxpad(:,1)=[0 saps]';
G1LSMuCDFmxpad(1,:)=[0 tails];
G1LSMuCDFp = G1LSMuCDF(1,:);
%--------------------------------------%

%--------------------------------------%
% G2BSMu - G2 BASE 
%--------------------------------------%
G2BSMuCDFmx = 1-cdf('norm',0,G2BSMuTSmx,1);
G2BSMuCDF = [linspace(0,1,numel(tails)); G2BSMuCDFmx];

G2BSMuCDFmxpad = padarray(G2BSMuCDFmx,[1 1],0,'pre');
G2BSMuCDFmxpad(:,1)=[0 saps]';
G2BSMuCDFmxpad(1,:)=[0 tails];
G2BSMuCDFp = G2BSMuCDF(1,:);
%--------------------------------------%
%--------------------------------------%
% G2BSMu - G2 LTP 
%--------------------------------------%
G2LSMuCDFmx = 1-cdf('norm',0,G2LSMuTSmx,1);
G2LSMuCDF = [linspace(0,1,numel(tails)); G2LSMuCDFmx];

G2LSMuCDFmxpad = padarray(G2LSMuCDFmx,[1 1],0,'pre');
G2LSMuCDFmxpad(:,1)=[0 saps]';
G2LSMuCDFmxpad(1,:)=[0 tails];
G2LSMuCDFp = G2LSMuCDF(1,:);
%--------------------------------------%


%=========================================================%
%					FIGURE SETUP
%---------------------------------------------------------%
fig97 = figure(97);figure(fig97)
set(97,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig97,'OuterPosition',pos1)
fig97 = figure(97);
figure(fig97)
set(gcf,'Color',[.9,.9,.9])
%=========================================================%
yt = [-1 saps];
xt = tails;
% G1Bcolor=[.9 .3 .5];
% G2Bcolor=[.2 .7 .8];
G1Bcolor=[.9 .1 .6];
G2Bcolor=[.2 .7 .8];


%======================================%
% G1BS Mu COLORMAP FIGURE
%--------------------------------------%
subplot('Position',[.06 .07 .44 .40]),...
imagesc(G1BSMuCDF)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
%set(gca,'YTick',0:4)
%set(gca,'YTickLabel', sprintf('%.1f|',yt),'FontSize',14)
set(gca,'XTickLabel', sprintf('%.1f|',xt),'FontSize',14)
set(gca,'YTickLabel', ['N' sprintf('%.0f|',yt(2)) sprintf('%.1f|',yt(2:end))],'FontSize',14)
set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(1,1,strcat('\color[rgb]{.7 .8 .9}G1 basal \phi = ',...
	num2str(G1BSMu)),'FontSize',16,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
%colorbar('location','North')
hold on
%-------------------------%
G1STBL = plot((G1STBASE*2+1),2:6,...
	'MarkerFaceColor',G1Bcolor,'MarkerEdgeColor',G1Bcolor,'MarkerSize',10,'Marker','o');
hold on
% G1STLL = plot((G1STLTP*2+1),2:6,'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',24,'Marker','+');
% hold on
% G2STBL = plot((G2STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.2 .7 .8],'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',8,'Marker','o');
% hold on
% G2STLL = plot((G2STLTP*2+1),2:6,'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',18,'Marker','x');
% hold on
% STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
% 	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
% 	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
% set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');
% SAPMAPTITLE=title(char('Gaussian: ', num2str(offsap)));
% set(SAPMAPTITLE, 'FontSize', 14);
%======================================%





%======================================%
% G1LS Mu COLORMAP FIGURE
%--------------------------------------%
subplot('Position',[.55 .07 .44 .40]),...
imagesc(G1LSMuCDF)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
%set(gca,'YTick',0:4)
% set(gca,'YTickLabel', sprintf('%.1f|',yt),'FontSize',14)
set(gca,'XTickLabel', sprintf('%.1f|',xt),'FontSize',14)
set(gca,'YTickLabel', ['N' sprintf('%.0f|',yt(2)) sprintf('%.1f|',yt(2:end))],'FontSize',14)
% set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
% set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(1,1,strcat('\color[rgb]{.7 .8 .9}G1 LTP \phi = ',...
	num2str(G1LSMu)),'FontSize',16,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
%colorbar('location','North')
hold on
%-------------------------%
% G1STBL = plot((G1STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.9 .3 .5],'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',8,'Marker','d');
% hold on
G1STLL = plot((G1STLTP*2+1),2:6,...
	'MarkerFaceColor',G1Bcolor,'MarkerEdgeColor',G1Bcolor,'MarkerSize',10,'Marker','o');
hold on
% G2STBL = plot((G2STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.2 .7 .8],'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',8,'Marker','o');
% hold on
% G2STLL = plot((G2STLTP*2+1),2:6,'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',18,'Marker','x');
% hold on
% STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
% 	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
% 	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
% set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');
% SAPMAPTITLE=title(char('Gaussian: ', num2str(offsap)));
% set(SAPMAPTITLE, 'FontSize', 14);
%======================================%





%======================================%
% G2BS Mu COLORMAP FIGURE
%--------------------------------------%
subplot('Position',[.06 .57 .44 .40]),...
imagesc(G2BSMuCDF)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
%set(gca,'YTick',0:4)
% set(gca,'YTickLabel', sprintf('%.1f|',yt),'FontSize',14)
set(gca,'XTickLabel', sprintf('%.1f|',xt),'FontSize',14)
set(gca,'YTickLabel', ['N' sprintf('%.0f|',yt(2)) sprintf('%.1f|',yt(2:end))],'FontSize',14)
% set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
% set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(1,1,strcat('\color[rgb]{.7 .8 .9}G2 basal \phi = ',...
	num2str(G2BSMu)),'FontSize',16,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
%colorbar('location','North')
hold on
%-------------------------%
% G1STBL = plot((G1STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.9 .3 .5],'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',8,'Marker','d');
% hold on
% G1STLL = plot((G1STLTP*2+1),2:6,'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',24,'Marker','+');
% hold on
G2STBL = plot((G2STBASE*2+1),2:6,...
	'MarkerFaceColor',G2Bcolor,'MarkerEdgeColor',G2Bcolor,'MarkerSize',10,'Marker','o');
hold on
% G2STLL = plot((G2STLTP*2+1),2:6,'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',18,'Marker','x');
% hold on
% STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
% 	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
% 	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
% set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');
% SAPMAPTITLE=title(char('Gaussian: ', num2str(offsap)));
% set(SAPMAPTITLE, 'FontSize', 14);
%======================================%





%======================================%
% G2LS Mu COLORMAP FIGURE
%--------------------------------------%
subplot('Position',[.55 .57 .44 .40]),...
imagesc(G2LSMuCDF)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
%set(gca,'YTick',0:4)
% set(gca,'YTickLabel', sprintf('%.1f|',yt),'FontSize',14)
set(gca,'XTickLabel', sprintf('%.1f|',xt),'FontSize',14)
set(gca,'YTickLabel', ['N' sprintf('%.0f|',yt(2)) sprintf('%.1f|',yt(2:end))],'FontSize',14)
% set(gca,'YTickLabel', '-')
% set(gca,'XTickLabel', '-')
%axis off

% set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
% set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(1,1,strcat('\color[rgb]{.7 .8 .9}G2 LTP \phi = ',...
	num2str(G2LSMu)),'FontSize',16,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
%colorbar('location','North')
hold on
%-------------------------%
% G1STBL = plot((G1STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.9 .3 .5],'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',8,'Marker','d');
% hold on
% G1STLL = plot((G1STLTP*2+1),2:6,'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',24,'Marker','+');
% hold on
% G2STBL = plot((G2STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.2 .7 .8],'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',8,'Marker','o');
% hold on
G2STLL = plot((G2STLTP*2+1),2:6,...
	'MarkerFaceColor',G2Bcolor,'MarkerEdgeColor',G2Bcolor,'MarkerSize',10,'Marker','d');
hold on
% STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
% 	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
% 	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
% set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');
% SAPMAPTITLE=title(char('Gaussian: ', num2str(offsap)));
% set(SAPMAPTITLE, 'FontSize', 14);
%======================================%


end


%-------------------------------------------%
% GAUSSIAN RECRUITMENT COLORMAP
%-------------------------------------------%
function [] = SAPRECFUN(G1STBASE,G1RTBASE,G1STLTP,G1RTLTP,G1BSMu,G1LSMu,...
G2STBASE,G2RTBASE,G2STLTP,G2RTLTP,G2BSMu,G2LSMu,...
S1,S1sz, sap, stepN,GT,GTab,...
G1P1SAPM, G1SP1, LTP1onG1, LTP1offG1, LTP2onG1, LTP2offG1)
 
%---------------------------------------------------------%
% INITIALIZE CONSTANTS
%---------------------------------------------------------%
% sap1 = sap(1); % SAPdotsPSD1
% sap2 = sap(2); % SAPdotsPSD2
% sap3 = sap(3); % SAPbetaPSD1
% sap4 = sap(4); % SAPtauPSD1
% sap5 = sap(5); % SAPL1PSD1
% sap6 = sap(6); % SAPbetaPSD2
% sap7 = sap(7); % SAPtauPSD2
% sap8 = sap(8); % SAPL1PSD2
% sap9 = sap(9); % SAPmuPSD1
% sap10 = sap(10); % SAPmuPSD2
% sap11 = sap(11); % SAPdTPSD1
% sap12 = sap(12); % SAPdTPSD2
% sap13 = sap(13); % SAPrhoPSD1
% sap14 = sap(14); % SAPrhoPSD2
% sap15 = sap(15); % SAPrPSD1
% sap16 = sap(16); % SAPrPSD2
% sap17 = sap(17); % doDynamicLeP1
% sap18 = sap(18); % doDynamicLeP2
% sap19 = sap(19); % doLTPS1
% sap20 = sap(20); % doLTPS2
 
 
 


%%
%---------------------------------------------------------%
%			SAP CLUSTER MATRIX PROBABILITIES
%---------------------------------------------------------%
S1Le = sap(5);
doDynamicLe = sap(17);
if doDynamicLe
L1x = linspace(1.0,3.2);
Tsz = numel(S1);			% Total matrix space
Ssz = numel(find(S1>0));	% Filled matrix space
Psz = round(Ssz/Tsz*100)+1;	% Percent filled space
S1Le = L1x(Psz);
end
%---------------------------------------------------------%
%	TEST MASK AGAINST CLUSTER - GET CONVOLUTION MX
%---------------------------------------------------------%

S = S1;
S(8,7) = 0; S(10,10) = 0; S(5,12) = 0; S(3,14) = 1;


%-------------------------------%
hkMask=[0 1 0; 1 0 1; 0 1 0];
dT = sap(11);
Le = S1Le;		% hi = shrink
B = sap(3);		% hi = shrink
rho = sap(13);	% hi = grow
r = sap(15);	% hi = grow
mu = sap(9);	% hi = shrink
%-------------------------------%
Pmx = rand(size(S));
Soc = (S>0);
Snoc = (S<1);
Sno = ~Soc;
%---
hk = convn(Soc,hkMask,'same');
Lhk = hk-Le;
%---
Pex = 1 ./ (1+exp((-1*B)*Lhk));
Pkex = (1-Soc) .* ( rho * r * dT * Pex );
Pen = (Soc) .* mu * dT;
%---
S0 = S;
Pex0 = Pex;
Pkex0 = Pkex;
Pen0 = Pen;

%=========================================================%
% figure(9);
% set(gcf,'OuterPosition',[50 50 800 400])
% subplot('Position',[.05 .06 .45 .90]),...
% imagesc(PhkSnoc);
% colormap('bone')
% colorbar('location','North')
% subplot('Position',[.50 .06 .45 .90]),...
% imagesc(Phk);
% caxis(caxis);
% cbr = colorbar('location','North');
% set(cbr(1), 'XColor', 'r');
%=========================================================%

stepN=2;
%---------------------------------------------%
% GluA RT SAP recruitment (basal and LTP)
%---------------------------------------------%
G1RT = G1RTBASE;
if (stepN >= LTP1onG1 && stepN <= LTP1offG1); G1RT=G1RTLTP; end;
if (stepN >= LTP2onG1 && stepN <= LTP2offG1); G1RT=G1RTLTP; end;

G1P1SAPM = [1:3;4:6;7:9];
G1P1SAPM(1,2) = 18;
G1P1SAPM(2,2) = 12;
G1P1SAPM(3,2) = 28;
% close all;
%---------------------------------------------------------%
%       Get Matrix Locations of SAPs near AMPARs
%---------------------------------------------------------%
G1P1poly0 = G1P1SAPM(:,2);
[~, ~, G1P1poly1] = find(G1P1poly0);
G1P1poly2 = G1P1poly1;

for nx = 1:(numel(G1P1poly1))
    if G1P1poly1(nx)<=8 && G1P1poly1(nx)>0
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*0)-1;
    end
    if G1P1poly1(nx)<=16 && G1P1poly1(nx)>8
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*1)-1;
    end
    if G1P1poly1(nx)<=24 && G1P1poly1(nx)>16
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*2)-1;
    end
    if G1P1poly1(nx)<=32 && G1P1poly1(nx)>24
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*3)-1;
    end
    if G1P1poly1(nx)<=40 && G1P1poly1(nx)>32
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*4)-1;
    end
    if G1P1poly1(nx)<=48 && G1P1poly1(nx)>40
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*5)-1;
    end
    if G1P1poly1(nx)<=56 && G1P1poly1(nx)>48
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*6)-1;
    end
    if G1P1poly1(nx)<=64 && G1P1poly1(nx)>56
        G1P1poly2(nx) = (G1P1poly1(nx)*2)+(S1sz*7)-1;
    end
end % for nx
%---------------------------------------------%

SG1oc = (S>2) .* 0.0;
G1P1n = numel(G1P1poly2);

%-------------------------------%
for polyid = 1:G1P1n
%-------------------------------%
G1P1IDo = G1P1poly2(polyid);

G1P1IDb = G1P1IDo+1;
G1P1IDr = G1P1IDo+16;
G1P1IDd = G1P1IDo+17;

if G1P1IDo>=16
	G1P1IDb = G1P1IDo+1;
	G1P1IDr = G1P1IDo+16+(round(rand)*-32); 
	G1P1IDd = G1P1IDr+1;
	% 50:50 WHETHER EAST OR WEST SAP
	% ALWAYS EAST SAP IF COMMENTED
	%G1P1IDd = G1P1IDo+17+(round(rand)*-34);
end
SG1oc(G1P1IDo)=1;
SG1oc(G1P1IDb)=1;
SG1oc(G1P1IDr)=1;
SG1oc(G1P1IDd)=1;
%-------------------------------%
end % for polyid = 1:G1P1n
%-------------------------------%

% Sno		SP non-occupied lattice locations
% SG1oc		G1 occupied lattice locations
% G1RT		G1 recruitment tail constant
% Pkex		SP P(exo) lattice

Gex = Sno .* SG1oc .* G1RT + Pkex;
Gp = 1 ./ (1+exp((-1*B)*Lhk)) + .0001;
Gpex = Gex .* Gp;

% Phkv	= 1-(1./(hk+G1RT+(1/(1+Le))));
% Gpex = Sno .* SG1oc .* Phkv;

%---
% Pkexo = Sno .* G1P1ex;
% Pkexo = Pkexo + Pkex0;

%---
Pkex = Gpex;
Pen = Pen0;
%---

%---------------------------------------------%
Sen = (Pen > Pmx);
Sex = (Pkex > Pmx);
S = (Soc-Sen) + Sex;
%---------------------------------------------%
%=========================================================%
figure(8) 
set(gcf,'OuterPosition',[500 55 1200 800])

%================================================%



% TOP LEFT (S0)
subplot('Position',[.02 .54 .3 .45]),...
imagesc(S0);
colormap('jet');
cbr = colorbar('location','South');
set(cbr(1), 'XColor', 'r');
text(1,1.5,'S0','FontSize',22,'Color',[.7 1 .7]);


% TOP MID (Pkex)
subplot('Position',[.35 .54 .3 .45]),...
imagesc(Pkex0);
caxis(caxis);
cbr = colorbar('location','South');
set(cbr(1), 'XColor', 'r');
text(1,1.5,'Pkex','FontSize',22,'Color',[.7 1 .7]);


% TOP RIGHT (Pen0)
subplot('Position',[.68 .54 .3 .45]),...
imagesc(Pen0);
colormap('jet');
cbr = colorbar('location','South');
set(cbr(1), 'XColor', 'r');
text(1,1.5,'Pen0','FontSize',22,'Color',[.7 1 .7]);


%=======================================%

hk = Sno .* hk;

% BOTTOM LEFT (hk)
subplot('Position',[.02 .05 .3 .45]),...
imagesc(hk);
caxis(caxis);
% cbr = colorbar('location','South');
% set(cbr(1), 'XColor', 'g');
% text(1,1.5,'hk','FontSize',22,'Color',[.7 1 .7]);


% BOTTOM MID (SG1oc)
subplot('Position',[.35 .05 .3 .45]),...
scfg1 = imagesc(SG1oc); hold on; 
colormap('jet');
cbr = colorbar('location','South');
set(cbr(1), 'XColor', 'g');
scfg2 = imagesc(Lhk);
alpha(.75);
hold off
text(1,1.5,'SG1oc + Lhk','FontSize',22,'Color',[.7 1 .7]);


% BOTTOM RIGHT (S)
subplot('Position',[.68 .05 .3 .45]),...
imagesc(S);
caxis(caxis);
cbr = colorbar('location','South');
set(cbr(1), 'XColor', 'g');
text(1,1.5,'S','FontSize',22,'Color',[.7 1 .7]);



%================================================%
% % TOP LEFT (S0)
% subplot('Position',[.02 .54 .3 .45]),...
% imagesc(S0);
% colormap('jet');
% cbr = colorbar('location','South');
% set(cbr(1), 'XColor', 'r');
% % xt = (get(gca,'XTick'))-.5;
% % set(gca,'XTickLabel', sprintf('%.0f|',xt))
% % yt = (get(gca,'YTick'))-.5;
% % set(gca,'YTickLabel', sprintf('%.0f|',yt))
% % haxes=axis;
% % ylim([haxes(3)+.5 haxes(4)-.5 ]);
% % xlim([haxes(1)+.5 haxes(2)-.5 ]);
% % grid on
% text(1,1.5,'S0','FontSize',22,'Color',[.7 1 .7]);
% text(XYLTpr1(1),XYLTpr1(2),...
% strcat(num2str(XYLTpr1(1)), '\bullet',num2str(XYLTpr1(2)),'\rightarrow'),...
% 'FontSize',12,'HorizontalAlignment','right');
% alpha('color'); alphamap('rampdown'); alphamap('increase',.1)
%=========================================================%
%%
%---------------------------%
% if stepN>500; keyboard; end

%---------------------------%
%======================================================================%
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, ['outputfigs/CLUSTREC.png']);
%======================================================================%
%%



%=========================================================%
figure(9);
set(gcf,'OuterPosition',[50 50 800 400])
subplot('Position',[.05 .06 .45 .90]),...
imagesc(PhkSnoc);
% colormap('gray');
colormap('jet');
colorbar('location','North')
subplot('Position',[.50 .06 .45 .90]),...
imagesc(Phk);
caxis(caxis);
cbr = colorbar('location','North');
set(cbr(1), 'XColor', 'r');
%=========================================================%












%=========================================================%
%=========================================================%
saps = 0:4;
tails = 0:.5:4;
offsap = G1BSMu;

param = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G1BSMu;

cdfMx = 1-cdf('norm',0,param,1);
cdfMx2 = [linspace(0,1,numel(tails)); cdfMx];

cdfMxpad = padarray(cdfMx,[1 1],0,'pre');
cdfMxpad(:,1)=[0 saps]';
cdfMxpad(1,:)=[0 tails];

ProbCDF = cdfMx2(1,:);

	
%=========================================================%
% FIGURE SETUP
%---------------------------------------------------------%
fig99 = figure(99);
figure(fig99)
set(99,'Units','pixels') 
scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/1.8  scnsize(4)/7  scnsize(3)/2.3  scnsize(4)/1.2];
% OuterPosition(left, bottom, width, height)
set(fig99,'OuterPosition',pos1)
fig99 = figure(99);
figure(fig99)
set(gcf,'Color',[.9,.9,.9])
%=========================================================%
yt = [-1 saps];
xt = tails; 

subplot(5,5,[16 22]), subplot('Position',[.075 .052 .7 .3]),...
imagesc(cdfMx2)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
set(gca,'YTickLabel', sprintf('%.1f|',yt))
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(.8,1.8,strcat('\color[rgb]{.7 .7 .5}Gaussian \phi = ',...
	num2str(offsap)),'FontSize',18,'HorizontalAlignment','left');
text(6,1,strcat('REFERENCE ROW'),'FontSize',14.5,'HorizontalAlignment','center');
% text(0,1,strcat('P(on)= ', num2str(ProbCDF(:)')),'FontSize',14.5,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
colorbar('location','North')
text(.5,1,strcat('-------------------------------------------------------'),...
	'FontSize',28,'HorizontalAlignment','left',...
	'BackgroundColor',[1 1 1]);
hold on

G1STBL = plot((G1STBASE*2+1),2:6,...
	'MarkerFaceColor',[.9 .3 .5],'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',8,'Marker','d');
hold on
G1STLL = plot((G1STLTP*2+1),2:6,'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',24,'Marker','+');
hold on
G2STBL = plot((G2STBASE*2+1),2:6,...
	'MarkerFaceColor',[.2 .7 .8],'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',8,'Marker','o');
hold on
G2STLL = plot((G2STLTP*2+1),2:6,'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',18,'Marker','x');
hold on
STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');







%=========================================================%
% CDF SETUP		% P(G|S) = P(x??):{x?N(?+?,1)}
%---------------------------------------------------------%
saps = 0:4; tails = 0:.5:4;

G1BSMuTSmx = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G1BSMu;

G1LSMuTSmx = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G1LSMu;
	 
G2BSMuTSmx = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G2BSMu;

G2LSMuTSmx = [tails+saps(1); tails+saps(2); tails+saps(3); tails+saps(4); tails+saps(5)]...
		-G2LSMu;
	 
%--------------------------------------%
% G1BSMu - G1 BASE 
%--------------------------------------%
G1BSMuCDFmx = 1-cdf('norm',0,G1BSMuTSmx,1);
G1BSMuCDF = [linspace(0,1,numel(tails)); G1BSMuCDFmx];

G1BSMuCDFmxpad = padarray(G1BSMuCDFmx,[1 1],0,'pre');
G1BSMuCDFmxpad(:,1)=[0 saps]';
G1BSMuCDFmxpad(1,:)=[0 tails];
G1BSMuCDFp = G1BSMuCDF(1,:);
%--------------------------------------%
%--------------------------------------%
% G1LSMu - G1 LTP 
%--------------------------------------%
G1LSMuCDFmx = 1-cdf('norm',0,G1LSMuTSmx,1);
G1LSMuCDF = [linspace(0,1,numel(tails)); G1LSMuCDFmx];

G1LSMuCDFmxpad = padarray(G1LSMuCDFmx,[1 1],0,'pre');
G1LSMuCDFmxpad(:,1)=[0 saps]';
G1LSMuCDFmxpad(1,:)=[0 tails];
G1LSMuCDFp = G1LSMuCDF(1,:);
%--------------------------------------%

%--------------------------------------%
% G2BSMu - G2 BASE 
%--------------------------------------%
G2BSMuCDFmx = 1-cdf('norm',0,G2BSMuTSmx,1);
G2BSMuCDF = [linspace(0,1,numel(tails)); G2BSMuCDFmx];

G2BSMuCDFmxpad = padarray(G2BSMuCDFmx,[1 1],0,'pre');
G2BSMuCDFmxpad(:,1)=[0 saps]';
G2BSMuCDFmxpad(1,:)=[0 tails];
G2BSMuCDFp = G2BSMuCDF(1,:);
%--------------------------------------%
%--------------------------------------%
% G2BSMu - G2 LTP 
%--------------------------------------%
G2LSMuCDFmx = 1-cdf('norm',0,G2LSMuTSmx,1);
G2LSMuCDF = [linspace(0,1,numel(tails)); G2LSMuCDFmx];

G2LSMuCDFmxpad = padarray(G2LSMuCDFmx,[1 1],0,'pre');
G2LSMuCDFmxpad(:,1)=[0 saps]';
G2LSMuCDFmxpad(1,:)=[0 tails];
G2LSMuCDFp = G2LSMuCDF(1,:);
%--------------------------------------%


%=========================================================%
%					FIGURE SETUP
%---------------------------------------------------------%
fig97 = figure(97);figure(fig97)
set(97,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig97,'OuterPosition',pos1)
fig97 = figure(97);
figure(fig97)
set(gcf,'Color',[.9,.9,.9])
%=========================================================%
yt = [-1 saps];
xt = tails;
% G1Bcolor=[.9 .3 .5];
% G2Bcolor=[.2 .7 .8];
G1Bcolor=[.9 .1 .6];
G2Bcolor=[.2 .7 .8];


%======================================%
% G1BS Mu COLORMAP FIGURE
%--------------------------------------%
subplot('Position',[.06 .07 .44 .40]),...
imagesc(G1BSMuCDF)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
%set(gca,'YTick',0:4)
%set(gca,'YTickLabel', sprintf('%.1f|',yt),'FontSize',14)
set(gca,'XTickLabel', sprintf('%.1f|',xt),'FontSize',14)
set(gca,'YTickLabel', ['N' sprintf('%.0f|',yt(2)) sprintf('%.1f|',yt(2:end))],'FontSize',14)
set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(1,1,strcat('\color[rgb]{.7 .8 .9}G1 basal \phi = ',...
	num2str(G1BSMu)),'FontSize',16,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
%colorbar('location','North')
hold on
%-------------------------%
G1STBL = plot((G1STBASE*2+1),2:6,...
	'MarkerFaceColor',G1Bcolor,'MarkerEdgeColor',G1Bcolor,'MarkerSize',10,'Marker','o');
hold on
% G1STLL = plot((G1STLTP*2+1),2:6,'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',24,'Marker','+');
% hold on
% G2STBL = plot((G2STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.2 .7 .8],'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',8,'Marker','o');
% hold on
% G2STLL = plot((G2STLTP*2+1),2:6,'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',18,'Marker','x');
% hold on
% STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
% 	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
% 	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
% set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');
% SAPMAPTITLE=title(char('Gaussian: ', num2str(offsap)));
% set(SAPMAPTITLE, 'FontSize', 14);
%======================================%





%======================================%
% G1LS Mu COLORMAP FIGURE
%--------------------------------------%
subplot('Position',[.55 .07 .44 .40]),...
imagesc(G1LSMuCDF)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
%set(gca,'YTick',0:4)
% set(gca,'YTickLabel', sprintf('%.1f|',yt),'FontSize',14)
set(gca,'XTickLabel', sprintf('%.1f|',xt),'FontSize',14)
set(gca,'YTickLabel', ['N' sprintf('%.0f|',yt(2)) sprintf('%.1f|',yt(2:end))],'FontSize',14)
% set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
% set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(1,1,strcat('\color[rgb]{.7 .8 .9}G1 LTP \phi = ',...
	num2str(G1LSMu)),'FontSize',16,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
%colorbar('location','North')
hold on
%-------------------------%
% G1STBL = plot((G1STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.9 .3 .5],'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',8,'Marker','d');
% hold on
G1STLL = plot((G1STLTP*2+1),2:6,...
	'MarkerFaceColor',G1Bcolor,'MarkerEdgeColor',G1Bcolor,'MarkerSize',10,'Marker','o');
hold on
% G2STBL = plot((G2STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.2 .7 .8],'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',8,'Marker','o');
% hold on
% G2STLL = plot((G2STLTP*2+1),2:6,'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',18,'Marker','x');
% hold on
% STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
% 	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
% 	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
% set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');
% SAPMAPTITLE=title(char('Gaussian: ', num2str(offsap)));
% set(SAPMAPTITLE, 'FontSize', 14);
%======================================%





%======================================%
% G2BS Mu COLORMAP FIGURE
%--------------------------------------%
subplot('Position',[.06 .57 .44 .40]),...
imagesc(G2BSMuCDF)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
%set(gca,'YTick',0:4)
% set(gca,'YTickLabel', sprintf('%.1f|',yt),'FontSize',14)
set(gca,'XTickLabel', sprintf('%.1f|',xt),'FontSize',14)
set(gca,'YTickLabel', ['N' sprintf('%.0f|',yt(2)) sprintf('%.1f|',yt(2:end))],'FontSize',14)
% set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
% set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(1,1,strcat('\color[rgb]{.7 .8 .9}G2 basal \phi = ',...
	num2str(G2BSMu)),'FontSize',16,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
%colorbar('location','North')
hold on
%-------------------------%
% G1STBL = plot((G1STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.9 .3 .5],'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',8,'Marker','d');
% hold on
% G1STLL = plot((G1STLTP*2+1),2:6,'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',24,'Marker','+');
% hold on
G2STBL = plot((G2STBASE*2+1),2:6,...
	'MarkerFaceColor',G2Bcolor,'MarkerEdgeColor',G2Bcolor,'MarkerSize',10,'Marker','o');
hold on
% G2STLL = plot((G2STLTP*2+1),2:6,'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',18,'Marker','x');
% hold on
% STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
% 	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
% 	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
% set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');
% SAPMAPTITLE=title(char('Gaussian: ', num2str(offsap)));
% set(SAPMAPTITLE, 'FontSize', 14);
%======================================%





%======================================%
% G2LS Mu COLORMAP FIGURE
%--------------------------------------%
subplot('Position',[.55 .57 .44 .40]),...
imagesc(G2LSMuCDF)
%axis image
colormap('bone')
%set(gca,'XTick',0:.5:4)
%set(gca,'YTick',0:4)
% set(gca,'YTickLabel', sprintf('%.1f|',yt),'FontSize',14)
set(gca,'XTickLabel', sprintf('%.1f|',xt),'FontSize',14)
set(gca,'YTickLabel', ['N' sprintf('%.0f|',yt(2)) sprintf('%.1f|',yt(2:end))],'FontSize',14)
% set(gca,'YTickLabel', '-')
% set(gca,'XTickLabel', '-')
%axis off

% set(get(gca,'YLabel'),'String','SAPS IN VACINITY \Theta','FontSize',18)
% set(get(gca,'XLabel'),'String','TETHERING CONSTANT \Gamma','FontSize',18)
caxis(caxis)
hold on
text(1,1,strcat('\color[rgb]{.7 .8 .9}G2 LTP \phi = ',...
	num2str(G2LSMu)),'FontSize',16,'HorizontalAlignment','left');
hold on
colorbar('YTickLabel',{num2str(linspace(0,1,11)')})
%colorbar('location','North')
hold on
%-------------------------%
% G1STBL = plot((G1STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.9 .3 .5],'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',8,'Marker','d');
% hold on
% G1STLL = plot((G1STLTP*2+1),2:6,'MarkerEdgeColor',[.7 .3 .5],'MarkerSize',24,'Marker','+');
% hold on
% G2STBL = plot((G2STBASE*2+1),2:6,...
% 	'MarkerFaceColor',[.2 .7 .8],'MarkerEdgeColor',[.2 .7 .8],'MarkerSize',8,'Marker','o');
% hold on
G2STLL = plot((G2STLTP*2+1),2:6,...
	'MarkerFaceColor',G2Bcolor,'MarkerEdgeColor',G2Bcolor,'MarkerSize',10,'Marker','d');
hold on
% STLEG=legend([G1STBL(1),G1STLL(1),G2STBL(1),G2STLL(1)],...
% 	'\color[rgb]{.7 .3 .5}G1 BASAL','\color[rgb]{.7 .3 .5}G1 LTP',...
% 	'\color[rgb]{0 .5 .5}G2 BASAL','\color[rgb]{0 .5 .5}G2 LTP');
% set(STLEG,'Location','SouthWest','Box', 'off','Color', 'none');
% SAPMAPTITLE=title(char('Gaussian: ', num2str(offsap)));
% set(SAPMAPTITLE, 'FontSize', 14);
%======================================%


end








for mm = 1:6
xdataw(:,mm) = Qr(mm)*sin(t);
ydataw(:,mm) = Qr(mm)*cos(t);
end
for mm = 1:5
xdatw(:,mm) = [xdataw(:,mm+1); flipud(xdataw(:,mm))];
ydatw(:,mm) = [ydataw(:,mm+1); flipud(ydataw(:,mm))];
end
patch(xdatw,ydatw,w);
hold on



clear all
clear xdata ydata cdata fdata vdata t xdat ydat cdat

t = linspace(0,2*pi,10);
for mm = 1:6
xdata(:,mm) = (mm-1)*sin(t);
ydata(:,mm) = (mm-1)*cos(t);
end

for mm = 1:5
xdat(:,mm) = [xdata(:,mm+1); flipud(xdata(:,mm))];
ydat(:,mm) = [ydata(:,mm+1); flipud(ydata(:,mm))];
end

hsz = numel(t);
sz = size(xdat);

cdat = repmat((1:sz(2))',1,sz(1))';
cdat(1:hsz,:) = cdat(1:hsz,:)+1;

figure(43)
p = patch(xdat,ydat,'w');
set(p,'FaceColor','interp',...
'CDataMapping','scaled',...
'EdgeColor','interp',...
'LineWidth',5)







clear all
clear xdata ydata cdata fdata vdata t xdat ydat cdat

t = linspace(0,2*pi,10);
t(end+1) = 0;

for mm = 1:6
xdata(:,mm) = (mm-1)*sin(t) + 10;
ydata(:,mm) = (mm-1)*cos(t) + 10;
end

for mm = 2
xdat(:,1) = [xdata(:,mm+1); flipud(xdata(:,mm))];
ydat(:,1) = [ydata(:,mm+1); flipud(ydata(:,mm))];
end

hsz = numel(t);
sz = size(xdat);

cdat = repmat((1:sz(2))',1,sz(1))';
cdat(1:hsz,:) = cdat(1:hsz,:)+1;

figure(43)
p = patch(xdat,ydat,'w');

figure(43)
p = patch(xdat,ydat,cdat);
set(p,'FaceColor','interp',...
'FaceVertexCData',cdat,...
'EdgeColor','interp',...
'LineWidth',5)



verts = [2 4; ...
        2 8; ...
        8 4; ...
        8 0; ...
        0 4; ...
        2 6; ...
        2 2; ...
        4 2; ...
        4 0; ...
        5 2; ...
        5 0 ];
faces = [1  2  3; ...
         1  3  4; ...
         5  6  1; ...
         7  8  9; ...
         11 10 4];
p = patch('Faces',faces,'Vertices',verts,'FaceColor','b');

clear cdata 
cdata = [0 0 1; 
         0 1 0;
         0 1 1;
         1 0 0;
         1 0 1;
         1 1 0;
         0 0 0;
         0.2 0.2 0.2;
         0.4 0.4 0.4;
         0.6 0.6 0.6;
         0.8 0.8 0.8];
set(p,'FaceColor','interp',...
'FaceVertexCData',cdata,...
'EdgeColor','interp',...
'LineWidth',5)








xverts = [0 0 2]';
yverts = [0 2 0]';
verts = [xverts yverts];
faces = [1  2  3];
cdata = [0 0 1; 0 1 0; 0 1 1];
p = patch('Faces',faces,'Vertices',verts,'FaceColor','b');
set(p,'FaceColor','interp',...
'FaceVertexCData',cdata,...
'EdgeColor','interp',...
'LineWidth',5)




clear all
clear xdata ydata cdata fdata vdata t xdat ydat cdat

t = linspace(0,2*pi,8);
t(end+1) = 0;

for mm = 1:3
xdata(:,mm) = (mm-1)*sin(t) + 10;
ydata(:,mm) = (mm-1)*cos(t) + 10;
end

for mm = 2
xdat(:,1) = [xdata(:,mm+1); flipud(xdata(:,mm))];
ydat(:,1) = [ydata(:,mm+1); flipud(ydata(:,mm))];
end

xverts = xdat;
yverts = ydat;
verts = [xverts yverts];
faces = [1:numel(xverts)];
Adatc = repmat([0 0 1], size(xdata,1), 1);
Bdatc = repmat([1 0 1], size(xdata,1), 1);
cdata = [Adatc; Bdatc];

p = patch('Faces',faces,'Vertices',verts,'FaceColor','b');
set(p,'FaceColor','interp',...
'FaceVertexCData',cdata,...
'EdgeColor','interp',...
'LineWidth',5)















clear all
clear xdata ydata cdata fdata vdata t xdat ydat cdat

t = linspace(0,2*pi,8);
% t(end+1) = 0;

for mm = 1:3
xdata(:,mm) = (mm-1)*sin(t);
ydata(:,mm) = (mm-1)*cos(t);
end
xdata(xdata<0.000001 & xdata>-0.000001)=0;
ydata(ydata<0.000001 & ydata>-0.000001)=0;
for mm = 2
xdat(:,1) = [xdata(:,mm+1); flipud(xdata(:,mm))];
ydat(:,1) = [ydata(:,mm+1); flipud(ydata(:,mm))];
end

xverts = xdat;
yverts = ydat;
verts = [xverts yverts];

faces = [1:numel(xverts)];
Adatc = repmat([0 0 1], size(xdata,1), 1);
Bdatc = repmat([1 0 1], size(xdata,1), 1);
cdata = [Adatc; Bdatc];

p = patch('Faces',faces,'Vertices',verts,'FaceColor','b');
set(p,'FaceColor','interp',...
'FaceVertexCData',cdata,...
'EdgeColor','interp',...
'LineWidth',5,...
'EdgeLighting','flat',...
'FaceLighting','flat',...
'BackFaceLighting','lit')
hold on
material shiny
axis vis3d off
%light('Position',[.5 0 1]);






%--------------------------
[x,y,z] = sphere;
figure(1)
axis square
daspect([1 1 1])
hold on
surf(x.*.8,y.*.8,z.*.1+.5)
shading interp
hold on
surf(x.*.5,y.*.5,z.*.2+1.1)
shading interp
hold on
surf(x.*.01,y.*.01,z.*.01+2)
shading interp
hold on
%--------------------------







clear all
clear xdata ydata cdata fdata vdata t xdat ydat cdat

nverts = 50;
t = linspace(0,2*pi,nverts);

for mm = 1:3
xdata(:,mm) = (mm-1)*sin(t);
ydata(:,mm) = (mm-1)*cos(t);
end
xdata(xdata<0.000001 & xdata>-0.000001)=0;
ydata(ydata<0.000001 & ydata>-0.000001)=0;

for mm = 2
xdat(:,1) = [xdata(:,mm+1); (xdata(:,mm))];
ydat(:,1) = [ydata(:,mm+1); (ydata(:,mm))];
end

xverts = xdat; yverts = ydat;
verts = [xverts yverts];

fvals = [1 2 (nverts+2) (nverts+1)];
fvals = repmat(fvals,nverts-1,1);
for mm = 1:(nverts-2); fvals(mm+1,:) = fvals(mm+1,:)+mm; end
faces = fvals;

%faces = [1:numel(xverts)];
Adatc = repmat([0 0 1], size(xdata,1), 1);
Bdatc = repmat([1 0 1], size(xdata,1), 1);
cdata = [Adatc; Bdatc];

p = patch('Faces',faces,'Vertices',verts,'FaceColor','b');
set(p,'FaceColor','interp','FaceVertexCData',cdata,...
'EdgeColor','interp','LineWidth',5,...
'EdgeLighting','flat','FaceLighting','flat',...
'BackFaceLighting','lit')
hold on; material shiny; axis vis3d off;











%{
clear all
clear xdata ydata cdata fdata vdata t xdat ydat cdat

nverts = 50;
nvert = nverts-1;
t = linspace(0,2*pi,nverts);

for mm = 1:5
xdata(:,mm) = (2*mm-1)*sin(t);
ydata(:,mm) = (2*mm-1)*cos(t);
end
xdata(xdata<0.000001 & xdata>-0.000001)=0;
ydata(ydata<0.000001 & ydata>-0.000001)=0;

mnms = 4;
for mm = 1:mnms
xdat(:,mm) = [xdata(:,mm+1); (xdata(:,mm))];
ydat(:,mm) = [ydata(:,mm+1); (ydata(:,mm))];
end

xverts = xdat(:); yverts = ydat(:);
verts = [xverts yverts];

fvals = [1 2 (nverts+2) (nverts+1)];
fvals = repmat(fvals,nverts,1);
for mm = 1:(nverts-1); fvals(mm+1,:) = fvals(mm+1,:)+mm; end
fvs = repmat(fvals,mnms*2,1);

for mm = 1:mnms*2-1
	fvs(mm*nvert:mm*nvert+nvert,:) = fvs(mm*nvert:mm*nvert+nvert,:) + mm*nverts;
end

faces = fvs';

%faces = [1:numel(xverts)];
Adatc = repmat([0 0 1], size(xdata,1), 1);
Bdatc = repmat([1 0 1], size(xdata,1), 1);
cdata = [Adatc; Bdatc];

p = patch('Faces',faces,'Vertices',verts,'FaceColor','b');
set(p,'FaceColor','interp','FaceVertexCData',cdata,...
'EdgeColor','interp','LineWidth',5,...
'EdgeLighting','flat','FaceLighting','flat',...
'BackFaceLighting','lit')
hold on; material shiny; axis vis3d off;
%}

clear all

% Rank 1 color
Rank1c = [1 0 0];
Rank2c = [1 1 0];
Rank3c = [0 1 0];
Rank4c = [0 1 1];
Rank5c = [0 0 1];
CRank = {Rank1c, Rank2c, Rank3c, Rank4c, Rank5c};
Qr = [0   22.3607   31.6228   38.7298   44.7214   50.0000];
QsRank = [1 2 3 4 5];
QcRank = [1 2 3 4 5];

for Nx = 1:4
	QsR = QsRank(Nx);
	QcR = QcRank(Nx);
	clear xverts yverts verts faces Adatc Bdatc fvals xdata ydata cdata t xdat ydat cdat
	nverts = 50;
	t = linspace(0,2*pi,nverts);	
	xdata(:,1) = Qr(QsR)*sin(t);
	ydata(:,1) = Qr(QsR)*cos(t);
	xdata(:,2) = Qr(QsR+1)*sin(t);
	ydata(:,2) = Qr(QsR+1)*cos(t);
	
	xdata(xdata<0.000001 & xdata>-0.000001)=0;
	ydata(ydata<0.000001 & ydata>-0.000001)=0;
	xdat(:,1) = [xdata(:,2); (xdata(:,1))];
	ydat(:,1) = [ydata(:,2); (ydata(:,1))];
	xverts = xdat; yverts = ydat;
	verts = [xverts yverts];
	fvals = [1 2 (nverts+2) (nverts+1)];
	fvals = repmat(fvals,nverts-1,1);
	for mm = 1:(nverts-2); fvals(mm+1,:) = fvals(mm+1,:)+mm; end
	faces = fvals;
	Adatc = repmat(CRank{QcR+1}, size(xdata,1), 1);
	Bdatc = repmat(CRank{QcR}, size(xdata,1), 1);
	cdata = [Adatc; Bdatc];

	Qfaces(Nx) = {faces};
	Qverts(Nx) = {verts};
	Qcdata(Nx) = {cdata};
end


for Ntx = 1:4
figure(700)
p = patch('Faces',Qfaces{Ntx},'Vertices',Qverts{Ntx},'FaceColor','b');
set(p,'FaceColor','interp','FaceVertexCData',Qcdata{Ntx},...
'EdgeColor','interp','LineWidth',1,...
'EdgeLighting','flat','FaceLighting','flat',...
'BackFaceLighting','lit')
hold on; material shiny; axis vis3d off;
end





















%%
clear all

[clrmap1 CRank] = ClrMap1;
Qr = [0   22.3607   31.6228   38.7298   44.7214   50.0000];
QsRank = [1 2 3 4 5];

QcRankA = [4 1 3 2 5];
QcRankB =  circshift(QcRankA,[0 -1]);

for Nx = 1:4
	QsR = QsRank(Nx);
	QcA = QcRankA(Nx);
	QcB = QcRankB(Nx);
	clear xverts yverts verts faces Adatc Bdatc fvals xdata ydata cdata t xdat ydat cdat
	nverts = 50;
	t = linspace(0,2*pi,nverts);	
	xdata(:,1) = Qr(QsR)*sin(t);
	ydata(:,1) = Qr(QsR)*cos(t);
	xdata(:,2) = Qr(QsR+1)*sin(t);
	ydata(:,2) = Qr(QsR+1)*cos(t);
	
	xdata(xdata<0.000001 & xdata>-0.000001)=0;
	ydata(ydata<0.000001 & ydata>-0.000001)=0;
	xdat(:,1) = [xdata(:,2); (xdata(:,1))];
	ydat(:,1) = [ydata(:,2); (ydata(:,1))];
	xverts = xdat; yverts = ydat;
	verts = [xverts yverts];
	fvals = [1 2 (nverts+2) (nverts+1)];
	fvals = repmat(fvals,nverts-1,1);
	for mm = 1:(nverts-2); fvals(mm+1,:) = fvals(mm+1,:)+mm; end
	faces = fvals;
	Adatc = repmat(CRank{QcB}, size(xdata,1), 1);
	Bdatc = repmat(CRank{QcA}, size(xdata,1), 1);
	cdata = [Adatc; Bdatc];

	Qfaces(Nx) = {faces};
	Qverts(Nx) = {verts};
	Qcdata(Nx) = {cdata};
end


figure1 = figure('Colormap',clrmap1);
axes1 = axes('Visible','off');
hold(axes1,'all');

for Ntx = 1:4
figure(figure1)
p = patch('Faces',Qfaces{Ntx},'Vertices',Qverts{Ntx},'FaceColor','b');
set(p,'FaceColor','interp','FaceVertexCData',Qcdata{Ntx},...
'EdgeColor','interp','LineWidth',1,...
'EdgeLighting','flat','FaceLighting','flat',...
'BackFaceLighting','lit')
hold on; material shiny; axis vis3d off;
end

colorbar('peer',axes1,'YTickLabel',{'Low','-','-','-','-','Mid','-','-','-','-','Hi'},'LineWidth',1);


%%





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					SPINE RADIAL DISTANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{S1G1CvL, S2G1CvL, S1G2CvL, S2G2CvL,G1P1r, G1P2r, G2P1r, G2P2r};

%--------------------------------------
% Make concentric circle annulus
%--------------------------------------
S1r = round(S1rad);
Qhr = S1r/2;
Qx = linspace(0,S1r,6);
QxL = Qx - Qhr;
yL = sqrt(Qhr.^2 - QxL.^2);
Qr = sqrt(Qx.^2 + yL.^2);	% Concentric circles

%--------------------------------------
% Tally GluRs in each spine annulus
%--------------------------------------
for qq = 1:4
Q0=0;Q1=0;Q2=0;Q3=0;Q4=0;Q5=0;Qnan=0;
for mm = 1:numel(SGr)
	QD1 = (SGr{mm}{qq}(SGr{mm}{qq+4}(:)));
	for Nin=1:numel(QD1)
		if	   QD1(Nin) <= Qr(1)
			Q0 = Q0+1;
		elseif QD1(Nin) <= Qr(2)
			Q1 = Q1+1;
		elseif QD1(Nin) <= Qr(3)
			Q2 = Q2+1;
		elseif QD1(Nin) <= Qr(4)
			Q3 = Q3+1;
		elseif QD1(Nin) <= Qr(5)
			Q4 = Q4+1;
		elseif QD1(Nin) <= Qr(6)
			Q5 = Q5+1;
		elseif QD1(Nin) > Qr(7)
			Qnan = Qnan+1;
		end
	end
end
Qqs(qq,:) = [Q1 Q2 Q3 Q4 Q5] ./(stepN/DATARATE);
end

%--------------------------------------
% Parse out G1 G2 in SPY1 SPY2
%--------------------------------------
% Qqs =   
%  anullus 1  2  3  4  5
% G1-SPY1: N  N  N  N  N
% G1-SPY2: N  N  N  N  N
% G2-SPY1: N  N  N  N  N 
% G2-SPY2: N  N  N  N  N

% Raw
G12SPY12 = Qqs;
G12SPY1 = [Qqs(1,:) ; Qqs(3,:)];
G12SPY2 = [Qqs(2,:) ; Qqs(4,:)];
G1SPY12 = [Qqs(1,:) ; Qqs(2,:)];
G2SPY12 = [Qqs(3,:) ; Qqs(4,:)];

% Means
G12SPY12m = mean(G12SPY12);
G12SPY1m = mean(G12SPY1);
G12SPY2m = mean(G12SPY2);
G1SPY12m = mean(G1SPY12);
G2SPY12m = mean(G2SPY12);

% Sums
G12SPY12s = sum(G12SPY12);
G12SPY1s = sum(G12SPY1);
G12SPY2s = sum(G12SPY2);
G1SPY12s = sum(G1SPY12);
G2SPY12s = sum(G2SPY12);

Qs = G12SPY1s;


%--------------------------------------
% Interpolate Annulus Patch Densities
%--------------------------------------
[clrmap1 CRank] = ClrMap2;
[QsortVal, QsortInd] = sort(Qs, 2);
QcRankA = [fliplr(QsortInd) 6];
QcRankB =  circshift(QcRankA,[0 -1]);

for Nx = 1:5
clear xverts yverts verts faces Adatc Bdatc fvals xdata ydata cdata t xdat ydat cdat
	QcA = QcRankA(Nx);
	QcB = QcRankB(Nx);
	nverts = 50;
	t = linspace(0,2*pi,nverts);	
	%xdata(:,1) = Qr(Nx)*sin(t);
	%ydata(:,1) = Qr(Nx)*cos(t);
	%xdata(:,2) = Qr(Nx+1)*sin(t);
	%ydata(:,2) = Qr(Nx+1)*cos(t);
	
	xdat = [(Qr(Nx+1)*sin(t)) (Qr(Nx)*sin(t))]';
	ydat = [(Qr(Nx+1)*cos(t)) (Qr(Nx)*cos(t))]';
	
	
	%xdat(:,1) = [xdata(:,2); (xdata(:,1))];
	%ydat(:,1) = [ydata(:,2); (ydata(:,1))];
	xverts = xdat; yverts = ydat;
	verts = [xverts yverts];
	fvals = [1 2 (nverts+2) (nverts+1)];
	fvals = repmat(fvals,nverts-1,1);
	for mm = 1:(nverts-2); fvals(mm+1,:) = fvals(mm+1,:)+mm; end
	faces = fvals;
	Adatc = repmat(CRank{QcB}, nverts, 1);
	Bdatc = repmat(CRank{QcA}, nverts, 1);
	cdata = [Adatc; Bdatc];

	Qfaces(Nx) = {faces};
	Qverts(Nx) = {verts};
	Qcdata(Nx) = {cdata};
end


figure1 = figure('Colormap',clrmap1);
axes1 = axes('Visible','off');
hold(axes1,'all');

for Ntx = 1:5
figure(figure1)
p = patch('Faces',Qfaces{Ntx},'Vertices',Qverts{Ntx},'FaceColor','b');
set(p,'FaceColor','interp','FaceVertexCData',Qcdata{Ntx},...
'EdgeColor','interp','LineWidth',.1,...
'EdgeLighting','flat','FaceLighting','flat',...
'BackFaceLighting','lit')
hold on; material shiny; axis vis3d off;
end

colorbar('peer',axes1,'YTickLabel',{num2str(QsortVal(1), '% 10.2f'),...
'-','-','-','-',num2str(QsortVal(3), '% 10.2f'),...
'-','-','-','-',num2str(QsortVal(end), '% 10.2f')},...
'LineWidth',1);
hold on;			  


%-------------------------------------
CQRall = [];
figure(figure1)
for mm = 200:2:numel(SGr)
QPh2 = scatter(SGr{mm}{9}(1,:),SGr{mm}{9}(2,:));
set(QPh2,'Marker','o','SizeData',20,'LineWidth',.5,...
	'MarkerFaceColor',[.2 1 1],'MarkerEdgeColor',[0 .4 .4])
%pause(.2)
CQRall = cat(2,CQRall,[SGr{mm}{9}(1,:);SGr{mm}{9}(2,:)]);
end
hold on

[PSDY,PSDX] = find(ACTINp1);
S1Ph2 = scatter(PSDX-50,PSDY-50);
set(S1Ph2,'Marker','o','SizeData',60,'LineWidth',.5,...
	'MarkerFaceColor',[.8 .2 .2],'MarkerEdgeColor',[.4 0 0])
hold on


%--------
for mm = 2:5
rectangle('Position',[-Qr(mm),-Qr(mm),Qr(mm)*2,Qr(mm)*2],'Curvature',[1,1])
hold on;
end
%--------



%===========================================================
% Alternative Anullus Plots:
%{
TQs = [PSDY PSDX]'-50;
TQRz = sqrt(sum(TQs)).*0+1;
CQRz = sqrt(sum(CQRall.^2)).*0+2;
figure2 = figure('XVisual','');
axes2 = axes('Parent',figure2,'PlotBoxAspectRatio',[0.555555555555556 0.555555555555556 1],...
	'OuterPosition',[0 0 1 1],...
	'DataAspectRatio',[10 10 1],...
	'CameraViewAngle',4.36580191767115);
view(axes2,[50.5 40]);
grid(axes2,'on');
hold(axes2,'all');
PhS3a = scatter3(TQs(1,:),TQs(2,:),TQRz,50, [.8 .2 .2],'filled');
set(PhS3a,'Marker','s')
hold on
PhS3b = scatter3(CQRall(1,:),CQRall(2,:),CQRz,10, [.2 .7 .0],'filled');
set(PhS3b,'Marker','o')
hold off
%}
%{

QallTG = cat(2,CQRall,TQs);
QallTGz = cat(2,CQRz,TQRz);

Qallz = cat(1,QallTG,QallTGz);
spinesurf2(Qallz(1,:)', Qallz(2,:)', Qallz(3,:)')
sftool(Qallz(1,:)', Qallz(2,:)', Qallz(3,:)')

%--------------------------------------
CQRz = sqrt(sum(CQRall.^2));
sja=0;
sjb=max(CQRz);
sjc=2;
sjd=1;
linsc = @(jx) (sjc*(1-(jx-sja)/(sjb-sja)) + sjd*((jx-sja)/(sjb-sja)));
for mm = 1:numel(CQRz)
	CQRz(mm) = linsc(CQRz(mm));
end
Qallz = cat(1,CQRall,CQRz);
spinesurf2(Qallz(1,:)', Qallz(2,:)', Qallz(3,:)')
sftool(Qallz(1,:)', Qallz(2,:)', Qallz(3,:)')
%}
%{

%-------------------------------------
%Generate values from the uniform distribution on the interval [a, b]:
% r = a + (b-a).*rand(100,1);
clear pxr pyr Qpr
for mm = 1:(numel(Qr)-1)
clear pxr pyr
	Qr1 = Qr(mm+1); Qr2 = Qr(mm); 
	for nn = 1:Qs(mm)*8e3 % <------- multiplier
		theta = 2*pi*rand;
		r = sqrt(rand*(Qr1^2-Qr2^2)+Qr2^2);
		pxr(nn) = r*cos(theta); 
		pyr(nn) = r*sin(theta);
	end
	Qpr{mm} = [pxr; pyr];
end

% figure
% scatter(Qpr{1}(1,:), Qpr{1}(2,:),3,[1 0 0],'filled')
% hold on
% scatter(Qpr{2}(1,:), Qpr{2}(2,:),3,[.8 .9 .1],'filled')
% hold on
% scatter(Qpr{3}(1,:), Qpr{3}(2,:),3,[.6 .5 0],'filled')
% hold on
% scatter(Qpr{4}(1,:), Qpr{4}(2,:),3,[.4 .9 .9],'filled')
% hold on
% scatter(Qpr{5}(1,:), Qpr{5}(2,:),3,[.2 0 1],'filled')
% hold on

%-------------------------------------
QA1 = round([Qpr{1}(1,:); Qpr{1}(2,:)]+51);
QA2 = round([Qpr{2}(1,:); Qpr{2}(2,:)]+51);
QA3 = round([Qpr{3}(1,:); Qpr{3}(2,:)]+51);
QA4 = round([Qpr{4}(1,:); Qpr{4}(2,:)]+51);
QA5 = round([Qpr{5}(1,:); Qpr{5}(2,:)]+51);
G1xy = zeros(100,100);
for xy = 1:numel(QA1(1,:));G1xy(QA1(1,xy),QA1(2,xy)) = Qs(1);end
for xy = 1:numel(QA2(1,:));G1xy(QA2(1,xy),QA2(2,xy)) = Qs(2);end
for xy = 1:numel(QA3(1,:));G1xy(QA3(1,xy),QA3(2,xy)) = Qs(3);end
for xy = 1:numel(QA4(1,:));G1xy(QA4(1,xy),QA4(2,xy)) = Qs(4);end
for xy = 1:numel(QA5(1,:));G1xy(QA5(1,xy),QA5(2,xy)) = Qs(5);end

G1xy(101:end,:) = [];
G1xy(:,101:end) = [];

GN=[4, .18, 9, .1];
ccMask = CircMask(GN);
ccMask = ccMask -.2;
G1xy(G1xy<min(Qs))=min(Qs)*.3;
G1Cv = convn(G1xy,ccMask,'same');

figure
imagesc(G1xy)
colormap(jet)
colorbar

figure
imagesc(G1Cv)
colormap(jet)
colorbar

%-------------------------------------
pxr=0; pyr=0; Qpr={0};
for mm = 1:(numel(Qr)-1)
clear pxr pyr
	Qr1 = Qr(mm+1); Qr2 = Qr(mm); 
	for nn = 1:round(20^(sqrt(Qs(mm))))
		theta = 2*pi*rand;
		r = sqrt(rand*(Qr1^2-Qr2^2)+Qr2^2);
		pxr(nn) = r*cos(theta); 
		pyr(nn) = r*sin(theta);
	end
	Qpr{mm} = [pxr; pyr];
end

Qall = [Qpr{1}(1,:) Qpr{2}(1,:) Qpr{3}(1,:) Qpr{4}(1,:) Qpr{5}(1,:);
		Qpr{1}(2,:) Qpr{2}(2,:) Qpr{3}(2,:) Qpr{4}(2,:) Qpr{5}(2,:)];

Qzz =	[(.8 + (.3).*rand(1,numel(Qpr{1}(1,:)))) ...
		(.6 + (.3).*rand(1,numel(Qpr{2}(1,:)))) ...
		(.4 + (.3).*rand(1,numel(Qpr{3}(1,:)))) ...
		(.2 + (.3).*rand(1,numel(Qpr{4}(1,:)))) ...
		(0 + (.3).*rand(1,numel(Qpr{5}(1,:))))];
Qallz = cat(1,Qall,Qzz);
Qallz(isinf(Qallz)) = 0;
Qallz(isnan(Qallz)) = 0;

spinesurf2(Qallz(1,:)', Qallz(2,:)', Qallz(3,:)')
%sftool(Qallz(1,:)', Qallz(2,:)', Qallz(3,:)')
%-------------------------------------
%}
%===========================================================



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







