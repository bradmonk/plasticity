clc, close all, clear all

% QRpoints: number of points to test in range specified below
QRpoints = 6;
NPtested = 3;
assignin('base','QRpoints',QRpoints);



% RATE CONSTANTS
Lon = 2.0;
Loff = 2.0;
%-- INDY (low=fast)
%Qon = [10 15];
Qon = [1 20];
Qoff = [.1 4];
%-- DEP (high=fast)
Ron = [.1 10];
Roff = [.1 3];

%-- Ratio of Ron:Roff
Rrat = [1.5 8];
Rbase = 1.5;



% profile on
%==========================================================%
% Sfull{25,50} = [];
% Lon = 1.8;  % On Energy (lower = more on events)
% Qon = 32;   % On Neighbor-Independant Rate (new growth) (lower = more on)
% Ron = 3;    % On Neighbor-Dependant Rate (cluster fill-in) (higher = more on)
%  
% Loff = 1.8; % Off Energy (higher = more off events)
% Qoff = 1.2; % Off Neighbor-Inependant Rate (uniform off)  (lower = more off)
% Roff = 3;   % Off Neighbor-Dependant Rate (edge off) (higher = more off)
%==========================================================%

PQon  = linspace(Qon(1),Qon(2),QRpoints);
PQoff = linspace(Qoff(1),Qoff(2),QRpoints);

PRon  = linspace(Ron(1),Ron(2),QRpoints);
PRoff = linspace(Roff(1),Roff(2),QRpoints);

PRrat = linspace(Rrat(1),Rrat(2),QRpoints);
PRRon = Rbase .* PRrat;
PRRoff = Rbase .* ones(1,numel(PRRon));



%==========================================================%
nn = 0;
for aa = 1:QRpoints; 
for bb = 1:QRpoints; 
for cc = 1:QRpoints;
disp(aa);disp(bb);disp(cc);
nn = nn+1;
%-------------------------

	Qon = PQon(aa);
	Qoff = PQoff(bb);
	
	Ron = PRRon(cc);
	Roff = PRRoff(cc);
	PPRrat = PRrat(cc);

	Spara = [Lon Loff Qon Qoff Ron Roff PPRrat];
	[Scell HKcell OKGO] = RobustSAP(Lon,Qon,Ron,Loff,Qoff,Roff);

	AOK(nn) = OKGO;
	LRPAR(nn,:) = Spara;
	Sfull{nn} = {Spara Scell OKGO};
	HKfull{nn} = {Spara HKcell OKGO};

%-------------------------
end;end;end;
%==========================================================%
% profile viewer
%==========================================================%
save('RobustData.mat', 'Sfull','HKfull')

%%
%==========================================================%
% clear
% load RobustData.mat
%==========================================================%

%-- [Qon Qoff Ron Roff PPRrat] --%
OKPARAMS = LRPAR((AOK>0),3:7);
NOKPARAMS = LRPAR((AOK<1),3:7);

%%
%======================================================================%
%							FIGURE SETUP
%----------------------------------------------------------------------%
%----------------------------------------------------------------------%
fig10 = figure(10);
set(10,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/2.5  scnsize(4)/2];
set(fig10,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(10); hold on
X = [OKPARAMS(:,1) OKPARAMS(:,2)];
V = OKPARAMS(:,5);
% plot3(X(:,1),X(:,2),(OKPARAMS(:,1).*0), '.b')
plot3(X(:,1),X(:,2),(OKPARAMS(:,5)), '.b')
grid
stem3(X(:,1),X(:,2),V,':ob','fill')

X = [NOKPARAMS(:,1) NOKPARAMS(:,2)];
V = NOKPARAMS(:,5);
% plot3(X(:,1),X(:,2),(NOKPARAMS(:,1).*0), '.r')
%plot3(X(:,1),X(:,2),(NOKPARAMS(:,5)), '.r')
stem3(X(:,1),X(:,2),V,':or','fill')
grid on
view(-10, 11);
hXLabel = xlabel('Qon');
hYLabel = ylabel('Qoff');
hZLabel = zlabel(num2str(['Ron:Roff (@ ' num2str(Rbase) ')']) );
set(gca,'FontName','Helvetica');
set([hXLabel, hYLabel, hZLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel, hZLabel],'FontSize',12);
%=============================%
% figure(10); hold on
% plot3(OKPARAMS(:,1),OKPARAMS(:,2),OKPARAMS(:,5), '.b')
% grid; hold on
% plot3(NOKPARAMS(:,1),NOKPARAMS(:,2),NOKPARAMS(:,5), '.r')
% view(-36, 60);
%=============================%

%=============================%
x = OKPARAMS(:,1);
y = OKPARAMS(:,2);
z = OKPARAMS(:,5);
xyz = [x y z];

% plot3(x,y,z,'.-')
tri = delaunay(x,y);
% plot(x,y,'.')
[r,c] = size(tri);
disp(r)
%=============================%
tri3 = delaunay(x,y,z);
figure(fig10)
h = trisurf(tri3, x, y, z);
axis vis3d
l = light('Position',[-50 -15 29]);
lighting phong
shading interp
hXLabel = xlabel('Qon');
hYLabel = ylabel('Qoff');
hZLabel = zlabel(num2str(['Ron:Roff (@ ' num2str(Rbase) ')']) );
%=============================%
% lighting flat
% lighting gouraud
% lighting phong
% lighting none
% shading flat 
% shading faceted 
% shading interp 
% Create light
%======================================================================%
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,'SPACEPLOT','png');
% saveas(gcf, ['outputfigs/SPACEPLOT.png']);
%======================================================================%

%%

keyboard

%=============================%
%{

Xok = OKPARAMS(:,1);
Yok = OKPARAMS(:,2);
Zok = OKPARAMS(:,5);
XYZok = [Xok Yok Zok];

%-- [Qon Qoff Ron Roff PPRrat] --%
PQon  = linspace(Qon(1),Qon(2),QRpoints);
PQoff = linspace(Qoff(1),Qoff(2),QRpoints);
PRon  = linspace(Ron(1),Ron(2),QRpoints);
PRoff = linspace(Roff(1),Roff(2),QRpoints);
PRrat = linspace(Rrat(1),Rrat(2),QRpoints);
PRRon = Rbase .* PRrat;
PRRoff = Rbase .* ones(1,numel(PRRon));
%--
QKon = linspace(min(Xok),max(Xok),numel(unique(Xok)));
QKoff = linspace(min(Yok),max(Yok),numel(unique(Yok)));
RKrat = linspace(min(Zok),max(Zok),numel(unique(Zok)));
%--




PQon  = linspace(Qon(1),Qon(2),QRpoints);
PQoff = linspace(Qoff(1),Qoff(2),QRpoints);
PRrat = linspace(Rrat(1),Rrat(2),QRpoints);





REPx = repmat(PQon,QRpoints,1);
REPy = repmat(PQoff',1,QRpoints);
Rz = zeros(size(REPx));

rr = 2
cc = 4

for rr = 1:numel(REPx(:,1))
for cc = 1:numel(REPx(1,:))
	
	Xnn = REPx(rr,cc); Ynn = REPy(rr,cc);
	
	Xfnd = XYZok(:,1) == Xnn;
	Yfnd = XYZok(:,2) == Ynn;
	XYfnd = [Xfnd Yfnd];
	XYZfnd = [Xfnd Yfnd XYZok(:,3)];
	
	


end
end





[XXok,XYok] = meshgrid(min(Xok):1:max(Xok));
XRok = sqrt(XXok.^2 + XYok.^2) + eps;
XZok = sin(XRok)./XRok;
mesh(XXok,XYok,XZok)

surf(Z) 
[PQonX,PQonY] = meshgrid(PQon);
[PQoffX,PQoffY] = meshgrid(PQoff);
[PRonX,PRonY] = meshgrid(PRon);
[PRoffX,PRoffY] = meshgrid(PRoff);


TOEPMx = toeplitz(PQoff,PRon);
mesh(PQonX,PQonY,TOEPMx)

% SLICE 3D VOLUME
[x,y,z] = meshgrid(PQon,PQoff,PRon);
FBFBT = ones(4,4,4);
FBFBT(:,:,1) = (PRoffX.*FBFBT(:,:,1));
FBFBT(:,:,2) = (PRoffX.*FBFBT(:,:,2));
FBFBT(:,:,3) = (PRoffX.*FBFBT(:,:,3));
FBFBT(:,:,4) = (PRoffX.*FBFBT(:,:,4));
v = FBFBT;
xslice = [PQon(2),PQon(3)]; yslice = PQoff(2); zslice = [PRon(2),PRon(3)];
slice(x,y,z,v,xslice,yslice,zslice)
colormap hsv



%=============================%
[x,y,z] = meshgrid(-2:.2:2,-2:.25:2,-2:.16:2);
v = x.*exp(-x.^2-y.^2-z.^2);
xslice = [-1.2,.8,2]; yslice = 2; zslice = [-2,0];
slice(x,y,z,v,xslice,yslice,zslice)
colormap hsv
for i = -2:.5:2
 hsp = surf(linspace(-2,2,20),linspace(-2,2,20),zeros(20)+i);
 rotate(hsp,[1,-1,1],30)
 xd = get(hsp,'XData');
 yd = get(hsp,'YData');
 zd = get(hsp,'ZData');
 delete(hsp)
 slice(x,y,z,v,[-2,2],2,-2) % Draw some volume boundaries
 hold on
 slice(x,y,z,v,xd,yd,zd)
 hold off
 axis tight 
 view(-5,10) 
 drawnow
 pause(.5)
end
%-----------------------------%
[xsp,ysp,zsp] = sphere;
slice(x,y,z,v,[-2,2],2,-2)  % Draw some volume boundaries
for i = -3:.2:3
 hsp = surface(xsp+i,ysp,zsp);
 rotate(hsp,[1 0 0],90)
 xd = get(hsp,'XData');
 yd = get(hsp,'YData');
 zd = get(hsp,'ZData');
 delete(hsp)
 hold on
 hslicer = slice(x,y,z,v,xd,yd,zd);
 axis tight 
 xlim([-3,3])
 view(-10,35)
 drawnow
 pause(.5)
 delete(hslicer)
 hold off
end
%-----------------------------%
PQonScat = [PQon' PQoff';PQon' PRon';PQon' PRoff';...
	PQoff' PRon'; PQoff' PRoff'; PRon' PRoff'];
X = PQonScat;
V = PQonScat(:,1);
hold on
plot3(X(:,1),X(:,2),zeros(numel(PQonScat(:,1)),1), '*r')
grid
stem3(X(:,1),X(:,2),V,'^','fill')
hold off
view(322.5, 30);




%=============================%
% Cell{ID}{data}{step}
%-----------------------------%
% ID = 1:QRpoints^QRpoints
%---
% data = 1 = param vals
% data = 2 = SAP Mx
% data = 3 = Robust yes/no indicator
%---
% step = 1 = step 10
% step = 2 = step 100
% step = 3 = step 1000
% step = 4 = step 5000
%---
% cellplot(Sfull{1})
%=============================%

Sfull{1}{3}




sum(Sfull{1}{2}{1})
sum(HKfull{1}{2}{1})

imagesc(Sfull{1}{2}{4})
imagesc(HKfull{1}{2}{4})

Sf = Sfull{1}{2}{4};
Sf1 = numel(find(Sf(1,1:end)));
Sf2 = numel(find(Sf(1:end,1)));

Sfull{1}{2}{4}
HKfull{1}{2}{4}

%}
%=============================%


%%
% cellplot(Sfull{1})
%==========================================================%
doMakeMovie = 0;
if doMakeMovie
%==========================================================%
fig1 = figure(1);
set(gcf,'OuterPosition',[800,300,500,500])

for kk = 1:QRpoints^NPtested
figure(1);
subplot('Position',[.05 .55 .40 .40]),imagesc(Sfull{kk}{2}{1});
set(gca,'XTick',[],'YTick',[])
subplot('Position',[.55 .55 .40 .40]),imagesc(Sfull{kk}{2}{2});
set(gca,'XTick',[],'YTick',[])
subplot('Position',[.05 .05 .40 .40]),imagesc(Sfull{kk}{2}{3});
set(gca,'XTick',[],'YTick',[])
subplot('Position',[.55 .05 .40 .40]),imagesc(Sfull{kk}{2}{4});
set(gca,'XTick',[],'YTick',[])
colormap('bone')
text(1,-2,['Lon       Loff       Qon       Qoff       Ron       Roff'],...
	'FontSize',16,'HorizontalAlignment','center')
text(1,-.5,['' num2str(Sfull{kk}{1}) ''],...
	'FontSize',16,'HorizontalAlignment','center')
Vid(kk) = getframe(gcf,[0 0 500 500]);
% writeVideo(ClusVid,Vid);
end
%==========================================================%
% PLAY MOVIE
Vidframes = numel(Vid);

if Vidframes >= 150
fig2 = figure(2);
set(gcf,'OuterPosition',[800,300,500,500])
movie(fig2,Vid(1:150),1,5)
movie(fig2,Vid(151:end),1,5)
else
fig2 = figure(2);
set(gcf,'OuterPosition',[800,300,500,500])
movie(fig2,Vid(1:end),1,5)
end

%==========================================================%
end % if doMakeMovie
%==========================================================%










%%
%{
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
ylim([0 haxes(4)*1.2 ]);
xlim([0 (haxes(2)-3)]);
%======================================================================%


%======================================================================%
%					SAVE FIGURE 1 PLOTS
%======================================================================%
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,'STARShiP2','png');
saveas(gcf, ['outputfigs/FIGURE1.png']);
%======================================================================%

%}
















