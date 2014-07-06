
%=================================================%
%			BRAD_MONK_MATLAB_NOTES    
%=================================================%
%{

SECTIONTEMPLATE


%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%



%=================================================%
%					TRIG
%=================================================%
%{

HYPOTENUSE
side opposite the right-angle

OPPOSITE
side opposite of interest-angle

ADJACENT
side connecting interest-angle with right-angle

sin(?) = Opposite / Hypotenuse
cos(?) = Adjacent / Hypotenuse
tan(?) = Opposite / Adjacent

http://www.mathsisfun.com/flash.php?path=/geometry/images/circle-triangle.swf



EXPONENTIAL FUNCTION e = 2.7183
MatLab: exp()
The slope of (assumed line) tangent to graph (x-axis) at each point 
is equal to its y coordinate at that point.

exp(-3)	= 0.0498
exp(-2)	= 0.1353
exp(-1)	= 0.3679
exp(0)	= 1
exp(1)	= 2.7183
exp(2)	= 7.3891
exp(3)	= 20.0855

xp = -3:3
yp = exp(xp)
plot( xp , yp )

%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%




%=================================================%
%       AXIS AXIS AXIS FUCKING AXIS  
%=================================================%
%{
get(gcf)
set(gcf,'PropName','val')


XT_Mu = xx_Mu';
		YT_Mu = yy_Mu';
		ET_Mu = ee_Mu';
		subplot('Position',sbpos),...
        [hl, hp] = boundedline(XT_Mu,YT_Mu, ET_Mu,...
		'cmap',cOLOR(itemN,:),'alpha','transparency', 0.4);
		ph = hl;
		po = hp;
		hax = gca;

h = findobj(gca,'Type','line')
findobj('Color','r')
h = gco(2)

sphere(36); 
h = findobj('Type','surface');
set(h,'FaceLighting','phong',...
      'FaceColor','interp',...
      'EdgeColor',[.4 .4 .4],...
      'BackFaceLighting','lit')
hold on
patch('faces',fac,'vertices',vert,'FaceColor','y');
light('Position',[1 3 2]);
light('Position',[-3 -1 3]);
material shiny
axis vis3d off
hold off


figure(1)
subplot(5,5,[3 25]), A = bar(10*rand(1,40))
Ah = get(A,'child')
set(Ah,'facea',.3)
hold on
figure(1)
subplot(5,5,[3 25]), B=imagesc(map)
Bh = get(B)
set(Bh,'facea',.3)


rectangle('Position',[10 10 6 6], 'LineWidth',2, 'EdgeColor','b');


% just get the handles from those axes 
% and use them as first argument in the plot
figure
hax1=axes
figure
hax2=axes
plot (hax1,t(:,1),t(:,2),'-r+')
plot (hax2,t(:,1),t(:,3),'-r+')



figure(1);
subplot(5,5,[3 25]), pfield = imagesc(PSDfield);
set(gca, 'nextplot', 'add')
hold on


grid on;
%set(scat2d,'marker','.', 'EraseMode','xor')
%set(scat2d,'marker','.', 'EraseMode','xor')
axis([xlim, ylim, zlim]);
set(gca, 'Box', 'off', 'GridLineStyle','-');
shading interp;
light;
lighting phong;
title('lighting phong', 'FontName', 'Courier', 'FontSize', 14);












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig12 = figure(12);
set(12,'Units','pixels');scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/3  scnsize(4)/5  scnsize(3)/1.5  scnsize(4)/1.5];
set(fig12,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%
c1= [.9 .2 .2]; c2= [.2 .4 .6]; c3= [.4 .8 .4]; c4= [.6 .6 .6];
c11=[.9 .3 .3]; c22=[.3 .5 .7]; c33=[.5 .9 .5]; c44=[.7 .7 .7];
cOLOR = [c1; c2; c3; c4; c1; c2; c3; c4];
%===========================================================%
% FIG2 TOP LEFT: Synaptic AMPARs
%===========================================================%
sbpos = [.055 .57 .4 .38]; ptype = 5;
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
set(ph2,'LineStyle','-','Color',c2,'LineWidth',1,...
'Marker','o','MarkerSize',MS1,'MarkerEdgeColor',c2,'MarkerFaceColor',c22);
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
ylim([haxes(3)/1.2 haxes(4)*1.2]);
xlim([0 (haxes(2)-3)]);
%======================================================================%
set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,'STARShiP2','png');
saveas(gcf, ['outputfigs/STARShiP2.png']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%=================================================%
figure(3)
cdfexit = gcf;
figure(cdfexit)
subplot(5,4,[1 12]),GluR2plot = cdfplot(GluR2exT30);
axis([xlim,ylim]);
set(GluR2plot,'color',[1 .3 1])
hold on
subplot(5,4,[1 12]),cdfplot(GluR1exT30)
CDFtitle = title(['CDF of Fraction Exited    D =' num2str(D) '   GluR2 N =' int2str(Sdots)...
	'    GluR1 N =' int2str(Sdots2)]);
set(CDFtitle, 'FontSize', 16);
leg1=legend('  GluR2', 'GluR1');
set(leg1,'Location','SouthEast');
figure(cdfexit)
subplot(5,4,[13 14]),cdfplot(exittime)	
subplot(5,4,[15 16]),ecdfhist(coECDF_G2, coX_G2)
subplot(5,4,[17 18]),cdfplot(exittime2)	
subplot(5,4,[19 20]),ecdfhist(coX_G1, coX_G1)



% set(0,'DefaultFigureWindowStyle','normal')
% set(gcf,'OuterPosition',[xorigin,yorigin,width,height])
% subplot('Position',[left bottom width height]);

figure(3)
subplot(2,2,4), subplot('Position',[.55 .05 .4 .4]),...
	plot([(Ddata(:,1)) (Ddata(:,2))]);
yt = get(gca,'YTick');
set(gca,'YTickLabel', sprintf('%.4f|',yt))
set(get(gca,'YLabel'),'String','scaled diffusion rate')
xt = (get(gca,'XTick'))*10;
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time-Step')
title('PSD DIFFUSION RATE'); leg4=legend('D-PSD1', 'D-PSD2');
set(leg4,'Location','SouthWest');


% overlay alpha
close all
clear all
fPSD1 = ones(6);                        % PSD1 SIZE
fPSD2 = ones(6);                        % PSD2 SIZE
pfPSD1 = padarray(fPSD1,[8 8], 0);      % PAD PSD1
pfPSD2 = padarray(fPSD2,[8 8], 0);      % PAD PSD2
PSDfield = cat(1, pfPSD1, pfPSD2);      % CONCAT PSD FIELDS
map = PSDfield

figure(1)
subplot(5,5,[3 25]), A = bar(10*rand(1,40))
Ah = get(A,'child')
set(Ah,'facea',.3)
hold on
figure(1)
subplot(5,5,[3 25]), B=imagesc(map)
Bh = get(B)
set(Bh,'facea',.3)


rectangle('Position',[10 10 6 6], 'LineWidth',2, 'EdgeColor','b');


just get the handles from those axes 
and use them as first argument in the plot
figure
hax1=axes
figure
hax2=axes
plot (hax1,t(:,1),t(:,2),'-r+')
plot (hax2,t(:,1),t(:,3),'-r+')


figure(1);
subplot(5,5,[3 25]), pfield = imagesc(PSDfield);
set(gca, 'nextplot', 'add')
hold on


rectangle('Position',[x,y,w,h]) 
draws the rectangle from the point x,y and having a width of w and 
a height of h. Specify values in axes data units.

Note that, to display a rectangle in the specified proportions, 
you need to set the axes data aspect ratio so that one unit is of 
equal length along both the x and y axes. You can do this with the 
command 'axis equal' or daspect([1,1,1]).

%---
figure(1);
subplot(5,5,[3 25]), 
scat2d = gscatter(Caxyl(1,:),Caxyl(2,:)); view(0, 80); % <--FIG-----######
axis image;
grid on;
set(scat2d,'marker','.','markersize',[2])
%set(scat2d,'marker','.', 'EraseMode','xor')
axis([xlim, ylim, zlim]);
set(gca, 'Box', 'off', 'GridLineStyle','-');
%---
keyboard
get(gcf)



-------------
LIGHTING
surf(X, Y, Z); view(30, 30);
shading interp;
light;
lighting phong;
title('lighting phong', 'FontName', 'Courier', 'FontSize', 14);
-------------


figure(4)
subplot(2,2,3), plot([SAPdata(:,1) SAPdata(:,2)]);
yt = get(gca,'YTick');
set(gca,'YTickLabel', sprintf('%.4f|',yt))
set(get(gca,'YLabel'),'String','scaled diffusion rate')
xt = (get(gca,'XTick'))*10;
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time-Step')
title('PSD DIFFUSION RATE'); leg4=legend('D-PSD1', 'D-PSD2');
set(leg4,'Location','SouthWest');
plotinset = get(gca,'TightInset');
subplot('Position',[plotinset(1) plotinset(2) .4 .4]),...
	plot([SAPdata(:,1) SAPdata(:,2)]);
set(get(gca,'YLabel'),'String','scaled diffusion rate')
set(get(gca,'XLabel'),'String','Time-Step')


% h1 = zoom;
% h2 = pan;
% figure(2)
% subplot(3,2,6), plot([(Ddata(:,1)) (Ddata(:,2))]);
% set(gca,'XTickLabelMode','auto')
% set(gca,'XTickLabel',num2str(get(gca,'XTick').'))
% set(h1,'ActionPostCallback',@mypostcallbackX);
% set(h2,'ActionPostCallback',@mypostcallbackX);
% title('PSD DIFFUSION RATE'); leg4=legend('D-PSD1', 'D-PSD2');
% set(leg4,'Location','SouthWest');
% function mypostcallbackX(obj,evd)
% set(gca,'XTickLabelMode','auto')
% set(gca,'XTickLabel',num2str(get(gca,'XTick').'))
% end



%-- Scatter Plot from Matrix Data --%
Mx = [1:5; 1:5; 1:5; 1:5; 1:5;]
My = Mx+5
gscatter(Mx, My)
load seamount
figure
hs(1) = subplot(2,1,1);
hs(2) = subplot(2,1,2);
scatter3(hs(1),x,y,z,'MarkerFaceColor',[0 .75 .75])
scatter3(hs(2),x,y,z,'*')


%-- Axis --%
axis([xmin xmax ymin ymax])
axis([xmin xmax ymin ymax zmin zmax cmin cmax])
v = axis
axis auto
axis manual
axis tight      % sets the axis limits to the range of the data.
axis fill
axis ij
axis xy
axis equal
axis image
axis square
axis vis3d
axis normal
axis off
axis on

axis(axes_handles,...)      % applies the axis command to the specified axes. 
% For example, the following statements:
h1 = subplot(221);
h2 = subplot(222);
axis([h1 h2],'square')

axis(subplot(2), off)

axis(axes_handles,...)
[mode,visibility,direction] = axis('state')

axis manual and axis(axis) 
freezes the scaling at the current limits, 
so that if hold is on, subsequent plots use the same limits. 
This sets the XLimMode, YLimMode, and ZLimMode properties to manual.





gplotmatrix(ratings(:,1:2),ratings(:,[4 7]),group,... 
            'br','.o',[],'on','',categories(1:2,:),... 
             categories([4 7],:))


%}
%=================================================%
%=================================================%


%=================================================%
%			help codetools  
%=================================================%
%{

help codetools

Editing and publishing
    edit      - Edit or create a file
    grabcode  - Copy MATLAB code from published HTML
    mlint     - Check files for possible problems
    notebook  - Open MATLAB Notebook in Microsoft Word
    publish   - Publish file containing cells to output file
    snapnow   - Force snapshot of image for published document


%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%


%=================================================%
%			generate help file
%=================================================%
%{

function [xp yp varargout] = circus(r,varargin)
% CIRCUS	circle boundary values
%
% %% Syntax
%
%	[xp yp] = circus(r)
%	[xp yp] = circus(r,xc,yc)
%	[xp yp zp] = circus(r,xc,yc,zc)
% 
% %% Description
% 
%	[xp yp] = circus(r)  takes a radius 'r' and returns vectors xp and yp
%	corresponding to the x,y boundary points of a circle centered at 0
% 
%	[xp yp] = circus(r,xc,yc) takes a radius 'r' and the x,y coordinate-center 
%	positions 'xc' and 'yc' and returns vectors xp and yp
%	corresponding to the x,y boundary points of a circle centered at xc,yc
% 
%	[xp yp zp] = circus(r,xc,yc,zc) same as circus(r,xc,yc) and includes
%	z-axis coordinates
%
% %% Examples
% 
%		r = 5;
%		xc = 10;
%		yc = 10;
%		zc = 4;
% 
%	% 2D Circle
%
%		[xp yp] = circus(r,xc,yc);
%		figure
%		plot(xp, yp);
%		axis square;
%
%	% 3D Circle
%
%		[xp yp zp] = circus(r,xc,yc,zc);
%		figure
%		plot3(xp, yp, zp);
%		axis square;
%
% See also INPOLYGON, RECTANGLE, PLOT3
%=====================================================================%

ang=0:0.01:2*pi; 

if nargin == 1
	xc=0; yc=0;
	xp = r*cos(ang) + xc;
	yp = r*sin(ang) + yc;
elseif nargin == 2
	xc=varargin{1}; yc=0;
	xp = r*cos(ang) + xc;
	yp = r*sin(ang) + yc;
elseif nargin == 3
	xc=varargin{1}; yc=varargin{2};
	xp = r*cos(ang) + xc;
	yp = r*sin(ang) + yc;
elseif nargin == 4
	xc=varargin{1}; yc=varargin{2}; zc=varargin{3};
	xp = r*cos(ang) + xc;
	yp = r*sin(ang) + yc;
	zp = zeros(numel(xp)) + zc;
	varargout = {zp};
else
	disp('bad job');
end

end


%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%

  


%=================================================%
%		Generate Variable Names from Workspace    
%=================================================%
%{

for i=1:10
  eval(sprintf('A%d = [1:i]', i));
end

%=====================================================%
% 	Generate Variable Names from Workspace
%=====================================================%
vnames = genvarname(vars{1}, who);
for m=1:11; eval([vnames{m} ' = vars{m+1};']); end;
%=====================================================%


%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%




%=================================================%
%			keyboard Debugging      
%=================================================%
%{

% debugging functions with:
% "keyboard"


keyboard
dbstack
dbquit





%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%



%=================================================%
%			PRINTING & SAVING      
%=================================================%
%{


set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 -tiff CIplots1
printpreview



%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%




%=================================================%
%			DRAWING POLYGONS      
%=================================================%
%{

%# draw a rectangle
rectangle('Position',[355 220 100 100], 'LineWidth',2, 'EdgeColor','b');


% pderect(xy) draws a rectangle with corner coordinates defined by xy=[xmin xmax ymin ymax].
% draw graph objects
pderect([-1 0 -1 0]) 
pdecirc(xc,yc,radius) 
pdepoly([-1 0 0 1 1 -1],[0 0 1 1 -1 -1]); 

%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%




%=================================================%
%			COOL FIGURES AND ANIMATIONS      
%=================================================%
%{



%=================================================%
% Generate The Klein Bottle
%===============================%
n = 12;
a = .2;                         % the diameter of the small tube
c = .6;                         % the diameter of the bulb
t1 = pi/4 : pi/n : 5*pi/4;      % parameter along the tube
t2 = 5*pi/4 : pi/n : 9*pi/4;    % angle around the tube
u  = pi/2 : pi/n : 5*pi/2;
[X,Z1] = meshgrid(t1,u);
[Y,Z2] = meshgrid(t2,u);

% The handle
len = sqrt(sin(X).^2 + cos(2*X).^2);
x1 = c*ones(size(X)).*(cos(X).*sin(X) ...
    - 0.5*ones(size(X))+a*sin(Z1).*sin(X)./len);
y1 = a*c*cos(Z1).*ones(size(X));
z1 = ones(size(X)).*cos(X) + a*c*sin(Z1).*cos(2*X)./len;
handleHndl=surf(x1,y1,z1,X);
set(handleHndl,'EdgeColor',[.5 .5 .5]);
hold on;

% The bulb
r = sin(Y) .* cos(Y) - (a + 1/2) * ones(size(Y));
x2 = c * sin(Z2) .* r;
y2 = - c * cos(Z2) .* r;
z2 = ones(size(Y)) .* cos(Y);
bulbHndl=surf(x2,y2,z2,Y);
set(bulbHndl,'EdgeColor',[.5 .5 .5])

colormap(hsv);
axis vis3d
view(-37,30);
axis off
light('Position',[2 -4 5])
light
hold off





% Half the bottle
%===============================%
shading interp
c = X;
[row col] = size(c);
c(1:floor(row/2),:) = NaN*ones(floor(row/2),col);
set(handleHndl,'CData',c);

c = Y;
[row col] = size(c);
c(1:floor(row/2),:) = NaN*ones(floor(row/2),col);
set(bulbHndl,'CData',c);
set([handleHndl bulbHndl],'FaceAlpha',1);


% Transparent Bottle
%===============================%
shading faceted;
set(handleHndl,'CData',X);
set(bulbHndl,'CData',Y);
set([handleHndl bulbHndl], ...
    'EdgeColor',[.5 .5 .5], ...
    'FaceAlpha',.5);
%=================================================%




%=================================================%
cla reset;
load topo;
[x y z] = sphere(45);
s = surface(x,y,z,'facecolor','texturemap','cdata',topo);
set(s,'edgecolor','none','facealpha','texture','alphadata',topo);
set(s,'backfacelighting','unlit');
colormap(topomap1);
alpha('direct');
alphamap([.1;1])
axis off vis3d;
campos([2 13 10]);
camlight;
lighting gouraud;
%=================================================%



%=================================================%
%# some sample images
I = imread('coins.png');
I_transp = imread('peppers.png');

%# create a gaussian mask for transparency
[r,c,~] = size(I_transp);
M = fspecial('gaussian', [r c], mean([r c]./5));
M = (M-min(M(:)))./range(M(:));

%# show overlayed images
figure, imshow(I, 'XData',[1 c], 'YData',[1 r]), hold on
hImg = imshow(I_transp);
set(hImg, 'AlphaData',M);
%=================================================%


%=================================================%
%# STREAM WIND (STAND ALONE)
figure
load wind
[sx sy sz] = meshgrid(80,20:10:50,0:5:15);
verts = stream3(x,y,z,u,v,w,sx,sy,sz);
div = divergence(x,y,z,u,v,w);
streamtube(verts,x,y,z,-div);
% Define viewing and lighting
view(3)
axis tight
shading interp
camlight; lighting gouraud
%=================================================%



%=================================================%
%# STREAM AWESOMENESS (STAND ALONE)
% Example ? Displaying Curl with Stream Ribbons
figure
% Define 3-D arrays x, y, z, u, v, w
xmin = -7; xmax = 7;
ymin = -7; ymax = 7; 
zmin = -7; zmax = 7; 
x = linspace(xmin,xmax,30);
y = linspace(ymin,ymax,20);
z = linspace(zmin,zmax,20);
[x y z] = meshgrid(x,y,z);
u = y; v = -x; w = 0*x+1;
[cx cy cz] = meshgrid(linspace(xmin,xmax,30),...
   linspace(ymin,ymax,30),[-3 4]);
h = coneplot(x,y,z,u,v,w,cx,cy,cz,'quiver');
set(h,'Color','k');

% Plot two sets of streamribbons
[sx sy sz] = meshgrid([-1 0 1],[-1 0 1],-6);
streamribbon(x,y,z,u,v,w,sx,sy,sz);
[sx sy sz] = meshgrid([1:6],[0],-6);
streamribbon(x,y,z,u,v,w,sx,sy,sz);

% Define viewing and lighting
shading interp
view(-30,20) ; axis off tight
camproj perspective; camva(66); camlookat; 
camdolly(0,0,.5,'fixtarget')
camlight
%=================================================%




%==========================================================%
MU1 = [1 2];
SIGMA1 = [2 0; 0 .5];
MU2 = [-3 -5];
SIGMA2 = [1 0; 0 1];
X = [mvnrnd(MU1,SIGMA1,1000);
mvnrnd(MU2,SIGMA2,1000)];
scatter(X(:,1),X(:,2),10,'.')
options = statset('Display','final');
obj = gmdistribution.fit(X,2,'Options',options);
hold on
h = ezcontour(@(x,y)pdf(obj,[x y]),[-8 6],[-8 6]);
hold off
%==========================================================%



%==========================================================%
% COOL PROBABILITY DENSITIES
mu = @(x) -1.9+.23*x;
x = 5:.1:15;
yhat = mu(x);
dy = -3.5:.1:3.5; sz = size(dy); k = (length(dy)+1)/2;
x1 =  7*ones(sz); y1 = mu(x1)+dy; z1 = normpdf(y1,mu(x1),1);
x2 = 10*ones(sz); y2 = mu(x2)+dy; z2 = normpdf(y2,mu(x2),1);
x3 = 13*ones(sz); y3 = mu(x3)+dy; z3 = normpdf(y3,mu(x3),1);
plot3(x,yhat,zeros(size(x)),'b-', ...
      x1,y1,z1,'r-', x1([k k]),y1([k k]),[0 z1(k)],'r:', ...
      x2,y2,z2,'r-', x2([k k]),y2([k k]),[0 z2(k)],'r:', ...
      x3,y3,z3,'r-', x3([k k]),y3([k k]),[0 z3(k)],'r:');
zlim([0 1]);
xlabel('X'); ylabel('Y'); zlabel('Probability density');
% grid on; view([-45 45]);
grid on; view([-45 25]);
%==========================================================%
% STEM PLOT DENSITY
mu = @(x) exp(-1.9+.23*x);
x = 5:.1:15;
yhat = mu(x);
x1 =  7*ones(1,5);  y1 = 0:4; z1 = poisspdf(y1,mu(x1));
x2 = 10*ones(1,7); y2 = 0:6; z2 = poisspdf(y2,mu(x2));
x3 = 13*ones(1,9); y3 = 0:8; z3 = poisspdf(y3,mu(x3));
plot3(x,yhat,zeros(size(x)),'b-', ...
      [x1; x1],[y1; y1],[z1; zeros(size(y1))],'r-', x1,y1,z1,'r.', ...
      [x2; x2],[y2; y2],[z2; zeros(size(y2))],'r-', x2,y2,z2,'r.', ...
      [x3; x3],[y3; y3],[z3; zeros(size(y3))],'r-', x3,y3,z3,'r.');
zlim([0 1]);
xlabel('X'); ylabel('Y'); zlabel('Probability');
grid on; view([-45 45]);
%==========================================================%
set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 -tiff TestImage
printpreview

%==========================================================%
x = [ 0 30  30   0  0]
y = [35 35 -35 -35 35]
z = ones(5)

Grid = unigrid(0,1,z)
mesh(Grid)

hplot = plot3(x,y,z)
grid on
box on
set(gca,'XGrid','on','YGrid','on','ZGrid','on')


%------------------------------------------%

[X,Y] = meshgrid(0:30,0:70)
% Z = zeros(length(Y),length(X))
Z = X.*0
Z=randn(size(Z))./20
surf(X,Y,Z)
mesh(X,Y,Z)
axis equal
axis([0 30 0 70 -10 10])
set(gca,'ZTick',[])
colormap(bone)

view([-25 12]);
rotate3d on

hTitle  = title ('Simulation Field');
hXLabel = xlabel('X (dendrite width)');
hYLabel = ylabel('Y (segment legth)');
set(gca,'FontName','Helvetica');
set([hTitle, hXLabel, hYLabel],'FontName','AvantGarde');
set([hXLabel, hYLabel],'FontSize',10);
set( hTitle,'FontSize',12,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
'XMinorTick','on','YMinorTick','on','YGrid','on', ...
'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'LineWidth',1);
set(gca,'ZTick',zeros(1,0),'YMinorTick','on',...
	'YColor',[0.3 0.3 0.3],...
	'XMinorTick','on',...
	'XColor',[0.3 0.3 0.3],...
	'TickDir','out',...
	'TickLength',[0.02 0.02],...
	'Projection','perspective',...
	'PlotBoxAspectRatio',[1.5 3.5 1],...
	'LineWidth',1,...
	'DataAspectRatio',[1 1 1],...
	'CameraViewAngle',8.85939248000315,...
	'CameraUpVector',[0.145721284742063 0.30974496228629 0.939586805735046],...
	'CameraTarget',[15.8888506772839 34.3644592005569 -0.545828345438358],...
	'CameraPosition',[-95.9099952128365 -203.275013750716 95.1335708267175]);

set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 -tiff TestImage
%======================================================================%

%==========================================================%


%==========================================================%
latlim = [-80 80];
lonlim = [100 -120];
figure('Color','white')
axesm('lambertstd','MapLatLimit',latlim,'MapLonLimit',lonlim, ...
   'Frame','on','Grid','on','MeridianLabel','off', ...
   'ParallelLabel','off')
axis off
%==========================================================%

%==========================================================%
x = 0:4; y=-2:2; s2 = 1/sqrt(2);
z(3,:,:) = [0 1 s2 0 -s2 -1 0].'*[1 1 1 1 1];
z(2,:,:) = [1 0 s2 1 s2 0 -1].'*[0 1 0 -1 0];
z(1,:,:) = [1 0 s2 1 s2 0 -1].'*[1 0 -1 0 1];
sph = csape({x,y},z,{'clamped','periodic'});
fnplt(sph), axis equal, axis off
figure
fnplt(fncmb(sph,[1 0 0; 0 0 1])), axis equal, axis off
%==========================================================%




%=================================================%
% Ball Rolling through Box
[x,y,z] = meshgrid(-2:.2:2,-2:.25:2,-2:.16:2);
v = x.*exp(-x.^2-y.^2-z.^2); % Create volume data
[xi,yi,zi] = sphere; % Plane to contour
contourslice(x,y,z,v,xi,yi,zi)
view(3)

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
 delete(hslicer)
 hold off
end
%=================================================%



%=================================================%
% 2D gradient with arrows and polygons
v = -2:0.2:2;
[x,y] = meshgrid(v);
z = x .* exp(-x.^2 - y.^2);
[px,py] = gradient(z,.2,.2);
contour(v,v,z), hold on, quiver(v,v,px,py), hold off

F(:,:,1) = magic(3); F(:,:,2) = pascal(3);
gradient(F)
[PX,PY,PZ] = gradient(F,0.2,0.1,0.2) 
%=================================================%


%=================================================%
% My 3D particle matrix
z = ones(4)
z = padarray(z,[2 2],0)
z = vertcat(z,z)
[FX,FY] = gradient(z)
h = streamslice(-FX,-FY)
set(h,'color','k')
for i=1:length(h); 
 zi = interp2(z,get(h(i),'xdata'),get(h(i),'ydata'));
 set(h(i),'zdata',zi);
end
view(30,50); axis tight
%=================================================%



%=================================================%
% 3D particle internalization
load wind
[sx sy sz] = meshgrid(80,20:1:55,5);
verts = stream3(x,y,z,u,v,w,sx,sy,sz);
sl = streamline(verts);
iverts = interpstreamspeed(x,y,z,u,v,w,verts,.025);
axis tight; view(30,30); daspect([1 1 .125])
camproj perspective; camva(8)
set(gca,'DrawMode','fast')
box on
streamparticles(iverts,35,'animate',10,'ParticleAlignment','on')
%=================================================%


%=================================================%
% 2D particle stream

load wind
daspect([1 1 1]); view(2)
[verts averts] = streamslice(x,y,z,u,v,w,[],[],[5]); 
sl = streamline([verts averts]);
axis tight off;
set(sl,'Visible','off')
iverts = interpstreamspeed(x,y,z,u,v,w,verts,.05);
set(gca,'DrawMode','fast','Position',[0 0 1 1],'ZLim',[4.9 5.1])
set(gcf,'Color','black')
streamparticles(iverts, 200, ...
    'Animate',100,'FrameRate',40, ...
    'MarkerSize',10,'MarkerFaceColor','yellow')
%=================================================%




%=================================================%
% Meshgrid with contour plot
figure
z = peaks;
surf(z)
shading interp
hold on
[c ch] = contour3(z,20); set(ch,'edgecolor','b')
[u v] = gradient(z); 
h = streamslice(-u,-v); 
set(h,'color','k')
for i=1:length(h); 
 zi = interp2(z,get(h(i),'xdata'),get(h(i),'ydata'));
 set(h(i),'zdata',zi);
end
view(30,50); axis tight 
%=================================================%



%=================================================%
% Meshgrid with images
load clown
surface(peaks,flipud(X),...
   'FaceColor','texturemap',...
   'EdgeColor','none',...
   'CDataMapping','direct')
colormap(map)
view(-35,45)
%=================================================%



%=================================================%
% Simple Meshgrids

[X,Y] = meshgrid(-2:.2:2, -2:.2:2)
Z = X .* exp(-X.^2 - Y.^2)
surf(X,Y)

[X,Y] = meshgrid(1:40, 1:40)
X = [1:40 1:40]
Y = [1:40 1:40]
Z = [ 1:40; ones(1,40); 21:60]
surf(Z)
%=================================================%





%=================================================%
% SOLVE EQUATIONS
a=2;
b=5.2;
eq1=sprintf('%d-(%d+1)*x+x^2*y',a,b);
eq2=sprintf('%d*x-x^2*y',b);
sol=solve(eq1,eq2,'x','y');
sol.x
sol.y

eq1=sprintf('2/(1+Y^4)-X=0');
eq2=sprintf('2/(1+X^4)-Y=0');
sol=solve(eq1,eq2,'X','Y');
sol.X
sol.Y
x=eval(sol.X)
y=eval(sol.Y)
k=find(imag(x)==0);
sol=[x(k) y(k)]
%=================================================%


%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%





%=================================================%
%		TriRep Superclass    DelaunayTri Class      
%=================================================%
%{



%=================================================%
%   TriRep Superclass    DelaunayTri Class
%=================================================%
%{
% use the PDE tool to draw a polygon and then export the polygon
pdepoly([-1 0 0 1 1 -1],[0 0 1 1 -1 -1]); 
pdetool

%------
poly1 = gd(3:end).*10;
polyX = poly1(1:8); polyX(9) = polyX(1);
polyY = poly1(9:16); polyY(9) = polyY(1);
polyXY = [polyX polyY]
XYL = ([1:9 6; 1:9 4].*1.0)';
PXm = mean(polyX)
PYm = mean(polyY)
XYL = ([PXm.*(randn(size(polyX)))*4 , PYm.*(randn(size(polyX)))*1.5])
%------


%------
IN = inpolygon(XYL(:,1),XYL(:,2),polyX,polyY)
dt = DelaunayTri(polyX,polyY)
e = edges(dt)
k = convexHull(dt)

figure(1)
hold off
plot(dt.X(:,1),dt.X(:,2), '.', 'markersize',10); hold on;
plot(dt.X(k,1),dt.X(k,2), 'r'); hold off; hold on;
triplot(dt);
%------


%------
[SI,BC] = pointLocation(dt,XYL)
[PI,D]  = nearestNeighbor(dt,XYL)
PDdt = dt.X(PI,:)

figure(2)
hold off
triplot(dt); hold on;
plot(PDdt(:,1),PDdt(:,2), '.', 'markersize',10); hold on;
scatter(XYL(:,1),XYL(:,2))
%------


%------
NPN = dsearchn(polyXY,XYL)
NPNmx = polyXY(NPN,1:2)

figure(3)
hold off
plot(polyX,polyY,'r')
hold on
scatter(NPNmx(:,1),NPNmx(:,2))
hold on
scatter(XYL(:,1),XYL(:,2))
%------
%}

%=================================================%
% delaunay triangulation
load seamount
plot(x,y,'.','markersize',12)
xlabel('Longitude'), ylabel('Latitude')
grid on
tri = delaunay(x,y);
hold on, triplot(tri,x,y), hold off
figure
hidden on
trimesh(tri,x,y,z)
xlabel('Longitude'),ylabel('Latitude'),zlabel('Depth in Feet');
%=================================================%


%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%





%=================================================%
%			nearestNeighbor      
%=================================================%
%{


%=================================================%
%			nearestNeighbor
%=================================================%
%{

[PI,D] = nearestNeighbor(dt,QX,QY)
PI = nearestNeighbor(DT,QX,QY,QZ)
[PI,D] = nearestNeighbor(DT,QX,...) 
% allow the query points to be specified in column vector format when working 
% in 2-D and 3-D. [PI,D] returns the index of the nearest point in DT.X for 
% each query point location in QX. The corresponding Euclidean distances between the 
% query points and their nearest neighbors are returned in D.


SI = pointLocation(DT,QX,QY)
SI = pointLocation(DT,QX,QY,QZ) 
% allow the query point locations to be specified in alternative column vector format 
% when working in 2-D and 3-D.
% [SI, BC] returns the barycentric coordinates BC.




k = dsearchn(X,T,XI) 
% returns the indices k of the closest points in X for each 
% point in XI. X is an m-by-n matrix representing m points in n-dimensional space. 
% XI is a p-by-n matrix, representing p points in n-dimensional space. 
% T is a numt-by-n+1 matrix, a triangulation of the data X generated by delaunayn. 
% The output k is a column vector of length p. 

k = dsearchn(X,T,XI,outval) 
% returns the indices k of the closest points in X for 
% each point in XI, unless a point is outside the convex hull. If XI(J,:) is outside 
% the convex hull, then K(J) is assigned outval, a scalar double. Inf is often used 
% for outval. If outval is [], then k is the same as in the case k = dsearchn(X,T,XI).

k = dsearchn(X,XI) 
% performs the search without using a triangulation. 
% With large X and small XI, this approach is faster and uses much less memory. 

[k,d] = dsearchn 
% also returns the distances d to the closest points. 
% d is a column vector of length p. 

%}






%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%




%=================================================%
%			USEFUL TRICKS      
%=================================================%
%{


%=================================================%
%               USEFUL TRICKS        ALT+Space
%=================================================%
%{

if 1-exist('in1','var')
	inval=Ndots;
end


----------------
linspace(a,b,n) 
generates a row vector of n column-values linearly spaced between a and b.
linspace(-50,50,101)'
----------------


----------------
if t > 100
    keyboard
end
----------------



B = reshape(A,2,[])
          
B =
    1    3    5    7    9   11
    2    4    6    8   10   12



----------------
profile off
format bank
format compact;
format short;
clear all;
close all;
----------------


----------------
padarray
A = [2 2 2; 2 1 2; 2 2 2];
B = padarray(A,[3 2],'replicate','post')

C = padarray(B, [3 3], 0)
----------------



----------------
circshift
Circularly shift first dimension values down by 1 and second dimension values to the left by 1.
B = circshift(A,[1 -1]);
B = 
    8     9     7
    2     3     1
    5     6     4
----------------


----------------
find
Find all non-negative ints and put their row,col index into two arrays
A = [2 2 2; 2 1 2; 2 2 2];
B = padarray(A, [3 3], 0)
C = (B==1)

[xylrow xylcol] = find(C)
----------------



DelLoc = (xylFRAP==99);
[xylrow xylcol] = find(DelLoc);
DelID = xylcol';
xylF = xylFRAP;
xylF(:,DelID) = [];


----------------
pause
Puts a pause(s) of 's' seconds into a loop. Usful right before 'draw now'
pause(.1) 
drawnow
hold off;
----------------



----------------
Dees=[.003:.003:3];
for i=1:1000
D=Dees(i);
Hight = abs((.1+D)/(D-.1));
Width = 1/((.1+D)/(D-.1));

%---
xodom = [-1 1]; yodom = [1 2.4];
plot([0 Width],[1 Hight]);
axis manual; axis([xodom, yodom]);
pause(.001)
end
----------------



----------------

----------------


----------------

----------------


%}




%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%





%=================================================%
%			KEYBOARD SHORTCUTS            
%=================================================%
%{


CMD+.         Abort 
CMD+K         Clear Command Window
CMD+SHFT+K    Clear Workspace Variables
CMD+ALT+L>    Highlight Line Left
CMD+ALT+R>    Highlight Line Right
CMD+ALT+U>    Highlight Line Up
CMD+ALT+D>    Highlight Line Down
CMD+ALT+ENT   Evaluate Highlighted Selection
ALT+SBAR      Comment Line
ALT+CTR+SBAR  Unomment Line



%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%



%=================================================%
%			CONVENTIONS OF MATLAB       
%=================================================%
%{

%{

m|n or r|c -- m-by-n 'm' rows by 'n' columns
i,j -- row i column j
Vc = vector
Mx = 



[ ]

Brackets are used to form vectors and matrices. [6.9 9.64 sqrt(-1)] is a
vector with three elements separated by blanks. [6.9, 9.64, i] is the same
thing. [1+j 2-j 3] and [1 +j 2 -j 3] are not the same. The first has three
elements, the second has five.

[11 12 13; 21 22 23] is a 2-by-3 matrix. The semicolon ends the first row.

Vectors and matrices can be used inside [ ] brackets. [A B;C] is allowed if
the number of rows of A equals the number of rows of B and the number of
columns of A plus the number of columns of B equals the number of columns
of C. This rule generalizes in a hopefully obvious way to allow fairly
complicated constructions.

A = [ ] stores an empty matrix in A. A(m,:) = [ ] deletes row m of A.
A(:,n) = [ ] deletes column n of A. A(n) = [ ] reshapes A into a column
vector and deletes the nth element.

[A1,A2,A3...] = function assigns function output to multiple variables.

For the use of [ and ] on the left of an "=" in multiple assignment
statements, see lu, eig, svd, and so on.

{ }

Curly braces are used in cell array assignment statements. For example,
A(2,1) = {[1 2 3; 4 5 6]}, or A{2,2} = ('str'). See help paren for more
information about { }.

( )

Parentheses are used to indicate precedence in arithmetic expressions in
the usual way. They are used to enclose arguments of functions in the usual
way. They are also used to enclose subscripts of vectors and matrices in a
manner somewhat more general than usual. If X and V are vectors, then X(V)
is [X(V(1)), X(V(2)), ..., X(V(n))]. The components of V must be integers
to be used as subscripts. An error occurs if any such subscript is less
than 1 or greater than the size of X. Some examples are

X(3) is the third element of X.

X([1 2 3]) is the first three elements of X.

See help paren for more information about ( ).



If X has n components, X(n:?1:1) reverses them. The same indirect subscripting works in matrices. If V has m components and W has n components, then A(V,W) is the m-by-n matrix formed from the elements of A whose subscripts are the elements of V and W. For example, A([1,5],:) = A([5,1],:) interchanges rows 1 and 5 of A.


%{
%=================================================%
%               TEMPLATE_BOX     
%=================================================%
%{.
Mx = randi([-10 10],3,5)
%}
%{%}
%}






%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%




%=================================================%
%			RANDOM NUMBER GENERATORS           
%=================================================%
%{


%{

randi([imin,imax],m,n)        % m-by-n Mx of 
Mx = randi([-10 10],4,3)

randi(10) % returns a scalar



r = randn(m,n)  Normally distributed pseudorandom numbers



Generate values from a normal distribution with 
mean = 1
sd = 2

r = 1 + 2.*randn(2,20)

r = 1 + .2.*randn(2,20)



r = randn(size(A)) returns an array the same size as A.

%{
%=================================================%
%               TEMPLATE_BOX     
%=================================================%
%{.
Mx = randi([-10 10],3,5)
%}
%}





%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%




%=================================================%
%			RANDOM NUMBERS DISTRIBUTIONS      
%=================================================%
%{


%{

Mx = randi([-10 10])

if you have the stats toolbox..
Random Number Generators.
    betarnd - Beta random numbers.
    binornd - Binomial random numbers.
    chi2rnd - Chi square random numbers.
    exprnd - Exponential random numbers.
    frnd - F random numbers.
    gamrnd - Gamma random numbers.
    geornd - Geometric random numbers.
    hygernd - Hypergeometric random numbers.
    lognrnd - Lognormal random numbers.
    mvnrnd - Multivariate normal random numbers.
    mvtrnd - Multivariate t random numbers.
    nbinrnd - Negative binomial random numbers.
    ncfrnd - Noncentral F random numbers.
    nctrnd - Noncentral t random numbers.
    ncx2rnd - Noncentral Chi-square random numbers.
    normrnd - Normal (Gaussian) random numbers.
    poissrnd - Poisson random numbers.
    random - Random numbers from specified distribution.
    raylrnd - Rayleigh random numbers.
    trnd - T random numbers.
    unidrnd - Discrete uniform random numbers.
    unifrnd - Uniform random numbers.
    weibrnd - Weibull random numbers.

%{
%=================================================%
%               TEMPLATE_BOX     
%=================================================%
%{.
Mx = randi([-10 10],3,5)
%}
%}



%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%



%=================================================%
%			ROUNDING      
%=================================================%
%{



x = fix(A)  % rounds towards zero
x = int8(A) % rounds to interger towards zero 
B = floor(A) % rounds toward negative 
x = fix(10*rand(1,10))  % rounds toward zero
x = int8(x)



floor
Rounds each element of the input signal to the nearest integer value towards minus infinity.

ceil
Rounds each element of the input signal to the nearest integer towards positive infinity.

round
Rounds each element of the input signal to the nearest integer.

fix
Rounds each element of the input signal to the nearest integer towards zero.




%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%


%=================================================%
%			LOOPS      
%=================================================%
%{

%{

for i = 1:10
    for j = 1:10
        HvecI(i,j) = PSD(i,(j-1)) + PSD(i,(j+1));
        HvecJ(i,j) = PSD((i-1),j) + PSD((i+1),j);
    end
end


%{
%=================================================%
%               TEMPLATE_BOX     
%=================================================%
%{.
Mx = randi([-10 10],3,5)
%}
%{%}
%}




%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%



%=================================================%
%			FORCE SAVE VARIABLE TO MAIN WORKSPACE           
%=================================================%
%{


%{

putvar(var1, var2, var3, var4)
assignin('base', 'VarName', VarValue)

%{
%=================================================%
%               TEMPLATE_BOX     
%=================================================%
%{.
Mx = randi([-10 10],3,5)
%}
%{%}
%}



%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%



%=================================================%
%			STRINGS      
%=================================================%
%{


%{

head = ['    Time     Dendrite     PSD1      PSD2'];       
       data = [stepN PSD0n PSD1n PSD2n;stepN PSD0 PSD1 PSD2];
       disp(' ')
       disp(head)
       disp(data)

%}



%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%




%=================================================%
%			ANIMATION      
%=================================================%
%{

%{


-----------------------
linkdata on
-----------------------
%{
x = [1:20];
y = rand(20,3);
area(x,y)
linkdata on
----------------------------
%} 



-----------------------
refreshdata(hf,'caller')
-----------------------
%{
hold off
x = rand(100,1)*4-2
y = rand(100,1)*4-2
z = x.*exp(-x.^2-y.^2)

F = TriScatteredInterp(x,y,z);

fig1 = figure(1)
ti = -2:.25:2;
[qx,qy] = meshgrid(ti,ti)
qz = F(qx,qy);
hm = mesh(qx,qy,qz)
set(hm,'XDataSource','qx')
set(hm,'YDataSource','qy')
set(hm,'ZDataSource','qz')
hold on;
hp = plot3(x,y,z,'o');
set(hp,'XDataSource','x')
set(hp,'YDataSource','y')
set(hp,'ZDataSource','z')


for i = 1:20
z = i.*exp(-x.^2-y.^2);
qz = qz./(.8);

pause(.1)
refreshdata(hm,'caller')
refreshdata(hp,'caller')
drawnow
end
----------------------------
%} 




%}




%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%



%=================================================%
%			XXXXX      
%=================================================%
%{






%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%


%=================================================%
%        MAKE MOVIE     
%=================================================%
%{

eye(10)

-----------------------
getframe(gcf)
-----------------------
%{
In this example, the movie frame contains the entire figure. 
To play the movie so that it looks like the original figure, 
make the playback axes fill the figure window.

h = uicontrol('style','slider','position',...
	[10 50 20 300],'Min',1,'Max',16,'Value',1)
for k = 1:16
	plot(fft(eye(k+16)))
	axis equal
	set(h,'Value',k)
	M(k) = getframe(gcf);
end


clf
axes('Position',[0 0 1 1])
movie(M,4)
----------------------------
%} 






%--------------------------%
% SCATTER PLOT TO MOVIE
%--------------------------%
%{
Vid = VideoWriter('vidtest');
open(Vid);

Mx = [1:5; 1:5; 16:20; 1:5; 16:20]
My = [1:5; 16:20; 16:20; 11:15; 1:5]

figure('Renderer','zbuffer');
gscatter(Mx, My);
axis tight
set(gca,'nextplot','replacechildren');
currFrame = getframe(gcf);
writeVideo(Vid,currFrame);

for k = 1:20 
   Mx = randi(20,5);
   My = randi(20,5);
   gscatter(Mx, My);
   currFrame = getframe(gcf);
   writeVideo(Vid,currFrame);
end
close(Vid);
%}




%--------------------------%
% EXAMPLE CODE
%--------------------------%
%{

Vid = VideoWriter('vidtest');
open(Vid);

Mx = randi(20,5);
surf(Mx) 
axis tight
set(gca,'nextplot','replacechildren');

for k = 1:20 
   surf(sin(2*pi*k/20)*Mx, Mx)
   
   % Write each frame to the file.
   currFrame = getframe;
   writeVideo(Vid,currFrame);
end

close(Vid);

%}




%--------------------------%
% Autonumber and evaluate
%--------------------------%
%{



for d=1:10
   s = ['MovieFrame' num2str(d) '=Var' num2str(d) ';']
   eval(s)
end

s =
   MovieFrame1=Var1;
s =
   MovieFrame2=Var2;
s =
   MovieFrame3=Var3;
%}




%{
%=================================================%
%               TEMPLATE_BOX     
%=================================================%
%{.
Mx = randi([-10 10],3,5)
%}
%{%}
%}
%=================================================%
%=================================================%







%=================================================%
%        COLOR MAP FROM MATRIX     
%=================================================%
%{


% ----------------------------------------- %
% Turn a matrix into a color map
% ----------------------------------------- %
Map(6:12, 6:12) = 1
LogicMap = (Map>0)

mask=[
    0 1 0; 
    1 0 1; 
    0 1 0];

LogicMask = convn(LogicMap,mask,'same')

figure(1)
imagesc(LogicMask)
drawnow
% ----------------------------------------- %


%{
%=================================================%
%               TEMPLATE_BOX     
%=================================================%
%{.
Mx = randi([-10 10],3,5)
%}
%{%}
%}
%=================================================%
%=================================================%




%=================================================%
%       CDF & PDF     
%=================================================%
%{

%-- pdf of X in a normal distribution --%
normpdf(x,mu,sigma)
normpdf(1,0,1) % .24 == pdf('norm',1,0,1)


%{
%=================================================%
%               TEMPLATE_BOX     
%=================================================%
%{.
Mx = randi([-10 10],3,5)
%}
%{%}
%}
%=================================================%
%=================================================%





%=================================================%
%       GRAPHING    
%=================================================%
%{
get(gcf)
set(gcf,'PropName','val')


figure(3)
cdfexit = gcf;
figure(cdfexit)
subplot(5,4,[1 12]),GluR2plot = cdfplot(GluR2exT30);
axis([xlim,ylim]);
set(GluR2plot,'color',[1 .3 1])
hold on
subplot(5,4,[1 12]),cdfplot(GluR1exT30)
CDFtitle = title(['CDF of Fraction Exited    D =' num2str(D) '   GluR2 N =' int2str(Sdots)...
	'    GluR1 N =' int2str(Sdots2)]);
set(CDFtitle, 'FontSize', 16);
leg1=legend('  GluR2', 'GluR1');
set(leg1,'Location','SouthEast');
figure(cdfexit)
subplot(5,4,[13 14]),cdfplot(exittime)	
subplot(5,4,[15 16]),ecdfhist(coECDF_G2, coX_G2)
subplot(5,4,[17 18]),cdfplot(exittime2)	
subplot(5,4,[19 20]),ecdfhist(coX_G1, coX_G1)



% set(0,'DefaultFigureWindowStyle','normal')
% set(gcf,'OuterPosition',[xorigin,yorigin,width,height])
% subplot('Position',[left bottom width height]);

figure(3)
subplot(2,2,4), subplot('Position',[.55 .05 .4 .4]),...
	plot([(Ddata(:,1)) (Ddata(:,2))]);
yt = get(gca,'YTick');
set(gca,'YTickLabel', sprintf('%.4f|',yt))
set(get(gca,'YLabel'),'String','scaled diffusion rate')
xt = (get(gca,'XTick'))*10;
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time-Step')
title('PSD DIFFUSION RATE'); leg4=legend('D-PSD1', 'D-PSD2');
set(leg4,'Location','SouthWest');


% overlay alpha
close all
clear all
fPSD1 = ones(6);                        % PSD1 SIZE
fPSD2 = ones(6);                        % PSD2 SIZE
pfPSD1 = padarray(fPSD1,[8 8], 0);      % PAD PSD1
pfPSD2 = padarray(fPSD2,[8 8], 0);      % PAD PSD2
PSDfield = cat(1, pfPSD1, pfPSD2);      % CONCAT PSD FIELDS
map = PSDfield

figure(1)
subplot(5,5,[3 25]), A = bar(10*rand(1,40))
Ah = get(A,'child')
set(Ah,'facea',.3)
hold on
figure(1)
subplot(5,5,[3 25]), B=imagesc(map)
Bh = get(B)
set(Bh,'facea',.3)


rectangle('Position',[10 10 6 6], 'LineWidth',2, 'EdgeColor','b');


just get the handles from those axes 
and use them as first argument in the plot
figure
hax1=axes
figure
hax2=axes
plot (hax1,t(:,1),t(:,2),'-r+')
plot (hax2,t(:,1),t(:,3),'-r+')


figure(1);
subplot(5,5,[3 25]), pfield = imagesc(PSDfield);
set(gca, 'nextplot', 'add')
hold on


rectangle('Position',[x,y,w,h]) 
draws the rectangle from the point x,y and having a width of w and 
a height of h. Specify values in axes data units.

Note that, to display a rectangle in the specified proportions, 
you need to set the axes data aspect ratio so that one unit is of 
equal length along both the x and y axes. You can do this with the 
command 'axis equal' or daspect([1,1,1]).

%---
figure(1);
subplot(5,5,[3 25]), 
scat2d = gscatter(Caxyl(1,:),Caxyl(2,:)); view(0, 80); % <--FIG-----######
axis image;
grid on;
set(scat2d,'marker','.','markersize',[2])
%set(scat2d,'marker','.', 'EraseMode','xor')
axis([xlim, ylim, zlim]);
set(gca, 'Box', 'off', 'GridLineStyle','-');
%---
keyboard
get(gcf)



-------------
LIGHTING
surf(X, Y, Z); view(30, 30);
shading interp;
light;
lighting phong;
title('lighting phong', 'FontName', 'Courier', 'FontSize', 14);
-------------


figure(4)
subplot(2,2,3), plot([SAPdata(:,1) SAPdata(:,2)]);
yt = get(gca,'YTick');
set(gca,'YTickLabel', sprintf('%.4f|',yt))
set(get(gca,'YLabel'),'String','scaled diffusion rate')
xt = (get(gca,'XTick'))*10;
set(gca,'XTickLabel', sprintf('%.1f|',xt))
set(get(gca,'XLabel'),'String','Time-Step')
title('PSD DIFFUSION RATE'); leg4=legend('D-PSD1', 'D-PSD2');
set(leg4,'Location','SouthWest');
plotinset = get(gca,'TightInset');
subplot('Position',[plotinset(1) plotinset(2) .4 .4]),...
	plot([SAPdata(:,1) SAPdata(:,2)]);
set(get(gca,'YLabel'),'String','scaled diffusion rate')
set(get(gca,'XLabel'),'String','Time-Step')


% h1 = zoom;
% h2 = pan;
% figure(2)
% subplot(3,2,6), plot([(Ddata(:,1)) (Ddata(:,2))]);
% set(gca,'XTickLabelMode','auto')
% set(gca,'XTickLabel',num2str(get(gca,'XTick').'))
% set(h1,'ActionPostCallback',@mypostcallbackX);
% set(h2,'ActionPostCallback',@mypostcallbackX);
% title('PSD DIFFUSION RATE'); leg4=legend('D-PSD1', 'D-PSD2');
% set(leg4,'Location','SouthWest');
% function mypostcallbackX(obj,evd)
% set(gca,'XTickLabelMode','auto')
% set(gca,'XTickLabel',num2str(get(gca,'XTick').'))
% end



%-- Scatter Plot from Matrix Data --%
Mx = [1:5; 1:5; 1:5; 1:5; 1:5;]
My = Mx+5
gscatter(Mx, My)
load seamount
figure
hs(1) = subplot(2,1,1);
hs(2) = subplot(2,1,2);
scatter3(hs(1),x,y,z,'MarkerFaceColor',[0 .75 .75])
scatter3(hs(2),x,y,z,'*')


%-- Axis --%
axis([xmin xmax ymin ymax])
axis([xmin xmax ymin ymax zmin zmax cmin cmax])
v = axis
axis auto
axis manual
axis tight      % sets the axis limits to the range of the data.
axis fill
axis ij
axis xy
axis equal
axis image
axis square
axis vis3d
axis normal
axis off
axis on

axis(axes_handles,...)      % applies the axis command to the specified axes. 
% For example, the following statements:
h1 = subplot(221);
h2 = subplot(222);
axis([h1 h2],'square')

axis(subplot(2), off)

axis(axes_handles,...)
[mode,visibility,direction] = axis('state')

axis manual and axis(axis) 
freezes the scaling at the current limits, 
so that if hold is on, subsequent plots use the same limits. 
This sets the XLimMode, YLimMode, and ZLimMode properties to manual.





gplotmatrix(ratings(:,1:2),ratings(:,[4 7]),group,... 
            'br','.o',[],'on','',categories(1:2,:),... 
             categories([4 7],:))


%}
%=================================================%
%=================================================%



%=================================================%
%			GUI Interface          
%=================================================%
%{


%{
%=================================================%
%               GUI Interface    
%=================================================%
%{.



function ReDiClusMenu(varargin)
% function SpikingNeurons(varargin)
%
% This Graphical User Interface is generated with user-supplied menu definition data
% buildMainGUI.m version 1.0 has been used to generate this GUI
%
fh = figure('Visible','off', ...
            'NumberTitle','off', ...  % turns off MATLAB GUI window heading
            'Name','SpikingNeurons', ...    % now, define my own
            'Units','normalized', ...
            'Units','normalized', ...
            'Position',[1.000000e-001 1.000000e-001 6.000000e-001 6.000000e-001], ...
            'Color','white', ...   % match bg of MODE_new.jpg
            'Resize','off');
set(fh,'MenuBar','none'); % removes MATLAB default menu bar
% create custom menu bar items as defined by GID loaded above:
% "File", "Model", "Articles", "Tutorials", "Examples", "Run", "Code", "Help"
%
% The menu item handles used here are all local, the naming convention is
% mh_ijk where i, j, and k represent the menu item and sub menu items
% First, set up all the menu bars (but without any callbacks . . .
% The menu bars are: File, Model, Articles, Tutorial, Examples, Run, Code, Help
% For menu bar item: 'File'
mh1 = uimenu(fh,'Label', 'File');
mh(1,1,1,1) = uimenu(mh1,'Label', 'Exit');


% For menu bar item: 'Model'
mh2 = uimenu(fh,'Label', 'Model');
mh(2,1,1,1) = uimenu(mh2,'Label', 'System diagram');


% For menu bar item: 'Articles'
mh3 = uimenu(fh,'Label', 'Articles');
mh(3,1,1,1) = uimenu(mh3,'Label', 'Izhikevich (2003)');


% For menu bar item: 'Tutorial'
mh4 = uimenu(fh,'Label', 'Tutorial');
mh(4,1,1,1) = uimenu(mh4,'Label', 'Abstract');
mh(4,2,1,1) = uimenu(mh4,'Label', 'Tutorial');


% For menu bar item: 'Examples'
mh5 = uimenu(fh,'Label', 'Examples');
mh(5,1,1,1) = uimenu(mh5,'Label', 'Preset image 1: Gail');
mh(5,2,1,1) = uimenu(mh5,'Label', 'Preset image 2: Wedding');
mh(5,3,1,1) = uimenu(mh5,'Label', 'Preset image 3: Dog');
mh(5,4,1,1) = uimenu(mh5,'Label', 'Preset image 4: Optical illusion1');
mh(5,5,1,1) = uimenu(mh5,'Label', 'Preset image 5: Optical illusion2');
mh(5,6,1,1) = uimenu(mh5,'Label', 'Preset image 6: Robert');


% For menu bar item: 'Run'
mh6 = uimenu(fh,'Label', 'Run');
mh(6,1,1,1) = uimenu(mh6,'Label', 'Load an image and spike!');


% For menu bar item: 'Code'
mh7 = uimenu(fh,'Label', 'Code');
mh(7,1,1,1) = uimenu(mh7,'Label', 'Source code');


% For menu bar item: 'Help'
mh8 = uimenu(fh,'Label', 'Help');
mh(8,1,1,1) = uimenu(mh8,'Label', 'Contact');
mh(8,2,1,1) = uimenu(mh8,'Label', 'Credit');
mh(8,3,1,1) = uimenu(mh8,'Label', 'License');




% Next, setup the callbacks for all relevant menus

set(mh(1,1,1,1),'callback', {@exit_Callback});
set(mh(2,1,1,1),'callback', {@opendoc1_Callback 'System_diagram.pdf'});
set(mh(3,1,1,1),'callback', {@opendoc1_Callback 'Izh_spikes.pdf'});
set(mh(4,1,1,1),'callback', {@opendoc1_Callback 'Abstract.pdf'});
set(mh(4,2,1,1),'callback', {@opendoc1_Callback 'Tutorial.pdf'});
set(mh(5,1,1,1),'callback', {@run1_Callback 'reducedSpikingDemo(1)'});
set(mh(5,2,1,1),'callback', {@run1_Callback 'reducedSpikingDemo(2)'});
set(mh(5,3,1,1),'callback', {@run1_Callback 'reducedSpikingDemo(3)'});
set(mh(5,4,1,1),'callback', {@run1_Callback 'reducedSpikingDemo(4)'});
set(mh(5,5,1,1),'callback', {@run1_Callback 'reducedSpikingDemo(5)'});
set(mh(5,6,1,1),'callback', {@run1_Callback 'reducedSpikingDemo(6)'});
set(mh(6,1,1,1),'callback', {@run1_Callback 'reducedSpikingDemo(0)'});
set(mh(7,1,1,1),'callback', {@opendoc1_Callback 'HTML/Index.html'});
set(mh(8,1,1,1),'callback', {@opendoc1_Callback 'Contact.html'});
set(mh(8,2,1,1),'callback', {@opendoc1_Callback 'Credit.html'});
set(mh(8,3,1,1),'callback', {@opendoc1_Callback 'License.html'});
% Displays a user-provided image on this GUI
[X,map] = imread('izhik.gif');
imshow(X,map); % display image on front page


set(fh,'Visible','on');
%}



%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%


%=================================================%
%       VECTORS     
%=================================================%
%{


% ----------------------------------------- %
% Vectors
% ----------------------------------------- %

vcol(3:5,1) = vrow(1)

% row    elements use blanks or commas
% column elements use semicolons
% Surround the entire list of elements with square brackets, [ ]

% vectors
vrow = [1 2 3 4 5 6 7 8 9];  % 1 row 9 col
vcol = [1 2 3 4 5 6 7 8 9]'; % 9 row 1 col 

v(5:end)     % fifth through the last elements
v(2:end-1)   % second through the next-to-last elements
v(1:2:end)   % Extract all the odd elements
v(end:-1:1)  % Reverse the order of elements

v([2 3 4]) = [10 15 20]   % Replace some elements of v
v([2 3]) = 30             % Replace 2nd and 3rd elements with 30

vrow, vcol      % row and column vectors

vrow(3)         % 3rd col in vrow
vrow([3 7])     % 3rd and 7th col in vrow
vcol([3 7])     % 3rd and 7th row in vcol
        % In this way vectors are not as picky as matrixes 
vswap = vrow([5:9 1:4])

         % Notice how these two are not the same
vrow*vcol   % ans=285
vcol*vrow   % ans=<9x9>
        % ... Bringing us to Matrixes and Linear Algebra






%{
%=================================================%
%               TEMPLATE_BOX     
%=================================================%
%{.
Mx = randi([-10 10],3,5)
%}
%{%}
%}
%=================================================%
%=================================================%



%=================================================%
%           Matrixes
%=================================================%
%{

%--------------------------%
B = reshape(A,m,n) 
returns the m-by-n matrix B whose elements are taken 
column-wise from A. An error results if A does not have m*n elements.
%--------------------------%


%--------------------------%
%    Updating the Matrix
%--------------------------%

A(A > 0.5) = 0






%--------------------------%
%    Matrix Basics
%--------------------------%

A = magic(4);
A(2,4)          % Extract the element in row 2, column 4
A(2:4,1:2)      % In a matrix, both subscripts are vectors
A(3,:)          % Extract third row
A(:,end)        % Extract last column

% The element in row i and column j of A is denoted by A(i,j). 
% For example, A(4,2) is the number in the fourth row and second column. 

rand(1,10)          % 10 numbers from 0 to 1
10*rand(1,10)       % 10 numbers from 1 to 10
fix(10*rand(1,10))  % rounded to int <X






m = [1:5;6:10;11:15;16:20;21:25];
n = [1:5;6:10;11:15;16:20;21:25]';
% m =
%           1.00          2.00          3.00          4.00          5.00
%           6.00          7.00          8.00          9.00         10.00
%          11.00         12.00         13.00         14.00         15.00
%          16.00         17.00         18.00         19.00         20.00
%          21.00         22.00         23.00         24.00         25.00

sum(m)      % Sums the columns

r = [1:3;4:6;7:9]
s = [1:3;7:9;2:4]' 

% equivelant
r+s
s+r

% not equivelant
r*s
s*r

% equivelant
r.*s
s.*r




%{
%=================================================%
%               TEMPLATE_BOX     
%=================================================%
%{.
Mx = randi([-10 10],3,5)
%}
%{%}
%}
%=================================================%
%=================================================%




%=================================================%
%       LOGICAL INDEXING     
%=================================================%
%{

%--------------------------%
%       Logical indexing
%--------------------------%

A = magic(4)
    16     2     3
     5    11    10    
     9     7     6

B = isprime(A)
     0     1     1     
     1     1     0     
     0     1     0     

% Logical indexing
A(~B) = 0;                       

A =
     0     2     3 
     5    11     0     
     0     7     0     

find(B)





%{
%=================================================%
%               TEMPLATE_BOX     
%=================================================%
%{.
Mx = randi([-10 10],3,5)
%}
%{%}
%}
%=================================================%
%=================================================%



%=================================================%
%       WRITING DATA OUTPUT TO FILE     
%=================================================%
%{

time1 = datestr(clock, 0);
time2 = num2str(time(1,[1:2 16:17 19:20]));
time3 = strcat(['Cluster1data' time2 '.txt']);
time4 = strcat(['Cluster2data' time2 '.txt']);
XLSdata1 = dataset({data1, 'Time' 'SAP' 'ES' 'PSD1' 'PSD2'});
export(XLSdata1,'file',eval('time3'));
XLSdata2 = dataset({data2, 'Time' 'SAP' 'ES' 'PSD1' 'PSD2'});
export(XLSdata2,'file','time4');





%{
%=================================================%
%               TEMPLATE_BOX     
%=================================================%
%{.
Mx = randi([-10 10],3,5)
%}
%{%}
%}
%=================================================%
%=================================================%



%=================================================%
%			POLYGONS AND POLYBOOL      
%=================================================%
%{

%=================================================%
%       [IN ON] = inpolygon(X,Y,xv,yv)     
%=================================================%
%{

Description

IN = inpolygon(X,Y,xv,yv) returns a matrix IN the same size as X and Y. 
Each element of IN is assigned the value 1 or 0 depending on whether the 
point (X(p,q),Y(p,q)) is inside the polygonal region whose vertices are 
specified by the vectors xv and yv.

[IN ON] = inpolygon(X,Y,xv,yv) returns a second matrix ON the 
same size as X and Y. Each element of ON is assigned the value 1 or 0 
depending on whether the point (X(p,q),Y(p,q)) is on the boundary of 
the polygonal region whose vertices are specified by the vectors xv and yv.

L = linspace(0,2.*pi,6); xv = cos(L)';yv = sin(L)';
xv = [xv ; xv(1)]; yv = [yv ; yv(1)];
x = randn(250,1); y = randn(250,1);
in = inpolygon(x,y,xv,yv);
plot(xv,yv,x(in),y(in),'r+',x(~in),y(~in),'bo')

lat2 = [0 1 -1 0]';
lon2 = [0 2  2 0]';
axesm miller
plotm(lat2,lon2,'r')

[latb,lonb] = bufferm(madagascar.Lat, madagascar.Lon, .75, 'in');
geoshow(latb, lonb, 'DisplayType', 'polygon', 'FaceColor', 'green')

%}

lat2 = [0 1 -1 0]';
lon2 = [0 2  2 0]';
axesm miller
plotm(lat2,lon2,'r')

xv = [XYLTpr1(1) XYLBpr1(1) XYRBpr1(1) XYRTpr1(1) XYLTpr1(1)]';
yv = [XYLTpr1(2) XYLBpr1(2) XYRBpr1(2) XYRTpr1(2) XYLTpr1(2)]';
plot(xv,yv);
IN = inpolygon(GluR2xyl(1,:)',GluR2xyl(2,:)',xv,yv);
plot(xv,yv,GluR2xyl(1,IN),GluR2xyl(2,IN),'r+',GluR2xyl(1,~IN),GluR2xyl(2,~IN),'bo');

xv = [XYLTpr2(1) XYLBpr2(1) XYRBpr2(1) XYRTpr2(1) XYLTpr2(1)]';
yv = [XYLTpr2(2) XYLBpr2(2) XYRBpr2(2) XYRTpr2(2) XYLTpr2(2)]';
plot(xv,yv);
IN = inpolygon(GluR2xyl(1,:)',GluR2xyl(2,:)',xv,yv);
plot(xv,yv,GluR2xyl(1,IN),GluR2xyl(2,IN),'r+',GluR2xyl(1,~IN),GluR2xyl(2,~IN),'bo');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%=================================================%
%       [x,y] = polybool(flag,x1,y1,x2,y2)     
%=================================================%
%{


Description


[x,y] = polybool(flag,x1,y1,x2,y2) performs the polygon set operation 
identified by flag. A valid flag string is any one of the following alternatives:

Region intersection: 'intersection', 'and', '&'

Region union: 'union', 'or', '|', '+', 'plus'

Region subtraction: 'subtraction', 'minus', '-'

Region exclusive or: 'exclusiveor' , 'xor'


Ax = {[1 1 6 6 1], [2 5 5 2 2], [2 5 5 2 2]};
Ay = {[1 6 6 1 1], [2 2 3 3 2], [4 4 5 5 4]};
Bx = {[0 0 7 7 0], [1 3 3 1 1], [4 6 6 4 4]};
By = {[0 7 7 0 0], [1 1 6 6 1], [1 1 6 6 1]};

%}

Ax = {[XYLTp1(1) XYRTp1(1) XYRBp1(1) XYLBp1(1) XYLTp1(1)],...
	  [XYLTp2(1) XYRTp2(1) XYRBp2(1) XYLBp2(1) XYLTp2(1)]};
Ay = {[XYLTp1(2) XYRTp1(2) XYRBp1(2) XYLBp1(2) XYLTp1(2)],...
	  [XYLTp2(2) XYRTp2(2) XYRBp2(2) XYLBp2(2) XYLTp2(2)]};

Bx = {[XYLTpr1(1) XYRTpr1(1) XYRBpr1(1) XYLBpr1(1) XYLTpr1(1)],...
      [XYLTpr2(1) XYRTpr2(1) XYRBpr2(1) XYLBpr2(1) XYLTpr2(1)]};
By = {[XYLTpr1(2) XYRTpr1(2) XYRBpr1(2) XYLBpr1(2) XYLTpr1(2)],...
      [XYLTpr2(2) XYRTpr2(2) XYRBpr2(2) XYLBpr2(2) XYLTpr2(2)]};

subplot(2, 3, 1)
[f, v] = poly2fv(Ax, Ay);
patch('Faces', f, 'Vertices', v, 'FaceColor', 'r', ...
  'EdgeColor', 'none')
axis equal, hold on
for k = 1:numel(Ax), plot(Ax{k}, Ay{k}, 'Color', 'k'), end
title('A')

subplot(2, 3, 4);
[f, v] = poly2fv(Bx, By);
patch('Faces', f, 'Vertices', v, 'FaceColor', 'r', ...
  'EdgeColor', 'none')
axis equal, hold on
for k = 1:numel(Bx), plot(Bx{k}, By{k}, 'Color', 'k'), end
title('B')

subplot(2, 3, 2)
[Cx, Cy] = polybool('union', Ax, Ay, Bx, By);
[f, v] = poly2fv(Cx, Cy);
patch('Faces', f, 'Vertices', v, 'FaceColor', 'r', ...
  'EdgeColor', 'none')
axis equal, hold on
for k = 1:numel(Cx), plot(Cx{k}, Cy{k}, 'Color', 'k'), end
title('A \cup B')

subplot(2, 3, 3)
[Dx, Dy] = polybool('intersection', Ax, Ay, Bx, By);
[f, v] = poly2fv(Dx, Dy);
patch('Faces', f, 'Vertices', v, 'FaceColor', 'r', ...
  'EdgeColor', 'none')
axis equal, hold on
for k = 1:numel(Dx), plot(Dx{k}, Dy{k}, 'Color', 'k'), end
title('A \cap B')

subplot(2, 3, 5)
[Ex, Ey] = polybool('subtraction', Ax, Ay, Bx, By);
[f, v] = poly2fv(Ex, Ey);
patch('Faces', f, 'Vertices', v, 'FaceColor', 'r', ...
  'EdgeColor', 'none')
axis equal, hold on
for k = 1:numel(Ex), plot(Ex{k}, Ey{k}, 'Color', 'k'), end
title('A - B')

subplot(2, 3, 6)
[Fx, Fy] = polybool('xor', Ax, Ay, Bx, By);
[f, v] = poly2fv(Fx, Fy);
patch('Faces', f, 'Vertices', v, 'FaceColor', 'r', ...
  'EdgeColor', 'none')
axis equal, hold on
for k = 1:numel(Fx), plot(Fx{k}, Fy{k}, 'Color', 'k'), end
title('XOR(A, B)')






PSD1xv = [XYLTp1(1) XYRTp1(1) XYRBp1(1) XYLBp1(1) XYLTp1(1)]';
PSD1yv = [XYLTp1(2) XYRTp1(2) XYRBp1(2) XYLBp1(2) XYLTp1(2)]';
PSD2xv = [XYLTp2(1) XYRTp2(1) XYRBp2(1) XYLBp2(1) XYLTp2(1)]';
PSD2yv = [XYLTp2(2) XYRTp2(2) XYRBp2(2) XYLBp2(2) XYLTp2(2)]';

PERI1xv = [XYLTpr1(1) XYRTpr1(1) XYRBpr1(1) XYLBpr1(1) XYLTpr1(1)]';
PERI1yv = [XYLTpr1(2) XYRTpr1(2) XYRBpr1(2) XYLBpr1(2) XYLTpr1(2)]';
PERI2xv = [XYLTpr2(1) XYRTpr2(1) XYRBpr2(1) XYLBpr2(1) XYLTpr2(1)]';
PERI2yv = [XYLTpr2(2) XYRTpr2(2) XYRBpr2(2) XYLBpr2(2) XYLTpr2(2)]';

[xa, ya] = polybool('union', PSD1xv, PSD1yv, PERI1xv, PERI1yv);
[xb, yb] = polybool('intersection', PSD1xv, PSD1yv, PERI1xv, PERI1yv);
[xc, yc] = polybool('xor', PSD1xv, PSD1yv, PERI1xv, PERI1yv);
[xd, yd] = polybool('subtraction', PSD1xv, PSD1yv, PERI1xv, PERI1yv);


subplot(2, 2, 1)
patch(xa, ya, 1, 'FaceColor', 'r')
axis equal, axis off, hold on
plot(PSD1xv, PSD1yv, PERI1xv, PERI1yv, 'Color', 'k')
title('union')

subplot(2, 2, 2)
patch(xb, yb, 1, 'FaceColor', 'r')
axis equal, axis off, hold on
plot(PSD1xv, PSD1yv, PERI1xv, PERI1yv, 'Color', 'k')
title('intersection')

subplot(2, 2, 3)
[f, v] = poly2fv(xc, yc);
patch('Faces', f, 'Vertices', v, 'FaceColor', 'r', ...
  'EdgeColor', 'none')
axis equal, axis off, hold on
plot(PSD1xv, PSD1yv, PERI1xv, PERI1yv, 'Color', 'k')
title('Exclusive Or')

subplot(2, 2, 4)
patch(xd, yd, 1, 'FaceColor', 'r')
axis equal, axis off, hold on
plot(PSD1xv, PSD1yv, PERI1xv, PERI1yv, 'Color', 'k')
title('subtraction')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%




%=================================================%
% Surface Plots of Nonuniformly Sampled Data    
%=================================================%
%{

Surface Plots of Nonuniformly Sampled Data


You can use meshgrid to create a grid of uniformly sampled data points at 
which to evaluate and graph the sinc function. MATLAB then constructs the 
surface plot by connecting neighboring matrix elements to form a mesh of quadrilaterals. 

To produce a surface plot from nonuniformly sampled data, use TriScatteredInterp 
to interpolate the values at uniformly spaced points, and then use mesh and surf in the usual way.
Example ? Displaying Nonuniform Data on a Surface


This example evaluates the sinc function at random points within a specific range 
and then generates uniformly sampled data for display as a surface plot. The process involves these tasks:

Use linspace to generate evenly spaced values over the range of your unevenly sampled data.
Use meshgrid to generate the plotting grid with the output of linspace.
Use TriScatteredInterp to interpolate the irregularly sampled data to the regularly spaced grid returned by meshgrid.
Use a plotting function to display the data.

Generate unevenly sampled data within the range [-8, 8] and use it to evaluate the function:
x = rand(100,1)*8;
y = rand(100,1)*8;
r = sqrt(x.^2 + y.^2) + eps;
z = sin(r)./r;
% z = [zeros(1,40) ones(1,20) zeros(1,40)]'

x = [zeros(1,40) ones(1,20) zeros(1,40)]'
y = [ones(1,40) zeros(1,20) ones(1,40)]'
z = rand(100,1)*8;

Zval = ones(5)*3		
Zval = padarray(Zval,[3 3], 2)
Zval = padarray(Zval,[3 3], 0)
mesh(Zval)
surf(Zval)

x = rand(100,1)*8;
y = rand(100,1)*8;
figure
mesh(x,y,Zval) %interpolated
axis tight; hold on
plot3(x,y,Zval,'.','MarkerSize',15) %nonuniform


The linspace function provides a convenient way to create uniformly spaced data with 
the desired number of elements. The following statements produce vectors over the range 
of the random data with the same resolution as that generated by the -8:.5:8 statement 
in the previous sinc example:

xlin = linspace(min(x),max(x),33);
ylin = linspace(min(y),max(y),33);



Now use these points to generate a uniformly spaced grid:
[X,Y] = meshgrid(xlin,ylin);



The key to this process is to use TriScatteredInterp to interpolate the values of the 
function at the uniformly spaced points, based on the values of the function at the 
original data points (which are random in this example). This statement uses the default 
linear interpolation to generate the new data:
f = TriScatteredInterp(z,x,z);;
Z = f(X,Y);


Plot the interpolated and the nonuniform data to produce:
figure
mesh(Z,X,Y) %interpolated
axis tight; hold on
plot3(x,y,z,'.','MarkerSize',15) %nonuniform

%}
%=================================================%
%=================================================%



%=================================================%
%       COOL STUFF TO TRY     
%=================================================%



%=================================================%
%			DRAW STUFF      
%=================================================%
%{


% Draw a Sphere
 r=20;               % resolution
 a=2;                % size of the sphere
 xc=1;yc=2;zc=3;     % coordinates of the center
 alpha=0.5;          % transparency 
[x,y,z]=sphere(r);
figure(1)
clf;
surf(xc+x*a,yc+y*a,zc+z*a,'facecolor','blue','edgecolor','k','facealpha',alpha)
%colormap('jet')   % colormap can be used if facecolor not used 
xlabel('x');
ylabel('y');
zlabel('z');


% Another sphere (particles)
figure
[x,y,z] = sphere(16);
X = [x(:)*.5 x(:)*.75 x(:)];
Y = [y(:)*.5 y(:)*.75 y(:)];
Z = [z(:)*.5 z(:)*.75 z(:)];
S = repmat([1 .75 .5]*10,numel(x),1);
C = repmat([1 2 3],numel(x),1);
scatter3(X(:),Y(:),Z(:),S(:),C(:),'filled'), view(-60,60)
view(40,35)
%=================================================%


%=================================================%
% DRAW A CUBE
xc=1; yc=1; zc=1;    % coordinated of the center
L=1;                 % cube size (length of an edge)
alpha=0.8;           % transparency (max=1=opaque)

X = [0 0 0 0 0 1; 1 0 1 1 1 1; 1 0 1 1 1 1; 0 0 0 0 0 1];
Y = [0 0 0 0 1 0; 0 1 0 0 1 1; 0 1 1 1 1 1; 0 0 1 1 1 0];
Z = [0 0 1 0 0 0; 0 0 1 0 0 0; 1 1 1 0 1 1; 1 1 1 0 1 1];

% C='blue';                                                   % unicolor
C= [0.1 0.5 0.9 0.9 0.1 0.5];                                 % color/face
% C = [0.1 0.8 1.1 1.1 0.1 0.8 ; 0.2 0.9 1.2 1.2 0.2 0.8 ;  
%      0.3 0.9 1.3 1.3 0.3 0.9 ; 0.4 0.9 1.4 1.4 0.4 0.9 ];   % color scale/face

X = L*(X-0.5) + xc;
Y = L*(Y-0.5) + yc;
Z = L*(Z-0.5) + zc; 

fill3(X,Y,Z,C,'FaceAlpha',alpha);    % draw cube
axis equal;

AZ=-20;         % azimuth
EL=25;          % elevation
view(AZ,EL);    % orientation of the axes 
%=================================================%



%=================================================%
% Draw Brownian Motion

d=4;      % distance
t=200;    % time
figure(1); clf;
x=[]; y=[];
x(1)=50; y(1)=50;

for i=2:t
    pause(0.02)
    r=rand();				% random radius
    dx=d*sin(2*pi*r);		% Ly = D * x:D(circ)
    dy=d*cos(2*pi*r);		% Lx = D * x:D(circ)
    xp=x(i-1)+dx;
    yp=y(i-1)+dy;
    xp=max(0,min(100,xp));
    yp=max(0,min(100,yp));    
    x(i)=xp;
    y(i)=yp;
    hold off;
    plot(x,y,'k');
    hold on;
    plot(x(1),y(1),'o','LineWidth',1,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5)
    hold on;
    plot(x(end),y(end),'o','LineWidth',1,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10);
    xlim([0 100]);
    ylim([0 100]);
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    axis square;
    drawnow;
end

%=================================================%






%=================================================%
%=================================================%
%}
%=================================================%
%=================================================%











