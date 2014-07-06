function [varargout] = FIGROT(Fh,varargin)

Fh = gcf;

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
	do2D = 0;
	pauseT1=1;
	pauseT2=.001;
end


ax = axis;
%----------------------
figure(Fh)
xlabel('X');ylabel('Y');zlabel('Z');
getme = {'CameraPosition', 'CameraTarget', 'CameraViewAngle'};
gotme = get(gca,getme);
set(gca,'CameraPosition',[20, -20, 20],'CameraTarget',[0, 0, 0],'CameraViewAngle',[60]);
%----------------------


magX = 1;
nfram = 100;
%----------------------
figure(Fh)
xax = ax(2);
xxp = linspace((-xax*magX),(xax*magX),nfram*magX);
for fram = 1:nfram*magX
    campos([-xxp(fram),ax(3),ax(6)])
    drawnow
end
%----------------------
figure(Fh)
ax2 = axis;
xax = ax2(2);
xxp = fliplr(linspace((-xax*magX),(xax*magX),nfram*magX));
for fram = 1:nfram*magX
    campos([-xxp(fram),ax2(3),ax2(6)])
    drawnow
end
%----------------------



%----------------------
if do2D
figure(Fh)
camT = camT1;
cTx = 15;
cTy = -15;
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

%{
for i=360:-5:0
    view(i,23.5);     %  Rotate axis tilts by 23.5 degrees
    drawnow
end


set(gca,'CameraPosition',[x, y, z],'CameraTarget',[x, y, z],'CameraViewAngle',[x, y, z])
get(gca,'CameraPosition','CameraTarget',CameraViewAngle')
get(gcf,'CurrentAxes')

camzoom(zoom_factor)
camdolly(dx,dy,dz)
camdolly(dx,dy,dz,'targetmode') 
camdolly(dx,dy,dz,'targetmode','coordsys')
camdolly(axes_handle,...)
surf(peaks)
axis vis3d
t = 0:pi/20:2*pi;
dx = sin(t)./40;
dy = cos(t)./40;
for i = 1:length(t);
    camdolly(dx(i),dy(i),0)
    drawnow
end

camorbit(dtheta,dphi)
camorbit(dtheta,dphi,'coordsys')
camorbit(dtheta,dphi,'coordsys','direction')
camorbit(axes_handle,...)
surf(peaks)
axis vis3d
for i=1:36
   camorbit(10,0,'data',[0 1 0])
   drawnow
end

%}






%----------------------
set(gca,'CameraPosition',gotme{1},'CameraTarget',gotme{2},'CameraViewAngle',gotme{3});
varargout = {campos, camtarget};
%----------------------
end








%{

plot(HSshape,HSline,'LineWidth',3)

set & get
	set(H,'PropertyName',PropertyValue,...)
	set(gca,'PropertyName',PropertyValue,...)
	a = get(h)
	a = get(0)
	a = get(h,'Default')
	h = gca
	get(gcf,'CurrentAxes')
	get(0,'DefaultLineLineWidth')

set(gca,'CameraPosition',[x, y, z],'CameraTarget',[x, y, z],'CameraViewAngle',[x, y, z])
get(gca,'CameraPosition','CameraTarget',CameraViewAngle')
get(gcf,'CurrentAxes')



CameraPosition
    [x, y, z] axes coordinates
Location of the camera. Position from which the camera views the scene. 
Specify the point in axes coordinates.




CameraTarget
    [x, y, z] axes coordinates
Camera aiming point. Specifies the location in the axes that the camera points to. 




CameraUpVector
    [x, y, z] axes coordinates
Camera rotation. Specifies the rotation of the camera around the viewing axis 
defined by the CameraTarget and the CameraPosition properties. Specify CameraUpVector 
as a three-element array containing the x, y, and z components of the vector. 
For example, [0 1 0] specifies the positive y-axis as the up direction. 
The default CameraUpVector is [0 0 1], which defines the positive z-axis as the up direction.




CameraViewAngle
    [0<A<180] scalar greater than 0 and less than or equal to 180 (angle in degrees)
Field of view. Determines the camera field of view. Changing this value affects the 
size of graphics objects displayed in the axes, but does not affect the degree of 
perspective distortion. The greater the angle, the larger the field of view, and 
the smaller objects appear in the scene. 




PlotBoxAspectRatio
    [px py pz]
	Relative scaling of axes plot box. Controls the relative scaling of the plot box 
	in the x, y, and z directions. The plot box is a box enclosing the axes data region 
	as defined by the x-, y-, and z-axis limits.

daspect
	The data aspect ratio determines the relative scaling of the data units 
	along the x-, y-, and z-axes. 
daspect() 
	with no arguments returns the data aspect ratio of the current axes. 
daspect([aspect_ratio]) 
	sets the data aspect ratio in the current axes to the specified value. Specify the 
	aspect ratio as three relative values representing the ratio of the x-, y-, 
	and z-axis scaling (e.g., [1 1 3] means one unit in x is equal in length to one 
	unit in y and three units in z).
daspect([aspect_ratio])
daspect('mode')
daspect('auto')
daspect('manual')
daspect(axes_handle,...)



LineWidth
    line width in points
	Width of axis lines. Specifies the width, in points, of the x-, y-, and z-axis lines. 
	The default line width is 0.5 points 1 point



OuterPosition
	[left bottom width height]
    four-element vector
	Position of axes including labels, title, and a margin. Specifies a rectangle that 
	locates the outer bounds of the axes, including axis labels, the title, and a margin. 
	The vector is as follows: [left bottom width height]




Title
    handle of text object
set(get(gca,'Title'),'Color','r')
	To create a new title, set this property to the handle of the text object you want to use:
set(gca,'Title',text('String','New Title','Color','r'))
	However, it is simpler to use the title command to create or replace an axes title:
title('New Title','Color','r') % Make text color red
title({'This title','has 2 lines'}) % Two line title





XGrid, YGrid, ZGrid
    on | {off}
set(gca,'XGrid','on')



XLabel, YLabel, ZLabel
	set(get(gca,'XLabel'),'String','axis label')
	xlabel('string')
	xlabel({'first line';'second line'})
	ylabel('George''s Popularity','fontsize',12,'fontweight','b')



XMinorGrid, YMinorGrid, ZMinorGrid
    on | {off}

XMinorTick, YMinorTick, ZMinorTick
    on | {off}


XScale, YScale, ZScale
    {linear} | log


XTick, YTick, ZTick
	set(gca,'XTickLabel',{'One';'Two';'Three';'Four'})
	set(gca,'XTickLabel',{'1';'10';'100'})
	set(gca,'XTickLabel','1|10|100')
	set(gca,'XTickLabel',[1;10;100])
	set(gca,'XTickLabel',['1  ';'10 ';'100'])





Zoom in using aspect ratio and limits:
sphere
set(gca,'DataAspectRatio',[1 1 1],...
        'PlotBoxAspectRatio',[1 1 1],'ZLim',[-0.6 0.6])



Zoom in and out using the CameraViewAngle:
sphere
set(gca,'CameraViewAngle',get(gca,'CameraViewAngle')-5)
set(gca,'CameraViewAngle',get(gca,'CameraViewAngle')+5)







Querying the plot box aspect ratio shows that the plot box is square.
pbaspect
ans = 
     1  1  1

It is also interesting to look at the data aspect ratio selected by MATLAB.
daspect
ans = 
     4  4  1

To illustrate the interaction between the plot box and data aspect ratios, 
set the data aspect ratio to [1 1 1] and again query the plot box aspect ratio.
daspect([1 1 1])
pbaspect([1 1 1])





sphere(36); 
h = findobj('Type','surface');
set(h,'FaceLighting','phong',...
      'FaceColor','interp',...
      'EdgeColor',[.4 .4 .4],...
      'BackFaceLighting','lit')
hold on

vert = [
1 1 1;
1 2 1;
2 2 1;
2 1 1;
1 1 2;
1 2 2;
2 2 2;
2 1 2]

fac = [
1 2 3 4;
2 6 7 3;
4 3 7 8;
1 5 8 4;
1 2 6 5;
5 6 7 8]



patch('faces',fac,'vertices',vert,'FaceColor','y');
light('Position',[1 3 2]);
light('Position',[-3 -1 3]);
material shiny
axis vis3d off
hold off




camzoom(zoom_factor)


campos
campos([camera_position])
campos('mode')
campos('auto')
campos('manual')
campos(axes_handle,...)
	surf(peaks)
	axis vis3d off
	for x = -200:5:200
		campos([x,5,10])
		drawnow
	end


camdolly(dx,dy,dz)
camdolly(dx,dy,dz,'targetmode') 
camdolly(dx,dy,dz,'targetmode','coordsys')
camdolly(axes_handle,...)
	surf(peaks)
	axis vis3d
	t = 0:pi/20:2*pi;
	dx = sin(t)./40;
	dy = cos(t)./40;
	for i = 1:length(t);
		camdolly(dx(i),dy(i),0)
		drawnow
	end

camorbit(dtheta,dphi)
camorbit(dtheta,dphi,'coordsys')
camorbit(dtheta,dphi,'coordsys','direction')
camorbit(axes_handle,...)
	surf(peaks)
	axis vis3d
	for i=1:36
	   camorbit(10,0,'data',[0 1 0])
	   drawnow
	end


camtarget
camtarget([camera_target])
camtarget('mode')
camtarget('auto')
camtarget('manual')
camtarget(axes_handle,...)
	surf(peaks); 
	axis vis3d
	xp = linspace(-150,40,50);
	xt = linspace(25,50,50);
	for i=1:50
		 campos([xp(i),25,5]);
		 camtarget([xt(i),30,0])
		 drawnow
	end


camva
camva(view_angle)
camva('mode')
camva('auto')
camva('manual')
camva(axes_handle,...)
	% Set the range checking in the callback statements to keep 
	% the values for the camera view angle in the range greater 
	% than zero and less than 180.
	uicontrol('Style','pushbutton',...
	  'String','Zoom In',...
	  'Position',[20 20 60 20],...
	  'Callback','if camva <= 1;return;else;camva(camva-1);end');
	uicontrol('Style','pushbutton',...
	  'String','Zoom Out',...
	  'Position',[100 20 60 20],...
	  'Callback','if camva >= 179;return;else;camva(camva+1);end');
	% Now create a graph to zoom in and out on:
	surf(peaks);

camup
camup([up_vector])
camup('mode')
camup('auto')
camup('manual')
camup(axes_handle,...)





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

set(gca, 'Box', 'on'),
set(gca,'xticklabel',[]),
set(gca,'yticklabel',[]),
set(gca,'zticklabel',[]),


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
%---
text(XYLTpr1(1),XYLTpr1(2),...
strcat(num2str(XYLTpr1(1)), '\bullet',num2str(XYLTpr1(2)),'\rightarrow'),...
'FontSize',12,'HorizontalAlignment','right');

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







