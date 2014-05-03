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
% %% Full Example (and formatting help)
%
% %------------------------------------------%
% r = 5;
% xc = 10;
% yc = 10;
% zc = 4;
% %------------------------------------------%
% % 2D Circle
% [xp yp] = circus(r,xc,yc);
% figure
% plot(xp, yp);
% axis square;
% %------------------------------------------%
% % 3D Circle
% [xp yp zp] = circus(r,xc,yc,zc);
% figure
% [ph1] = plot3(xp, yp, zp);
% axis square;
% %------------------------------------------%
% c1= [.9 .2 .2]; c11=[.9 .3 .3]; MS1 = 5;
% %------------------------------------------%
% set(ph1,'LineStyle','-','Color',c1,'LineWidth',2);
% hTitle  = title ('MAKE CIRCLES WITH CIRCUS()');
% hXLabel = xlabel('X-AXIS'); hYLabel = ylabel('Y-AXIS'); hZLabel = zlabel('Z-AXIS');
% set(gca,'FontName','Helvetica');
% set([hXLabel, hYLabel, hZLabel],'FontSize',12);
% set( hTitle,'FontSize',14);
% set(gca,'Box','off','TickDir','out','TickLength',[.02 .02], ...
% 'XMinorTick','on','YMinorTick','on','YGrid','on','XGrid','on', ...
% 'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'ZColor',[.3 .3 .3],'LineWidth',1);
% haxes=axis;
% ylim([0 haxes(4)]);
% xlim([0 haxes(2)]);
% %------------------------------------------%
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf, ['outputfigs/CIRCUS.png']);
% %------------------------------------------%
%
% See also PI, CIRCLE, SPLINE, PATCH,  RECTANGLE, PLOT3
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