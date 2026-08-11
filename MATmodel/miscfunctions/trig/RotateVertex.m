function [varargout] = RotateVertex(varargin)


%---------------
if nargin == 13
% LineRota(x,y,z,a,b,c,u,v,w,tta)
Pr = [varargin{1};varargin{2};varargin{3}];
Pt = [varargin{4};varargin{5};varargin{6}];
Po = [varargin{7};varargin{8};varargin{9}];
Dv = Pt-Po;

x=Pr(1);y=Pr(2);z=Pr(3);
a=Pt(1);b=Pt(2);c=Pt(3);
u=Dv(1);v=Dv(2);w=Dv(3);

tta = varargin{13};
doPoint = 0;
doLine = 1;


%---------------
elseif nargin == 7
% PointRota(x,y,z,u,v,w,tta)
x = varargin{1};
y = varargin{2};
z = varargin{3};
u = varargin{4};
v = varargin{5};
w = varargin{6};
tta = varargin{7};
doPoint = 1;
doLine = 0;

%---------------
elseif nargin == 4
% LineRota(x,y,z,a,b,c,u,v,w,tta)
Pr = varargin{1};
Pt = varargin{2};
Po = varargin{3};
Dv = Pt-Po;

x=Pr(1);y=Pr(2);z=Pr(3);
a=Pt(1);b=Pt(2);c=Pt(3);
u=Dv(1);v=Dv(2);w=Dv(3);

tta = varargin{4};
doPoint = 0;
doLine = 1;

%---------------
elseif nargin == 3
% PointRota(x,y,z,u,v,w,tta)
Pr = varargin{1};
Dv = varargin{2};

x=Pr(1);y=Pr(2);z=Pr(3);
u=Dv(1);v=Dv(2);w=Dv(3);

tta = varargin{3};
doPoint = 1;
doLine = 0;

%---------------
else

Pr = [5;1;1];
Pt = [5;4;5];
Po = [.5;.5;.5];
Dv = Pt-Po;

x=Pr(1);y=Pr(2);z=Pr(3);
a=Pt(1);b=Pt(2);c=Pt(3);
u=Dv(1);v=Dv(2);w=Dv(3);

tta = rand*2*pi;

doPoint = 1;
doLine = 1;
end
%---------------



if doPoint
PointRota = @(x,y,z,u,v,w,tta) ...
	([(u*(u*x+v*y+w*z)*(1-cos(tta))+(u^2+v^2+w^2)*x*cos(tta)+...
		sqrt(u^2+v^2+w^2)*(-w*y+v*z)*sin(tta))/(u^2 + v^2 + w^2);
	
	(v*(u*x+v*y+w*z)*(1-cos(tta))+(u^2+v^2+w^2)*y*cos(tta)+...
		sqrt(u^2+v^2+w^2)*(w*x-u*z)*sin(tta))/(u^2+v^2+w^2);
	
	(w*(u*x+v*y+w*z)*(1-cos(tta))+(u^2+v^2+w^2)*z*cos(tta)+...
		sqrt(u^2+v^2+w^2)*(-v*x+u*y)*sin(tta))/(u^2+v^2+w^2)]);

RotDot = PointRota(x,y,z,u,v,w,tta);
end


if doLine
LineRota = @(x,y,z,a,b,c,u,v,w,tta) ... 
	([(((a*(v^2+w^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(tta))+(u^2+v^2+w^2)*x*cos(tta)+...
		sqrt((u^2+v^2+w^2))*(-c*v+b*w-w*y+v*z)*sin(tta))/(u^2+v^2+w^2));
	
	(((b*(u^2+w^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(tta))+(u^2+v^2+w^2)*y*cos(tta)+...
		sqrt((u^2+v^2+w^2))*(c*u-a*w+w*x-u*z)*sin(tta))/(u^2+v^2+w^2));
	
	(((c*(u^2+v^2)-w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(tta))+(u^2+v^2+w^2)*z*cos(tta)+...
		sqrt((u^2+v^2+w^2))*(-b*u+a*v-v*x+u*y)*sin(tta))/(u^2+v^2+w^2))]);
	
RotDot = LineRota(x,y,z,a,b,c,u,v,w,tta);
end



varargout = {RotDot};






% x1=Po(1);y1=Po(2);z1=Po(3);
% x2=Pt(1);y2=Pt(2);z2=Pt(3);
% x0=Pr(1);y0=Pr(2);z0=Pr(3);
% dx=Dv(1);dy=Dv(2);dz=Dv(3);
% 
% % given a line that goes through Po(x1,y1,z1) (aka x1) and Pt(x2,y2,z2) (aka x2)
% % t gives a value that is proportional to the vector length |PoPt|
% % v gives the location of this value
% t1=-1;	
% v1 = [x1+(x2-x1)*t1;y1+(y2-y1)*t1;z1+(z2-z1)*t1];
% 
% % The squared distance between a point on the line with parameter t (aka v1)
% % and a point Pr(x0,y0,z0) (aka x0) is therefore
% d1sq = [(x1-x0)+(x2-x1)*t1]^2 + [(y1-y0)+(y2-y1)*t1]^2 + [(z1-z0)+(z2-z1)*t1]^2;
% d1 = sqrt(d1sq);
% 
% % To minimize the distance, set d(d^2)/dt=0 and solve for t  (via dot product) to obtain
% t0 = LV((dot((Po-Pr),(Pt-Po)))) / LV((Pt-Po))^2;
% % And again, the location of that point
% v0 = [x1+(x2-x1)*t0;y1+(y2-y1)*t0;z1+(z2-z1)*t0];
% 
% % The minimum distance can then be found by using the min t to obtain:
% d1sq_min = ((LV(Po-Pr)^2 * LV(Pt-Po)^2) - (dot((Po-Pr),(Pt-Po))^2)) / LV(Pt-Po)^2;
% d1_min = sqrt(d1sq_min);
% 
% % Here is the quickest way to find the min distance from a point 
% d2 = LV(cross((Pt-Po),(Po-Pr))) / LV(Pt-Po);
% d2 = LV(cross((Pr-Po),(Pr-Pt))) / LV(Pt-Po);
% 






% %---------------------------
% % Pr = [x;y;z]
% % Po = [d;e;f];
% % Pt = [a;b;c];
% % Dv = [u;v;w] = [d-a; e-b; f-c]
% % Ar = [theta]
% %---------------------------
% 
% Pr = [5;1;1];
% Pt = [5;4;5];
% Po = [.5;.5;.5];
% Dv = Pt-Po;
% 
% 
% x=Pr(1);y=Pr(2);z=Pr(3);
% a=Pt(1);b=Pt(2);c=Pt(3);
% u=Dv(1);v=Dv(2);w=Dv(3);
% 
% ArMx = linspace(.01,2*pi,25);
% for ang = 1:25
% 	% Ar = rand*2*pi;
% 	Ar=ArMx(ang);
% 	RMX1(:,ang) = PointRota(x,y,z,u,v,w,Ar);
% end
% 
% for ang = 1:25
% 	Ar=ArMx(ang);
% 	RMX2(:,ang) = LineRota(x,y,z,a,b,c,u,v,w,Ar);
% end



% %---------------------------------------------%
% % a function of of seven variables that yields the result of
% % rotating the point (x,y,z) about the axis ?u,v,w? by the angle ?.
% PointRota = @(x,y,z,u,v,w,tta) ...
% 	([(u*(u*x+v*y+w*z)*(1-cos(tta))+(u^2+v^2+w^2)*...
% 	x*cos(tta)+sqrt(u^2+v^2+w^2)*(-w*y+v*z)*sin(tta))/... 
% 	(u^2 + v^2 + w^2);
% 	(v*(u*x+v*y+w*z)*(1-cos(tta))+(u^2+v^2+w^2)*...
% 	y*cos(tta)+sqrt(u^2+v^2+w^2)*(w*x-u*z)*sin(tta))/... 
% 	(u^2+v^2+w^2);
% 	(w*(u*x+v*y+w*z)*(1-cos(tta))+(u^2+v^2+w^2)*...
% 	z*cos(tta)+sqrt(u^2+v^2+w^2)*(-v*x+u*y)*sin(tta))/... 
% 	(u^2+v^2+w^2)]);
% % MFxut = Fxut(1,2,3,4,5,6,.8)
% %---------------------------------------------%



% %---------------------------------------------%
% %{
% An axis of rotation can be formulated from a line between two points Po(d,e,f) -> Pt(a,b,c)
% Or we can define an arbitrary line by a single point and a direction vector Dv(u,v,w). 
% If we multiply a transformation and point-rotation matrix a point to rotate Pr(x,y,z)
% we obtain a 10-variable function that yields the result of rotating this point Pr(x,y,z)
% around the axis defined by Pt(a,b,c)-Po(d,e,f)=Dv(u,v,w) by an angle (theta).
% 
% % Point to rotate:
% Pr = [x;y;z]
% 
% % Line defining orthonormal axis:
% Po = [d;e;f];
% Pt = [a;b;c];
% 
% % Direction vector:
% Dv = [u;v;w] = [d-a; e-b; f-c]
% 
% % Angle of rotation
% Ar = [theta]
% 
% %}
% 
% 
% LineRota = @(x,y,z,a,b,c,u,v,w,tta) ... 
% 	([(((a*(v^2+w^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(tta))+(u^2+v^2+w^2)*x*cos(tta)+...
% 	sqrt((u^2+v^2+w^2))*(-c*v+v*w-w*y+v*z)*sin(tta))/(u^2+v^2+w^2));
% 	(((b*(u^2+w^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(tta))+(u^2+v^2+w^2)*y*cos(tta)+...
% 	sqrt((u^2+v^2+w^2))*(c*u-a*w+w*x-u*z)*sin(tta))/(u^2+v^2+w^2));
% 	(((c*(u^2+v^2)-w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(tta))+(u^2+v^2+w^2)*z*cos(tta)+...
% 	sqrt((u^2+v^2+w^2))*(-b*u+a*v-v*x+u*y)*sin(tta))/(u^2+v^2+w^2))]);
% % MFxaut = Fxaut(5,0,0,1,1,1,1,1,1,.8)
% %---------------------------------------------%


end

