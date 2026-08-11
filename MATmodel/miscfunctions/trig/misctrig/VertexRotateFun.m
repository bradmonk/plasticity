function [] = VertexRotateFun()
%----------------------------------------
%	  FIGURE SETUP
%----------------------------------------
sz = [2.5e-3 2.5e-3 1.6 1.2];
Fh25 = FigSetup(25,sz);


%---------------------------------------------%
% a function of of seven variables that yields the result of
% rotating the point (x,y,z) about the axis ?u,v,w? by the angle ?.
PointRota = @(x,y,z,u,v,w,tta) ...
	([(u*(u*x+v*y+w*z)*(1-cos(tta))+(u^2+v^2+w^2)*...
	x*cos(tta)+sqrt(u^2+v^2+w^2)*(-w*y+v*z)*sin(tta))/... 
	(u^2 + v^2 + w^2);
	(v*(u*x+v*y+w*z)*(1-cos(tta))+(u^2+v^2+w^2)*...
	y*cos(tta)+sqrt(u^2+v^2+w^2)*(w*x-u*z)*sin(tta))/... 
	(u^2+v^2+w^2);
	(w*(u*x+v*y+w*z)*(1-cos(tta))+(u^2+v^2+w^2)*...
	z*cos(tta)+sqrt(u^2+v^2+w^2)*(-v*x+u*y)*sin(tta))/... 
	(u^2+v^2+w^2)]);
% MFxut = Fxut(1,2,3,4,5,6,.8)
%---------------------------------------------%



%---------------------------------------------%
%{
An axis of rotation can be formulated from a line between two points Po(d,e,f) -> Pt(a,b,c)
Or we can define an arbitrary line by a single point and a direction vector Dv(u,v,w). 
If we multiply a transformation and point-rotation matrix a point to rotate Pr(x,y,z)
we obtain a 10-variable function that yields the result of rotating this point Pr(x,y,z)
around the axis defined by Pt(a,b,c)-Po(d,e,f)=Dv(u,v,w) by an angle (theta).

% Point to rotate:
Pr = [x;y;z]

% Line defining orthonormal axis:
Po = [d;e;f];
Pt = [a;b;c];

% Direction vector:
Dv = [u;v;w] = [d-a; e-b; f-c]

% Angle of rotation
Ar = [theta]

%}


LineRota = @(x,y,z,a,b,c,u,v,w,tta) ... 
	([(((a*(v^2+w^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(tta))+(u^2+v^2+w^2)*x*cos(tta)+...
	sqrt((u^2+v^2+w^2))*(-c*v+v*w-w*y+v*z)*sin(tta))/(u^2+v^2+w^2));
	(((b*(u^2+w^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(tta))+(u^2+v^2+w^2)*y*cos(tta)+...
	sqrt((u^2+v^2+w^2))*(c*u-a*w+w*x-u*z)*sin(tta))/(u^2+v^2+w^2));
	(((c*(u^2+v^2)-w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(tta))+(u^2+v^2+w^2)*z*cos(tta)+...
	sqrt((u^2+v^2+w^2))*(-b*u+a*v-v*x+u*y)*sin(tta))/(u^2+v^2+w^2))]);
% MFxaut = Fxaut(5,0,0,1,1,1,1,1,1,.8)
%---------------------------------------------%


%---------------------------
% Pr = [x;y;z]
% Po = [d;e;f];
% Pt = [a;b;c];
% Dv = [u;v;w] = [d-a; e-b; f-c]
% Ar = [theta]
%---------------------------

Pr = [5;1;1];
Pt = [5;4;5];
Po = [.5;.5;.5];
Dv = Pt-Po;


x=Pr(1);y=Pr(2);z=Pr(3);
a=Pt(1);b=Pt(2);c=Pt(3);
u=Dv(1);v=Dv(2);w=Dv(3);

ArMx = linspace(.01,2*pi,25);
for ang = 1:25
	% Ar = rand*2*pi;
	Ar=ArMx(ang);
	RMX1(:,ang) = PointRota(x,y,z,u,v,w,Ar);
end

for ang = 1:25
	Ar=ArMx(ang);
	RMX2(:,ang) = LineRota(x,y,z,a,b,c,u,v,w,Ar);
end

%---------------------------------------------%
figure(Fh25)
% Point Rotatoin (Red)
scatter3(RMX1(1,:),RMX1(2,:),RMX1(3,:),'.r')
hold on
% Line Rotation (Blue)
scatter3(RMX2(1,:),RMX2(2,:),RMX2(3,:),'.b')
hold on

% Pr(o) Pt(d) Dv(x)
scatter3(Pr(1),Pr(2),Pr(3),'o')
hold on
scatter3(Pt(1),Pt(2),Pt(3),'d')
hold on
scatter3(Dv(1),Dv(2),Dv(3),'x')
axis([-15 15 -15 15 -15 15])
hold on
%---------------------------------------------%

%---------------------------------------------%
% Plot Standard ijk Orthonormal Axis
sFilO = [-15 0 0; 0 -15 0; 0 0 -5];
sFilT = [15 0 0; 0 15 0; 0 0 15];
sFilO = cat(2,Po,sFilO);
sFilT = cat(2,Pt,sFilT);
Xop0 = [sFilO(1,:); sFilT(1,:)];
Yop0 = [sFilO(2,:); sFilT(2,:)];
Zop0 = [sFilO(3,:); sFilT(3,:)];
plot3(Xop0,Yop0,Zop0)
% axis vis3d
hold on
%---------------------------------------------%



%---------------------------------------------%
% Plot Line Connections

RMPr = repmat(Pr,1,25);
RMPt = repmat(Pt,1,25);
RMPo = repmat(Po,1,25);
RMDv = repmat(Dv,1,25);

RMx1PrX = [RMPr(1,:); RMX1(1,:)];
RMx1PrY = [RMPr(2,:); RMX1(2,:)];
RMx1PrZ = [RMPr(3,:); RMX1(3,:)];
RMx2PrX = [RMPr(1,:); RMX2(1,:)];
RMx2PrY = [RMPr(2,:); RMX2(2,:)];
RMx2PrZ = [RMPr(3,:); RMX2(3,:)];

RMx1PtX = [RMPt(1,:); RMX1(1,:)];
RMx1PtY = [RMPt(2,:); RMX1(2,:)];
RMx1PtZ = [RMPt(3,:); RMX1(3,:)];
RMx2PtX = [RMPt(1,:); RMX2(1,:)];
RMx2PtY = [RMPt(2,:); RMX2(2,:)];
RMx2PtZ = [RMPt(3,:); RMX2(3,:)];

RMx1PoX = [RMPo(1,:); RMX1(1,:)];
RMx1PoY = [RMPo(2,:); RMX1(2,:)];
RMx1PoZ = [RMPo(3,:); RMX1(3,:)];
RMx2PoX = [RMPo(1,:); RMX2(1,:)];
RMx2PoY = [RMPo(2,:); RMX2(2,:)];
RMx2PoZ = [RMPo(3,:); RMX2(3,:)];

RMx1DvX = [RMDv(1,:); RMX1(1,:)];
RMx1DvY = [RMDv(2,:); RMX1(2,:)];
RMx1DvZ = [RMDv(3,:); RMX1(3,:)];
RMx2DvX = [RMDv(1,:); RMX2(1,:)];
RMx2DvY = [RMDv(2,:); RMX2(2,:)];
RMx2DvZ = [RMDv(3,:); RMX2(3,:)];


% Connect lines from:
doPDv = 0;	doLDv = 0; 
doPPt = 0;	doLPt = 1; 
doPPo = 1;	doLPo = 0; 
doPPr = 0;	doLPr = 0;
%---------
if doPDv
plot3(RMx1DvX,RMx1DvY,RMx1DvZ,'r')
hold on
end
if doLDv
plot3(RMx2DvX,RMx2DvY,RMx2DvZ,'b')
hold on
end
%---------
if doPPt
plot3(RMx1PtX,RMx1PtY,RMx1PtZ,'r')
hold on
end
if doLPt
plot3(RMx2PtX,RMx2PtY,RMx2PtZ,'b')
hold on
end
%---------
if doPPo
plot3(RMx1PoX,RMx1PoY,RMx1PoZ,'r')
hold on
end
if doLPo
plot3(RMx2PoX,RMx2PoY,RMx2PoZ,'b')
hold on
end
%---------
if doPPr
plot3(RMx1PrX,RMx1PrY,RMx1PrZ,'r')
hold on
end
if doLPr
plot3(RMx2PrX,RMx2PrY,RMx2PrZ,'b')
hold on
end
%---------------------------------------------%


%---------------------------------------------%
% RotaFig(0,.5,.01)
FIGROT(Fh25,0,.5,.01)
%---------------------------------------------%

end

