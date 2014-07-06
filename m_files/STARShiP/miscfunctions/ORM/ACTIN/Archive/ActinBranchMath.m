function [] = ActinBranchMath(varargin)
clc; close all; clear all;


%---------------------------------------------%
%				ANGLE SETUP
%---------------------------------------------%
% unitsratio('rad', 'deg')
rad2dg = 180/pi;
dg2rad = 1/(180/pi);
Si = .00001;
Vn = [.66446;.66446;.342];
Veq = [.58;.58;.58];
Ovec = [1.0 28.7 56.4 84.1 111.8 139.5 167.2 194.8 222.5 250.2 277.9 305.6 333.3];
Ovec = Ovec * dg2rad;
%---------------------------------------------%



%---------------------------------------------%
%			STARTING VECTOR LOCATIONS
%---------------------------------------------%
%{
nSV = 1;

PoX = 0;
PoY = 0;
PoZ = 0;
PtX = .66446;
PtY = .66446;
PtZ = .342;
PrX = 0;
PrY = 0;
PrZ = 1;

Po = [PoX;PoY;PoZ];			
Pt = [PtX;PtY;PtZ];		
Pr = [PrX;PrY;PrZ];		
Dv = Pt-Po;
%}


Po = [0;0;0];			
Pt = [.66446;.66446;.342];		
Pr = [0;0;1];		
Dv = Pt-Po;



Po = [.5;.3;1.5];
Pt = [-3;2;1];
Pr = Pt+(Veq.*sign(Pt+Si));
Dv = Pt-Po;

%{
for fID = 1:5
act(fID).n = 270;
act(fID).ax = 70;
act(fID).ay = 20;
act(fID).az = 70;
act(fID).ox = 0;
act(fID).oy = 0;
act(fID).oz = 0;
act(fID).tx = 0;
act(fID).ty = 0;
act(fID).tz = 0;
act(fID).or = 0;
act(fID).ov = Ovec;
act(fID).r = act(fID).n * GaSize;
act(fID).Oxyz = {act(fID).ox, act(fID).oy, act(fID).oz};
act(fID).Txyz = {act(fID).tx, act(fID).ty, act(fID).tz};
act(fID).Amx = {zeros(2,3)};	% Angle Mx			{?,?,?; ?,?,?} 
act(fID).Rmx = {zeros(4,3)};	% Euler Rotation Mx	{x0,y0,z0; x1,y1,z1; x2,y2,z2; X3,Y3,Z3} 
end

fN = numel(actin); % gives number of filaments
fID = 1;
act(fID).n(:)
act(fID).Amx{:}
act(fID).Amx{:}(1,1)

for fID = 1:fN
	act(fID).Amx
end
%}
%---------------------------------------------%



%---------------------------------------------%
%		 VECTOR INTERACTION MATH
%---------------------------------------------%
%{

t1 : value proportional to the vector length |PoPt|
v1 : coordinate location of t1
d1 : distance between the point v1 and Pr
t0 : vector magnitude to minimize distance between Pr and line PoPt
v0 : line location of vector magnitude point
d2 : distance between v0 and Pr

%}


t1=1; 
v1 = [Po(1)+(Pt(1)-Po(1))*t1;Po(2)+(Pt(2)-Po(2))*t1;Po(3)+(Pt(3)-Po(3))*t1];
d1 = sqrt(((Po(1)-Pr(1))+(Pt(1)-Po(1))*t1)^2 + ((Po(2)-Pr(2))+(Pt(2)-Po(2))*t1)^2 + ((Po(3)-Pr(3))+(Pt(3)-Po(3))*t1)^2);
t0 = LV((dot((Po-Pr),(Pt-Po)))) / LV((Pt-Po))^2;
v0 = [Po(1)+(Pt(1)-Po(1))*t0;Po(2)+(Pt(2)-Po(2))*t0;Po(3)+(Pt(3)-Po(3))*t0];
d0 = LV(cross((Pr-Po),(Pr-Pt))) / LV(Pt-Po);
%---------------------------------------------%
%{
t2=1+Vn; 
v2 = [Po(1)+(Pt(1)-Po(1))*t2(1);...
	  Po(2)+(Pt(2)-Po(2))*t2(2);...
	  Po(3)+(Pt(3)-Po(3))*t2(3)];
d2 = sqrt(((Po(1)-Pr(1))+(Pt(1)-Po(1))*t2(1))^2 + ...
		  ((Po(2)-Pr(2))+(Pt(2)-Po(2))*t2(2))^2 + ...
		  ((Po(3)-Pr(3))+(Pt(3)-Po(3))*t2(3))^2);
%}


%====================================================%
Fdots = 14;
%---------------------------------------------%
%	  RotateVertex AROUND LINE LOOP
%---------------------------------------------%
for ang = 1:Fdots
	Rmono = ceil(200 * rand);	% Random monomer along segment
	Rmono13 = mod(ang,13)+1;	% Get monomer repeat among the 13 rotational axis angles
	tta = Ovec(Rmono13);		% Rotational angle of new branch
	RMX0(:,ang) = RotateVertex(Pr(1),Pr(2),Pr(3),Pt(1),Pt(2),Pt(3),Po(1),Po(2),Po(3),...
		Dv(1),Dv(2),Dv(3),tta);	
end
%====================================================%
%{
%====================================================%
% Fdots = 35;
ArMx = linspace(.01,2*pi,Fdots);
%---------------------------------------------%
%	  RotateVertex AROUND POINT LOOP
%---------------------------------------------%
for ang = 1:Fdots
	tta=ArMx(ang);
	RMX1(:,ang) = RotateVertex(Pr(1),Pr(2),Pr(3),Dv(1),Dv(2),Dv(3),tta);
end
%---------------------------------------------%
%	  RotateVertex AROUND LINE LOOP
%---------------------------------------------%
for ang = 1:Fdots
	
	Rmono = ceil(200 * rand);	% Random monomer along segment
	Rmono13 = mod(ang,13)+1;	% Get monomer repeat among the 13 rotational axis angles
	tta = Ovec(Rmono13);		% Rotational angle of new branch
	RMX2(:,ang) = RotateVertex(Pr(1),Pr(2),Pr(3),Pt(1),Pt(2),Pt(3),Po(1),Po(2),Po(3),...
		Dv(1),Dv(2),Dv(3),tta);
	
	% Update Rotation Point
	Pr = Pr+v00;
	v0 = v0 + v00;
	
	% Update Line Origin Point
	PoMX(:,ang)  = PoV00;
	PoV00 = PoV00+v00;
	
end
%====================================================%
%}
%====================================================%



%---------------------------------------------%
%				FIGURE SETUP
%---------------------------------------------%
sz = [2.5e-3 2.5e-3 1.6 1.2];
Fh25 = FigSetup(25,sz);
%---------------------------------------------%
%		PLOT ORIGIN->TIP LINE
%---------------------------------------------%
figure(Fh25);
[XMx YMx ZMx] = plot3prep({Pt v1 Pr},{Po Po Pt});
hA1 = plot3(XMx,YMx,ZMx,'LineWidth',2,'MarkerSize',10);
grid on;
xlabel('X');ylabel('Y');zlabel('Z');
axis([-15 15 -15 15 -5 10])
axis vis3d
hold on;
addaxis;

%---------------------------------------------%


%{
keyboard

set(gca,'CameraTarget',Po);
zoom(4)

Raxes = [0 0 1];
rotate(hA1(3),Raxes,20,Po); figure(Fh25);
Raxes = [1 0 0];
rotate(hA1(3),Raxes,70,Po); figure(Fh25);
%}




%---------------------------------------------%
% DOT: Pr ROTATED AROUND LINE (BLUE)
%---------------------------------------------%
% figure(Fh25)
scatter3(RMX0(1,:),RMX0(2,:),RMX0(3,:),'.b');
hold on
%---------------------------------------------%


%---------------------------------------------%
%	LINE: Pr-Pt Circle to Point
%---------------------------------------------%
PtMx = repmat(Pt,1,Fdots);
PtPrMxX = [PtMx(1,:); RMX0(1,:)];
PtPrMxY = [PtMx(2,:); RMX0(2,:)];
PtPrMxZ = [PtMx(3,:); RMX0(3,:)];

hA2 = plot3(PtPrMxX,PtPrMxY,PtPrMxZ,'b');
FMx.Fdots = {1.5};
set(hA2,{'LineWidth'},repmat(FMx.Fdots,1,Fdots)');
hold on;
%---------------------------------------------%


%---------------------------------------------%
% DOT: {Po Po Pt Dv Pr v0 v1} {'d';'d';'o';'x';'s';'v';'v'};
%---------------------------------------------%
[XMx YMx ZMx] = scatter3prep({Po Pt Pr Dv v0 v1});
hA3 = scatter3(XMx,YMx,ZMx,'.','Tag','TMxH');
%----------------
hA3c = get(hA3,'Children');
mrkr = {'o';'d';'x';'s';'v';'v'};
set(hA3c,{'Marker'},mrkr);
hold on
colr = {[0 0 0]; [1 0 0]; [0 0 0]; [0 1 0]; [0 .8 .8]; [1 .5 .1]};
set(hA3c,{'MarkerFaceColor'},colr);
hold on
set(hA3c,'LineWidth',1,'MarkerSize',11);
% legend(hAc99,fliplr({'Po','Pt','Pr','Dv','v0','v1'}))
%---------------------------------------------%






%---------------------------------------------%
% PLOTSUGAR: ANOTATE ANIMATE
%---------------------------------------------%

%---------------------------------------------%
simlegend({'Po','Pt','Dv','Pr','v0','v1'},hA3);
%---------------------------------------------%

%---------------------------------------------%
% anotext({Po, 'Po'},{Pt, 'Pt'},{Pr, 'Pr'},{v0, 'v0'},{v1, 'v1'});
%---------------------------------------------%

keyboard

%---------------------------------------------%
FIGROT(gca,0,.5,.01);
%---------------------------------------------%

%---------------------------------------------%
set(gca,'CameraTarget',Po);
for cz = 1:.005:1.05
zoom(cz);
drawnow;
pause(.01);
end
%---------------------------------------------%









%======================================================================%
end
%======================================================================%














%---------------------------------------------%
%	ADD STANDARD IJK ORTHONORMAL AXIS
%---------------------------------------------%
function [varargout] = addaxis(varargin)

if nargin == 1
	axvals = varargin{1};
XAxO = axvals{1}; XAxT = axvals{2};
YAxO = axvals{3}; YAxT = axvals{4};
ZAxO = axvals{5}; ZAxT = axvals{6};
else
XAxO = [-15 0 0]; XAxT = [15 0 0];
YAxO = [0 -15 0]; YAxT = [0 15 0];
ZAxO = [0 0 -5]; ZAxT = [0 0 10];
end
[XAx YAx ZAx] = plot3prep({XAxT YAxT ZAxT},{XAxO YAxO ZAxO});
%-----------------
hAax = plot3(XAx,YAx,ZAx,'k','LineWidth',1,'MarkerSize',10);
set(hAax,{'LineWidth'},{1,1,1}'); hold on
varargout={hAax};
end
%---------------------------------------------%

%==================================================%
%	3D VECTOR LENGTH FUNCTION
%--------------------------------------------------%
function [L] = LV(V)

L = sqrt(sum(V.^2));

end
%==================================================%


%==================================================%
%	FIGROT FUNCTION
%--------------------------------------------------%
% function [varargout] = FIGROT(Fh,varargin)
%--------------------------------------------------%


%==================================================%
%	plot3prep & scatter3prep FUNCTIONS
%--------------------------------------------------%
% function [] = plot3prep()
% function [TMx] = scatter3prep(varargin)
%--------------------------------------------------%







%==================================================%
%				v1 t1 d1 Math
%--------------------------------------------------%
%{
x1=Po(1);y1=Po(2);z1=Po(3);
x2=Pt(1);y2=Pt(2);z2=Pt(3);
x0=Pr(1);y0=Pr(2);z0=Pr(3);
dx=Dv(1);dy=Dv(2);dz=Dv(3);

% given a line that goes through Po(x1,y1,z1) (aka x1) and Pt(x2,y2,z2) (aka x2)
% t gives a value that is proportional to the vector length |PoPt|
% v gives the location of this value
t1=15;	
v1 = [x1+(x2-x1)*t1;y1+(y2-y1)*t1;z1+(z2-z1)*t1];

% The squared distance between a point on the line with parameter t (aka v1)
% and a point Pr(x0,y0,z0) (aka x0) is therefore
d1sq = [(x1-x0)+(x2-x1)*t1]^2 + [(y1-y0)+(y2-y1)*t1]^2 + [(z1-z0)+(z2-z1)*t1]^2;
d1 = sqrt(d1sq);

% To minimize the distance, set d(d^2)/dt=0 and solve for t  (via dot product) to obtain
t0 = LV((dot((Po-Pr),(Pt-Po)))) / LV((Pt-Po))^2;
% And again, the location of that point
v0 = [x1+(x2-x1)*t0;y1+(y2-y1)*t0;z1+(z2-z1)*t0];

% The minimum distance can then be found by using the min t to obtain:
d1sq_min = ((LV(Po-Pr)^2 * LV(Pt-Po)^2) - (dot((Po-Pr),(Pt-Po))^2)) / LV(Pt-Po)^2;
d1_min = sqrt(d1sq_min);

% Here is the quickest way to find the min distance from a point 
d2 = LV(cross((Pt-Po),(Po-Pr))) / LV(Pt-Po);
d2 = LV(cross((Pr-Po),(Pr-Pt))) / LV(Pt-Po);




% t1 : value proportional to the vector length |PoPt|
t1=20; 
% v1 : coordinate location of t1
v1 = [Po(1)+(Pt(1)-Po(1))*t1;Po(2)+(Pt(2)-Po(2))*t1;Po(3)+(Pt(3)-Po(3))*t1];
 
% d1 : distance between the point v1 and Pr
d1sq = [(Po(1)-Pr(1))+(Pt(1)-Po(1))*t1]^2 + [(Po(2)-Pr(2))+(Pt(2)-Po(2))*t1]^2 + [(Po(3)-Pr(3))+(Pt(3)-Po(3))*t1]^2;
d1 = sqrt(d1sq);
 
% t0 : vector magnitude to minimize distance between Pr and line PoPr
t0 = LV((dot((Po-Pr),(Pt-Po)))) / LV((Pt-Po))^2;
% v0 : line location of vector magnitude point
v0 = [Po(1)+(Pt(1)-Po(1))*t0;Po(2)+(Pt(2)-Po(2))*t0;Po(3)+(Pt(3)-Po(3))*t0];
  
% d2 : distance between v0 and Pr
d2 = LV(cross((Pt-Po),(Po-Pr))) / LV(Pt-Po);
d2 = LV(cross((Pr-Po),(Pr-Pt))) / LV(Pt-Po);





% Po = Po00;
% Pt = Pt00;
% Pr = Pr00;
% Dv = Dv00;
% t0 = t00;
% v0 = v00;
% d2 = d200;
% 
% 
% t1=20; 
% v1 = [Po(1)+(Pt(1)-Po(1))*t1;Po(2)+(Pt(2)-Po(2))*t1;Po(3)+(Pt(3)-Po(3))*t1];
% %---------------------------------------------%
% 
% %---------------------------------------------%
% % [XMx YMx ZMx] = plot3prep({T1 T2 T3 T4},{O1 O2 O3 O4});
% [XMx YMx ZMx] = plot3prep({Pt v1},{Po Po});
% %--------------
% % Plot
% phAx = plot3(XMx,YMx,ZMx,'LineWidth',2,'MarkerSize',10);
% % set(phAx,{'LineWidth'},{2,2,2,1,1}')
% hold on
% %---------------------------------------------%
%}
%--------------------------------------------------%

%==================================================%
%			MAIN VERTEX ROTATION FUNCTION
%--------------------------------------------------%
%{
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
%}
%--------------------------------------------------%
















%==================================================%
%				UNUSED FUNCTIONS
%--------------------------------------------------%
%{

%==================================================%
%				FIGURE SETUP FUNCTION
%--------------------------------------------------%
function [varargout] = FigSetup(varargin)

scsz = get(0,'ScreenSize');


if nargin == 1 
	Fnum=varargin{1};
	pos = scsz./[2.5e-3 2.5e-3 1.5 2];
elseif nargin == 2
	Fnum=varargin{1};
	pos=scsz./varargin{2};
else
	Fnum=1;
	pos = scsz./[2.5e-3 2.5e-3 1.5 2];
end

Fh = figure(Fnum);
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])
varargout = {Fh};

end



%==================================================%
%				ROTATION FUNCTION
%--------------------------------------------------%
function [fBPsy] = rotfun(xyz, phi, tta, psy)

fDPhi = ([cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1] * xyz);
fCTta = ([1 0 0; 0 cos(tta) sin(tta); 0 -sin(tta) cos(tta)] * fDPhi);
fBPsy = ([cos(psy) sin(psy) 0; -sin(psy) cos(psy) 0; 0 0 1] * fCTta);

end



%==================================================%
%				FIGURE SETUP FUNCTION
%--------------------------------------------------%
function [varargout] = FigAnime(OTx, OTy, OTz, varargin)


%----------------------------------------
%		FIGURE SETUP
%----------------------------------------
sz = [2.5e-3 2.5e-3 1.6 1.2];
Fh1 = FigSetup(1,sz);
%------------
figure(Fh1)
plot3(OTx,OTy,OTz)
axis([-15 15 -15 15 -5 25])
xlabel('X');ylabel('Y');zlabel('Z');
grid on;


% figure(Fh1)
camP0 = campos;
camT0 = camtarget;

nfram = 50;
cxp = linspace(-10,25,nfram);
for cf = 1:nfram
    campos([cxp(cf),-25,25])
    drawnow
	pause(.02)
end


% figure(Fh1)
camP1 = campos;
camT1 = camtarget;



% camT = camT1+5;
camT = camT1;
cTx = 5;
cTy = -5;
cTz = 10;
nfram = 30;
cxt = fliplr(linspace(camT(1)-cTx,camT(1),nfram));
cyt = fliplr(linspace(camT(2)-cTy,camT(2),nfram));
czt = fliplr(linspace(camT(3)-cTz,camT(3),nfram));
cxt = [cxt fliplr(cxt)];
cyt = [cyt fliplr(cyt)];
czt = [czt fliplr(czt)];

for cf = 1:nfram*2
	camtarget([cxt(cf),cyt(cf),czt(cf)])
    drawnow
	pause(.01)
end


% figure(Fh1)
% campos(camP1)
% camtarget(camT1)

% [cpxyz] = campos;
% nfram = 50;
% cxp = linspace(0,2,nfram);
% for cf = 1:nfram
%     campos([cxp(cf),-25,25])
%     drawnow
% 	pause(.03)
% end


figure(Fh1)
view(0,0)
pause(1)
view(0,90)
pause(1)
view(90,0)
pause(1)
view(18,28)

varargout = {campos, camtarget};

end



%==================================================%
%				ROTATION FUNCTION
%--------------------------------------------------%
function [Rxyz varargout] = ARxFun(xyz, phi, tta, psy)

%---------------------------
Dphi = ([cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1] * xyz);
Ctta = ([1 0 0; 0 cos(tta) sin(tta); 0 -sin(tta) cos(tta)] * Dphi);
Bpsy = ([cos(psy) sin(psy) 0; -sin(psy) cos(psy) 0; 0 0 1] * Ctta);
%---
ERx = Bpsy;
%---------------------------


%---------------------------
Q = [
(cos(psy)*cos(tta)*cos(phi) - sin(psy)*sin(phi)),...
(sin(psy)*cos(tta)*cos(phi) + cos(psy)*sin(phi)), (-sin(tta)*cos(phi));
(-cos(psy)*cos(tta)*sin(phi) - sin(psy)*cos(tta)),...
(-sin(psy)*cos(tta)*sin(phi) + cos(psy)*cos(phi)), (sin(tta)*sin(phi));
(cos(psy)*sin(tta)), (sin(psy)*sin(tta)), (cos(tta))];
%---
QRx = Q*xyz;
%---------------------------


%---------------------------
a11	= cos(psy)*cos(phi) - cos(tta)*sin(phi)*sin(psy);
a12	= cos(psy)*sin(phi) + cos(tta)*cos(phi)*sin(psy);
a13	= sin(psy)*sin(tta);
a21	= -sin(psy)*cos(phi) - cos(tta)*sin(phi)*cos(psy);
a22	= -sin(psy)*sin(phi) + cos(tta)*cos(phi)*cos(psy);
a23	= cos(psy)*sin(tta);
a31	= sin(tta)*sin(phi);
a32	= -sin(tta)*cos(phi);
a33	= cos(tta);
AMx  = [a11 a12 a13;
		a21 a22 a23; 
		a31 a32 a33];
%---
ARx = AMx * xyz;
%---------------------------


%---------------------------
Rxyz = ERx;
Rxyz = QRx;
Rxyz = ARx;
%---------------------------
varargout = {ARx, QRx, ERx};

end


%==================================================%
%				ROTATION FUNCTION
%--------------------------------------------------%
function [t1 t2 t3] = EulerAng(A,rotSet,space)
%INPUT
% A 3x3 Direction Cosine Matrix
% rotSet 1x3 (or 3x1) array of rotation axes to use
% space logical, true for Space rotation, false for Body
%OUTPUT
% t1,t2,t3 Euler angles (in radians)
% assumes Body rotations unless told otherwise

% n = length(unique(rotSet));
n = numel(rotSet);


if ~exist('space','var')
	space = false; 
end


ax2neginds = [8 ,3 ,4]; %negative elements for 2 axis rotations
i = rotSet(1); j = rotSet(2);

switch n
	case 3
		A = A .* [1 -1 1; 1 1 -1; -1 1 1];
		if space, A = A.';end
		k = rotSet (3);
		c2 = sqrt(A(i,i)^2 + A(i,j)^2);
		t1 = atan2(A(j,k) / c2,A(k,k) / c2);
		t2 = atan2(A(i,k),c2);
		t3 = atan2(A(i,j) / c2,A(i,i) / c2);
	case 2
		A(ax2neginds(j)) = -A(ax2neginds(j));
		if space, A = A.';end
		p = 6 - (i+j); %element missing from rotSet
		s2 = sqrt(A(i,p)^2 + A(i,j)^2);
		t1 = atan2(A(j,i) / s2,A(p,i) / s2);
		t2 = atan2(s2,A(i,i));
		t3 = atan2(A(i,j) / s2,A(i,p) / s2);
end

end



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

%}

%==================================================%
%					NOTES
%--------------------------------------------------%
%{

% F1
%------------------
Oxyz1 = Origin;
Txyz1 = Tip;
% Oxyz1 = rotfun(Origin, phi, tta, psy);
% Txyz1 = rotfun(Tip, phi, tta, psy);
%------------------


% F2
%------------------
Oxyz2 = rotfun(Oxyz1, phi, tta, psy);
Txyz2 = rotfun(Txyz1, phi, tta, psy);
%---
% OTcor = (Oxyz2 - Oxyz1);
% Oxyz2 = Oxyz2 - OTcor;
% Txyz2 = Txyz2 - OTcor;
%------------------


% F3
%------------------
% Oxyz3 = rotfun(Txyz2, phi, tta, psy);
% Txyz3 = rotfun(Txyz2+1, phi, tta, psy);
% %---
% OTcor = (Oxyz3- Oxyz2);
% Oxyz3 = Oxyz3 - OTcor;
% Txyz3 = Txyz3 - OTcor;
%------------------
[Oxyz3 OARx OQ OQRx] = ARxFun(Oxyz2, phi, tta, psy);
[Txyz3 TARx TQ TQRx] = ARxFun(Txyz2, phi, tta, psy);
%------------------


% F4
%------------------
[Oxyz4 OARx OQ OQRx] = ARxFun(Oxyz3, phi, tta, psy);
[Txyz4 TARx TQ TQRx] = ARxFun(Txyz3, phi, tta, psy);
%------------------
%---
% R = Length;
% Tx = R * sind(tta) * cosd(phi) + Ox;
% Ty = R * sind(tta) * sind(phi) + Ox;
% Tz = R * cosd(tta) + Oz;
%------------------




% fDPhi = @(xyz,phi) ([cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1] * xyz);
% fCTta = @(xyz,tta) ([1 0 0; 0 cos(tta) sin(tta); 0 -sin(tta) cos(tta)] * xyz);
% fBPsy = @(xyz,psy) ([cos(psy) sin(psy) 0; -sin(psy) cos(psy) 0; 0 0 1] * xyz);

Tx = R * sind(tta) * cosd(phi) + Ox;
Ty = R * sind(tta) * sind(phi) + Ox;
Tz = R * cosd(tta) + Oz;

theta = (Ovec-180)*pi/180;
r = 2*ones(size(theta));
[u,v] = pol2cart(theta,r);
feather(u,v);

OTx = {OT(:).x}
OTy = {OT(:).y}
OTz = {OT(:).z}

[OTx{1}(:)]
[OTx{:}]

OTx = [OT(:).x]

% for rlen = 1:numel(xp(:,1))
% flen(rlen) = sqrt(xp(rlen,2)^2 + yp(rlen,2)^2 + zp(rlen,2)^2);
% end
% disp(flen)


view([25 25])
campos([25,-25,25])

% for rotfig = 1:10
% 	figure(Fh1)
% 	view(azel+rot)
% 	rot = rot+rot1;
% 	pause(.3)
% end

% figure(Fh1)
% axis vis3d
% for i=1:36
%    % camorbit(10,0,'camera','z')
%    camorbit(10,0,'camera','z')
%    % pause(.2)
%    drawnow
% end


a11	= cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi)
a12	= cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi)
a13	= sin(psi)*sin(theta)
a21	= -sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi)
a22	= -sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi)
a23	= cos(psi)*sin(theta)
a31	= sin(theta)*sin(phi)
a32	= -sin(theta)*cos(phi)
a33	= cos(theta)
ARx  = [a11 a12 a13; 
		a21 a22 a23; 
		a31 a32 a33];

RxA = [1 0 0; 0 cos(a) sin(a); 0 -sin(a) cos(a)];
RyB = [cos(b) 0 -sin(b); 0 1 0; sin(b) 0 cos(b)];
RzG = [cos(g) sin(g) 0; -sin(g) cos(g) 0; 0 0 1];


RdPhi = [1			0			0;
		 0			cos(Phi)	sin(Phi);
		 0			-sin(Phi)	cos(Phi)]
	 
RcTta = [cos(Tta)	0			-sin(Tta); 
		 0			1			0; 
		 sin(Tta)	0			cos(Tta)]
	 
RbPsi = [cos(Psi)	sin(Psi)	0; 
		 -sin(Psi)	cos(Psi)	0; 
		 0			0			1]

RaMx  = [a11 a12 a13; 
		 a21 a22 a23; 
		 a31 a32 a33];


      |
      |
      |     /|
      |    / |
      | ? /  |
      |  /   |
      | /    |
      |/_____|_________
     /  \    |
    /     \  |
   /  ?     \|   
  /
 /
/

% radPdeg = unitsratio('radian', 'degrees')
% degPrad = unitsratio('degrees', 'radian')

% 2D Rotation Matrix
% R2dMx = [cos(Tta) -sin(Tta); sin(Tta) cos(Tta)];
% fR2dMx = @(ang) [cos(ang) -sin(ang); sin(ang) cos(ang)];
% fR2dMx(Tta)

r = vrrotvec([1 1 1],[1 2 3])

ax = 1
ay = 1
az = 1
t = 70

M = makehgtform('xrotate',t) 
M = makehgtform('yrotate',t) 
M = makehgtform('zrotate',t) 
M = makehgtform('axisrotate',[ax,ay,az],t) 

M = makehgtform('xrotate',t) returns a transform that rotates around the x-axis by t radians.
M = makehgtform('yrotate',t) returns a transform that rotates around the y-axis by t radians.
M = makehgtform('zrotate',t) returns a transform that rotates around the z-axis by t radians.
M = makehgtform('axisrotate',[ax,ay,az],t) Rotate around axis [ax ay az] by t radians.


? (or phi ? ?) alpha is the angle between the x axis and the N axis.
? (or theta ?) beta is the angle between the z axis and the Z axis.
? (or psi ?) gamma is the angle between the N axis and the X axis.

This definition implies that:

? represents a rotation around the z axis,
? represents a rotation around the N axis,
? represents a rotation around the Z axis.

for ? and ?, the range is defined modulo 2? radians. A valid range could be [??,??].
for ?, the range covers ? radians (but is not modulo ?)

? = atan2(z1, -z2)
? = arccos(Z3)
? = atan2(X3, Y3)

[x,y,z] = sph2cart(THETA,PHI,R)

xyzO = [0;
		0;
		0];

xyzT = [5;
		5;
		5];

xyzA = [45;
		20;
		70];

RxA = [1 0 0; 0 cos(a) sin(a); 0 -sin(a) cos(a)];
RyB = [cos(b) 0 -sin(b); 0 1 0; sin(b) 0 cos(b)];
RzG = [cos(g) sin(g) 0; -sin(g) cos(g) 0; 0 0 1];


% RdPhi = [1 0 0; 0 cos(Phi) sin(Phi); 0 -sin(Phi) cos(Phi)];
% RcTta = [cos(Tta) 0 -sin(Tta); 0 1 0; sin(Tta) 0 cos(Tta)];
% RbPsi = [cos(Psi) sin(Psi) 0; -sin(Psi) cos(Psi) 0; 0 0 1];

% RaMx  = [	a11 a12 a13 ; 
%			a21 a22 a23 ; 
%			a31 a32 a33 ];


% nO1 = fDPhi(Oxyz,zPhi)
% nO2 = fCTta(nO1,xTta)
% nO3 = fBPsy(nO2,zPsy)
% nT1 = fDPhi(Txyz,zPhi)
% nT2 = fCTta(nT1,xTta)
% nT3 = fBPsy(nT2,zPsy)

% nO4 = fDPhi(nO3,zPhi)
% nO5 = fCTta(nO4,xTta)
% nO6 = fBPsy(nO5,zPsy)
% nO6 = nT3;
% nT4 = fDPhi(nT3,zPhi)
% nT5 = fCTta(nT4,xTta)
% nT6 = fBPsy(nT5,zPsy)

% xp = [origin0 tip0; 
%       origin1 tip1]
% xp = [Oxyz(1) Txyz(1); nO3(1) nT3(1); nO6(1) nT6(1)];
% yp = [Oxyz(2) Txyz(2); nO3(2) nT3(2); nO6(2) nT6(2)];
% zp = [Oxyz(3) Txyz(3); nO3(3) nT3(3); nO6(3) nT6(3)];


% Set or get the value of the camera view angle
% camva
% to half or double the zoom use: camzoom(.5) or camzoom(2)
% camzoom(.5)

% vis3d
% campos([-10,-25,25])
% camzoom(.12)


%------------------
% Tx = R * sind(tta) * cosd(phi) + Ox;
% Ty = R * sind(tta) * sind(phi) + Ox;
% Tz = R * cosd(tta) + Oz;

Oxyz3 = Txyz2;
R = 14.1421;

Tx = R * sind(phi) * cosd(tta) + Oxyz3(1);
Ty = R * sind(phi) * sind(tta) + Oxyz3(2);
Tz = R * cosd(phi) + Oxyz3(3);

Txyz3 = [Tx Ty Tz]';

% Txyz3 = rotfun([Tx Ty Tz]', phi, tta, psy);
%------------------


%------------------
% Tx = R * sin(tta) * cos(phi) + Ox;
% Ty = R * sin(tta) * sin(phi) + Ox;
% Tz = R * cos(tta) + Oz;
keyboard
n = ones(3);
p = ones(3);
nxyz = cos(phi) * eye(3) + (1 - cos(p)) * (n.*n) + sin(p) *... 
[0 n(3) -n(2); -n(3) 0 n(1); n(2) -n(1) 0]
%------------------

%}

%{

% Plot Line Connections
RMPr = repmat(Pr,1,Fdots);
RMPt = repmat(Pt,1,Fdots);
RMPo = repmat(Po,1,Fdots); % replace this
RMDv = repmat(Dv,1,Fdots);


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
doPPt = 0;	doLPt = 0; 
doPPo = 0;	doLPo = 0; 
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


disp('t1 v1 d1 t0 v0 d2')
disp({t1, v1, d1, t0, v0, d2})
disp('v1 v0')
disp([v1 v0])


% text(v1(1),v1(2),v1(3),...
% strcat('\leftarrow  v1(',sprintf('%.2f|',v1'),')'),...
% 'FontSize',12,'HorizontalAlignment','left');
% 

% text([t1, v1, d1, t0, v0, d2],...
% strcat(num2str(XYLTpr1(1)), '\bullet',num2str(XYLTpr1(2)),'\rightarrow'),...
% 'FontSize',12,'HorizontalAlignment','right');
% text(XYRTpr1(1),XYRTpr1(2),...
% strcat('\leftarrow ', num2str(XYRTpr1(1)), '\bullet',num2str(XYRTpr1(2))),...
% 'FontSize',12,'HorizontalAlignment','left');
%==================================================%
% OTx = [OT(:).x];
% OTy = [OT(:).y];
% OTz = [OT(:).z];
% OTflen = [OT(:).flen];
% varargout = {OTx, OTy, OTz, OTflen};

% XP = cross(Tip,TipN);
% TipN = XP;

% XP = cross(uv,uvRTip);
% ttaCP = acosd(XP);
% TipN = XP .* Length + Tip;

% Po = [1;1;1];
% Pt = [4;3;5];
% Pr = [5;1;3];
% Dv = Pt-Po; %[3;2;4]
% Tx = r * sind(tta) * cosd(phi);
% Ty = r * sind(tta) * sind(phi);
% Tz = r * cosd(tta);
% tta = atan(z/(x^2+y^2)^.5)
% phi = atan(x/y)
% LDv = sqrt(u^2+v^2+w^2);
% tta2 = atan2(nw,sqrt(nu^2+nv^2));
% phi2 = atan2(nu,nv);
% tta2 = acos(nw/nLDv);
% phi2 = atan(nDv(2)/nDv(1))
% (deal({Dv./LDv}));
% d1sq = [(x1-x0)+(x2-x1)*t]^2 + [(y1-y0)+(y2-y1)*t]^2 + [(z1-z0)+(z2-z1)*t]^2
% d1 = sqrt(d)
% t = |dot((x1-x0),(x2-x1))| / |(x2-x1)|
% d1sq_min = ((|(x1-x0)|^2 * |(x2-x1)|^2) - (dot((x1-x0),(x2-x1))^2)) / |(x2-x1)|^2
% d2 = |cross((x2-x1),(x1-x0))| / |(x2-x1)|
% d2 = |cross((x0-x1),(x0-x2))| / |(x2-x1)|
% L = sqrt(sum(Dv.^2))

%==================================================%
Pr = [5;1;1];
Pt = [5;4;5];
Po = [.5;.5;.5];
Dv = Pt-Po;

x=Pr(1);y=Pr(2);z=Pr(3);
a=Pt(1);b=Pt(2);c=Pt(3);
d=Po(1);e=Po(2);f=Po(3);
u=Dv(1);v=Dv(2);w=Dv(3);

% tta = rand*2*pi;
% doPoint = 1;
% doLine = 1;


ArMx = linspace(.01,2*pi,25);

for ang = 1:25
	tta=ArMx(ang);
	RMX1(:,ang) = RotateVertex(x,y,z,u,v,w,tta,1,0);
end

for ang = 1:25
	tta=ArMx(ang); % Ar = rand*2*pi;
	RMX2(:,ang) = RotateVertex(x,y,z,a,b,c,d,e,f,u,v,w,tta,0,1);
end

% Point = RotateVertex(x,y,z,a,b,c,u,v,w,tta,doPoint,doLine);
% Line = RotateVertex(x,y,z,a,b,c,u,v,w,tta,doPoint,doLine);



%----------------------------------------
%	  FIGURE SETUP
%----------------------------------------
sz = [2.5e-3 2.5e-3 1.6 1.2];
Fh25 = FigSetup(25,sz);
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
return
%==================================================%

%}

%{
%==================================================%
%					SETUP
%--------------------------------------------------%

if nargin == 2
	rot=varargin{1};
	azel=varargin{2};
	do2D=0; do3D=1;
elseif nargin == 1 
	rot=varargin{1};
	azel = [-60 25];
	do2D=0; do3D=1;
else
	rot = [10 0];
	rot1 = rot;
	azel = [-60 25];
	do2D=0; do3D=1;
	sz = [2.5e-3 2.5e-3 1.6 1.2];
end


Ovec = [1.0 28.7 56.4 84.1 111.8 139.5 167.2 194.8 222.5 250.2 277.9 305.6 333.3];
Ovec = Ovec - 180;

%==================================================%
%					3D ROTATION
%--------------------------------------------------%
if do3D
%------------



% STARTING ACTIN FILAMENT SETUP
%---------------------------------------------%
% unitsratio('rad', 'deg')
rad2dg = 180/pi;
dg2rad = 1/(180/pi);
rad360 = 2*pi;

% PHI THETA PSY EULER ANGLES
%---------------------------
phi = 0 *dg2rad;			% Phi Z Rotation
tta = 0 *dg2rad;			% Theta X Rotation
psy = 0 *dg2rad;			% Psy Z Rotation
%---
 
% XYZ VECTOR LOCATION
%---------------------------
Origin = [1; 1; 3];	% XYZ Origin 
Tip = [5; 4; 8];	% XYZ Tip



Length = sqrt(sum((Tip - Origin).^2));
uv = (Tip - Origin) ./ Length;
uvRTip = ARxFun(uv, phi, tta, psy);
TipN = uvRTip .* Length + Tip;



OriginO = Origin;
TipO = TipN;


%---------------------------------------------%

%{
%---------------------------------------------%
nFils = 13;
for iFil = 1:nFils
	
	FilID = mod(iFil,13)+1;
	phi = Ovec(FilID) *dg2rad;	% Phi Z Rotation
	tta = 70 *dg2rad;			% Theta X Rotation
	psy = 0 *dg2rad;			% Psy Z Rotation

	
	Length = sqrt(sum((Tip - Origin).^2));
	uv = (Tip - Origin) ./ Length;
	uvRTip = ARxFun(uv, phi, tta, psy);
	TipN = uvRTip .* Length + Tip;
	
	Oxyz(:,iFil) = Tip;
	Txyz(:,iFil) = TipN;

end
%---------------------------------------------%
%}

%{
%---------------------------------------------%
phi = 0;			% Phi Z Rotation
tta = -70;			% Theta X Rotation
psy = 0;			% Psy Z Rotation
phi=phi*dg2rad;
tta=tta*dg2rad;
psy=psy*dg2rad;
%---

for iFil = 1:nFils
	
	[Oxyz1] = ARxFun(Oxyz0, phi, tta, psy);
	[Txyz1] = ARxFun(Txyz0, phi, tta, psy);
	
	%------------------
	% OTm = (Oxyz1 - Oxyz0);
	% Oxyz1 = Oxyz1 - OTm;
	% Txyz1 = Txyz1 - OTm;
	%------------------
	
	Oxyz(:,iFil) = ARxFun(Oxyz(:,iFil), phi, tta, psy);
	Txyz(:,iFil) = ARxFun(Txyz(:,iFil), phi, tta, psy);
	
	% Txyz(3,iFil) = abs(Txyz(3,iFil));
	
end
%---------------------------------------------%
%}




sFilO = [-15 0 0; 0 -15 0; 0 0 -5];
sFilT = [15 0 0; 0 15 0; 0 0 15];
Oxyz = cat(2,sFilO,Oxyz);
Txyz = cat(2,sFilT,Txyz);

%---------------------------------------------%
for nfil = 1:numel(Oxyz(1,:))
OT(nfil).x = [Oxyz(1,nfil); Txyz(1,nfil)];
OT(nfil).y = [Oxyz(2,nfil); Txyz(2,nfil)];
OT(nfil).z = [Oxyz(3,nfil); Txyz(3,nfil)];
end
OTx = [OT(:).x];
OTy = [OT(:).y];
OTz = [OT(:).z];
%------------------
OTxR = abs(OTx(2,:) - OTx(1,:));
OTyR = abs(OTy(2,:) - OTy(1,:));
OTzR = abs(OTz(2,:) - OTz(1,:));
for rlen = 1:numel(OTx(1,:))
flen(rlen) = sqrt(OTxR(rlen)^2 + OTyR(rlen)^2 + OTzR(rlen)^2);
OT(rlen).flen = flen(rlen);
end
OTflen = [OT(:).flen];
disp([Length OTflen])
%------------------

OTx = cat(2,[OriginO(1) TipO(1)]',OTx);
OTy = cat(2,[OriginO(2) TipO(2)]',OTy);
OTz = cat(2,[OriginO(3) TipO(3)]',OTz);

%------------------
FigAnime(OTx, OTy, OTz);
%------------------


%----------------------------------------
end
%==================================================%






varargout = {OT};
%==================================================%

%}
