function [] = Trig()
%==================================================%

% Trig Images
%--------------
[euler,map,alpha] = imread('euler.png','BackgroundColor',[1 1 1]);
figure('color','w'); image(euler)
colormap(map); axis off; axis image 
%---------------------------------------------%


% unitsratio('rad', 'deg')
%---------------------------
rad2dg = 180/pi;
dg2rad = 1/(180/pi);
rad360 = 2*pi;
%---------------------------------------------%
% sin(?) =	opp/h	= tan(?)/sec(?) = 1/csc(?)
% cos(?) =	adj/h	= cot(?)/csc(?) = 1/sec(?)
% tan(?) =	opp/adj	= sin(?)/cos(?) = 1/cot(?)
% sec(?) =	h/adj	= csc(?)/cot(?) = 1/cos(?)
% cot(?) =	adj/opp	= cos(?)/sin(?) = 1/tan(?)
% csc(?) =	h/opp	= sec(?)/tan(?) = 1/sin(?)
%---------------------------------------------%


% XYZ POINT VECTOR
%---------------------------
O = [0; 0; 0];	% XYZ Origin 
A = [0; 0; 5];	% XYZ Tip
%---------------------------------------------%


% Length in R3
%--------------
L = sqrt(sum((A - O).^2));
L = sqrt((A(1)-O(1))^2 + (A(2)-O(2))^2 + (A(3)-O(3))^2);
%---------------------------------------------%


% Standard Basis Vectors i,j,k
%--------------
i = [1;0;0];
j = [0;1;0];
k = [0;0;1];

% if a = [a1,a2,a3] then a = a1*i + a2*j + a3*k
% let a = [3;4;5]
a1=3; a2=4; a3=5;
a = a1*i + a2*j +a3*k
%---------------------------------------------%


% Unit Vectors u 
%--------------
% unit vectors have a length = 1 (i,j,k are all unit vectors)
% let a = [3;4;5] and length of a = |a| and c = 1/|a|

% a = [3;4;5];
a = [4;3;0];
aL = sqrt(sum(a.^2));
c = 1/aL;

u = a / aL
uL = c*aL
%---------------------------------------------%



% Dot Product
%--------------
% The dot product of two vectors is defined as:
% a dot b = |a||b|cos(theta)
% where theta is the angle between two vectors
% it can also be defined in component form as:
% a dot b = a1b1 + a2b2 + a3b3

a = [4;3;0];
b = [4;0;0];

aL = sqrt(sum(a.^2));
bL = sqrt(sum(b.^2));
tta = 0.6435; % radians

% 3 ways to get the dot product:

adb = aL*bL*cos(tta)

adb = sum(a.*b)

adb = dot(a,b)

% acos(.8) * (180/pi)
%---------------------------------------------%






%---------------------------------------------%
% sin(?) =	opp/h	= tan(?)/sec(?) = 1/csc(?)
% cos(?) =	adj/h	= cot(?)/csc(?) = 1/sec(?)
% tan(?) =	opp/adj	= sin(?)/cos(?) = 1/cot(?)
% sec(?) =	h/adj	= csc(?)/cot(?) = 1/cos(?)
% cot(?) =	adj/opp	= cos(?)/sin(?) = 1/tan(?)
% csc(?) =	h/opp	= sec(?)/tan(?) = 1/sin(?)
%---------------------------------------------%

% Defining a point in a rotated axis
%--------------
% if a vector point is defined in the standard i,j,k axis base
% it may be useful to define it in another orthonormal axis
% for the purpose of rotating a point around a non-standard axis.

% First lets define a few things:

% the standard bases unit vectors i,j,k [e1;e2;e3]:
e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];
e = [e1 e2 e3];


% the point vector in R3 at standard axis 'a'
% a = av.*e1 + av.*e2 + av.*e3;
% and the new axis direction vector 'n'
o = [0;0;0];
a = [4;11;0];
n = [4;0;0]; 

% we now need to convert the vectors into unit vectors:
aL = sqrt(sum(a.^2));
ua = a / aL;

nL = sqrt(sum(n.^2));
un = n / nL;

% the dot product between the two vectors provides the cosine of the angle 
% between the two vectors, which can be converted to the angle (theta)
% by taking the inverse cosine (acosd) of the dot product:

dotP = dot(ua,un)
dotP = (ua(1)*un(1) + ua(2)*un(2) + ua(3)*un(3))
ttaDP = acosd(dotP)

crossP = cross(un,ua)
ttaCP = acosd(crossP)

OrthV = crossP .* ((nL+aL)/2)


OTx = [0 0 0; 4 4 OrthV(1)];
OTy = [0 0 0; 11 0 OrthV(2)];
OTz = [0 0 0; 0 0 OrthV(3)];
plot3(OTx,OTy,OTz)
axis([-10 10 -10 10 0 10])
xlabel('X');ylabel('Y');zlabel('Z');
grid on; % view(0,90)


















% EULER ANGLES
%---------------------------
phi = 0 *dg2rad;			% Phi Z Rotation
tta = 0 *dg2rad;			% Theta X Rotation
psy = 0 *dg2rad;			% Psy Z Rotation
%---




%---------------------------------------------%



% Rotation About an Arbitrary Axis in 3 Dimensions
%--------------
%{
Step-1: Translate space so O is at the orthonormal axis origin

Step-2: Rotate space so the rotation axis lies on the xz plane
		That is, rotate 'a' around the z-axis

Step-3: Rotate space around the y-axis so the rotation axis
		lies along the z-axis. 'a' moves to the z-axis

Step-4: Perform the desired rotation by theta and phi along the
		z- and y- axes respectively

Step-5: Apply inverse of Step-3, then Step-2, then Step-1
%}


% Translation Matrix:
TMx = @(ai,bi,ci) ([1 0 0 0; 0 1 0 0; 0 0 1 0; -ai -bi -ci 1]');
TP1 = TMx(1,2,3)

%---------------------------------------------%
% Rotation matrix Rzyx (rotates about x, then y, then z);
Rzyx = @(alp,bet,gam) ...
	([cos(bet)*cos(gam) cos(gam)*sin(alp)*sin(bet)-cos(alp)*sin(gam)...
	cos(alp)*cos(gam)*sin(bet)+sin(alp)*sin(gam) 0;...
	cos(bet)*sin(gam) cos(alp)*cos(gam)+sin(alp)*sin(bet)*sin(gam)...
	-cos(gam)*sin(alp)+cos(alp)*sin(bet)*sin(gam) 0;...
	-sin(bet) cos(bet)*sin(alp) cos(alp)*cos(bet) 0; 0 0 0 1]);
MRzyx = Rzyx(5,10,15)

% Rotation matrix Rxyz (rotates about z, then y, then x);
Rxyz = @(alp,bet,gam) ...
	([cos(bet)*cos(gam) -cos(bet)*sin(gam) sin(bet) 0;
	cos(alp)*sin(gam)+sin(alp)*sin(bet)*cos(gam) ...
	cos(alp)*cos(gam)-sin(alp)*sin(bet)*sin(gam) ...
	-sin(alp)*cos(bet) 0; ...
	sin(alp)*sin(gam)-cos(alp)*sin(bet)*cos(gam) ...
	sin(alp)*cos(gam)+cos(alp)*sin(bet)*sin(gam) ...
	cos(alp)*cos(bet) 0; 0 0 0 1;]);
MRxyz = Rxyz(5,10,15)


%---------------------------------------------%
% Rotation matrix to rotate a vector about the z-axis to the xz-plane
Txz = @(u,v) ...
	([u/sqrt(u^2+v^2) v/sqrt(u^2+v^2) 0 0;
	-v/sqrt(u^2+v^2) u/sqrt(u^2+v^2) 0 0;
	0 0 1 0; 0 0 0 1;])
MTxz = Txz(1,2)

% Rotation matrix to rotate the vector in the xz-plane to the z-axis
Tz = @(u,v,w) ...
	([w/sqrt(u^2+v^2+w^2) 0 -sqrt(u^2+v^2)/sqrt(u^2+v^2+w^2) 0; 0 1 0 0;
	sqrt(u^2+v^2)/sqrt(u^2+v^2+w^2) 0 w/sqrt(u^2+v^2+w^2) 0; 0 0 0 1]);
MTz = Tz(1,2,3)

%---------------------------------------------%
% a function of of seven variables that yields the result of
% rotating the point (x,y,z) about the axis ?u,v,w? by the angle ?.
Fxut = @(x,y,z,u,v,w,tta) ...
	([(u*(u*x+v*y+w*z)*(1-cos(tta))+(u^2+v^2+w^2)*...
	x*cos(tta)+sqrt(u^2+v^2+w^2)*(-w*y+v*z)*sin(tta))/... 
	(u^2 + v^2 + w^2);
	(v*(u*x+v*y+w*z)*(1-cos(tta))+(u^2+v^2+w^2)*...
	y*cos(tta)+sqrt(u^2+v^2+w^2)*(w*x-u*z)*sin(tta))/... 
	(u^2+v^2+w^2);
	(w*(u*x+v*y+w*z)*(1-cos(tta))+(u^2+v^2+w^2)*...
	z*cos(tta)+sqrt(u^2+v^2+w^2)*(-v*x+u*y)*sin(tta))/... 
	(u^2+v^2+w^2)]);
MFxut = Fxut(1,2,3,4,5,6,.8)

for ang = 1:20
	tta = rand*2*pi;
	RMX(:,ang) = Fxut(1,2,3,4,5,6,tta);
end

xyz = [1;2;3];
uvw = [4;5;6];
RMXO = repmat(uvw,1,20);

Xop = [RMXO(1,:); RMX(1,:)];
Yop = [RMXO(2,:); RMX(2,:)];
Zop = [RMXO(3,:); RMX(3,:)];

plot3(Xop,Yop,Zop)
hold on
% scatter3(Xop(2,:),Yop(2,:),Zop(2,:))


%---------------------------------------------%
% The matrix for rotation about an arbitrary line
% If we multiply this times ?x,y,z? we can obtain a function of of ten variables 
% that yields the result of rotating the point (x,y,z) about the 
% line through (a,b,c) with direction vector ?u,v,w? by the angle ?.
Fxaut = @(x,y,z,a,b,c,u,v,w,tta) ... 
	([(((a*(v^2+w^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(tta))+(u^2+v^2+w^2)*x*cos(tta)+...
	sqrt((u^2+v^2+w^2))*(-c*v+v*w-w*y+v*z)*sin(tta))/(u^2+v^2+w^2));
	(((b*(u^2+w^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(tta))+(u^2+v^2+w^2)*y*cos(tta)+...
	sqrt((u^2+v^2+w^2))*(c*u-a*w+w*x-u*z)*sin(tta))/(u^2+v^2+w^2));
	(((c*(u^2+v^2)-w*(a*u+c*w-u*x-v*y-w*z))*(1-cos(tta))+(u^2+v^2+w^2)*z*cos(tta)+...
	sqrt((u^2+v^2+w^2))*(-b*u+a*v-v*x+u*y)*sin(tta))/(u^2+v^2+w^2))]);
MFxaut = Fxaut(5,0,0,1,1,1,1,1,1,.8)

for ang = 1:20
	tta = rand*2*pi;
	RMX(:,ang) = Fxaut(5,0,0,1,1,1,1,1,1,tta);
end

xyz = [5;0;0];
abc = [1;1;1];
uvw = [0;0;0];
RMXO = repmat(uvw,1,20);

Xop = [RMXO(1,:); RMX(1,:)];
Yop = [RMXO(2,:); RMX(2,:)];
Zop = [RMXO(3,:); RMX(3,:)];

plot3(Xop,Yop,Zop)
% scatter3(Xop(2,:),Yop(2,:),Zop(2,:))
%---------------------------------------------%




xyz = [5;0;0];
abc = [1;1;1];
uvw = [0;0;0];

x=xyz(1);y=xyz(2);z=xyz(3);
a=abc(1);b=abc(2);c=abc(3);
u=uvw(1);v=uvw(2);w=uvw(3);
UVW0 = repmat(uvw,1,50);

for ang = 1:50
	tta = rand*2*pi;
	RMX1(:,ang) = Fxut(x,y,z,u,v,w,tta);
end

Xop1 = [UVW0(1,:); RMX1(1,:)];
Yop1 = [UVW0(2,:); RMX1(2,:)];
Zop1 = [UVW0(3,:); RMX1(3,:)];

for ang = 1:50
	tta = rand*2*pi;
	RMX2(:,ang) = Fxaut(x,y,z,a,b,c,u,v,w,tta);
end

Xop2 = [UVW0(1,:); RMX2(1,:)];
Yop2 = [UVW0(2,:); RMX2(2,:)];
Zop2 = [UVW0(3,:); RMX2(3,:)];


figure(1)
plot3(Xop1,Yop1,Zop1)
hold on
plot3(Xop2,Yop2,Zop2)

figure(2)
scatter3(Xop1(2,:),Yop1(2,:),Zop1(2,:),'.r')
hold on
scatter3(Xop2(2,:),Yop2(2,:),Zop2(2,:),'.b')
hold on
scatter3(xyz(1),xyz(2),xyz(3),'o')
hold on
scatter3(abc(1),abc(2),abc(3),'d')
hold on
scatter3(uvw(1),uvw(2),uvw(3),'x')
hold on
axis([-15 15 -15 15 -15 15])

%---------------------------------------------%



% Dot Product
%--------------
%
%---------------------------------------------%



% Dot Product
%--------------
%
%---------------------------------------------%
















%==================================================%
end
%==================================================%





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



function [] = ExampleMainFunction()
clc; close all; clear all;


%---------------------------------------------%
% a function of of seven variables that yields the result of
% rotating the point (x,y,z) about the axis ?u,v,w? by the angle ?.
Fxut = @(x,y,z,u,v,w,tta) ...
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
% We will define an arbitrary line by a point the line goes through and a 
% direction vector. If the axis of rotation is given by two points P1 = (a,b,c) 
% and P2 = (d,e,f), then a direction vector can be obtained by 
% ?u,v,w? = ?d-a,e-b,f-c?
% We can now write a transformation for the rotation of a point about this line.
% 
% The matrix for rotation about an arbitrary line
% If we multiply this times ?x,y,z? we can obtain a function of of ten variables 
% that yields the result of rotating the point (x,y,z) about the 
% line through (a,b,c) with direction vector ?u,v,w? by the angle ?.
Fxaut = @(x,y,z,a,b,c,u,v,w,tta) ... 
	([(((a*(v^2+w^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(tta))+(u^2+v^2+w^2)*x*cos(tta)+...
	sqrt((u^2+v^2+w^2))*(-c*v+v*w-w*y+v*z)*sin(tta))/(u^2+v^2+w^2));
	(((b*(u^2+w^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(tta))+(u^2+v^2+w^2)*y*cos(tta)+...
	sqrt((u^2+v^2+w^2))*(c*u-a*w+w*x-u*z)*sin(tta))/(u^2+v^2+w^2));
	(((c*(u^2+v^2)-w*(a*u+c*w-u*x-v*y-w*z))*(1-cos(tta))+(u^2+v^2+w^2)*z*cos(tta)+...
	sqrt((u^2+v^2+w^2))*(-b*u+a*v-v*x+u*y)*sin(tta))/(u^2+v^2+w^2))]);
% MFxaut = Fxaut(5,0,0,1,1,1,1,1,1,.8)
%---------------------------------------------%




xyz = [3;3;3];
abc = [5;0;0];
def = [1;1;1];
uvw = abc-def;
% uvw = uvw ./ sqrt(sum(uvw.^2));

abcT = abc;
abcO = def;

x=xyz(1);y=xyz(2);z=xyz(3);
a=abc(1);b=abc(2);c=abc(3);
u=uvw(1);v=uvw(2);w=uvw(3);
UVW0 = repmat(uvw,1,50);

for ang = 1:50
	tta = rand*2*pi;
	RMX1(:,ang) = Fxut(x,y,z,u,v,w,tta);
end

Xop1 = [UVW0(1,:); RMX1(1,:)];
Yop1 = [UVW0(2,:); RMX1(2,:)];
Zop1 = [UVW0(3,:); RMX1(3,:)];

for ang = 1:50
	tta = rand*2*pi;
	RMX2(:,ang) = Fxaut(x,y,z,a,b,c,u,v,w,tta);
end

Xop2 = [UVW0(1,:); RMX2(1,:)];
Yop2 = [UVW0(2,:); RMX2(2,:)];
Zop2 = [UVW0(3,:); RMX2(3,:)];


% figure(1)
% plot3(Xop1,Yop1,Zop1)
% hold on
% plot3(Xop2,Yop2,Zop2)

figure(2)
% scatter3(Xop1(2,:),Yop1(2,:),Zop1(2,:),'.r')
% hold on
scatter3(Xop2(2,:),Yop2(2,:),Zop2(2,:),'.b')
hold on
scatter3(xyz(1),xyz(2),xyz(3),'o')
hold on
scatter3(abc(1),abc(2),abc(3),'d')
hold on
scatter3(uvw(1),uvw(2),uvw(3),'x')
axis([-15 15 -15 15 -15 15])
hold on

sFilO = [-15 0 0; 0 -15 0; 0 0 -5];
sFilT = [15 0 0; 0 15 0; 0 0 15];
sFilO = cat(2,abcO,sFilO);
sFilT = cat(2,abcT,sFilT);
Xop0 = [sFilO(1,:); sFilT(1,:)];
Yop0 = [sFilO(2,:); sFilT(2,:)];
Zop0 = [sFilO(3,:); sFilT(3,:)];
plot3(Xop0,Yop0,Zop0)
axis vis3d
%---------------------------------------------%

keyboard

end







