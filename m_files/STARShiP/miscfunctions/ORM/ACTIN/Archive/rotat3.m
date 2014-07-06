function [varargout] = rotat3(varargin)

clc; close all; clear all;

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
phi = 180;			% Phi Z Rotation
tta = 360;			% Theta X Rotation
psy = 180;			% Psy Z Rotation

% phi0 = 20;			% Phi Z Rotation
% tta0 = 70;			% Theta X Rotation
% psy0 = 0;			% Psy Z Rotation
%---

% phi=phi*dg2rad;
% tta=tta*dg2rad;
% psy=psy*dg2rad;

% XYZ VECTOR LOCATION
%---------------------------
% Origin = [0; 0; 0];	% XYZ Origin 
% Tip = [0; 0; 5];	% XYZ Tip
Origin = [1; 1; 3];	% XYZ Origin 
Tip = [5; 4; 8];	% XYZ Tip

	Length = sqrt(sum((Tip - Origin).^2));
	uv = (Tip - Origin) ./ Length;
	uvRTip = ARxFun(uv, phi, tta, psy);
	TipON = uvRTip .* Length + Origin;
	
OriginO = Origin;
TipO = TipON;

%---
Length = sqrt(sum((Tip - Origin).^2));
nFils = 13;
Oxyz = Origin;
Txyz = Tip;
%---------------------------------------------%

%{.
%---------------------------------------------%
for iFil = 1:nFils
	
	FilID = mod(iFil,13)+1;
	phi = Ovec(FilID);	% Phi Z Rotation
	tta = 70;			% Theta X Rotation
	psy = 0;			% Psy Z Rotation
 	phi=phi*dg2rad;
 	tta=tta*dg2rad;
 	psy=psy*dg2rad;
	
	uv = (Tip - Origin) ./ Length;
	uvRTip = ARxFun(uv, phi, tta, psy);
	TipN = uvRTip .* Length + Origin;
	
	Oxyz(:,iFil) = Origin;
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
disp([Length OTflen]')
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




%==================================================%
% OTx = [OT(:).x];
% OTy = [OT(:).y];
% OTz = [OT(:).z];
% OTflen = [OT(:).flen];
% varargout = {OTx, OTy, OTz, OTflen};

varargout = {OT};
%==================================================%

end



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









