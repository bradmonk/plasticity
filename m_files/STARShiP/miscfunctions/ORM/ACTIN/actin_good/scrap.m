function [] = scrap(varargin)
clc; close all; clear all;


%---------------------------------------------%
%				ANGLE SETUP
%---------------------------------------------%
% unitsratio('rad', 'deg')
rad2dg = 180/pi;
dg2rad = 1/(180/pi);
Si = .00001;
Vn = [.66446;.66446;.342];
Veq = [.5774;.5774;.5774];
% Veq = [1;1;1];
Ov = [1.0 28.7 56.4 84.1 111.8 139.5 167.2 194.8 222.5 250.2 277.9 305.6 333.3];
Ov = Ov * dg2rad;
%---------------------------------------------%

% (90 deg = pi/2) (70: 7*pi/18) (60: pi/3) (45: pi/4) (30: pi/6)
%TPi = (pi/2);
TPi = (7*pi/18);
PPi = 0;

Lv = @(V) sqrt(sum(V.^2));
Lv(Veq);

Po = [-1;-1;-1];
Pt = [-2;-2;-2];
Pv = Pt-Po;

tL = sqrt(sum((Pv).^2));	% Length of vector PoPt (aka Pv)
Pu = (Pv) ./ tL;			% Unit vector of Pv

tTheta = acos(Pu(3)) + TPi;	% angle theta
tPhi = atan(Pu(2)/Pu(1)) + PPi;	% angle phi

x = sin(tTheta) * cos(tPhi);
y = sin(tTheta) * sin(tPhi);
z = cos(tTheta);
Pr = [x;y;z]+Pt;


for fN = 1:13
rOv = mod(fN,13)+1;			% Get monomer repeat among the 13 rotational axis angles
tta = Ov(rOv);				% Rotational angle of new branch
Pn(:,fN) = RotateVertex(Pr(1),Pr(2),Pr(3),Pt(1),Pt(2),Pt(3),Po(1),Po(2),Po(3),...
						 Pv(1),Pv(2),Pv(3),tta);

end



%---------------------------------------------%
% DOT: {Po Po Pt Dv Pr v0 v1} {'d';'d';'o';'x';'s';'v';'v'};
%---------------------------------------------%
addaxis(3); hold on;
[XMx YMx ZMx] = scatter3prep({Po Pt Pr Pv Pu});
Ph1 = scatter3(XMx,YMx,ZMx,'.','Tag','TMxH');
grid on;
%----------------
Ph1c = get(Ph1,'Children');
mrkr = {'d';'o';'o';'o';'o'};
set(Ph1c,{'Marker'},mrkr);
hold on
colr = {[1 0 0]; [1 .3 .7]; [0 1 0]; [0 0 0]; [1 1 1]};
set(Ph1c,{'MarkerFaceColor'},colr);
hold on
set(Ph1c,'LineWidth',1,'MarkerSize',9);
simlegend({'Po','Pt','Pr','Pv','Pu'},Ph1);
hold on
Ph2 = scatter3(Pn(1,:)',Pn(2,:)',Pn(3,:)','.');
view(0,0)
%---------------------------------------------%

end