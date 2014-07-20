function [] = RotVertDemo(varargin)
clc; close all; clear all;


%---------------------------------------------%
%				ANGLE SETUP
%---------------------------------------------%
% unitsratio('rad', 'deg')
dg2rad = 1/(180/pi);
Ov = [2.0 29.7 57.4 85.1 112.8 140.5 168.2 195.8 223.5 251.2 278.9 306.6 334.3];
Ov = Ov * dg2rad;
%---------------------------------------------%
% (90 deg = pi/2) (70: 7*pi/18) (60: pi/3) (45: pi/4) (30: pi/6)
TPi = (7*pi/18);
PPi = 0;

Lv = @(V) sqrt(sum(V.^2));

Po = [1;1;-1];
Pt = [1;1;2];
Pv = Pt-Po;

tL = sqrt(sum((Pv).^2));	% Length of vector PoPt (aka Pv)
Pu = (Pv) ./ tL;			% Unit vector of Pv

tTheta = acos(Pu(3)) + TPi;	% angle theta
tPhi = atan2(Pu(2),Pu(1)) + PPi;	% angle phi

x = sin(tTheta) * cos(tPhi);
y = sin(tTheta) * sin(tPhi);
z = cos(tTheta);
Pr = [x;y;z]+Pt;


for fN = 1:13
rOv = mod(fN,13)+1;			% Get monomer repeat among the 13 rotational axis angles
Otta = Ov(rOv);				% Rotational angle of new branch
Pn(:,fN) = RotateVertex(Pr(1),Pr(2),Pr(3),Pt(1),Pt(2),Pt(3),Po(1),Po(2),Po(3),...
						 Pv(1),Pv(2),Pv(3),Otta);

end


%----------------------------
Po2 = Pt;
Pt2 = Pn(:,4);
Pv2 = Pt2-Po2;
tL2 = sqrt(sum((Pv2).^2));		% Length of vector PoPt (aka Pv)
Pu2 = (Pv2) ./ tL2;				% Unit vector of Pv
tTheta2 = acos(Pu2(3));			% angle theta	+TPi;
tPhi2 = atan2(Pu2(2),Pu2(1));		% angle phi		+PPi;

% New branch Angles
TPi2 = tTheta2;
PPi2 = tPhi2;
%----------------------------

% TIP of current branch		(not actual tip, just rotational point) 
Pt3x = 2 * sin(tTheta2) * cos(tPhi2) + Po2(1);
Pt3y = 2 * sin(tTheta2) * sin(tPhi2) + Po2(2);
Pt3z = 2 * cos(tTheta2) + Po2(3);



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


%---------------------------------------------%
%	LINE: Pr-Pt Circle to Point
%---------------------------------------------%
PtPrMxX = [Pt3x; Po2(1)];
PtPrMxY = [Pt3y; Po2(2)];
PtPrMxZ = [Pt3z; Po2(3)];
hA2 = plot3(PtPrMxX,PtPrMxY,PtPrMxZ,'b');
FMx.Fdots = {1.5}; Fdots=1;
set(hA2,{'LineWidth'},repmat(FMx.Fdots,1,Fdots)');
hold on;
%---------------------------------------------%





end