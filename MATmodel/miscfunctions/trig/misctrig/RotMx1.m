function [varargout] = RotMx1(varargin)
clc, close all, clear all

Po = [1.5;.5;.5];	% LinePoint Origin
Pt = [2;-1.5;1.5];	% LinePoint Tip
Dv = Pt-Po;			% Vector Length

% Various Branching Angles
rad2dg = 180/pi;
dg2rad = 1/(180/pi);

Ovec = linspace(0,360,20) * dg2rad;

% Selected Branching Angle
TPi = 70 * dg2rad;
PPi = 0;

TL0 = sqrt(sum((Pt - Po).^2));
Txyz = (Pt - Po) ./ TL0;
TL = sqrt(sum((Txyz).^2));
Ttta = acos(Txyz(3)/TL) + TPi;
Tphi = atan(Txyz(2)/Txyz(1)) + PPi;
TL=3;
LTx = (TL) * sin(Ttta) * cos(Tphi);
LTy = (TL) * sin(Ttta) * sin(Tphi);
LTz = (TL) * cos(Ttta);
Pr = [LTx;LTy;LTz]+Pt;

fPo = Po;
fPt = Pt;
fPr = Pr;
fDv = Dv;
fDr = [];



Novec = numel(Ovec);
for Lnums = 1:3
%---------------
	for fN = 1:Novec
	%---------------
	PMX.Po(:,fN) = fPo;
	PMX.Pt(:,fN) = fPt;
	PMX.Pr(:,fN) = fPr;
	PMX.Dv(:,fN) = fDv;
	
 	Rmono = mod(fN,numel(Ovec))+1;			% Get monomer repeat among the 13 rotational axis angles
	tta = Ovec(Rmono);						% Rotational angle of new branch
	fDr(:,fN) = RotateVertex(fPr(1),fPr(2),fPr(3),fPt(1),fPt(2),fPt(3),fPo(1),fPo(2),fPo(3),...
		fDv(1),fDv(2),fDv(3),tta);
		
	PMX.Dr(:,fN) = fDr(:,fN);
	%---------------
	end


	CMX.Po(:,Lnums) = fPo;
	CMX.Pt(:,Lnums) = fPt;
	CMX.Pr(:,Lnums) = fPr;
	CMX.Dv(:,Lnums) = fDv;
	CMX.Dr{Lnums} = PMX.Dr;
	CMX.PMX{Lnums} = PMX;
	
	
	Rfil = ceil(200 * rand);
	Pfil = mod(Rfil,Novec)+1;
	%-----
	fPo = fPt;
	fPt = fDr(:,Pfil);
	fDv = fPt-fPo;
	%-----
	TL0 = sqrt(sum((fPt - fPo).^2));
	Txyz = (fPt - fPo) ./ TL0;
	TL = sqrt(sum((Txyz).^2));
	Ttta = acos(Txyz(3)/TL) + TPi;
	Tphi = atan2(Txyz(2),Txyz(1)) + PPi;
	TL=3;
	LTx = (TL) * sin(Ttta) * cos(Tphi);
	LTy = (TL) * sin(Ttta) * sin(Tphi);
	LTz = (TL) * cos(Ttta);
	Pe = [LTx;LTy;LTz]+fPt;
	%-----
	fPr = Pe;
	%-----
	
%---------------
end



for PLoop = 1:Lnums
%---------------
	Po = CMX.Po(:,PLoop);
	Pt = CMX.Pt(:,PLoop);
	Pr = CMX.Pr(:,PLoop);
	Dv = CMX.Dv(:,PLoop);
	Dr = CMX.Dr{PLoop};

	%---------------
	% FIGURE SETUP
	%---------------
	scsz = get(0,'ScreenSize'); sz = [2.5e-3 2.5e-3 1.6 1.2];
	pos=scsz./sz; Fh25 = figure(25);
	set(Fh25,'OuterPosition',pos,'Color',[.9,.9,.9])
	%---------------
	% PLOT ORIGIN->TIP LINE
	%---------------
	[XMx YMx ZMx] = plot3prep({Pt Pr},{Po Pt});
	figure(Fh25);
	hA1 = plot3(XMx,YMx,ZMx,'LineWidth',2,'MarkerSize',10);
	xlabel('X');ylabel('Y');zlabel('Z');
	axis equal; grid on;
	lne = {'-';':'};
	set(hA1,{'LineStyle'},lne);
	hold on;
	%---------------
	axvals = {[-10 10 -10 10 -10 10]};
	addaxis(axvals);
	axis(axvals{1})
	%---------------
	% DOT: Pr ROTATED AROUND LINE (BLUE)
	%---------------
	scatter3(Dr(1,:),Dr(2,:),Dr(3,:),'.b');
	hold on
%---------------
end


%---------------
end % END RotMx FUNCTION
%======================================================================%


%---------------------------------------------%
%	ADD STANDARD IJK ORTHONORMAL AXIS
%---------------------------------------------%
function [varargout] = addaxis(varargin)

if nargin == 1
	axvals = varargin{1};
XAxO = [axvals{1}(1) 0 0]; XAxT = [axvals{1}(2) 0 0];
YAxO = [0 axvals{1}(3) 0]; YAxT = [0 axvals{1}(4) 0];
ZAxO = [0 0 axvals{1}(5)]; ZAxT = [0 0 axvals{1}(6)];




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









