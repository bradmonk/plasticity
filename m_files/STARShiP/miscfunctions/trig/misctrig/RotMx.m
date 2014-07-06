function [varargout] = RotMx(varargin)
clc

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

%---------------------------------------------%

fPo = Po;
fPt = Pt;
fPr = Pr;
fDv = Dv;
fDr = [];
%---------------------------------------------%


Novec = numel(Ovec);

%====================================================%
Fdots = Novec;
%---------------------------------------------%
%	  RotateVertex AROUND LINE LOOP
%---------------------------------------------%
for Lnums = 1:4
	
for fN = 1:Fdots
	
	PMX.Po(:,fN) = fPo;
	PMX.Pt(:,fN) = fPt;
	PMX.Pr(:,fN) = fPr;
	PMX.Dv(:,fN) = fDv;
	
 	Rmono = mod(fN,numel(Ovec))+1;			% Get monomer repeat among the 13 rotational axis angles
	tta = Ovec(Rmono);						% Rotational angle of new branch
	fDr(:,fN) = RotateVertex(fPr(1),fPr(2),fPr(3),fPt(1),fPt(2),fPt(3),fPo(1),fPo(2),fPo(3),...
		fDv(1),fDv(2),fDv(3),tta);
		
	PMX.Dr(:,fN) = fDr(:,fN);
	
end


	CMX.Po(:,Lnums) = fPo;
	CMX.Pt(:,Lnums) = fPt;
	CMX.Pr(:,Lnums) = fPr;
	CMX.Dv(:,Lnums) = fDv;
	CMX.Dr{Lnums} = PMX.Dr;
	CMX.PMX{Lnums} = PMX;
	
	%---------------------------------------------%
	Rfil = ceil(200 * rand);
	Pfil = mod(Rfil,Fdots)+1;
	%-----
	fPo = fPt;
	fPt = fDr(:,Pfil); % fPt = fPr;
	fDv = fPt-fPo;
	%---------------------------------------------%
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
	%---------------------------------------------%
	
	DMX.Po(:,Lnums) = fPo;
	DMX.Pt(:,Lnums) = fPt;
	DMX.Pr(:,Lnums) = fPr;
	DMX.Dv(:,Lnums) = fDv;
	DMX.Dr{Lnums} = PMX.Dr;

%---------------------------
PMXCELL{Lnums} = PMX;
DMXCELL{Lnums} = DMX;
CMXCELL{Lnums} = CMX;
end


%====================================================%
for PLoop = 1:Lnums
%====================================================%
Po = CMX.Po(:,PLoop);
Pt = CMX.Pt(:,PLoop);
Pr = CMX.Pr(:,PLoop);
Dv = CMX.Dv(:,PLoop);
Dr = CMX.Dr{PLoop};

%---------------------------------------------%
%				FIGURE SETUP
%---------------------------------------------%
sz = [2.5e-3 2.5e-3 1.6 1.2];
fpos = [.05 .05 .95 .95];
Fh25 = FigSetup(25,sz,fpos);
%---------------------------------------------%
%		PLOT ORIGIN->TIP LINE
%---------------------------------------------%
figure(Fh25);
[XMx YMx ZMx] = plot3prep({Pt Pr},{Po Pt});
hA1 = plot3(XMx,YMx,ZMx,'LineWidth',2,'MarkerSize',10);
grid on;
xlabel('X');ylabel('Y');zlabel('Z');

axis equal
hold on;

axvals = {[-8 8 -8 8 -8 8]};
addaxis(axvals);
axis(axvals{1})

%---------------------------------------------%


%---------------------------------------------%
% DOT: Pr ROTATED AROUND LINE (BLUE)
%---------------------------------------------%
scatter3(Dr(1,:),Dr(2,:),Dr(3,:),'.b');
hold on
%---------------------------------------------%

%---------------------------------------------%
%	LINE: Pr-Pt Circle to Point
%---------------------------------------------%
PtMx = repmat(Pt,1,Fdots);
PtPrMxX = [Dr(1,:); PtMx(1,:)];
PtPrMxY = [Dr(2,:); PtMx(2,:)];
PtPrMxZ = [Dr(3,:); PtMx(3,:)];

hA2 = plot3(PtPrMxX,PtPrMxY,PtPrMxZ,'b');
FMx.Fdots = {1.5};
set(hA2,{'LineWidth'},repmat(FMx.Fdots,1,Fdots)');
hold on;
%---------------------------------------------%


% %---------------------------------------------%
% %	LINE: Pr-Pt Circle to Point
% %---------------------------------------------%
% PtMx = repmat(Pt,1,Fdots);
% PtPrMxX = [PtMx(1,:); Dr(1,:)];
% PtPrMxY = [PtMx(2,:); Dr(2,:)];
% PtPrMxZ = [PtMx(3,:); Dr(3,:)];
% 
% hA2 = plot3(PtPrMxX,PtPrMxY,PtPrMxZ,'b');
% FMx.Fdots = {1.5};
% set(hA2,{'LineWidth'},repmat(FMx.Fdots,1,Fdots)');
% hold on;
% %---------------------------------------------%

%---------------------------------------------%
[XMx YMx ZMx] = scatter3prep({Po Pt Pr Dv});
hA3 = scatter3(XMx,YMx,ZMx,'.','Tag','TMxH');
%----------------
hA3c = get(hA3,'Children');
mrkr = {'x';'o';'d';'s'}; %markers are reversed
set(hA3c,{'Marker'},mrkr);
hold on
colr = {[0 0 0]; [1 0 0]; [0 0 0]; [0 1 0]};
set(hA3c,{'MarkerFaceColor'},colr);
hold on
set(hA3c,'LineWidth',1,'MarkerSize',11);
% legend(hAc99,fliplr({'Po','Pt','Pr','Dv','v0','v1'}))
%---------------------------------------------%



%====================================================%
end
%====================================================%



%---------------------------------------------%
% PLOTSUGAR: ANOTATE ANIMATE
%---------------------------------------------%
% simlege({'Po','Pt','Dv','Pr'},hA3);
% FIGRO();
% anotex({Po, 'Po'},{Pt, 'Pt'},{Pr, 'Pr'},{v0, 'v0'},{v1, 'v1'});
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

%==================================================%
%	3D VECTOR LENGTH FUNCTION
%--------------------------------------------------%
function [L] = LV(V)

L = sqrt(sum(V.^2));

end
%==================================================%


%==================================================%
%	3D VECTOR LENGTH FUNCTION
%--------------------------------------------------%
function [varargout] = FigSetup(varargin)

scsz = get(0,'ScreenSize');


if nargin == 1 
	Fnum=varargin{1};
	pos = scsz./[2.5e-3 2.5e-3 1.5 2];
elseif nargin >= 2
	Fnum=varargin{1};
	pos=scsz./varargin{2};
else
	Fnum=1;
	pos = scsz./[2.5e-3 2.5e-3 1.5 2];
end

Fh = figure(Fnum);
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])



if nargin >= 3
	fpos=varargin{3};
	set(gca,'Position',fpos)
end

varargout = {Fh};

end
%--------------------------------------------------%



%==================================================%
%	FIGROT
%--------------------------------------------------%
function [varargout] = FIGRO(varargin)


if nargin == 1
	Fh = varargin{1};
else
	Fh = gcf;
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
set(gca,'CameraPosition',gotme{1},'CameraTarget',gotme{2},'CameraViewAngle',gotme{3});
varargout = {campos, camtarget};
%----------------------
end
%--------------------------------------------------%


%==================================================%
%	simlege
%--------------------------------------------------%
function [varargout] = simlege(LABS,varargin)
add = numel(varargin);
LABS = fliplr(LABS);
if add
	hP = varargin{1};
	hAc = get(hP,'Children');
	hL = legend(hAc,LABS);
	[hL,hOb,hP,hTxt] = legend;
	hold on
end
varargout = {hL,hOb,hP,hTxt};
end
%--------------------------------------------------%


%==================================================%
%	simlege
%--------------------------------------------------%
function [] = anotex(varargin)

tbgc = [.91 .91 .91];
nt = nargin;

% inputname(1)

for an = 1:nt
val=[];pos=[];txt=[];

if numel(varargin{an}) >= 1
val = varargin{an}{1};
if numel(val) == 3
val=[val(1) val(2) val(3)];
else
val=[0 0 0]+an;
end
else
val=[0 0 0]+an;
end

if numel(varargin{an}) >= 2
txt = varargin{an}{2};
txt = strcat('\leftarrow ',txt);
else
txt='\bullet';
end

if numel(varargin{an}) >= 3
pos = varargin{an}{3};
	if numel(pos) == 3
	pos=[pos(1) pos(2) pos(3)];
	else
	pos=val;
	end
else
pos=val;
end

% hT = text(pos(1),pos(2),pos(3),...
% strcat('\bullet ', txt,'(', sprintf('%.2f|',val'),')'),...
% 'FontSize',12,'HorizontalAlignment','left');

hT = text(pos(1),pos(2),pos(3),...
strcat('\bullet ', txt,'(', sprintf('%.2f|',val'),')'),...
'FontSize',12,'HorizontalAlignment','left',...
'BackgroundColor',tbgc);

end

% hBGC = findobj('BackgroundColor',tbgc);
alpha(.05);

end
%--------------------------------------------------%




