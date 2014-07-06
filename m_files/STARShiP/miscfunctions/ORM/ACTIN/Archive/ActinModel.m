function [] = ActinModel()
clc; close all; clear all;
%{
r = 10;
za = 45;
xa = 70;
ya = 90 - xa;

x = r * sind(za) * cosd(xa);
y = r * sind(za) * sind(xa);
z = r * cosd(za);

h = sqrt(x^2 + y^2);


a = [0 0 0]; % origin
b = [x y z]; % tip
c = [x y 0];
d = [0 0 0];

xx= [a(1) b(1) c(1) d(1)];
yy= [a(2) b(2) c(2) d(2)];
zz= [a(3) b(3) c(3) d(3)];

figure(1)
plot3(xx,yy,zz)
axis([0 10 0 10 0 10])
xlabel('X');ylabel('Y');zlabel('Z');
az=60;el=42;
view([az el])
grid on
%}


%----------------------------------------------------------------------------%
							Actin = zeros(5,12);
%[Nact	Xang	Xorg	Xtip	Yang	Yorg	Ytip	Zang	Zorg	Ztip	Lact	OrO ]%
%[1		2		3		4		5		6		7		8		9		10		11		12	]%
%----------------------------------------------------------------------------%

GaSize = 5.1 / 2;	% Actin Size

% N monomers in 5 Starting Filaments
Actin(:,1) = 270;

% Length of 5 Starting Filaments
Actin(:,11) = Actin(:,1) .* GaSize; 

% Branching Angles
ArpX = 70;
ArpY = 20;
ArpZ = 20;
Actin(:,2) = ArpX;
Actin(:,5) = ArpY;
Actin(:,8) = ArpZ;

% Angle of 5 Starting Filaments
Actin(:,2) = 90;	% X angle
Actin(:,5) = 90;	% Y angle
Actin(:,8) = 0;		% Z angle

% Origin of 5 Starting Filaments
Actin(2,3) = 50;		% X origin
Actin(2,6) = 50;		% Y origin
Actin(3,3) = 50;		% X origin
Actin(3,6) = -50;		% Y origin
Actin(4,3) = -50;		% X origin
Actin(4,6) = 50;		% Y origin
Actin(5,3) = -50;		% X origin
Actin(5,6) = -50;		% Y origin

% Filament Tip Locations
x = Actin(:,11) * sind(ArpZ) * cosd(ArpX);
y = Actin(:,11) * sind(ArpZ) * sind(ArpX);
z = Actin(:,11) * cosd(ArpZ);
Actin(:,4) = x;
Actin(:,7) = y;
Actin(:,10) = z;


%----------------------------------------
% Spine Dimensions
SPYneckXY = 150;
SPYheadZN = 1000; % North Z dim
SPYheadZS = 700; % South Z dim
SPYheadX = 400; % X dim
SPYheadY = 400; % Y dim

PSDproxy = 100;
inPSD = SPYheadZN - PSDproxy;
%----------------------------------------


%----------------------------------------
% Rate Parameters
PolyRate = .15;
PolyRate0 = PolyRate;
PolyMol = 50 * 10^4;
PolyMol0 = PolyMol;
PolyDec = .0005;
PolyRev = .7;

DePolyRate = .0005;
DePolyN = 40;
DePolyStart = 200;

ArpRate = .005;
ArpRate0= ArpRate;
ArpMol = 10 * 10^4;
ArpMol0 = ArpMol;
ArpScalar = 120;
ArpDec = .001;
ArpRev = .7;

Nsteps = 8000;

% 14 possible rotational angles
Ovec = [25.7143   51.4286   77.1429  102.8572  128.5715  154.2858  180.0001...
		205.7144  231.4287  257.1430  282.8573  308.5716 334.2859  360.000];
	
Ovec = [0   27.6923   55.3846   83.0769  110.7692  138.4615  166.1538...
			193.8461  221.5384  249.2307  276.9230  304.6153  332.3076];
%----------------------------------------



%----------------------------------------
%			FIGURE SETUP
%----------------------------------------
Fh1 = FigSetup(1);
%----------------------------------------



%---------------------------------------------------------------
for nT = 1:Nsteps
%---------------------------------------------------------------

	% radial distance to spine shaft membrane
	ActinXYh = sqrt(Actin(:,4).^2 + Actin(:,7).^2);
	
	Xspi = ActinXYh > SPYneckXY;			% Logical array of filaments beyond SPYneckXY
	Yspi = ActinXYh > SPYneckXY;			% Logical array of filaments beyond SPYneckXY
	Zspi = (abs(Actin(:,10)) > SPYheadZN);	% Logical array of filaments beyond SPYheadZN

	Xhed = ActinXYh >= SPYheadX;
	Yhed = ActinXYh >= SPYheadY;
	Zhed = (abs(Actin(:,10)) <= SPYheadZN) & (abs(Actin(:,10)) >= SPYheadZS);
		
	ActMem = [Xspi Yspi Zspi Xhed Yhed Zhed];
	
 	Act0 = (Actin(:,1)>0);	% assures no negative actin values
	Actin(:,1) = Actin(:,1) .* (Actin(:,1)>0);
	% Actin(:,1) = Actin(:,1) .* (Actin(:,1)>0);
	Actin(:,1) = Actin(:,1) .* (Actin(:,10) > -20);


	NFact = numel(Actin(:,1));	% current number of actin filaments

	%-------------------------------------
	for aN=1:NFact
	%-------------------------------------
	
	rv=rand(6,1); % Generate a few random vaules from uniform{0:1}

		%---------------
		if PolyRate > rv(1)
		if ((~sum(ActMem(aN,:)) && Act0(aN)) || Zhed(aN))
		if (~Xhed(aN) && ~Yhed(aN))

			Actin(aN,1) = Actin(aN,1)+1;

			PolyRate = PolyRate - (PolyRate*PolyDec);
			% PolyMol = PolyMol - 2;
			% PolyRate = PolyRate*(PolyMol/PolyMol0);

		end
		end
		end
		%---------------

		
		%---------------
		if ( ArpRate * (Actin(aN,1) / ArpScalar ) ) > rv(2)
		if ((~sum(ActMem(aN,:)) && Act0(aN)) || Zhed(aN))
		if (~Xhed(aN) && ~Yhed(aN))

			Actin(NFact+1,1) = 4;		% create branch: add 2 subunits to new branch
			Actin(NFact+1,11) = Actin(NFact+1,1) .* GaSize;

			OrO = Actin(aN,12);			% Orientation of first monomer in mother filament

			Nmono = Actin(aN,1);		% N monomers in mother filament

			Rmono = ceil(Nmono * rand); % Random monomer along segment

			Rmono13 = mod(Rmono,13)+1;	% Get monomer repeat among the 13 rotational axis angles

			Rang = Ovec(Rmono13);		% Rotational angle of new branch
			Actin(NFact+1,12) = Rang;	% Store this rotational angle


			% New branch XYZ origin coordinates
			Ox = (Rmono*GaSize) * sind(Actin(aN,8)) * cosd(Actin(aN,2)) + Actin(aN,3);
			Oy = (Rmono*GaSize) * sind(Actin(aN,8)) * sind(Actin(aN,2)) + Actin(aN,6);
			Oz = (Rmono*GaSize) * cosd(Actin(aN,8)) + Actin(aN,9);
			Actin(NFact+1,3) = Ox;	% X origin
			Actin(NFact+1,6) = Oy;	% Y origin
			Actin(NFact+1,9) = Oz;	% Z origin
			
			% Store angle
			Actin(NFact+1,2) = Actin(aN,2) + ArpX;
			Actin(NFact+1,5) = Actin(aN,5) + ArpY;
			Actin(NFact+1,8) = Actin(NFact+1,12);

			% New branch XYZ tip coordinates
			Tx = Actin(NFact+1,11) * sind(Actin(NFact+1,8)) * cosd(Actin(NFact+1,2)) + Actin(NFact+1,3);
			Ty = Actin(NFact+1,11) * sind(Actin(NFact+1,8)) * sind(Actin(NFact+1,2)) + Actin(NFact+1,6);
			Tz = Actin(NFact+1,11) * cosd(Actin(NFact+1,8)) + Actin(NFact+1,9);
			Actin(NFact+1,4) = Tx;	% X tip
			Actin(NFact+1,7) = Ty;	% Y tip
			Actin(NFact+1,10) = Tz;	% Z tip

			ArpRate = ArpRate - (ArpRate*ArpDec);
			% ArpMol = ArpMol - 1;
			% ArpRate = ArpRate*(ArpMol/ArpMol0);

		end
		end
		end
		%---------------
		if nT == round(Nsteps*ArpRev); ArpRate=ArpRate0; end;
		if nT == round(Nsteps*PolyRev); PolyRate=PolyRate0; end;
		%---------------

		%{.
		%---------------
		if nT > DePolyStart
		if DePolyRate > rv(3)
			Actin(aN,1) = Actin(aN,1)-DePolyN;
			NullAct = (Actin(aN,1)>0);
			Actin(aN,1) = Actin(aN,1) .* NullAct;
			
			PolyRate = PolyRate + (DePolyRate*DePolyN*NullAct);
			% PolyMol = PolyMol + DePolyN*2;
			% PolyRate = PolyRate*(PolyMol/PolyMol0);
		end
		end
		%---------------
		%}



	%-------------------------------------
	end
	%-------------------------------------
	
	NoAct = find(Actin(:,1)<1);
	Actin(NoAct,:) = [];
	if numel(NoAct); 
		ArpRate = ArpRate + (ArpRate*ArpDec*numel(NoAct)); 
	end;
	
	Actin(:,11) = Actin(:,1).*GaSize;	% Length of Factin segments
	% branch XYZ tip coordinates
	Actin(:,4) = Actin(:,11) .* sind(Actin(:,8)) .* cosd(Actin(:,2)) + Actin(:,3);
	Actin(:,7) = Actin(:,11) .* sind(Actin(:,8)) .* sind(Actin(:,2)) + Actin(:,6);
	Actin(:,10) = Actin(:,11) .* cosd(Actin(:,8)) + Actin(:,9);
	
	
	
	
	%==================================================%
	%				LIVE PLOT
	%--------------------------------------------------%
	if mod(nT,500) == 0
		LivePlot(nT,Actin,inPSD,Fh1)
		disp(nT);
	end
	%--------------------------------------------------%



%---------------------------------------------------------------
end
%---------------------------------------------------------------



%{
%======================================================================%
%							FIGURE SETUP
%----------------------------------------------------------------------%
figure
scsz = get(0,'ScreenSize');
pos = [scsz(3)/3.5  scsz(4)/5.5  scsz(3)/2  scsz(4)/1.5];
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%
plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]')
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[1,1,1])
%----------------------------------------------------------------------%


%======================================================================%
%							FIGURE SETUP
%----------------------------------------------------------------------%
figure
pos = [scsz(3)/3  scsz(4)/5  scsz(3)/2  scsz(4)/1.5];
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%--------------------
ph11c = plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]');
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[1,1,1])
hold on;
ph11a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph11b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[1,1,1])
set(ph11a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph11b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
set(ph11c,'LineStyle','-','Color',[.7 .7 .7],'LineWidth',.1);
%----------------------------------------------------------------------%



%======================================================================%
%							FIGURE SETUP
%----------------------------------------------------------------------%
figure
scsz = get(0,'ScreenSize');
pos = [scsz(3)/2.8  scsz(4)/4.8  scsz(3)/2  scsz(4)/1.5];
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%
set(gca,'Color',[1,1,1])
hold on;
ph12a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph12b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[1,1,1])
set(ph12a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph12b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
%----------------------------------------------------------------------%



%======================================================================%
%							FIGURE SETUP
%----------------------------------------------------------------------%
figure
scsz = get(0,'ScreenSize');
pos = [scsz(3)/2.8  scsz(4)/4.8  scsz(3)/2  scsz(4)/1.5];
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%
set(gca,'Color',[1,1,1])
hold on;
ph12a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph12b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[1,1,1])
set(ph12a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph12b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
%----------------------------------------------------------------------%

PSDXYZ = [PSDTips(:,1) PSDTips(:,2) PSDTips(:,3)];
PSDXY = round([PSDTips(:,1) PSDTips(:,2)]);



PSDactMx = zeros(SPYheadY+100,SPYheadX+100);
for mxp = 1:numel(PSDXY(:,1))
PSDactMx(PSDXY(mxp,2)+SPYheadY+10, PSDXY(mxp,1)+SPYheadX+10) = 1;
end

ActMask=[1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1];
ActMx = convn(PSDactMx,ActMask,'same');
ActMx = (ActMx>0).*1.0;
figure(99)
scatter(PSDXY(:,1), PSDXY(:,2))
hold on
imagesc(ActMx)
colormap(bone)

%}


end
%=============================================================================%
%=============================================================================%
%							END MAIN FUNCTION
%=============================================================================%
%=============================================================================%










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
%					LIVE PLOT
%--------------------------------------------------%
function [] = LivePlot(nT,Actin,inPSD,Fh)

%--------------------
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%--------------------
figure(Fh)
subplot('Position',[.08 .05 .40 .90]), 
ph11c = plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]');
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[1,1,1])
hold on;
ph11a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph11b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[1,1,1])
set(ph11a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph11b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
set(ph11c,'LineStyle','-','Color',[.7 .7 .7],'LineWidth',.1);
hold off;
%--------------------
figure(Fh)
subplot('Position',[.55 .05 .40 .90]), 
ph12a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph12b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis([-500 500 -500 500 0 1000])
view([0 90])
grid off
set(gca,'Color',[1,1,1])
set(ph12a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph12b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
%--------------------

set(gca,'XTickLabel', sprintf('%.1f|',nT),'FontSize',14)
hold off;
%--------------------
end


















