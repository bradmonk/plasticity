clc; close all; clear all;

%----------------------------------------------------------------------------%
							Actin = zeros(5,11);
%[Nact	Xang	Xorg	Xtip	Yang	Yorg	Ytip	Zang	Zorg	Ztip	Lact]%
%[1		2		3		4		5		6		7		8		9		10		11	]%
%----------------------------------------------------------------------------%

% Actin Size
actsz = 5.1 / 2;

Actin(:,1) = 10;					% N monomers in Factin segments
Actin(:,11) = Actin(:,1).*actsz;	% Length of Factin segments

Actin(:,2) = 90;		% X angle
Actin(2,3) = 50;		% X origin
Actin(2,6) = 50;		% Y origin
Actin(3,3) = 50;		% X origin
Actin(3,6) = -50;		% Y origin
Actin(4,3) = -50;		% X origin
Actin(4,6) = 50;		% Y origin
Actin(5,3) = -50;		% X origin
Actin(5,6) = -50;		% Y origin

%-----------------------------------------
Actin(:,4) = cosd(Actin(:,2)) .* Actin(:,11) + Actin(:,3);
Actin(:,7) = sind(Actin(:,5)) .* Actin(:,11) + Actin(:,6);
Actin(:,10) = secd(Actin(:,8)) .* Actin(:,11) + Actin(:,9);
%-----------------------------------------

% Branching Angles
ArpX=70;
ArpY=70;
ArpZ=0;

% Spine dimensions
SPYneckXY = 150;
SPYheadZN = 1000; % North Z dim
SPYheadZS = 700; % South Z dim
SPYheadX = 400; % X dim
SPYheadY = 400; % Y dim

PSDproxy = 100;
inPSD = SPYheadZN - PSDproxy;


%----------------------------------------
PolyRate = .15;
PolyRate0 = PolyRate;
PolyDec = .0008;
PolyRev = .6;

DePolyRate = .001;
DePolyN = 50;
DePolyStart = 200;

ArpRate = .005;
ArpRate0= ArpRate;
ArpDec = .001;
ArpRev = .6;

Nsteps = 7000;
%----------------------------------------




%-----------------------------------------
for nT = 1:Nsteps
%-----------------------------------------

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
	Actin(:,1) = Actin(:,1) .* (Actin(:,1)>0);
	
	
	NFact = numel(Actin(:,1));	% current number of actin filaments
	
	%-------------------------------------
	for aN=1:NFact
	%-------------------------------------

		rv=rand(6,1);

		%---------------
		if PolyRate > rv(1)
			if ((~sum(ActMem(aN,:)) && Act0(aN)) || Zhed(aN))
			if (~Xhed(aN) && ~Yhed(aN))

			Actin(aN,1) = Actin(aN,1)+1;

			PolyRate = PolyRate - (PolyRate*PolyDec);

			end
			end
		end
		%---------------


		%---------------
		if ArpRate > rv(2)
			if ((~sum(ActMem(aN,:)) && Act0(aN)) || Zhed(aN))
			if (~Xhed(aN) && ~Yhed(aN))

			Actin(NFact+1,1) = 1;
			% X coordinates
			Actin(NFact+1,2) = Actin(aN,2)+ ArpX * (round(exp(round(rand)))-2);	% X angle
			Actin(NFact+1,3) = Actin(aN,4);					% X origin

			% Y coordinates
			Actin(NFact+1,5) = Actin(aN,5)+ ArpY * (round(exp(round(rand)))-2);	% Y angle
			Actin(NFact+1,6) = Actin(aN,7);					% Y origin

			% Z coordinates
			Actin(NFact+1,8) = Actin(aN,8)+ ArpZ;	% Z angle
			Actin(NFact+1,9) = Actin(aN,10);		% Z origin

			ArpRate = ArpRate - (ArpRate*ArpDec);

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
		end
		end
		%---------------
		%}
		
	%-------------------------------------
	end
	%-------------------------------------

Actin(:,11) = Actin(:,1).*actsz;	% Length of Factin segments
Actin(:,4) = cosd(Actin(:,2)) .* Actin(:,11) + Actin(:,3);
Actin(:,7) = sind(Actin(:,5)) .* Actin(:,11) + Actin(:,6);
Actin(:,10) = secd(Actin(:,8)) .* Actin(:,11) + Actin(:,9);

%-----------------------------------------
end
%-----------------------------------------



%-----------------------------------------
figure(10)
plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]')
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[1,1,1])
%-----------------------------------------


%-----------------------------------------
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%--------------------
figure(11)
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
%-----------------------------------------

%-----------------------------------------
figure(12)
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
%-----------------------------------------





%{
r = 10;
za = 45;
xa = 70;
ya = 90 - xa;

x = r * sind(za) * cosd(xa);
y = r * sind(za) * sind(xa);
z = r * cosd(za);

h = sqrt(x^2 + y^2);


a = [0 0 0];
b = [x y z];
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