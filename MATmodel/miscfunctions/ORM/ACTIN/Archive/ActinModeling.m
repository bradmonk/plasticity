%{
r = 10;
ArpZ = 45;
ArpX = 70;
ArpY = 90 - ArpX;

x = r * sind(ArpZ) * cosd(ArpX);
y = r * sind(ArpZ) * sind(ArpX);
z = r * cosd(ArpZ);

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
clc; close all; clear all;
%----------------------------------------------------------------------------%
							Actin = zeros(1,12);
%[Nact	Xang	Xorg	Xtip	Yang	Yorg	Ytip	Zang	Zorg	Ztip	Lact	OrO ]%
%[1		2		3		4		5		6		7		8		9		10		11		12	]%
%----------------------------------------------------------------------------%

Actin(:,1) = 35;	% N monomers in Factin segments
 

GaSize = 5.1 / 2;	% G-actin Monomer (dyad) Size

Length = Actin(:,1) .* GaSize; % length of segment (hypotenuse) 

Actin(:,11) = Length;	% Length of Factin segments

% Angles
ArpX = 90;
ArpY = 90;
ArpZ = 0;
Actin(:,2) = ArpX;
Actin(:,5) = ArpY;
Actin(:,8) = ArpZ;

Xtip = Length * sind(ArpZ) * cosd(ArpX);
Ytip = Length * sind(ArpZ) * sind(ArpX);
Ztip = Length * cosd(ArpZ);

Actin(:,4) = Xtip;
Actin(:,7) = Ytip;
Actin(:,10) = Ztip;


%----------------------------------------
% 13 possible rotational angles
Ovec = [0   27.6923   55.3846   83.0769  110.7692  138.4615  166.1538...
			193.8461  221.5384  249.2307  276.9230  304.6153  332.3076];

Nsteps = 4;
ArpX = 70;
ArpY = 20;
ArpZ = 70;

ax = repmat((Actin(1) * GaSize * 1.5),1,6); ax(1:2:6) = ax(1:2:6).*-1; ax(5) = -10; 
%======================================================================%
%							FIGURE SETUP
%----------------------------------------------------------------------%
% fig21 = figure(21);
% set(21,'Units','pixels');scsz = get(0,'ScreenSize');
% pos = [scsz(3)/3  scsz(4)/5  scsz(3)/2  scsz(4)/1.5];
% set(fig21,'OuterPosition',pos)
% set(gcf,'Color',[.9,.9,.9])
% %----------------------------------------
% xx= [Actin(:,3) Actin(:,4)];
% yy= [Actin(:,6) Actin(:,7)];
% zz= [Actin(:,9) Actin(:,10)];
% plot3(xx',yy',zz','LineWidth',3)
% axis([ax])
% xlabel('X');ylabel('Y');zlabel('Z');
% az=24;el=16;
% view([az el])
% grid on
%----------------------------------------------------------------------%



%-----------------------------------------
for nT = 1:Nsteps
%-----------------------------------------


NFact = numel(Actin(:,1));	% current number of actin filaments

	%-------------------------------------
	for aN=1:NFact
	%-------------------------------------




Actin(NFact+1,1) = 20;		% create branch: add X subunits to new branch
Actin(NFact+1,11) = Actin(NFact+1,1) .* GaSize;

OrO = Actin(aN,12);			% Orientation of first monomer in mother filament

Nmono = Actin(aN,1);		% N monomers in mother filament
Rmono = ceil(Nmono * rand); % Random monomer along segment
Rmono13 = mod(Rmono,13)+1;	% Get monomer repeat among the 13 rotational axis angles

Rang = Ovec(Rmono13);		% Rotational angle of new branch
Actin(NFact+1,12) = Rang;	% Store this rotational angle


% New branch XYZ origin coordinates
Ox = (Rmono*GaSize) * (sind(Actin(aN,8)) * cosd(Actin(aN,2))) + Actin(aN,3);
Oy = (Rmono*GaSize) * (sind(Actin(aN,8)) * sind(Actin(aN,2))) + Actin(aN,6);
Oz = (Rmono*GaSize) * (cosd(Actin(aN,8))) + Actin(aN,9);
Actin(NFact+1,3) = Ox;	% X origin
Actin(NFact+1,6) = Oy;	% Y origin
Actin(NFact+1,9) = Oz;	% Z origin

% store angle
Actin(NFact+1,2) = mod([Actin(NFact+1,12) + (Actin(aN,2))],360) - 180;
Actin(NFact+1,5) = mod([Actin(NFact+1,12) + (Actin(aN,5))],360) - 180;
Actin(NFact+1,8) = ArpZ; % Actin(aN,8) + ArpZ;
% Actin(NFact+1,2) = Actin(aN,2) + ArpX;
% Actin(NFact+1,5) = Actin(aN,5) + ArpY;
% Actin(NFact+1,8) = Actin(NFact+1,12);

% New branch XYZ tip coordinates
Tx = Actin(NFact+1,11) * sind(Actin(NFact+1,8)) * cosd(Actin(NFact+1,2)) + Actin(NFact+1,3);
Ty = Actin(NFact+1,11) * sind(Actin(NFact+1,8)) * sind(Actin(NFact+1,2)) + Actin(NFact+1,6);
Tz = Actin(NFact+1,11) * cosd(Actin(NFact+1,8)) + Actin(NFact+1,9);
Actin(NFact+1,4) = Tx;	% X origin
Actin(NFact+1,7) = Ty;	% Y origin
Actin(NFact+1,10) = Tz;	% Z origin

	end
	
end






%======================================================================%
%							FIGURE SETUP
%----------------------------------------------------------------------%
figure(100)
scsz = get(0,'ScreenSize');
pos = [scsz(3)/3.5  scsz(4)/5.5  scsz(3)/2  scsz(4)/1.5];
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])
%----------------------------------------
xx= [Actin(:,3) Actin(:,4)];
yy= [Actin(:,6) Actin(:,7)];
zz= [Actin(:,9) Actin(:,10)];

plot3(xx',yy',zz','LineWidth',3)
axis([ax])
xlabel('X');ylabel('Y');zlabel('Z');
az=45;el=8;
view([az el])
grid on
%----------------------------------------------------------------------%



%======================================================================%
%							FIGURE SETUP
%----------------------------------------------------------------------%
figure
pos = [scsz(3)/3  scsz(4)/5  scsz(3)/2  scsz(4)/1.5];
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])
%----------------------------------------------------------------------%

% Spine Dimensions
SPYneckXY = 150;
SPYheadZN = 1000; % North Z dim
SPYheadZS = 700; % South Z dim
SPYheadX = 400; % X dim
SPYheadY = 400; % Y dim

PSDproxy = 100;
inPSD = SPYheadZN - PSDproxy;

ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);

%--------------------

ph11c = plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]');
%--
axis(ax)
xlabel('X');ylabel('Y');zlabel('Z');
az=24;el=16; view([az el]); grid on;
set(gcf,'Color',[1,1,1]); set(gca,'Color',[1,1,1]);
hold on;


ph11a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',70,'ob');
hold on
ph11b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',70,'or');
%--
axis(ax)
xlabel('X');ylabel('Y');zlabel('Z');
view([24 16])
grid on
set(ph11c,'LineStyle','-','Color',[.7 .7 .7],'LineWidth',2);
% set(ph11a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
% set(ph11b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
%----------------------------------------------------------------------%

disp('    Nact	Xang	Xorg	Xtip	Yang	Yorg	Ytip    Zang	Zorg	Ztip	Lact	OrO')
Actin