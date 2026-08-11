function [] = ActinRateModel2()


clc; close all; clear all;
GaSize = 5.1 / 2;	% Actin Size
fID = 1; % filament ID
% 13 possible rotational angles
Ovec = [1.0 28.7 56.4 84.1 111.8 139.5 167.2 194.8 222.5 250.2 277.9 305.6 333.3];
clear actin

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

x0 = 9
y0 = 9
z0 = 9

x3 = 10
y3 = 10
z3 = 10

degPrad = 57.2958;
radPdeg = 0.0175;

xyzO = [0;
		0;
		0];

xyzT = [5;
		5;
		5];

RxA = [1 0 0; 0 cos(a) sin(a); 0 -sin(a) cos(a)];
RyB = [cos(b) 0 -sin(b); 0 1 0; sin(b) 0 cos(b)];
RzG = [cos(g) sin(g) 0; -sin(g) cos(g) 0; 0 0 1];


RdPhi = [1 0 0; 0 cos(Phi) sin(Phi); 0 -sin(Phi) cos(Phi)];
RcTta = [cos(Tta) 0 -sin(Tta); 0 1 0; sin(Tta) 0 cos(Tta)];
RbPsi = [cos(Psi) sin(Psi) 0; -sin(Psi) cos(Psi) 0; 0 0 1];

RaMx  = [a11 a12 a13; 
		 a21 a22 a23; 
		 a31 a32 a33];

% a11	= cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi)
% a12	= cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi)
% a13	= sin(psi)*sin(theta)
% a21	= -sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi)
% a22	= -sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi)
% a23	= cos(psi)*sin(theta)
% a31	= sin(theta)*sin(phi)
% a32	= -sin(theta)*cos(phi)
% a33	= cos(theta)


%{
a11	= cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi)
a12	= cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi)
a13	= sin(psi)*sin(theta)
a21	= -sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi)
a22	= -sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi)
a23	= cos(psi)*sin(theta)
a31	= sin(theta)*sin(phi)
a32	= -sin(theta)*cos(phi)
a33	= cos(theta)

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

%}


actin(1).n(:)
















%----------------------------------------------------------------------------%
							Actin = zeros(5,12);
%[Nact	Xang	Xorg	Xtip	Yang	Yorg	Ytip	Zang	Zorg	Ztip	Lact	OrO ]%
%[1		2		3		4		5		6		7		8		9		10		11		12	]%
%----------------------------------------------------------------------------%



Actin(:,1) = 270;	% N monomers in 5 Starting Filaments

Actin(:,11) = Actin(:,1) .* GaSize; % Length of 5 Starting Filaments

% Branching Angles
ArpX = 70;
ArpY = 20;
ArpZ = 70;
Actin(:,2) = ArpX;
Actin(:,5) = ArpY;
Actin(:,8) = ArpZ;

% Angle of 5 Starting Filaments
Actin(:,2) = 90;	% X angle
Actin(:,5) = 90;	% Y angle
Actin(:,8) = 0;		% Z angle

% 13 possible rotational angles
Ovec = [0   27.6923   55.3846   83.0769  110.7692  138.4615  166.1538...
			193.8461  221.5384  249.2307  276.9230  304.6153  332.3076];

% Origin of 5 Starting Filaments
Actin(2,3) = 50;		% X origin
Actin(2,6) = 50;		% Y origin
Actin(3,3) = 50;		% X origin
Actin(3,6) = -50;		% Y origin
Actin(4,3) = -50;		% X origin
Actin(4,6) = 50;		% Y origin
Actin(5,3) = -50;		% X origin
Actin(5,6) = -50;		% Y origin


% MATH - Filament Tip Locations
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

dims = [SPYneckXY SPYheadZN SPYheadZS SPYheadX SPYheadY PSDproxy inPSD];
%----------------------------------------


%----------------------------------------
% Quantitative Values (mol, volume, rate parameters)

% Units
mM = 1e-3;
uM = 1e-6;
upM = 1e3;
dnM = 1e-3;
% 1 cm^3 = 1 mL
% SI units conversion from m^3 to L are:
% 1 cm^3 = 1 mL
% To convert cubic volume to mL you may first need to convert
% by an order of magnitude; here's a reminder of cubic volume conversions:
% 1 m^3 = 1 cm^3 * .01^3;		% 1e-6
% 1 m^3 = 1 mm^3 * .001^3;		% 1e-9
% 1 cm^3 = 1 mm^3 * .1^3;		% 1e-3
% 1 cm^3 = 1 um^3 * .0001^3;	% 1e-12
% Thus, if dendritic spines have an average volume of 0.1 um^3
% that would be equivalent to X uL
%
% 0.1 um^3 * (1 cm^3 / 1e12 um^3) * (1 mL / 1 cm^3) * (1000 uL / 1 mL)
% 0.1*(1/1e12)*(1/1)*(1000/1)
% 0.1 um^3 = .1e-12 cm^3 = .1e-9 uL
%
% and 0.1 um^3 equivalent to X L
% 
% 0.1*(1/1e12)*(1/1)*(1/1000)
% 1e-16

% Constants
Nsteps = 20000;
dT = 1/1000;
mol = 6e23;		% in N

% Spine volume (0.1 um^3) in L and uL
SpyV = 1e-16;	% in L
SpyVu = .1e-9;	% in uL

% Actin Polymerization (12 N/µM*s)
Act_Na = 1e5;
% Act_Nb = 1e6;
% Act_N = Act_Na;
Act_PRnT = 1000;

% Actin Depolymerization (2 N/s)	
Act_DRnT = 200;
DePSum = 0;

% Cofilin Depoly
CofR = .0005;
CofN = 40;
CofS = 200;

% Actin Poly Math
% Cytosolic concentration of actin in cells ranges from .1 to .5 mM
% Given a spine volume of 'SpyV' we can find how many actin monomers
% are in an average spine:
%
% .1 mM (mmol/L) * (1 mol / 1000 mmol) * SpyV = 1e-17 mol
% 1e-17 mol * (6e23 units / 1 mol) = 6e3 monomer units

Act_N = .1 * (1/1000) * SpyV * 6e23; % 6e3 monomer units

% we can check our math starting with a set number of actin monomers
% and calculate the spine molarity (6e3 monomer units as an example):
% 
% 6e3 units/SpyV * (1 mol / 6e23 units) * (1000 mmol / 1 mol)
% 6e3/SpyV*(1/6e23) 

Act_mM = Act_N/SpyV * (1/6e23) * (1000/1);	% 1.6e-10 

% Act_mM = Act_N / SpyVu / mol;	% 1.6e-10 
% Act_N = Act_N / SpyVu / mol;
Act_N = Act_N;
Act_PR = 12 * Act_mM * dT;
Act_DR = 2 * dT;



% Arp Branching Rate
ArpRa = .0005;
ArpRb = .008;
ArpR = ArpRa;
ArpST = 3000;
ArpScalar = 100;
ArpDec = .001;

% Arp Branching Math
ArpN = 1e3;
ArpOn = 5;
ArpOff = 1;
Arp_uM = ArpN / SpyVu / mol;	% 1.6 - 16 uM
Arp_PR = ArpOn * Arp_uM * dT;
Arp_DR = ArpOff * dT;


%----------------------------------------


%{

% 5e8 in cells; 1e4 synapses per neuron; 5e8/1e4 = 5e4

Polymerization Rate
(+)end: .012 N/µM*ms
** (12 N/mM*ms) (12 N/µM*s)
** thus at 1 µM free ATP-actin, .012 subunits will be added to the (+)end per ms
** at .1 mM free ATP-actin, 1.2 subunits will be added to the (+)end per ms

Depolymerization Rate
(+)end: 1.4 N/s
(-)end: 0.8 N/s
** dissociation is independent of free actin concentration

To find the critical concentration (Cc) for growth we set the two rate 
equations equal to each other:

12.0/µM*s = 1.4/s
12/1.4 = 1/µM
(+)Cc = .12 µM
(-)Cc = .6 µM





Actin monomer size:
- 5.5 nm x 5.5 nm x 3.5 nm

Factin filaments
- ?-helix composed of 2 strands of subunits
- 28 subunits (14 in each strand) in 1 full 360? turn 
- 180? turn: 36 nm
- 360? turn: 72 nm
- 28 subunits per 72 nm
--------
* For 1 filament to span a 1000 nm spine would require:
-- (mean spine length: 1.0 µm or 1000 nm)
-- each monomer spans 5.1 nm
-- every 5.1 nm requires 2 monomers (cuz double helix)
-- 13.89 turns (~14 turns)
-- 388.89 actin monomers (~400 monomers)
-- 195 monomers per strand
--------

The ATP-binding cleft of actin faces the (+) end of the filament
(+)end grows
(-)end shrinks

Polymerization Rate
(+)end: ~12.0 subunits/µM*s
(-)end: ~1.3 subunits/µM*s
** thus if there is 1 µM of free ATP-Gactin then 12 subunits will be added 
to the (+)end per second and 1.3 subunits will be added to the (-)end every second

Depolymerization Rate
(+)end: ~1.4 subunits/s
(-)end: ~0.8 subunits/s
** dissociation is independent of free Gactin concentration

To find the critical concentration (Cc) for growth we set the two rate 
equations equal to each other:

12.0/µM*s = 1.4/s
12/1.4 = 1/µM
(+)Cc = .12 µM
(-)Cc = .6 µM

Thus when the free actin concentration >.12 µM filaments will grow at the (+) end 
and when >.6 µM filaments will grow at the (-) end too.


----
There are 2 isoforms of actin
- ?-actin & ?-actin

?-actin is enriched in dendritic spines and builds filament stress fibers

G-actin = monomeric "globular" actin
F-actin = filamentous actin

Each actin molecule contains a Mg2+ ion complexed with ATP or ADP

----
Cytosolic concentration of actin in cells ranges from .1 to .5 mM

Actin makes up 1-5% of all cellular protein
(this suggests there is 10 mM total proteins in cells)


A typical cell contains around 5e8 actin molecules

The average number of synapses per neuron was 1e4
http://www.ncbi.nlm.nih.gov/pubmed/2778101

5e8 / 1e4 = 5e4
A typical spine contains 5e4 actin molecules


Total spines volume averaged 0.09 ?m^3 and ranged from 0.01 to 0.38 ?m^3 (n = 133). 
The mode peak value was 0.06 ?m^3. 
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2518053/

0.09 ?m^3 = 9e-17 L
5e4 / 9e-17 L = 5.55e20
5.55 N/L * (1 L / 6e23 N) = 0.000925 M = .9 mM
.5 - .9 mM actin / spine


%}



%----------------------------------------
%		STATIC FINAL FIGURES SETUP
%----------------------------------------
Fh1 = FigSetup(1);


%----------------------------------------
%	  ANIMATED REAL-TIME FIGURE SETUP
%----------------------------------------
rot0 = [5 0];
rot1 = rot0;
azel0 = [-32 12];
sz = [2.5e-3 2.5e-3 1.6 1.2];
Fh2Live = FigSetup(2,sz);
%----------------------------------------





%===============================================================%
%						MAIN OUTER LOOP
for nT = 1:Nsteps
%---------------------------------------------------------------%

	% if nT == Act_PRnT; Act_N = Act_Nb; end;
	if nT == ArpST; ArpR = ArpRb; end;

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

	%===============================================================%
	%						MAIN INNER LOOP
	for aN=1:NFact
	%---------------------------------------------------------------%
	
	rv=rand(6,1); % Generate a few random vaules from uniform{0:1}

		%---------------
		% POLYMERIZATION
		if Act_PR > rv(1)
		if ((~sum(ActMem(aN,:)) && Act0(aN)) || Zhed(aN))
		if (~Xhed(aN) && ~Yhed(aN))

			Actin(aN,1) = Actin(aN,1) + ceil(Act_PR);
			
			Act_N = Act_N - ceil(Act_PR);

		end
		end
		end
		%---------------
		
		
		%---------------
		% BRANCHING
		if ( ArpR * (Actin(aN,1) / ArpScalar ) ) > rv(2)
		if ((~sum(ActMem(aN,:)) && Act0(aN)) || Zhed(aN))
		if (~Xhed(aN) && ~Yhed(aN))

			Actin(NFact+1,1) = 10;		% create branch: add 10 subunits to new branch
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
			
			
			% Actin(NFact+1,2) = 180 - Actin(NFact+1,12);
			% Actin(NFact+1,5) = 180 - Actin(NFact+1,12);
			% Actin(NFact+1,8) = Actin(aN,8) + cosd(ArpZ);
			
			
			% r = sqrt(x^2 + y ^2 + z^2)
			% r / sqrt(x^2 + y^2) = z
			% zang = arccos(z/r)
			% zang = arccos(r / sqrt(x^2 + y^2))
			% Actin(NFact+1,8) = arccos(Actin(NFact+1,11) / sqrt(Actin(NFact+1,2)^2 + Actin(NFact+1,5)^2))
			% 360° = 2? rad

			
			
			% New branch XYZ tip coordinates
			Tx = Actin(NFact+1,11) * sind(Actin(NFact+1,8)) * cosd(Actin(NFact+1,2)) + Actin(NFact+1,3);
			Ty = Actin(NFact+1,11) * sind(Actin(NFact+1,8)) * sind(Actin(NFact+1,2)) + Actin(NFact+1,6);
			Tz = Actin(NFact+1,11) * cosd(Actin(NFact+1,8)) + Actin(NFact+1,9);
			Actin(NFact+1,4) = Tx;	% X tip
			Actin(NFact+1,7) = Ty;	% Y tip
			Actin(NFact+1,10) = Tz;	% Z tip

			ArpR = ArpR - (ArpR*ArpDec);
			% ArpN = ArpN - 1;
			% ArpRate = ArpRate*(ArpMol/ArpMol0);

		end
		end
		end
		%---------------


		%---------------
		% DEPOLYMERIZATION
		if nT > Act_DRnT
		if Act_DR > rv(3);
			
			Actin(aN,1) = Actin(aN,1)-ceil(Act_DR) .* (Actin(aN,1)>0);
			
			DePSum = DePSum +1;
			
			Act_N = Act_N + ceil(Act_DR);
			
		end
		end
		%---------------
		
		%---------------
		if nT > CofS
		if CofR > rv(3)
			Actin(aN,1) = Actin(aN,1)-CofN;
			
			Actin(aN,1) = Actin(aN,1) .* (Actin(aN,1)>0);
			
			Act_N = Act_N + CofN;
			% Act_N = Act_N + (CofN*(Actin(aN,1)>0));
		end
		end
		%---------------
		
		
		
		%---------------
		% ADJUST RATE VALUES
		
		Act_mM = Act_N/SpyV * (1/6e23) * (1000/1);
		Act_PR = 12e3 * (Act_mM * dT);
		
		% Act_mM = Act_N / SpyVu / mol;	% .17 uM
		% Act_N = Act_mM * (1/1000) * SpyV * 6e23;
		% Act_PR = 12 * Act_uM * dT;
		
		if Act_N < 0; Act_N = 0; Act_PR = 0; end;
		%---------------
		


	%---------------------------------------------------------------%
	end
	%						MAIN INNER LOOP
	%===============================================================%
	
	
	NoAct = find(Actin(:,1)<1);
	Actin(NoAct,:) = [];
	if numel(NoAct); 
		ArpR = ArpR + (ArpR*ArpDec*numel(NoAct)); 
	end;
	
	Actin(:,11) = Actin(:,1).*GaSize;	% Length of Factin segments
	
	% MATH - branch XYZ tip coordinates
	Actin(:,4) = Actin(:,11) .* sind(Actin(:,8)) .* cosd(Actin(:,2)) + Actin(:,3);
	Actin(:,7) = Actin(:,11) .* sind(Actin(:,8)) .* sind(Actin(:,2)) + Actin(:,6);
	Actin(:,10) = Actin(:,11) .* cosd(Actin(:,8)) + Actin(:,9);
	
	
	
	%==================================================%
	%				LIVE PLOT
	%--------------------------------------------------%
	if mod(nT,500) == 0
		LivePlot(Fh2Live,nT,Actin,inPSD,rot1,azel0,dims);
		rot1 = rot1 + rot0;
		% disp(nT); disp(Act_N);disp(Act_PR);disp(Act_mM);
	end
	%--------------------------------------------------%

% 	if aN >= 12; 
% 		LivePlot(Fh2Live,nT,Actin,inPSD,rot1,azel0,dims)
% 		keyboard; 
% 	end;
	
	
	
	
ActData(nT,:) = [Act_mM, Act_PR, Act_N];
DePData(nT,:) = DePSum;
% DePSum = 0;

% if nT == 125; keyboard; end;

%---------------------------------------------------------------%
end
%						MAIN OUTER LOOP
%===============================================================%




%--------------------------------------------------%
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%--------------------------------------------------%



%============================================================%
%				OUTPUT FIGURES
%------------------------------------------------------------%

%{
%==================================================%
%			MATRIX SURFACE FIGURE
%--------------------------------------------------%
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%--------------------------------------------------%
PSDXYZ = [PSDTips(:,1) PSDTips(:,2) PSDTips(:,3)];
PSDXY = round([PSDTips(:,1) PSDTips(:,2)]);
PSDactMx = zeros(SPYheadY+100,SPYheadX+100);
for mxp = 1:numel(PSDXY(:,1))
PSDactMx(PSDXY(mxp,2)+SPYheadY+10, PSDXY(mxp,1)+SPYheadX+10) = 1;
end
ActMask=[1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1];
ActMx = convn(PSDactMx,ActMask,'same');
ActMx = (ActMx>0).*1.0;
%--------------------------------------------------%
figure
subplot('Position',[.08 .05 .40 .90]), 
imagesc(ActMx)
colormap(bone)
subplot('Position',[.55 .05 .40 .90]), 
scatter(PSDXY(:,1), PSDXY(:,2))
%--------------------------------------------------%
%==================================================%
%				DEPOLY FIGURE
%--------------------------------------------------%
sz = [2.5e-3 2.5e-3 1.2 2];
Fh2 = FigSetup(2,sz);
%--------------------------------------------------%
figure(Fh2)
subplot('Position',[.03 .05 .28 .90]),
plot(ActData(Act_PRnT:nT,1)); title('Act uM');
subplot('Position',[.35 .05 .28 .90]),
plot(ActData(Act_PRnT:nT,2)); title('Act PR');
subplot('Position',[.68 .05 .28 .90]),
plot(ActData(Act_PRnT:nT,3)); title('Act N');

figure
plot(DePData(Act_PRnT:nT)); title('Depoly Sum');
%--------------------------------------------------%



%==================================================%
%					FIGURE
%--------------------------------------------------%
figure
scsz = get(0,'ScreenSize');
pos = [scsz(3)/3.5  scsz(4)/5.5  scsz(3)/2  scsz(4)/1.5];
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])
%--------------------------------------------------%
plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]')
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[1,1,1])
%--------------------------------------------------%


%==================================================%
%					FIGURE
%--------------------------------------------------%
figure
pos = [scsz(3)/3  scsz(4)/5  scsz(3)/2  scsz(4)/1.5];
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])

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
%--------------------------------------------------%



%==================================================%
%					FIGURE
%--------------------------------------------------%
figure
scsz = get(0,'ScreenSize');
pos = [scsz(3)/2.8  scsz(4)/4.8  scsz(3)/2  scsz(4)/1.5];
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])
%--------------------------------------------------%
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
%--------------------------------------------------%



%==================================================%
%				FIGURE (Fh1)
%--------------------------------------------------%
figure
scsz = get(0,'ScreenSize');
pos = [scsz(3)/2.8  scsz(4)/4.8  scsz(3)/2  scsz(4)/1.5];
set(gcf,'OuterPosition',pos,'Color',[.9,.9,.9])
%--------------------------------------------------%
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
%--------------------------------------------------%
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
function [varargout] = LivePlot(Fh,nT,Actin,inPSD,varargin)

% keyboard

scsz = get(0,'ScreenSize');
basepos = scsz./[2.5e-3 2.5e-3 1.5 2];
baseazel = [-32 12];
baserot=[0 0];

if nargin >= 2
	rot=varargin{1};
	azel=varargin{2};
elseif nargin == 1 
	rot=varargin{1};
	azel = baseazel;
else
	rot=baserot;
	azel = baseazel;
end



%--------------------
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%--------------------
figure(Fh)
subplot('Position',[.08 .15 .45 .70]), 
ph11c = plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]');
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
view(azel)
grid off
set(gca,'Color',[1,1,1])
hold on;
ph11a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph11b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
view(azel+rot)
grid off
set(gca,'Color',[1,1,1])
set(ph11a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph11b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
set(ph11c,'LineStyle','-','Color',[.7 .7 .7],'LineWidth',.1);
hold off;
%--------------------
figure(Fh)
subplot('Position',[.6 .55 .38 .38]), 
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
set(gca,'XTickLabel', sprintf('%.1f|',nT),'FontSize',10)
hold off;
%--------------------



%==================================================%
if nargin >= 3
dims=varargin{3};
%-----------------------------------%
% Spine Dimensions
SPYneckXY = dims(1);
SPYheadZN = dims(2);
SPYheadZS = dims(3);
SPYheadX = dims(4);
SPYheadY = dims(5);
PSDproxy = dims(6);
inPSD = dims(7);
%-----------------------------------%
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%-----------------------------------%
PSDXYZ = [PSDTips(:,1) PSDTips(:,2) PSDTips(:,3)];
PSDXY = round([PSDTips(:,1) PSDTips(:,2)]);
PSDactMx = zeros(SPYheadY+100,SPYheadX+100);
for mxp = 1:numel(PSDXY(:,1))
PSDactMx(PSDXY(mxp,2)+SPYheadY+10, PSDXY(mxp,1)+SPYheadX+10) = 1;
end
ActMask=[1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1];
ActMx = convn(PSDactMx,ActMask,'same');
ActMx = (ActMx>0).*1.0;
%===================================%
%				FIGURE
%-----------------------------------%
figure(Fh)
subplot('Position',[.6 .1 .38 .38]), 
imagesc(ActMx)
colormap(bone)
% subplot('Position',[.55 .05 .40 .40]), 
% scatter(PSDXY(:,1), PSDXY(:,2),'r')
%-----------------------------------%

end
%==================================================%

varargout = {Fh};
end















