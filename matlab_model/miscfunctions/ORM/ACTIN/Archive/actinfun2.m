function [] = actinfun2()
clc; clear all; close all;

Nsteps=4000;

Actin = zeros(3,7); % [units, angle, Xloc, Yloc]
Actin(1,1) = 200;		% N monomers in Factin segment
Actin(1,2) = 90;	% Filament angle (starting angle = 90 degrees)
Actin(1,3) = 0;		% Filament start location X
Actin(1,4) = 0;		% Filament start location Y
Actin(1,5) = -50;	% Filament end location X
Actin(1,6) = 0;		% Filament end location Y
Actin(1,7) = 0;		% Filament end location Z
Actin(2,1) = 200;		% N monomers in Factin segment
Actin(2,2) = 90;	% Filament angle (starting angle = 90 degrees)
Actin(2,3) = 50;	% Filament start location X
Actin(2,4) = 0;		% Filament start location Y
Actin(2,5) = 0;		% Filament end location X
Actin(2,6) = 0;		% Filament end location Y
Actin(1,7) = 0;		% Filament end location Z
Actin(3,1) = 200;		% N monomers in Factin segment
Actin(3,2) = 90;	% Filament angle (starting angle = 90 degrees)
Actin(3,3) = -50;	% Filament start location X
Actin(3,4) = 0;		% Filament start location Y
Actin(3,5) = 0;		% Filament end location X
Actin(3,6) = 0;		% Filament end location Y
Actin(1,7) = 50;	% Filament start location Z
Actin(1,8) = 0;		% Filament end location Z

Arp = zeros(1,4); % [Actin ID, angle, Xloc, Yloc]
Arp(1,1) = 1;		% Initial Factin ID = 1
Arp(1,2) = 70;		% +70 or -70
Arp(1,3) = 0;		% Branching location X
Arp(1,4) = 0;		% Branching location Y



PolyRate = .2;
DePolyRate = .01;
ArpRate = .0018;

ActinSize = 1;
ArpAngle = 70;

rz=180-(rand(100).*360);

for nT = 1:Nsteps
	
	NFact = numel(Actin(:,1));
	
	% if nT == 1500; keyboard; end;
		
	Xmem = (abs(Actin(:,5)) > 150);
	Ymem = (abs(Actin(:,6)) > 800);
		
 	HeadY = (abs(Actin(:,6)) <= 800) & (abs(Actin(:,6)) >= 500);
	HeadX = (abs(Actin(:,5)) >= 400);
	HeadZ = (abs(Actin(:,5)) >= 400);
		
	ActMem = (Xmem + Ymem + HeadX);
	Act0 = (Actin(:,1)>0);
		
 	% Actin(:,1) = Actin(:,1) - Xmem - Ymem + HeadY - HeadX;
		
	Actin(:,1) = Actin(:,1) .* (Actin(:,1)>0);
	Actin(:,1) = Actin(:,1) .* (Actin(:,1)>0);

	for aN=1:NFact
	
		rv=rand(3,1);
	
		if PolyRate > rv(1)
			if (~ActMem(aN) && Act0(aN)) || HeadY(aN)
			if ~HeadX(aN)			
			Actin(aN,1) = Actin(aN,1)+ActinSize;
			Actin(aN,5) = cosd(Actin(aN,2)) * Actin(aN,1) + Actin(aN,3);
			Actin(aN,6) = sind(Actin(aN,2)) * Actin(aN,1) + Actin(aN,4);
% 			Actin(aN,8) = atan2(Actin(aN,5),Actin(aN,6)) * Actin(aN,1) + Actin(aN,7);
			Actin(aN,8) = atan2(Actin(aN,5),Actin(aN,6)) * rz(aN);
			PolyRate = PolyRate - (PolyRate*.00001);
			end
			end
		end
		
		if ArpRate > rv(2)
			if (~ActMem(aN) && Act0(aN)) || HeadY(aN)
			if ~HeadX(aN)
			Arp(aN+1,3) = Actin(aN,5);
			Arp(aN+1,4) = Actin(aN,6);
			Actin(NFact+1,1) = 1;							% N monomers in Factin segment
			Actin(NFact+1,2) = ((Actin(aN,2)) + ArpAngle * (round(exp(round(rand)))-2));	% Filament angle
			Actin(NFact+1,3) = Actin(aN,5);				% Filament start location X
			Actin(NFact+1,4) = Actin(aN,6);				% Filament start location Y
			Actin(NFact+1,7) = Actin(aN,8);				% Filament start location Z
			Actin(NFact+1,5) = Actin(aN,5);				% Filament end location X
			Actin(NFact+1,6) = Actin(aN,6);				% Filament end location Y
			Actin(NFact+1,8) = Actin(aN,8);				% Filament end location Z
			ArpRate = ArpRate-(ArpRate*.001);
			end
			end
		end
		
		if DePolyRate > rv(3)
			Actin(aN,1) = Actin(aN,1)-ActinSize;
			Actin(aN,5) = cosd(Actin(aN,2))*Actin(aN,1) + Actin(aN,3);
			Actin(aN,6) = sind(Actin(aN,2))*Actin(aN,1) + Actin(aN,4);
		end
		

	end
	
	
		
		
	
	doLivePlot=0;
	if doLivePlot
	figure(1)
	scatter([Actin(:,3); Actin(:,5)], [Actin(:,4); Actin(:,6)])
	axis([-300 300 0 1000])
	hold on;
	end

end




	figure(2)
	for pnum=1:numel(Actin(:,1))
	plot([Actin(pnum,3) Actin(pnum,5)], [Actin(pnum,4) Actin(pnum,6)])
	axis([-500 500 0 1000])
	hold on;
	end
	
	
	figure(3)
	for pnum=1:numel(Actin(:,1))
	plot3([Actin(pnum,3) Actin(pnum,5)], [Actin(pnum,4) Actin(pnum,6)], [Actin(pnum,7) Actin(pnum,8)])
	axis([-500 500 0 1000])
	hold on;
	end

	
	
% 	zvals = [zeros(NFact+1,1) randn(NFact+1,1)];
% 
% 	figure(2)
% 	for pnum=1:numel(Actin(:,1))
% 	plot3([Actin(pnum,3) Actin(pnum,5)], [Actin(pnum,4) Actin(pnum,6)], [zvals(pnum,1)  zvals(pnum,2)])
% 	axis([-1000 1000 -500 500 -3 3])
% 	hold on;
% 	end
	% drawnow;


end