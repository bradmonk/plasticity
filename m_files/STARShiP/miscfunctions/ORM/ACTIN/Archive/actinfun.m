function [] = actinfun()
clc; clear all; close all;

Nsteps=3000;

Actin = zeros(1,6); % [units, angle, Xloc, Yloc]
Actin(1,1) = 0;		% N monomers in Factin segment
Actin(1,2) = 90;	% Filament angle (starting angle = 90 degrees)
Actin(1,3) = 0;		% Filament start location X
Actin(1,4) = 0;		% Filament start location Y
Actin(1,5) = 0;		% Filament end location X
Actin(1,6) = 0;		% Filament end location Y


Arp = zeros(1,4); % [Actin ID, angle, Xloc, Yloc]
Arp(1,1) = 1;		% Initial Factin ID = 1
Arp(1,2) = 70;		% +70 or -70
Arp(1,3) = 0;		% Branching location X
Arp(1,4) = 0;		% Branching location Y



PolyRate = .2;
ArpRate = .001;

ActinSize = 1;
ArpAngle = 70;



for nT = 1:Nsteps
	
	NFact = numel(Actin(:,1));

	for aN=1:NFact
	
		rv=rand(2,1);
	
		if PolyRate > rv(1)
			Actin(aN,1) = Actin(aN,1)+1;
		end
	
		if PolyRate > rv(1)
			Actin(aN,5) = cosd(Actin(aN,2))*Actin(aN,1) + Actin(aN,3);
			Actin(aN,6) = sind(Actin(aN,2))*Actin(aN,1) + Actin(aN,4);
		end
		
		if ArpRate > rv(2)
			Arp(aN+1,3) = Actin(aN,5);
			Arp(aN+1,4) = Actin(aN,6);
			Actin(NFact+1,1) = 1;							% N monomers in Factin segment
			Actin(NFact+1,2) = (Actin(aN,2)) + ArpAngle * (round(exp(round(rand)))-2);	% Filament angle
			Actin(NFact+1,3) = Actin(aN,5);				% Filament start location X
			Actin(NFact+1,4) = Actin(aN,6);				% Filament start location Y
			Actin(NFact+1,5) = Actin(aN,5);				% Filament end location X
			Actin(NFact+1,6) = Actin(aN,6);				% Filament end location Y
		end
		

	end
	
	doLivePlot=0;
	if doLivePlot
	figure(1)
	scatter([Actin(:,3); Actin(:,5)], [Actin(:,4); Actin(:,6)])
	axis([-300 300 0 500])
	hold on;
	end

end

	

	figure(2)
	for pnum=1:numel(Actin(:,1))
	plot([Actin(pnum,3) Actin(pnum,5)], [Actin(pnum,4) Actin(pnum,6)])
	axis([-500 500 0 1000])
	hold on;
	end
	% drawnow;


end