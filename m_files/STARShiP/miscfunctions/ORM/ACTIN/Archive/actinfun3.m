function [] = actinfun3()
clc; clear all; close all;

Nsteps=4000;

Actin = zeros(1,10);
Actin(1,1) = 200;		% N monomers in Factin segment
Actin(1,2) = 90;		% Filament X angle
Actin(1,3) = 0;			% Filament Y angle
Actin(1,4) = 0;			% Filament Z angle
Actin(1,5) = 0;			% Filament X origin
Actin(1,6) = 0;			% Filament Y origin
Actin(1,7) = 50;		% Filament Z origin
Actin(1,8) = 0;			% Filament X tip
Actin(1,9) = 0;			% Filament Y tip
Actin(1,10) = 0;		% Filament Z tip

Arp = zeros(1,4); % [Actin ID, angle, Xloc, Yloc]
Arp(1,1) = 1;		% Initial Factin ID = 1
Arp(1,2) = 70;		% +70 or -70
Arp(1,3) = 0;		% Branching location X
Arp(1,4) = 0;		% Branching location Y
Arp(1,5) = 0;		% Branching location Z



PolyRate = .2;
DePolyRate = .01;
ArpRate = .0018;

ActinSize = 1;
ArpAngle = 70;

rz=180-(rand(100).*360);

for nT = 1:Nsteps
	
	NFact = numel(Actin(:,1));
	
	Xmem = (abs(Actin(:,8)) > 150);
	Ymem = (abs(Actin(:,9)) > 150);
	Zmem = (abs(Actin(:,10)) > 800);

	HeadX = (abs(Actin(:,8)) >= 400);
	HeadY = (abs(Actin(:,9)) >= 400);
	HeadZ = (abs(Actin(:,10)) <= 800) & (abs(Actin(:,10)) >= 500);
		
	ActMem = [Xmem Ymem Zmem HeadX HeadY HeadZ];
	
	keyboard
	
	
	Act0 = (Actin(:,1)>0);
	Actin(:,1) = Actin(:,1) .* (Actin(:,1)>0);
	Actin(:,1) = Actin(:,1) .* (Actin(:,1)>0);

	for aN=1:NFact
	
		rv=rand(6,1);
	
		if PolyRate > rv(1)
			if ((~sum(ActMem(aN,:)) && Act0(aN)) || HeadZ(aN))
			if (~HeadX(aN) && ~HeadY(aN))

			Actin(aN,1) = Actin(aN,1)+ActinSize;
			PolyRate = PolyRate - (PolyRate*.00001);
			
			end
			end
		end
		
			
		
		if ArpRate > rv(2)
			if ((~sum(ActMem(aN,:)) && Act0(aN)) || HeadZ(aN))
			if (~HeadX(aN) && ~HeadY(aN))
			Arp(aN+1,3) = Actin(aN,8); Arp(aN+1,4) = Actin(aN,9); Arp(aN+1,5) = Actin(aN,10);
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
	
	Actin(aN,5) = cosd(Actin(aN,2)) * Actin(aN,1) + Actin(aN,3);
	Actin(aN,6) = sind(Actin(aN,2)) * Actin(aN,1) + Actin(aN,4);
 	Actin(aN,8) = atan2(Actin(aN,5),Actin(aN,6)) * Actin(aN,1) + Actin(aN,7);
	
		
		
	% if nT == 1500; keyboard; end;
	doLivePlot=0;
	if doLivePlot
	figure(1)
	scatter([Actin(:,3); Actin(:,5)], [Actin(:,4); Actin(:,6)])
	axis([-300 300 0 1000])
	hold on;
	end

end



	% 2D FIGURE
	figure(2)
	for pnum=1:numel(Actin(:,1))
	plot([Actin(pnum,3) Actin(pnum,5)], [Actin(pnum,4) Actin(pnum,6)])
	axis([-500 500 0 1000])
	hold on;
	end
	
	%3D FIGURE
	figure(3)
	for pnum=1:numel(Actin(:,1))
	plot3([Actin(pnum,3) Actin(pnum,5)], [Actin(pnum,4) Actin(pnum,6)], [Actin(pnum,7) Actin(pnum,8)])
	axis([-500 500 0 1000])
	hold on;
	end


	
	
	
%{
% 	zvals = [zeros(NFact+1,1) randn(NFact+1,1)];
% 
% 	figure(2)
% 	for pnum=1:numel(Actin(:,1))
% 	plot3([Actin(pnum,3) Actin(pnum,5)], [Actin(pnum,4) Actin(pnum,6)], [zvals(pnum,1)  zvals(pnum,2)])
% 	axis([-1000 1000 -500 500 -3 3])
% 	hold on;
% 	end
	% drawnow;
%}
end