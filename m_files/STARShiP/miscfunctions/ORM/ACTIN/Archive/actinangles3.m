clc; close all; clear all;

ArpX=70;
ArpY=20;
ArpZ=70;

Actin = zeros(5,10);
% [Nactin Xang Xorig Xtip Yang Yorig Ytip Zang Zorig Ztip]

Actin(:,1) = 50;		% N monomers in Factin segment
Actin(:,2) = 45;		% X angle
Actin(:,5) = 45;		% Y angle
Actin(:,8) = 90;		% Z angle

Actin(2,3) = 50;		% X origin
Actin(2,6) = 50;		% Y origin
Actin(3,3) = 50;		% X origin
Actin(3,6) = -50;		% Y origin
Actin(4,3) = -50;		% X origin
Actin(4,6) = 50;		% Y origin
Actin(5,3) = -50;		% X origin
Actin(5,6) = -50;		% Y origin


Actin(:,4) = Actin(:,1) .* sind(Actin(:,8)) .* cosd(Actin(:,2));
Actin(:,7) = Actin(:,1) .* sind(Actin(:,8)) .* sind(Actin(:,2));
Actin(:,10) = Actin(:,1) .* cosd(Actin(:,8));


PolyRate = .2;
PolyRate0 = PolyRate;
PolyDec = .0005;
PolyRev = .6;

DePolyRate = .001;
DePolyN = 50;

ArpRate = .005;
ArpRate0= ArpRate;
ArpDec = .001;
ArpRev = .6;

Nsteps = 5000;

%{.
for nT = 1:Nsteps
		
	ActinXYh = sqrt(Actin(:,4).^2 + Actin(:,7).^2);
	Xspi = ActinXYh > 150;
	Yspi = ActinXYh > 150;
	Zspi = (abs(Actin(:,10)) > 800);

	Xhed = ActinXYh >= 400;
	Yhed = ActinXYh >= 400;
	Zhed = (abs(Actin(:,10)) <= 800) & (abs(Actin(:,10)) >= 500);
		
	ActMem = [Xspi Yspi Zspi Xhed Yhed Zhed];
	
	Act0 = (Actin(:,1)>0);
	Actin(:,1) = Actin(:,1) .* (Actin(:,1)>0);
	Actin(:,1) = Actin(:,1) .* (Actin(:,1)>0);
	
	
NFact = numel(Actin(:,1));

for aN=1:NFact
	
	rv=rand(6,1);
	
	if PolyRate > rv(1)
		if ((~sum(ActMem(aN,:)) && Act0(aN)) || Zhed(aN))
		if (~Xhed(aN) && ~Yhed(aN))
			
		Actin(aN,1) = Actin(aN,1)+1;
		
		PolyRate = PolyRate - (PolyRate*PolyDec);
		
		end
		end
	end
	
	
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
		Actin(NFact+1,9) = Actin(aN,10);				% Z origin
		
		ArpRate = ArpRate - (ArpRate*ArpDec);
		
		end
		end
	end
	if nT == round(Nsteps*ArpRev); ArpRate=ArpRate0; end;
	if nT == round(Nsteps*PolyRev); PolyRate=PolyRate0; end;
	
	%{.
	if DePolyRate > rv(3)
		Actin(aN,1) = Actin(aN,1)-DePolyN;
		NullAct = (Actin(aN,1)>0);
		Actin(aN,1) = Actin(aN,1) .* NullAct;
		PolyRate = PolyRate + (DePolyRate*DePolyN*NullAct);
	end
	%}

end

Actin(:,4) = Actin(:,1) .* sind(Actin(:,8)) .* cosd(Actin(:,2));
Actin(:,7) = Actin(:,1) .* sind(Actin(:,8)) .* sind(Actin(:,2));
Actin(:,10) = Actin(:,1) .* cosd(Actin(:,8));

end
%}


figure(10)
plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]')
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[.9,.9,.9])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[.9,.9,.9])









