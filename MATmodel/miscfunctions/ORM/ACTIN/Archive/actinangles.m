clc; close all; clear all;

ArpX=70;
ArpY=70;
ArpZ=0;

Actin = zeros(5,10);
% [Nactin Xang Xorig Xtip Yang Yorig Ytip Zang Zorig Ztip]

Actin(:,1) = 50;		% N monomers in Factin segment
Actin(:,2) = 90;		% X angle

Actin(2,3) = 50;		% X origin
Actin(2,6) = 50;		% Y origin
Actin(3,3) = 50;		% X origin
Actin(3,6) = -50;		% Y origin
Actin(4,3) = -50;		% X origin
Actin(4,6) = 50;		% Y origin
Actin(5,3) = -50;		% X origin
Actin(5,6) = -50;		% Y origin



Actin(:,4) = cosd(Actin(:,2)) .* Actin(:,1) + Actin(:,3);
Actin(:,7) = sind(Actin(:,5)) .* Actin(:,1) + Actin(:,6);
Actin(:,10) = secd(Actin(:,8)) .* Actin(:,1) + Actin(:,9);



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
	
	Xmem = (abs(Actin(:,4)) > 150);
	Ymem = (abs(Actin(:,7)) > 150);
	Zmem = (abs(Actin(:,10)) > 800);

	HeadX = (abs(Actin(:,4)) >= 400);
	HeadY = (abs(Actin(:,7)) >= 400);
	HeadZ = (abs(Actin(:,10)) <= 800) & (abs(Actin(:,10)) >= 500);
		
	ActMem = [Xmem Ymem Zmem HeadX HeadY HeadZ];
	
	Act0 = (Actin(:,1)>0);
	Actin(:,1) = Actin(:,1) .* (Actin(:,1)>0);
	Actin(:,1) = Actin(:,1) .* (Actin(:,1)>0);
	
	
NFact = numel(Actin(:,1));

for aN=1:NFact
	
	rv=rand(6,1);
	
	if PolyRate > rv(1)
		if ((~sum(ActMem(aN,:)) && Act0(aN)) || HeadZ(aN))
		if (~HeadX(aN) && ~HeadY(aN))
			
		Actin(aN,1) = Actin(aN,1)+1;
		
		PolyRate = PolyRate - (PolyRate*PolyDec);
		
		end
		end
	end
	
	
	if ArpRate > rv(2)
		if ((~sum(ActMem(aN,:)) && Act0(aN)) || HeadZ(aN))
		if (~HeadX(aN) && ~HeadY(aN))
	
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

% Actin(:,1) = Actin(:,1) + 2;
Actin(:,4) = cosd(Actin(:,2)) .* Actin(:,1) + Actin(:,3);
Actin(:,7) = sind(Actin(:,5)) .* Actin(:,1) + Actin(:,6);
Actin(:,10) = secd(Actin(:,8)) .* Actin(:,1) + Actin(:,9);

end
%}


figure(10)
plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]')
axis([-500 500 -500 500 0 1000])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid on










