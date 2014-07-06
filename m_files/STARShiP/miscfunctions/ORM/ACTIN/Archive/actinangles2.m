clc; close all; clear all;

%----------------------------------------------------------------------------%
							Actin = zeros(5,10);
%[Nact	Xang	Xorg	Xtip	Yang	Yorg	Ytip	Zang	Zorg	Ztip]%
%[1		2		3		4		5		6		7		8		9		10  ]%
%----------------------------------------------------------------------------%
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

%-----------------------------------------
Actin(:,4) = cosd(Actin(:,2)) .* Actin(:,1) + Actin(:,3);
Actin(:,7) = sind(Actin(:,5)) .* Actin(:,1) + Actin(:,6);
Actin(:,10) = secd(Actin(:,8)) .* Actin(:,1) + Actin(:,9);
%-----------------------------------------

ArpX=70;
ArpY=70;
ArpZ=0;

%-----------------------------------------
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
%-----------------------------------------




%-----------------------------------------
for nT = 1:Nsteps
%-----------------------------------------
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
			Actin(NFact+1,9) = Actin(aN,10);				% Z origin

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
		if DePolyRate > rv(3)
			Actin(aN,1) = Actin(aN,1)-DePolyN;
			NullAct = (Actin(aN,1)>0);
			Actin(aN,1) = Actin(aN,1) .* NullAct;
			PolyRate = PolyRate + (DePolyRate*DePolyN*NullAct);
		end
		%---------------
		%}
		
	%-------------------------------------
	end
	%-------------------------------------

Actin(:,4) = cosd(Actin(:,2)) .* Actin(:,1) + Actin(:,3);
Actin(:,7) = sind(Actin(:,5)) .* Actin(:,1) + Actin(:,6);
Actin(:,10) = secd(Actin(:,8)) .* Actin(:,1) + Actin(:,9);

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
figure(11)
scatter3([Actin(:,4)]', [Actin(:,7)]', [Actin(:,10)]','.')
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[1,1,1])
%-----------------------------------------




%-----------------------------------------
figure(12)
plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]')
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[1,1,1])
hold on;
scatter3([Actin(:,4)]', [Actin(:,7)]', [Actin(:,10)]','b')
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
[Zrow1,Zcol1] = find(ActinTips(:,3) > 700);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < 700);
SPYTips = ActinTips(Zrow2,:);
%--------------------
figure(13)
plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]')
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[1,1,1])
hold on;
scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]','b')
hold on;
scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]','r')
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
az=-32;el=12;
view([az el])
grid off
set(gca,'Color',[1,1,1])
%-----------------------------------------
