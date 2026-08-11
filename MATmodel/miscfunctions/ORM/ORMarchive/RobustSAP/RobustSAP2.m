function [stepN OKGO] = RobustSAP2(Lon,Bon,Ron,Loff,Boff,Roff,dT,vars)


%-------------------------------%
% Parameters
%-------------------------------%
NSteps = vars(1);
%---
doAMPARs = vars(4);
AMPARN = 2;
amparate = 50;
G1RT = .1;
%-------------------------------%
% Presets
%-------------------------------%
hkMask=[0 1 0; 1 0 1; 0 1 0];
S=padarray(ones(vars(5)),[vars(6) vars(6)], 0);
S0 = sum(nonzeros(S));

%===============================================%
for stepN = 1:NSteps
%-----------------------------------------------%

Pmx = rand(size(S));
Soc = (S>0);
Sno = ~Soc;
hk = convn(Soc,hkMask,'same');


Pon = 1 ./ (1+exp((hk-Lon).*(-Bon)));
Pkon = Sno .* ( Ron * dT * Pon );
Son = (Pkon>Pmx);

Poff = 1 ./ (1+exp(((-hk)+Loff).*(-Boff)));
Pkoff = Soc .* ( Roff * dT * Poff );
Soff = (Pkoff>Pmx);



%====================================%
if doAMPARs
if mod(stepN, amparate) == 0
%-------------------------------%
SG1oc = S .* 0.0;
%-------------------------------%
RPOS=randi([4 12],1,AMPARN);CPOS=randi([4 12],1,AMPARN);
SG1oc(RPOS,CPOS)=1;
SG1oc(RPOS,CPOS+1)=1;
SG1oc(RPOS+1,CPOS)=1;
SG1oc(RPOS+1,CPOS+1)=1;
%-------------------------------%
Gex = Sno .* SG1oc .* G1RT + Pkon;
Gp = 1 ./ (1+exp((-1*1)*hk)) + .001;
Gpex = Gex .* Gp;
%---

Son = (Gpex>Pmx);
%-------------------------------%
end; end;
%====================================%



%-----------------------------%
	S = (Soc-Soff) + Son;
%-----------------------------%


if mod(stepN,5)==0
Ssum = sum(sum(S));

	if (Ssum <= (S0*vars(2))) 
	OKGO = 0;
	return;
	end

	if (Ssum >= (S0*vars(3))) 
	OKGO = 2;
	return;
	end

end


%-----------------------------------------------%
end % end main loop
%===============================================%


OKGO = 1;
%---------------------------------------%
end % end main function
%---------------------------------------%


