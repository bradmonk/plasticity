function [stepN OKGO] = RobustSAP4(Lon,Bon,Ron,Loff,Boff,Roff,dT,vars,glu)


%-------------------------------%
% Parameters
%-------------------------------%
NSteps = vars(1);
%---
doAMPARs = glu(1);
AMPARN = glu(2);
amparate = glu(3);
G1RT = glu(4);
G1ST = glu(5);
%-------------------------------%
% Presets
%-------------------------------%
SYNsz = vars(4)+(vars(5)*2);
hkMask=[0 1 0; 1 0 1; 0 1 0];
GhkMask=[0 1 0; 1 1 1; 0 1 0];
S=padarray(ones(vars(4)),[vars(5) vars(5)], 0);
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
if doAMPARs; %if mod(stepN, amparate) == 0;
%-------------------------------%
SG1oc = zeros(SYNsz);
%-------------------------------%
GRPOS=randi([1 (SYNsz*SYNsz)],1,AMPARN);
SG1oc(GRPOS)=1;
GRhk = convn(SG1oc,GhkMask,'same');
GSk=(GRhk.*G1ST);GRk=(GRhk.*G1RT);
%-------------------------------%
Gon = Pkon+(GRk.*Pkon);
Goff = Pkoff+(GSk.*Pkoff);

Son = (Gon>Pmx);
Soff = (Goff>Pmx);

%-------------------------------%
end; %end;
%====================================%


%{
%====================================%
if doAMPARs; if mod(stepN, amparate) == 0;
%-------------------------------%
SG1oc = S .* 0.0;
%-------------------------------%
RPOS=randi([2 14],1,AMPARN);
CPOS=randi([2 14],1,1);
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
%}


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


