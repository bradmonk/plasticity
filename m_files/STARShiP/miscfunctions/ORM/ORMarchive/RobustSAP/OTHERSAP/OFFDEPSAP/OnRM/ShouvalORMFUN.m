function [varargout] = ShouvalORMFUN(QRvals,params)
clc, close all

%-------------------------------%
% Parameters
%-------------------------------%
NSteps = params(1);
dT = params(2);

PSDsz = 8;
PSAsz = 4;

hkMask=[0 1 0; 1 0 1; 0 1 0];
SYNsz = PSDsz+(PSAsz*2);
S=padarray(ones(PSDsz),[PSAsz PSAsz], 0);
S0=S;

%---------------------------
Lon = QRvals(5);	%(lo=fst)
Loff = QRvals(6);	%(hi=fst)
%-- INDY (low=fast)
Bon = QRvals(1);
Boff = QRvals(3);
%-- DEP (high=fast)
Ron = QRvals(2);
Roff = QRvals(4);
%---------------------------


%===============================================%
for stepN = 1:NSteps
%-----------------------------------------------%

Pmx = rand(size(S));
Soc = (S>0);
Sno = ~Soc;
hk = convn(Soc,hkMask,'same');

Lhon = (hk-Lon) .* (-Bon);
Pon = 1 ./ (1+exp(Lhon));
Pkon = Sno .* ( Ron * dT * Pon );
Son = (Pkon>Pmx);

Lhoff = ((-hk)+Loff) .* (-Boff);
Poff = 1 ./ (1+exp(Lhoff));
Pkoff = Soc .* ( Roff * dT * Poff );
Soff = (Pkoff>Pmx);

%--------------------%
S = (Soc-Soff) + Son;
%--------------------%

if sum(S) < 1
varargout = {NSteps,stepN,S0,S};
return;
end


%-----------------------------------------------%
end % end main loop
%===============================================%

varargout = {NSteps,stepN,S0,S};

%---------------------------------------%
end % end main function
%---------------------------------------%

