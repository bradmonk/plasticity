function [stepN OKGO] = RobustGR(Lon,Bon,Ron,Loff,Boff,Roff,dT,vars,glu,GTon,GToff)


%-------------------------------%
% Parameters
%-------------------------------%
NSteps = vars(1);
%---
doAMPARs = glu(1);
AMPARN = glu(2);
LTPv = 1e-6;
%-------------------------------%
% Presets
%-------------------------------%
SYNsz = vars(4)+(vars(5)*2);
GhkMask=[0 1 0; 1 1 1; 0 1 0];
S=padarray(ones(vars(4)),[vars(5) vars(5)], 0);
S0 = sum(nonzeros(S));

%----------------------------------------------------%
%				Mask Setup
%----------------------------------------------------%
hkMask=[0 1 0; 1 0 1; 0 1 0];

doGaussianMask = vars(6);
%-------------------------------%
if doGaussianMask
%-------------------------------%
A = 2;	x0=0; y0=0;	sx = .2; sy = .2; rx=sx; ry=sy;	res=2;

t = 0;
a = cos(t)^2/2/sx^2 + sin(t)^2/2/sy^2;
b = -sin(2*t)/4/sx^2 + sin(2*t)/4/sy^2 ;
c = sin(t)^2/2/sx^2 + cos(t)^2/2/sy^2;

[X, Y] = meshgrid((-sx*res):(rx):(sx*res), (-sy*res):(ry):(sy*res));
Z = A*exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;

hkMask=Z;
hk = convn(S,hkMask,'same');
hkor = hk(PSAsz+1,PSAsz+1);
LBR(1) = hkor-sqrt(A); LBR(2) = hkor+sqrt(A);
dT=dT*2;
end
%----------------------------------------------------%

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
GRk=(GRhk.*GTon); GSk=(GRhk.*GToff);
%-------------------------------%
Gon = Pkon+(GRk.*(Pkon+LTPv));
Goff = Pkoff+(GSk.*Pkoff);

Son = (Gon>Pmx);
Soff = (Goff>Pmx);

%-------------------------------%
end; %end;
%====================================%



%{

% scsz = get(0,'ScreenSize');
% Fh3 = figure(3);
% Ph3 = imagesc(S);
% set(Ph3,'CData',S);
% drawnow

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


