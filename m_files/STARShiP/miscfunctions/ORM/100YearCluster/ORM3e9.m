function [varargout] = ORM3e9()
clc, close all;

tic
%-------------------------------%
% Parameters
%-------------------------------%
% NSteps = 3e9;
NSteps = 1e7;
dT = .0014;
datarate = 100;

PSDsz = 10;
PSAsz = 4;

hkMask=[0 1 0; 1 0 1; 0 1 0];
S=padarray(ones(PSDsz),[PSAsz PSAsz], 0); S0=S;
Csave = zeros(1,(NSteps/datarate));

%-------------------------------%
%		RATE PARAMETERS
%-------------------------------%
Lon = 2;
Loff = 2;
Bon = 50;
Boff = 1;
Ron = 40;
Roff = 3;

%-------------------------------%
%		SPEED-UP PARAMETERS
%-------------------------------%
dTRon = dT*Ron;
dTRoff = dT*Roff;

%===============================================%
for stepN = 1:NSteps
%-----------------------------------------------%

Pmx = rand(size(S));
Soc = (S>0);
hk = convn(Soc,hkMask,'same');

Lhon = (hk-Lon) .* (-Bon);
Pkon = ~Soc .* ( dTRon ./ (1+exp(Lhon)) );

Lhoff = ((-hk)+Loff) .* (-Boff);
Pkoff = Soc .* ( dTRoff ./ (1+exp(Lhoff)) );

Son = (Pkon>Pmx);
Soff = (Pkoff>Pmx);

S = (Soc-Soff) + Son;

%-------
if mod(stepN, datarate) == 0
Csave(stepN/datarate) = sum(S(:));
end
%-------
if sum(S) < 1
PlotCsave(Csave,datarate,stepN,NSteps,Lon,Bon,Ron,Loff,Boff,Roff,dT)
varargout = {NSteps,stepN,S0,S,Csave};
return;
end
%-------


%-----------------------------------------------%
end % end main loop
%===============================================%

%--------------
% The original
%{
Pmx = rand(size(S));
Soc = (S>0);
Sno = ~Soc;
hk = convn(Soc,hkMask,'same');

Lhon = (hk-Lon) .* (-Bon);
Pon = 1 ./ (1+exp(Lhon));
Pkon = Sno .* ( dTRon * Pon );
Son = (Pkon>Pmx);

Lhoff = ((-hk)+Loff) .* (-Boff);
Poff = 1 ./ (1+exp(Lhoff));
Pkoff = Soc .* ( dTRoff * Poff );
Soff = (Pkoff>Pmx);

S = (Soc-Soff) + Son;
%}
%--------------

%--------------
% The fastest
%{
%-------------------------------%
%		SPEED-UP PARAMETERS
%-------------------------------%
dTRon = dT*Ron;
dTRoff = dT*Roff;

%===============================================%
for stepN = 1:NSteps
%-----------------------------------------------%

Pmx = rand(size(S));
Soc = (S>0);
hk = convn(Soc,hkMask,'same');

Lhon = (hk-Lon) .* (-Bon);
Pkon = ~Soc .* ( dTRon ./ (1+exp(Lhon)) );

Lhoff = ((-hk)+Loff) .* (-Boff);
Pkoff = Soc .* ( dTRoff ./ (1+exp(Lhoff)) );

Son = (Pkon>Pmx);
Soff = (Pkoff>Pmx);

S = (Soc-Soff) + Son;

%-------
if mod(stepN, datarate) == 0
Csave(stepN/datarate) = sum(S(:));
end
%-------
if sum(S) < 1
PlotCsave(Csave,datarate,stepN,NSteps,Lon,Bon,Ron,Loff,Boff,Roff,dT)
varargout = {NSteps,stepN,S0,S,Csave};
return;
end
%-------


%-----------------------------------------------%
end % end main loop
%===============================================%

%}
%--------------

%-----------------------------------------------%

PlotCsave(Csave,datarate,stepN,NSteps,Lon,Bon,Ron,Loff,Boff,Roff,dT)

varargout = {toc,NSteps,stepN,S0,S,Csave};

%===============================================%
end % end main function
%===============================================%


function [] = PlotCsave(Csave,datarate,stepN,NSteps,Lon,Bon,Ron,Loff,Boff,Roff,dT)
%-------------------------------------------%

disp(['dT: ' num2str(dT)])
disp(['Lon: ' num2str(Lon)])
disp(['Bon: ' num2str(Bon)])
disp(['Ron: ' num2str(Ron)])
disp(['Loff: ' num2str(Loff)])
disp(['Boff: ' num2str(Boff)])
disp(['Roff: ' num2str(Roff)])

s = 1;
m = s*60;
h = m*60;
d = h*24;
y = d*365;

sec = stepN/s;
min = stepN/m;
hour = stepN/h;
day = stepN/d;
year = stepN/y;

fsec = floor(sec);
fmin = floor(min);
fhour = floor(hour);
fday = floor(day);
fyear = floor(year);

rsec = sec - fsec;
rmin = min - fmin;
rhour = hour - fhour;
rday = day - fday;
ryear = year - fyear;

lyear = fyear;
lday = floor(ryear*y/d);
lhour = floor(rday*d/h);
lmin = floor(rhour*h/m);
lsec = floor(rmin*m/s);

disp([' '])
disp(['Cluster lifetime (of ' num2str(NSteps) ' steps requested)'])
disp(['seconds: ' num2str(sec)])
disp(['minutes: ' num2str(min)])
disp(['hous: ' num2str(hour)])
disp(['days: ' num2str(day)])
disp(['years: ' num2str(year)])

disp(['totaling: ' num2str(lyear) ' years ' num2str(lday) ' days '...
num2str(lhour) ' hours ' num2str(lmin) ' minutes ' num2str(lsec) ' seconds '])
disp([' '])


figure(2); hold on;
plot(Csave)
set(get(gca,'XLabel'),'String','Steps')
set(get(gca,'YLabel'),'String','SAP')
set(gca,'YLim',[0 (max(Csave)+10)])
xt = (get(gca,'XTick'))*datarate;
set(gca,'XTickLabel', sprintf('%.0f|',xt));
SAPtitle = title(['  Cluster lasted '...
num2str(stepN) ' of ' num2str(NSteps) ' steps. \bullet '...
num2str(lyear) ' yr ' num2str(lday) ' day '...
num2str(lhour) ' hr ' num2str(lmin) ' min '...
num2str(lsec) ' sec']);
set(SAPtitle, 'FontSize', 12);

end


