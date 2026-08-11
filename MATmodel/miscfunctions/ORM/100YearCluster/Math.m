%================================================%
close all; clear all; clc; scsz=get(0,'ScreenSize');
%------------------------------------------------%


%===============================================%
%				BASE CLUSTER STUFF
%===============================================%
Ssz=10;Psz=1;SYNsz=Ssz+Psz+Psz;
S=padarray(ones(Ssz),[Psz Psz], 0);
% S(5+(SYNsz*7))=1; S(2+(SYNsz*1))=1;
hkMask=[0 1 0; 1 0 1; 0 1 0];
%--------
dT = .01;
Lon = 2;
Loff = 2;
Bon = 40;
Boff = 1;
Ron = 30;
Roff = 4;
%--------

Soc = (S>0);
Sno = ~Soc;
hk = convn(Soc,hkMask,'same');
%---
Lhon = (hk-Lon) .* (-Bon);
Pon = 1 ./ (1+exp(Lhon));
Pkon = Sno .* ( Ron * dT * Pon );
%---
Lhoff = ((-hk)+Loff) .* (-Boff);
Poff = 1 ./ (1+exp(Lhoff));
Pkoff = Soc .* ( Roff * dT * Poff );

%--------
figure(10);
imagesc(Lhon); axis off; title('Lhon = (hk-Lon) .* (-Bon)'); 
colormap('bone'); colorbar;
%--------


%===============================================%
%				AMPAR STUFF
%===============================================%
GhkMask=[0 1 0; 1 1 1; 0 1 0]; 
AMPARN = 5;
GTon = 2;
GToff = -.5;
LTPv = 0;
%====================================%
SG1oc = zeros(SYNsz);
%-------------------------------%
% GRPOS=randi([1 (SYNsz*SYNsz)],1,AMPARN);
GRPOS=[(7+(SYNsz*2)) (4+(SYNsz*5))];
SG1oc(GRPOS)=1;
GRhk = convn(SG1oc,GhkMask,'same');
GRk=(GRhk.*GTon); GSk=(GRhk.*GToff);
%-------------------------------%

Pon = 1 ./ (1+exp(Lhon));
Pkon = Sno .* ( Ron * dT * Pon );
%---
Poff = 1 ./ (1+exp(Lhoff));
Pkoff = Soc .* ( Roff * dT * Poff );

Gon = Pkon+(GRk.*(Pkon+LTPv));
Goff = Pkoff+(GSk.*Pkoff);

%===============================================%
%					FIGURES
%===============================================%
fh1=figure(1); colormap('bone'); 
set(fh1,'OuterPosition',(scsz./[5e-3 5e-3 1.6 1.6]))
%--------
subplot('Position',[.03 .51 .29 .41])
  imagesc(hk); axis off; title('hk'); colorbar
subplot('Position',[.35 .51 .29 .41])
  imagesc((hk-Lon)); axis off; title('(hk-Lon)'); colorbar
subplot('Position',[.67 .51 .29 .41])
  imagesc(Lhon); axis off; title('Lhon = (hk-Lon).*(-Bon)'); colorbar
%--
subplot('Position',[.03 .02 .29 .41])
  imagesc(Pon); axis off; title('Pon = 1./(1+exp(Lhon))'); colorbar
subplot('Position',[.35 .02 .29 .41])
  imagesc((Ron*dT*Pon)); axis off; title('(Ron*dT*Pon)'); colorbar
subplot('Position',[.67 .02 .29 .41])
  imagesc(Pkon); axis off; title('Pkon = Sno.*(Ron*dT*Pon)'); colorbar
%===============================================%
fh3=figure(3); colormap('bone'); 
set(fh3,'OuterPosition',(scsz./[5e-3 5e-3 1.6 1.6]))
%--------
subplot('Position',[.03 .51 .29 .41])
  imagesc(-hk); axis off; title('-hk'); colorbar
subplot('Position',[.35 .51 .29 .41])
  imagesc((-hk+Loff)); axis off; title('(-hk+Lon)'); colorbar
subplot('Position',[.67 .51 .29 .41])
  imagesc(Lhoff); axis off; title('Lhoff=(-hk+Loff).*(-Boff)'); colorbar
%--
subplot('Position',[.03 .02 .29 .41])
  imagesc(Poff); axis off; title('Poff = 1./(1+exp(Lhoff))'); colorbar
subplot('Position',[.35 .02 .29 .41])
  imagesc((Roff*dT*Poff)); axis off; title('(Roff*dT*Poff)'); colorbar
subplot('Position',[.67 .02 .29 .41])
  imagesc(Pkoff); axis off; title('Pkoff = Soc.*(Roff*dT*Poff)'); colorbar
%================================================%
fh2=figure(2); colormap('bone');
set(fh2,'OuterPosition',(scsz./[3e-3 3e-3 2.1 1.6]))
%--------
subplot('Position',[.03 .53 .38 .42])
  imagesc(S.*(rand(size(S))+1)); axis off; title('S');
subplot('Position',[.51 .53 .43 .42])
  imagesc(hk); axis off; title('hk'); colorbar;
subplot('Position',[.03 .04 .43 .42])
  imagesc(GRk); axis off; title('GRk'); colorbar;
subplot('Position',[.51 .04 .43 .42])
  imagesc(S); hold on;
  imagesc(GSk); axis off; title('GSk'); colorbar;
  alpha(.2); 
%================================================%



clear all; close all; clc;
Ssz=10;Psz=3;SYNsz=Ssz+Psz+Psz;
S=padarray(ones(Ssz),[Psz Psz], 0);
% S(5+(SYNsz*7))=1; S(2+(SYNsz*1))=1;
hkMask=[0 1 0; 1 0 1; 0 1 0];
%--------
dT = .0014;
Lon = 2;
Loff = 2;
Bon = 10;
Boff = 1;
Ron = 15;
Roff = 4;
%--
TRn = dT*Ron;
TRf = dT*Roff;
%--------
Soc = (S>0);
Sno = ~Soc;
hk = convn(Soc,hkMask,'same');
%--

Ln = [-2:2] * -Bon;
Lf = [-2:2] * Boff;

Pn = TRn ./ (1+ exp(Ln)); % .* Sno
Pf = TRf ./ (1+ exp(Lf)); % .* Soc


sprintf('%0.3g \n',exp(Ln))
sprintf('%0.3g \n',exp(Lf))

sprintf('%0.3g \n',Pn)
sprintf('%0.3g \n',Pf)



%--------
clear all; close all; clc;
Nsteps = 1e4;
Sl = [1 1 1 1 1 1 1 1 1 1];
fP = .00151;
nP = .021;

scsz = get(0,'ScreenSize');
Fh3 = figure(3);
set(Fh3,'OuterPosition',(scsz./[2e-3 2e-3 3 8]))
Ph3 = imagesc(Sl); colormap('bone');

for ns = 1:Nsteps
Pmx = rand(size(Sl));
nS = (nP>Pmx);
fS = (fP>Pmx);
Soc = Sl>0;

Sl = (Soc-fS) + nS;

set(Ph3,'CData',Sl);
drawnow
end




%-------------------------------------------------------
clear all; close all; clc; scsz = get(0,'ScreenSize'); 
steps = []; step = []; hkMask = [1 0 1];
nP = .021; fP = .00151;

Nsteps = 5e4;
Nloops = 50000; 

nP = .2;
fP = .1;

S0 = [1 1 1 1 1 1 1 1 1 1];

% Fh3 = figure(3); set(Fh3,'OuterPosition',(scsz./[2e-3 2e-3 3 8]));
% Ph3 = imagesc([.1 S0 .9]); colormap('bone');

for nl = 1:Nloops
Sl = S0;
%---
for ns = 1:Nsteps
Pmx = rand(size(Sl)); Soc = Sl>0;
PnS = ((~Sl*nP)>Pmx); 
PfS = ((Sl*fP)>Pmx);

% Sl = (Soc-PfS) + PnS;
Sl = (Soc-PfS);

% set(Ph3,'CData',[.1 Sl .9]); drawnow
if sum(Sl)<10, break, end;
end
%---
steps(nl)=ns;
end
mean(steps)
%-------------------------------------------------------

% 1.53 = mean lifetime (50000 trials) for 1 of 10 being removed at fP=.1



























%-------------------------------------------------------
clear all; close all; clc; steps = []; step = [];

Nsteps = 5e4;
Nloops = 5e4; 

nP = .2;
fP = .1;

S0 = [1 1 1 1 1 1 1 1 1 1];

for nl = 1:Nloops
Sl = S0;
%---------------------------
for ns = 1:Nsteps
	
	Pmx = rand(size(Sl));
	
	PnS = ((~Sl*nP)>Pmx); 
	PfS = ((Sl*fP)>Pmx);

	Sl = ((Sl>0)-PfS);		%  + PnS;

	% set(Ph3,'CData',[.1 Sl .9]); drawnow
	if sum(Sl)<2, break, end;

end
%---------------------------
steps(nl)=ns;
end
mean(steps)
%-------------------------------------------------------
% mean lifetimes (5e5 trials) for 1:10 of 10 being removed at fP=.1
% 1: 1.5324
% 2: 2.5059
% 3: 3.699
% 4: 5.040
% 5: 6.63
% 6: 8.51
% 7: 10.89
% 8: 14.05
% 9: 18.75
%10: 28.28

log([1.5324 2.5059 3.699 5.040 6.63 8.51 10.89 14.05 18.8 28.2])

fP1 = [1.5324 2.5059 3.699 5.040 6.63 8.51 10.89 14.05 18.8 28.2]
exf = exp(.2:.1:1.1).^3



figure(88)
plot(fP1); hold on; plot(exf,'g')
subplot(1,2,1),plot(fP1);
subplot(1,2,2),plot(exf)

ss=[0.43 0.92 1.31 1.62 1.89 2.14 2.39 2.64 2.93 3.34];
exp(ss)


%-------------------------------------------------------
% [1.5324 2.5059 3.699 5.040 6.63 8.51 10.89 14.05 18.8 28.2]
% P10 = fliplr(1 - binocdf(0,1:10,.1));
close all; clear all; clc;

% probability something happens once in 1:10 trials when p(.1)

P1in10 = (1 - binocdf(0,10,.1));	% 0.6513

% The mean waiting time until the first event is 1/P(event) (SIM: 1.5324).

T1in10 = 1/P1in10;	% 1.5354

% The mean waiting time before 2 off events (SIM: 2.5059)

P1in10 = (1 - binocdf(0,10,.1));	% 0.6513
P1in9  = (1 - binocdf(0,9,.1));		% 0.6126

T1in10 = 1/P1in10;					% 3.7893
T1in9  = 1/P1in9;					% 1.6324

T2 = 1/(P1in10*P1in9);				% 2.5064

% The mean waiting time before 3 off events (SIM: 3.699) (jP=0.27)
% E(X) = p + (p*(1-p)^1) + (p*(1-p)^2) + (p*(1-p)^3) ...

P1i10 = (1 - binocdf(0,10,.1));		% 0.6513
P2i10 = (1 - binocdf(1,10,.1));		% 0.2639
P3i10 = (1 - binocdf(2,10,.1));		% 0.0702
P1i9  = (1 - binocdf(0,9,.1));		% 0.6126
P2i9  = (1 - binocdf(1,9,.1));		% 0.2252
P1i8  = (1 - binocdf(0,8,.1));		% 0.5695

% Every 3 steps, these scenarios have a chance of happening:
px1=(P1i10*P1i9*P1i8);				% 0.2272
px2=(P1i10*P2i9);					% 0.1467
px3=(P2i10*P1i8);					% 0.1503
px4=(P3i10);						% 0.0702

T3a = 1/px1;						% 4.4007
T3b = 1/px2;						% 6.6533
T3c = 1/px3;						% 6.8189
T3d = 1/px4;						% 14.2469



pdf1i10 = binopdf(1,10,.1);		% 0.6513
pdf2i10 = binopdf(2,10,.1);		% 0.2639
pdf3i10 = binopdf(3,10,.1);		% 0.0702
pdf1i9  = binopdf(1,9,.1);		% 0.6126
pdf2i9  = binopdf(2,9,.1);		% 0.2252
pdf1i8  = binopdf(1,8,.1);		% 0.5695

pdf1=(pdf1i10*pdf1i9*pdf1i8);		% 0.2272
pdf2=(pdf1i10*pdf2i9);				% 0.1467
pdf3=(pdf2i10*pdf1i8);				% 0.1503
pdf4=(pdf3i10);						% 0.0702

1/(pdf1+pdf2+pdf3+pdf4)



pf1i10 = P1i10 ./ (P1i10 + P2i10 + P3i10)
pf2i10 = P2i10 ./ (P1i10 + P2i10 + P3i10)
pf3i10 = P3i10 ./ (P1i10 + P2i10 + P3i10)
pf1i9 = P1i9 ./ (P1i9 + P2i9)
pf2i9 = P2i9 ./ (P1i9 + P2i9)




pall = P1i10 + (P1i10*(1-P1i10)^1) + (P1i10*(1-P1i10)^2) %equivalent: (1-binocdf(0,3,P1i10))



%-------------------------------------------------------
nP = .2;
fP = .1;

nPbino = fliplr(1 - binocdf(0,0:9,nP));
fPbino = 1 - binocdf(0,1:10,fP);

JP = [nPbino' fPbino' (1:10)']

pnP1st = (JP(1)) ./ (JP(1) + JP(2));
pfP1st = JP(2) ./ (JP(2) + JP(1));
%-------------------------------------------------------

fPbi = ones(1,10);
fPbi(1) = fPbino(1);

for nn=1:9

fPbi(nn+1) = fPbi(nn) * fPbino(nn+1);

end
fPbin = fPbi(10)

%-------------------------------------------------------
nP = .2;
fP = .1;

FR = fP;
R1 = 1-FR;

for nn=1:10
uN = 11-nn;

R1s = R1^uN;
F1s = -log(R1s);

FTx = 1/uN;
RTx = -log(FTx) / F1s;

RTall(nn) = RTx;
end
sum(RTall)
%-------------------------------------------------------

1/(-log(fPbino)/10)

%{


%===============================================%
%					FIGURES
%===============================================%
fh1=figure(1);
set(fh1,'OuterPosition',(scsz./[2e-3 2e-4 2.1 1.6]))
%--------
subplot('Position',[.03 .02 .45 .45])
imagesc(Pkon); axis off; colorbar
subplot('Position',[.5 .02 .45 .45])
imagesc(Pkoff); axis off
colormap('bone'); colorbar
subplot('Position',[.03 .51 .45 .45])
imagesc(Gon); axis off; colorbar
subplot('Position',[.5 .51 .45 .45])
imagesc(Goff); axis off
colormap('bone'); colorbar
%================================================%
fh2=figure(2); colormap('bone');
set(fh2,'OuterPosition',(scsz./[3e-3 3e-3 2.1 1.6]))
%--------
subplot('Position',[.03 .02 .45 .45])
imagesc(S); axis off; title('S'); colorbar;
%--
subplot('Position',[.5 .02 .45 .45])
imagesc(Sno); axis off; title('Sno'); colorbar;
%--
subplot('Position',[.03 .51 .45 .45])
imagesc(hk); axis off; title('hk'); colorbar;
%--
subplot('Position',[.5 .51 .45 .45])
imagesc(GRhk); axis off; title('GRhk'); colorbar;
%================================================%
fh3=figure(3); colormap('bone'); 
set(fh3,'OuterPosition',(scsz./[5e-3 5e-3 1.8 1.6]))
%--------
subplot('Position',[.03 .51 .30 .42])
  imagesc(hk); axis off; title('hk'); colorbar
subplot('Position',[.35 .51 .30 .42])
  imagesc((hk-Lon)); axis off; title('(hk-Lon)'); colorbar
subplot('Position',[.67 .51 .30 .42])
  imagesc(Lhon); axis off; title('Lhon=(hk-Lon).*(-Bon)'); colorbar
%--
subplot('Position',[.03 .02 .30 .42])
  imagesc(Pon); axis off; title('Pon'); colorbar
subplot('Position',[.35 .02 .30 .42])
  imagesc((Ron*dT*Pon)); axis off; title('(Ron*dT*Pon)'); colorbar
subplot('Position',[.67 .02 .30 .42])
  imagesc(Pkon); axis off; title('Pkon'); colorbar
%================================================%
%------------------------------------------------%
fh1=figure(1);
set(fh1,'OuterPosition',(scsz./[2e-3 2e-4 2.3 3]))
%--------
subplot('Position',[.03 .05 .45 .88])
imagesc(Pkon); axis off; colorbar
subplot('Position',[.5 .05 .45 .88])
imagesc(Pkoff); axis off
colormap('bone'); colorbar
%================================================%

%------------------------------------------------%
fh2=figure(2);
set(fh2,'OuterPosition',(scsz./[2e-3 3.2e-3 2.3 3]))
%--------
subplot('Position',[.03 .05 .45 .88])
imagesc(Gon); axis off; colorbar
subplot('Position',[.5 .05 .45 .88])
imagesc(Goff); axis off
colormap('bone'); colorbar
%================================================%
%}

