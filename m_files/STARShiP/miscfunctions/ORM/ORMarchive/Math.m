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

