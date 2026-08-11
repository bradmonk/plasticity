function [varargout] = BMSurfD()

Ndots = 40;					% Number of particles
Nsteps = 200;				% Number of steps
DiffRate = 0.1;				% µm²/s
DiffRate2 = 0.01;			% µm²/s
Scale = 1/10;				% µm scale of model
dt = 1000/1000;				% ms time step

SzX = 3;					% µm length of X-dim
SzY = 6;					% µm length of Y-dim
PSD1Sz = .4;				% µm² size of synapse 1
PSD2Sz = .4;				% µm² size of synapse 2
PeriPSD1Sz = .3;			% µm² size of perisynapse 1
PeriPSD2Sz = .3;			% µm² size of perisynapse 2


%---------------------------------------------------------%
% BASE DIFFUSION RATES EQUATIONS

sc = Scale;					% scale of model:life
dm = 2;                     % dimensions
D = DiffRate*dt/sc;			% Diffusion Rate ES (D = L² / 2d*t)
Dp = DiffRate2*dt/sc;		% Diffusion Rate PSD
Dr = D/Dp;					% Ratio of D:Dp (1/Ls)^2;
Dn = D/Dr;					% new D after scaling L
k = sqrt(dm*D);				% stdv of D's step size distrib
MSD = 2*dm*D;               % mean squared displacement
L = sqrt(2*dm*D);           % average diagonal (2D) step size
Lx = L/sqrt(2);             % average linear (1D) step size
Ls = 1/sqrt(Dr);			% scales Lx values for Dn
PSD1 = 1/sqrt(Dr);			% (PSD1) D Lx-Scalar (Ls scales Lx so Dn_PSD1 = Ds:PSD1Ds)
PSD2 = 1/sqrt(Dr);			% (PSD2) D Lx-Scalar (Ls scales Lx so Dn_PSD2 = Ds:PSD2Ds)


XYs = ones(2,Ndots);
XYl = ones(2,Ndots);



% DENDRITIC FIELD
fsizeX = fix(SzX/sc);
fsizeY = fix(SzY/sc);
PSD1size = fix(PSD1Sz/sc);
PSD2size = fix(PSD2Sz/sc);
periPSD1size = fix(PeriPSD1Sz/sc);
periPSD2size = fix(PeriPSD2Sz/sc);

[Frow Fcol row2L col2L row2F col2F row1L col1L row1F col1F PSDfield...
PSDLOC PSDSZE] = FIELDGEN(fsizeX, fsizeY, PSD1size, PSD2size, periPSD1size, periPSD2size);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stepN = 1;
for Nt = 1:Nsteps 
%%%%%%%%%%%%%%%%%%
    
  
    %-------------------------------%
    %     STEP DIRECTION & SIZE
    %-------------------------------%
    XYs = STEPGEN(Ndots, k);


	%-------------------------------%
    %     STEP DIRECTION & SIZE
    %-------------------------------%
	[XYl] = MOVEDOTS(Ndots,XYs,XYl,PSD1,PSD2,...
	Frow,Fcol,row1F,col1F,row1L,col1L,row2F,col2F,row2L,col2L);


	%-------------------------------%
    %     LIVE DIFFUSION PLOT
    %-------------------------------%
	DOTPLOT(XYl, Fcol, PSDLOC, PSDSZE);
	pause(.01)

%%%%%%%%%%%%%%%%%%
stepN = stepN+1;
end
%%%%%%%%%%%%%%%%%%


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%-------------------------------------------%
% STEP SIZE GENERATOR
%-------------------------------------------%
function XYs = STEPGEN(Ndots, k)
    XYs = (k * randn(2,Ndots));
end


%-------------------------------------------%
% MOVE PARTICLES MAIN FUNCTION
%-------------------------------------------%
function [XYl] = MOVEDOTS(Ndots,XYs,XYl,PSD1,PSD2,...
	Frow,Fcol,row1F,col1F,row1L,col1L,row2F,col2F,row2L,col2L)


    for j = 1:Ndots
        
       if XYl(1,j)>=col1F && XYl(1,j)<=col1L &&...
          XYl(2,j)>=row1F && XYl(2,j)<=row1L 

            XYs(:,j) = XYs(:,j)*(PSD1);
            
        elseif XYl(1,j)>=col1F && XYl(1,j)<=col1L &&...
          XYl(2,j)>=row2L && XYl(2,j)<=row2F

            XYs(:,j) = XYs(:,j)*(PSD2);
                
        elseif XYl(1,j)>Fcol
            XYl(1,j) = Fcol;
        elseif XYl(1,j)<0
            XYl(1,j) = 0;
        elseif XYl(2,j)>Frow
            XYl(2,j) = Frow/2;
        elseif XYl(2,j)<-Frow/2
            XYl(2,j) = -Frow/2;
	   end	   

        XYl(:,j) = XYl(:,j)+XYs(:,j);
    end
	%===========================%

end


%-------------------------------------------%
% FIELD MAP FOR DENDRITE AREA AND PSD AREAS
%-------------------------------------------%
function [Frow Fcol row2L col2L row2F col2F row1L col1L row1F col1F PSDfield...
PSDLOC PSDSZE] = FIELDGEN(fsizeX, fsizeY, PSD1size, PSD2size, periPSD1size, periPSD2size)

Syn1Size = PSD1size+(periPSD1size*2);
Syn2Size = PSD2size+(periPSD2size*2);

PSD1padX = fix((fsizeX - Syn1Size)/2);
PSD1padY = fix(((fsizeY - Syn1Size)/2)/2);
PSD2padX = fix((fsizeX - Syn2Size)/2);
PSD2padY = fix(((fsizeY - Syn2Size)/2)/2);

fPSD1 = ones(PSD1size);						% PSD1 SIZE
fPSD2 = ones(PSD2size);						% PSD2 SIZE
fPSD1 = padarray(fPSD1,[periPSD1size periPSD1size], 2);	% PSD1 SIZE
fPSD2 = padarray(fPSD2,[periPSD2size periPSD2size], 2);	% PSD2 SIZE
pfPSD1 = padarray(fPSD1,[PSD1padY PSD1padX], 0);	% PAD PSD1 [Y-rows X-cols]
pfPSD2 = padarray(fPSD2,[PSD2padY PSD2padX], 0);	% PAD PSD2 [Y-rows X-cols]

[pfPSD1R pfPSD1C] = size(pfPSD1);
[pfPSD2R pfPSD2C] = size(pfPSD2);


if pfPSD1C > pfPSD2C
		pfPSD2 = padarray(pfPSD2,[0 (pfPSD1C - pfPSD2C)], 0,'post');
elseif pfPSD1C < pfPSD2C
		pfPSD1 = padarray(pfPSD1,[0 (pfPSD2C - pfPSD1C)], 0,'post');
end



PSDfield = cat(1, pfPSD1, pfPSD2);					% CONCAT PSD FIELDS

% figure(5); imagesc(pfPSD1)
% figure(5); imagesc(pfPSD2)
% figure(5); imagesc(PSDfield)

[Frow Fcol] = size(PSDfield);           % Fcol: #cols | Frow: #cols x2
[row1F col1F] = find(pfPSD1,1,'first'); % PSD1 1st row & col
[row1L col1L] = find(pfPSD1,1,'last');  % PSD1 Lst row & col
[row2F col2F] = find(pfPSD2,1,'first'); % PSD2 1st row & col
[row2L col2L] = find(pfPSD2,1,'last');  % PSD2 Lst row & col
row2F = row2F*-1; col2F = col2F*-1;     % -PSD2 1st row & col
row2L = row2L*-1; col2L = col2L*-1;     % -PSD2 Lst row & col
P1TL = [row1F col1F]';					% PSD1 top R location
P1BR = [row1L col1L]';					% PSD1 btm L location
P2TL = [row2F col2F]';					% PSD1 top R location
P2BR = [row2L col2L]';					% PSD1 btm L location
PSDLOC = [P1TL P1BR P2TL P2BR];			% ALL PSD TL BR locations
PSDSZE = [PSD1size PSD2size; periPSD1size periPSD2size];



end


%-------------------------------------------%
% FIELD MAP FOR DENDRITE AREA AND PSD AREAS
%-------------------------------------------%
function [] = DOTPLOT(XYl, Fcol, PSDLOC, PSDSZE)
%-------------------------------------------%


%=================================%
%       SURFACE PARAMETERS
%---------------------------------%
PSDLOC1 = PSDLOC;
PSDLOC1 = [PSDLOC1(2,:); PSDLOC1(1,:)];

SYN1SZ = (PSDSZE(1,1) + (PSDSZE(2,1)*2))-1;
SYN2SZ = (PSDSZE(1,2) + (PSDSZE(2,2)*2))-1;
PSD1SZ = (PSDSZE(1,1))-1;
PSD2SZ = (PSDSZE(1,2))-1;

P1INs = PSDSZE(2,1);
P1INx = PSDLOC1(1,1) + P1INs;
P1INy = PSDLOC1(2,1) + P1INs;

P2INs = PSDSZE(2,2);
P2INx = (-1*PSDLOC(2,3)) + P2INs;
P2INy = PSDLOC(1,4) + P2INs;


P1OUTx=PSDLOC(2,1);
P1OUTy=PSDLOC(1,1);
P2OUTx=(-1*PSDLOC(2,3));
P2OUTy=PSDLOC(1,4);
%=================================%


%=================================%
%       MAIN 2D PLOT
%---------------------------------%
xlim = [0 Fcol];
ylim = [-Fcol Fcol];
zlim = [-2 2];
%---
figure(1)
subplot(5,5,[3 25]), 
AMPARPlot = gscatter(XYl(1,:),XYl(2,:));
axis([xlim, ylim]);
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(AMPARPlot,'marker','.','markersize',[6],'color',[1 0 0])
rectangle('Position',[P1OUTx,P1OUTy,SYN1SZ,SYN1SZ])
rectangle('Position',[P2OUTx,P2OUTy,SYN2SZ,SYN2SZ])
rectangle('Position',[P1INx,P1INy,PSD1SZ,PSD1SZ])
rectangle('Position',[P2INx,P2INy,PSD2SZ,PSD2SZ])

%=================================%
%           3D PLOTS
%---------------------------------%
figure(1);
subplot(5,5,[1 7]), 
gscatter(XYl(1,:),XYl(2,:)); view(20, 30);
axis normal;
grid off
axis([xlim, ylim, zlim]);
set(gca, 'Box', 'on');
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])

figure(1);
subplot(5,5,[11 22]), 
gscatter(XYl(1,:),XYl(2,:)); view(40, 70);
axis normal;
grid off
axis([xlim, ylim, zlim]);
set(gca, 'Box', 'on');
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
%=================================%


end