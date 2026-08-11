function [] = TrackMSD()


D = 0.3;		% Diffusion Rate
Scale = 1;		% scale of model
t = .1;			% time step


Dpsd = .1;		% Alternative Diffusion Rate (modifies k in special region)

trackMSD = 1;
t = 1;
GluR2Ndots = 100;
Nsteps = 100;
MSDdrop = 1;



% BASE DIFFUSION RATES EQUATIONS
d = 2;                      % dimensions
D = D*t/Scale;				% Diffusion Rate ES (D = L² / 2d*t)
Dp = Dpsd*t/Scale;			% Diffusion Rate PSD
Dr = D/Dp;					% Ratio of D:Ds (1/Ls)^2;
Dn = D/Dr;					% new D after scaling L
k = sqrt(d*D);	            % stdev of D's step size distribution
MSD = 2*d*D;                % mean squared displacement
L = sqrt(2*d*D);            % average diagonal (2D) step size
Lx = L/sqrt(2);             % average linear (1D) step size
Ls = 1/sqrt(Dr);			% scales Lx values for Dn

% PSD DIFFUSION RATES (DEFAULT: GLUR1-BASED)
Ds = D;						% ExtraSynaptic Model-Scaled Diffusion Rate 
PSD1Ds = Dpsd*t/Scale;		% (PSD1) Model-Scaled Diffusion Rate
PSD2Ds = Dpsd*t/Scale;		% (PSD2) Model-Scaled Diffusion Rate
Dr_PSD1 = Ds/PSD1Ds;		% (PSD1) Ratio Ds:PSD1Ds
Dr_PSD2 = Ds/PSD2Ds;		% (PSD2) Ratio Ds:PSD2Ds
Dn_PSD1 = Ds/Dr_PSD1;		% (PSD1) D value after scaling L
Dn_PSD2 = Ds/Dr_PSD2;		% (PSD2) D value after scaling L
PSD1 = 1/sqrt(Dr_PSD1);		% (PSD1) D Lx-Scalar (Ls scales Lx so Dn_PSD1 = Ds:PSD1Ds)
PSD2 = 1/sqrt(Dr_PSD2);		% (PSD2) D Lx-Scalar (Ls scales Lx so Dn_PSD2 = Ds:PSD2Ds)


% GLUR2
GluR2xyds = ones(2,GluR2Ndots);
GluR2xyl = ones(2,GluR2Ndots);
G2Z = ones(1,GluR2Ndots);


% TRACK MSD AND STEP SIZE MEANS (0:no 1:yes)
dT = t;
Nsteps = GluR2Ndots;
tracks = cell(GluR2Ndots, 1);
lims = ((D+1)^2)*10;
trackStepMeans = 0;
strdropES = 'D-ES';
strdropPSD1 = 'D-PSD1';
strdropPSD2 = 'D-PSD2';
testESMSD = strcmp(MSDdrop,strdropES);
testPSD1MSD = strcmp(MSDdrop,strdropPSD1);
testPSD2MSD = strcmp(MSDdrop,strdropPSD2);
MSDtest = [testESMSD testPSD1MSD testPSD2MSD];
MANstepsize = 0;


stepN = 1;
for Nt = 1:Nsteps 
%-------------------------------%
    %-------------------------------%
    %     STEP DIRECTION & SIZE
    %-------------------------------%
    GluR2xyds = STEPxyds(GluR2Ndots, k);
	%-------------------------------%
    %     DO MANUAL STEP SIZES
    %-------------------------------%
	if MANstepsize
		GluR2xyds = CIRCSTEPxyds(GluR2Ndots, k);
	end % xyd = DIRxyd(xyd);
	
		
	[GluR2xyl GluR2xyds] = MSDAMPARSTEP(GluR2Ndots, GluR2xyds, GluR2xyl,... 
	testPSD1MSD, testPSD2MSD, testESMSD, PSD1, PSD2, D, Dn_PSD1, Dn_PSD2);


    [tracks] = MSDfun(stepN, Nsteps, tracks, GluR2xyds);

    
    if trackStepMeans
     if mod(stepN, 10) == 0
        StepMeans = meanStep(GluR2xyds, k, D);
        head = ['MeanAbsStep CalcHalfNorm MeanPythStep CalcPythStep CalcHalfPythStep'];
       disp(head)
       disp(StepMeans)
     end
    end


 stepN = stepN+1;
end
output1 = tracks;
[rawMSDout] = MSDfunction(tracks, GluR2Ndots, Nsteps, D, Dn_PSD1, Dn_PSD2, dT, k, MSDtest);


%===========================================================%
%               LIVE PARTICLE DIFFUSION
%-----------------------------------------------------------%
for Nt = 1:Nsteps

    GluR2xyds = STEPxyds2(GluR2Ndots, k);
	[GluR2xyl] = AMPARSTEP2(GluR2Ndots, GluR2xyds, GluR2xyl);
	MAINPLOT2(GluR2xyl, lims);

end
%===========================================================%


%===========================================================%
%               MSD RANDOM STEPS ANALYSIS
%-----------------------------------------------------------%
tracks = cell(GluR2Ndots, 1);

stepN = 1;
for Nt = 1:Nsteps 

	GluR2xyds = STEPxyds2(GluR2Ndots, k);
	[GluR2xyl] = AMPARSTEP2(GluR2Ndots, GluR2xyds, GluR2xyl);
    [tracks] = MSDfun2(stepN, Nsteps, tracks, GluR2xyds);

stepN = stepN+1;
end
[ESMSDout] = MSDfunction2(tracks,GluR2Ndots,Nsteps,D,Dn,L,dT,k,Scale,MSDtest);

%===========================================================%
%               MSD UNIFORM STEPS ANALYSIS
%-----------------------------------------------------------%
stepN = 1;
for t = 1:Nsteps 

	GluR2xyds = UNIFORMSTEPS2(GluR2Ndots, Lx);
	[GluR2xyl GluR2xyds] = SCALEUNIFORMSTEPS2(GluR2Ndots, GluR2xyds, GluR2xyl, Ls, MSDtest);
    [tracks] = MSDfun2(stepN, Nsteps, tracks, GluR2xyds);

stepN = stepN+1;
end
[PSDMSDout] = MSDfunction2(tracks,GluR2Ndots,Nsteps,D,Dn,L,dT,k,Scale,MSDtest);

figure(80)
set(gcf,'OuterPosition',[600,400,300,300])
MSDfig = gcf;
figure(MSDfig)
subplot(2,1,1),text(0.5,0.5,strcat('MSD:',' \bullet  ', num2str(ESMSDout), '  µm²/s'),...
	'FontSize',24,'HorizontalAlignment','center','BackgroundColor',[.7 .9 .7]);
subplot(2,1,2),text(0.5,0.5,strcat('PSD:',' \bullet  ', num2str(PSDMSDout), '  µm²/s'),...
	'FontSize',24,'HorizontalAlignment','center','BackgroundColor',[.7 .9 .7]);

end % main function
%===========================================%





%===================================%
% STEP SIZE GENERATOR
%===================================%
function xyds = STEPxyds(Ndots, k)
    xyds = (k * randn(2,Ndots));
end
%===================================%



%===================================%
% MOVE PARTICLES MAIN FUNCTION
%===================================%
function [GluR1xyds GluR1xyl GluR2xyds GluR2xyl]...
    = MOVEGLUR(stepN,XWIDE,YHIGH,...
	GluR1Ndots,GluR1xyds,GluR1xyl,...
    GluR2Ndots,GluR2xyds,GluR2xyl)

%===================================%
GluR1xyl = GluR1xyl+GluR1xyds;
 
for j = 1:GluR1Ndots 
    if GluR1xyl(1,j)>(XWIDE) || GluR1xyl(1,j)<(0)
            GluR1xyl(1,j) = uint8(sign(GluR1xyl(1,j)))*(XWIDE);
    elseif GluR1xyl(2,j)>(YHIGH) || GluR1xyl(2,j)<(-YHIGH)
            GluR1xyl(2,j) = sign(GluR1xyl(2,j))*(YHIGH);
    end    
        
end
%===================================%
GluR2xyl = GluR2xyl+GluR2xyds;

for j = 1:GluR2Ndots 
	if GluR2xyl(1,j)>(XWIDE) || GluR2xyl(1,j)<(0)
            GluR2xyl(1,j) = uint8(sign(GluR2xyl(1,j)))*(XWIDE);
	elseif GluR2xyl(2,j)>(YHIGH) || GluR2xyl(2,j)<(-YHIGH)
            GluR2xyl(2,j) = sign(GluR2xyl(2,j))*(YHIGH);
	end	   
        
end
%===================================%
% if stepN == 500; keyboard; end
end
%===================================%


%===================================%
%		inboxfun
%===================================%
% Tests whether particles are in a box polygon
% and returns a logical vector
%===================================%
function [inbox] = inboxfun(LB,RT,xyl)

if LB(1)>RT(1)
	LBt=LB;
	RTt=RT;
	LB(1)=RTt(1);
	RT(1)=LBt(1);
end
if LB(2)>RT(2)
	LBt=LB;
	RTt=RT;
	LB(2)=RTt(2);
	RT(2)=LBt(2);
end

xylLB1 = xyl(1,:) > LB(1);
xylRT1 = xyl(1,:) < RT(1);
xylLB2 = xyl(2,:) > LB(2);
xylRT2 = xyl(2,:) < RT(2);

inbox = xylLB1 & xylRT1 & xylLB2 & xylRT2;

end
%===================================%


%-------------------------------------------%
% ONEDOTPLOT2
%-------------------------------------------%
function [] = ONEDOTPLOT(stepN,Nsteps,GluR2xyl,xyl2,XWIDE,YHIGH,...
				XYBOXpr1,XYBOXpr2,XYBOXp1,XYBOXp2)
%-------------------------------------------%
xlim = [0 XWIDE]; ylim = [-YHIGH YHIGH];
%---
 

LLg = 3;

if stepN > LLg
xp = xyl2(1,(stepN-LLg):stepN);
yp = xyl2(2,(stepN-LLg):stepN);
else
xp = xyl2(1,stepN);
yp = xyl2(2,stepN);
end


figure(1)
subplot(5,5,[3 25]), plot(xp, yp);
axis([xlim, ylim]);
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(get(gca,'XLabel'),'String','GluR2:Blue \bullet \bullet GluR1:Red')
rectangle('Position',[XYBOXpr1(1),XYBOXpr1(2),XYBOXpr1(3),XYBOXpr1(4)])
rectangle('Position',[XYBOXpr2(1),XYBOXpr2(2),XYBOXpr2(3),XYBOXpr2(4)])
rectangle('Position',[XYBOXp1(1),XYBOXp1(2),XYBOXp1(3),XYBOXp1(4)])
rectangle('Position',[XYBOXp2(1),XYBOXp2(2),XYBOXp2(3),XYBOXp2(4)])
hold on;
%---

% if mod(Nt,10)==0
% axis vis3d
% camorbit(10,0,'camera')
% drawnow
% end
 
end




%%		  MSD BROWNIAN MOTION ANALYSIS TOOLS
%-------------##########################------------------%
%			 DIFFUSION ANALYSIS FUNCTIONS
%-------------##########################------------------%
%-------------------------------------------%
% MSD PARTICLE MOVEMENT GENERATION
%-------------------------------------------%
function [GluR2xyl GluR2xyds] = MSDAMPARSTEP(GluR2Ndots, GluR2xyds, GluR2xyl,... 
	testPSD1MSD, testPSD2MSD, testESMSD, PSD1, PSD2, D, Dn_PSD1, Dn_PSD2)
	
	for j = 1:GluR2Ndots
        
		if testPSD1MSD
           GluR2xyds(:,j) = GluR2xyds(:,j)*PSD1;
        elseif testPSD2MSD
            GluR2xyds(:,j) = GluR2xyds(:,j)*PSD2;
		elseif testESMSD
            GluR2xyds(:,j) = GluR2xyds(:,j);
		end
		
	   
        GluR2xyl(:,j) = GluR2xyl(:,j)+GluR2xyds(:,j);
	end

end

%-------------------------------------------%
% MEAN SQUARED DISPLACEMENT FUNCTION
%-------------------------------------------%
function [tracks] = MSDfun(stepN, Nsteps, tracks, GluR2xyds)
    time = (0:Nsteps-1)';
    xymsd = GluR2xyds';
    xymsd = cumsum(xymsd,1);
    tracks{stepN} = [time xymsd];

end


%-------------------------------------------%
% MSDfunction - FINAL TRACKS ANALYSIS
%-------------------------------------------%
function [MSDout] = MSDfunction(tracks, GluR2Ndots, Nsteps, D, Dn_PSD1, Dn_PSD2, dT, k, MSDtest)

SPACE_UNITS = 'µm';
TIME_UNITS = 's';
N_PARTICLES = GluR2Ndots;
N_TIME_STEPS = Nsteps;
N_DIM = 2;

oDes = D;				% raw		µm^2/s
D  = D*.1;              % to-scale	µm^2/s

oDpsd1 = Dn_PSD1;		% raw		µm^2/s
Dpsd1 = Dn_PSD1*.1;		% to-scale	µm^2/s

oDpsd2 = Dn_PSD2;		% raw		µm^2/s
Dpsd2 = Dn_PSD2*.1;		% to-scale	µm^2/s

dTbase = dT;			% raw time-step 
dT = dT*1;				% to-scale time-step
k = k;					% stdv of step distribution

ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
ma = ma.addAll(tracks);
disp(ma)

figure
ma.plotTracks;
ma.labelPlotTracks;

ma = ma.computeMSD;
ma.msd;

t = (0 : N_TIME_STEPS)' * dT;
[T1, T2] = meshgrid(t, t);
all_delays = unique( abs(T1 - T2) );

figure
ma.plotMSD;


cla
ma.plotMeanMSD(gca, true)

mmsd = ma.getMeanMSD;
t = mmsd(:,1);
x = mmsd(:,2);
dx = mmsd(:,3) ./ sqrt(mmsd(:,4));
errorbar(t, x, dx, 'k')

[fo, gof] = ma.fitMeanMSD;
plot(fo)
ma.labelPlotMSD;
legend off


ma = ma.fitMSD;

good_enough_fit = ma.lfit.r2fit > 0.8;
Dmean = mean( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;
Dstd  =  std( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;

Dheader1 = ['Raw Unscaled Values'];
Dhead1 = ['     D.es	D.psd1	D.psd2'];
Ddat1 = [oDes oDpsd1 oDpsd2];
disp(' ')
disp(Dheader1)
disp(Dhead1)
disp(Ddat1)


yourtesthead = ['YOU ARE TESTING DIFFUSION FOR:'];
if MSDtest(1)
	yourtest = ['   Des:   extrasynaptic diffusion rate'];
elseif MSDtest(2)
	yourtest = ['   Dpsd1:  PSD-1 diffusion rate'];
elseif MSDtest(3)
	yourtest = ['   Dpsd2:  PSD-2 diffusion rate'];
else
	yourtest = ['   generic diffusion rate'];
end
disp(yourtesthead)
disp(yourtest)

disp(' ')
fprintf('Estimation of raw D coefficient from MSD:\n')
fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', ...
    Dmean, Dstd, sum(good_enough_fit));




% Retrieve instantaneous velocities, per track
 trackV = ma.getVelocities;

 % Pool track data together
 TV = vertcat( trackV{:} );

 % Velocities are returned in a N x (nDim+1) array: [ T Vx Vy ...]. So the
 % velocity vector in 2D is:
 V = TV(:, 2:3);

 % Compute diffusion coefficient
varV = var(V);
mVarV = mean(varV); % Take the mean of the two estimates
Dest = mVarV / 2 * dT;



Dheader2 = ['Scaling to model (10units = 1µm)...'];
Dhead2 = ['     D.es	D.psd1	D.psd2'];
Ddat2 = [D Dpsd1 Dpsd2];

disp(' ')
disp(Dheader2)
disp(Dhead2)
disp(Ddat2)
fprintf('Estimation from velocities histogram:\n')
fprintf('Tested D = %.3g %s, compare to scaled Des value of %.3g %s\n', ...
    Dest, [SPACE_UNITS '²/' TIME_UNITS], D, [SPACE_UNITS '²/' TIME_UNITS]);

% printf('D.psd target value was %.3g %s\n', ...
%     Dest, msdDpsd, [SPACE_UNITS '²/' TIME_UNITS]);

MSDout = D;

end



%-------------############################------------------%
%			  ##   MSD2 SUBFUNCTIONS    ##
%-------------############################------------------% 

%-------------------------------------------%
% STEP SIZE GENERATOR
%-------------------------------------------%
function GluR2xyds = STEPxyds2(GluR2Ndots, k)

    GluR2xyds = (k * randn(2,GluR2Ndots));

end


%-------------------------------------------%
% MOVE PARTICLES MAIN FUNCTION
%-------------------------------------------%
function [GluR2xyl] = AMPARSTEP2(GluR2Ndots, GluR2xyds, GluR2xyl)
	
	for j = 1:GluR2Ndots
        GluR2xyl(:,j) = GluR2xyl(:,j)+GluR2xyds(:,j);
	end
	
end


%-------------------------------------------%
% LIVE DIFFUSION PLOT
%-------------------------------------------%
function [] = MAINPLOT2(GluR2xyl, lims)
%-------------------------------------------%
xlim = [-lims lims];
ylim = [-lims lims];
zlim = [-5 5];

%=================================%
%       MAIN 2D PLOT
%---------------------------------%
figure(1)
subplot(2,1,1), 
AMPARPlot = gscatter(GluR2xyl(1,:),GluR2xyl(2,:));
axis([xlim, ylim]);
set(AMPARPlot,'marker','.','markersize',[6],'color',[1 0 0])


%=================================%
%           3D PLOT
%---------------------------------%
figure(1);
subplot(2,1,2), 
gscatter(GluR2xyl(1,:),GluR2xyl(2,:)); view(20, 30);
axis normal;
grid off
axis([xlim, ylim, zlim]);
set(gca, 'Box', 'on');

end


%-------------------------------------------%
% MANUAL STEP SIZE FUNCTION
%-------------------------------------------%
function GluR2xyds = UNIFORMSTEPS2(GluR2Ndots, Lx)
%-------------------------------------------%

   Lx(1:2,1:GluR2Ndots) = Lx;
   xyd = randi([0 1],GluR2Ndots,2)';
   xyd(xyd == 0) = -1;
   GluR2xyds = (Lx.*xyd);
   
end


%-------------------------------------------%
% MSD SCALED STEPS FUNCTION
%-------------------------------------------%
function [GluR2xyl GluR2xyds] = SCALEUNIFORMSTEPS2(GluR2Ndots, GluR2xyds, GluR2xyl, Ls, MSDtest)


if MSDtest(1)
	Ls = 1;	
end
	
	for j = 1:GluR2Ndots
        GluR2xyds(:,j) = GluR2xyds(:,j)*Ls;
        GluR2xyl(:,j) = GluR2xyl(:,j)+GluR2xyds(:,j);
	end	
	
end


%-------------------------------------------%
% MSD TRACKS GENERATOR
%-------------------------------------------%
function [tracks] = MSDfun2(stepN, Nsteps, tracks, GluR2xyds)
    time = (0:Nsteps-1)';
    xymsd = GluR2xyds';
    xymsd = cumsum(xymsd,1);
    tracks{stepN} = [time xymsd];

end


%-------------------------------------------%
% MSD TRACKS ANALYSIS
%-------------------------------------------%
function [Dest] = MSDfunction2(tracks,GluR2Ndots,Nsteps,D,Dn,L,dT,k,Scale,MSDtest)


% printf('D.psd target value was %.3g %s\n', ...
%     Dest, msdDpsd, [SPACE_UNITS '²/' TIME_UNITS]);

yourtesthead = ['-------TESTING DIFFUSION---------'];
if MSDtest(1)
	yourtest = ['   D:   original diffusion rate'];
elseif MSDtest(2)
	yourtest = ['   Dn:  PSD diffusion rate'];
elseif MSDtest(3)
	yourtest = ['   L:  step length'];
else
	yourtest = ['   generic diffusion rate'];
end
disp(yourtesthead)
disp(yourtest)



SPACE_UNITS = 'µm';
TIME_UNITS = 's';
N_PARTICLES = GluR2Ndots;
N_TIME_STEPS = Nsteps;
N_DIM = 2;

oD = D;				% raw		µm^2/s
D  = D*Scale;       % to-scale	µm^2/s

oDn = Dn;			% raw		µm^2/s
Dn = Dn*Scale;		% to-scale	µm^2/s

oL = L;				% raw		µm
L = L*Scale;		% to-scale	µm

dTbase = dT;		% raw time-step 
dT = dT*Scale;		% to-scale time-step
k = k;				% stdv of step distribution

ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
ma = ma.addAll(tracks);
disp(ma)

figure
ma.plotTracks;
ma.labelPlotTracks;

ma = ma.computeMSD;
ma.msd;

t = (0 : N_TIME_STEPS)' * dT;
[T1, T2] = meshgrid(t, t);
all_delays = unique( abs(T1 - T2) );

figure
ma.plotMSD;


cla
ma.plotMeanMSD(gca, true)

mmsd = ma.getMeanMSD;
t = mmsd(:,1);
x = mmsd(:,2);
dx = mmsd(:,3) ./ sqrt(mmsd(:,4));
errorbar(t, x, dx, 'k')

[fo, gof] = ma.fitMeanMSD;
plot(fo)
ma.labelPlotMSD;
legend off


ma = ma.fitMSD;

good_enough_fit = ma.lfit.r2fit > 0.8;
Dmean = mean( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;
Dstd  =  std( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;

Dheader1 = ['Raw Unscaled Values'];
Dhead1 = ['    D        Dn        L'];
Ddat1 = [oD oDn oL];
disp(' ')
disp(Dheader1)
disp(Dhead1)
disp(Ddat1)


yourtesthead = ['YOU ARE TESTING DIFFUSION FOR:'];
if MSDtest(1)
	yourtest = ['   D:   original diffusion rate'];
elseif MSDtest(2)
	yourtest = ['   Dn:  PSD diffusion rate'];
elseif MSDtest(3)
	yourtest = ['   L:  step length'];
else
	yourtest = ['   generic diffusion rate'];
end
disp(yourtesthead)
disp(yourtest)

disp(' ')
fprintf('Estimation of raw D coefficient from MSD:\n')
fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', ...
    Dmean, Dstd, sum(good_enough_fit));




% Retrieve instantaneous velocities, per track
 trackV = ma.getVelocities;

 % Pool track data together
 TV = vertcat( trackV{:} );

 % Velocities are returned in a N x (nDim+1) array: [ T Vx Vy ...]. So the
 % velocity vector in 2D is:
 V = TV(:, 2:3);

 % Compute diffusion coefficient
varV = var(V);
mVarV = mean(varV); % Take the mean of the two estimates
Dest = mVarV / 2 * dT;



Dheader2 = ['Scaling to model...'];
Dhead2 = ['    D        Dn        L'];
Ddat2 = [D Dn L];

disp(' ')
disp(Dheader2)
disp(Dhead2)
disp(Ddat2)
fprintf('Estimation from velocities histogram:\n')
fprintf('Computed scaled velocity D = %.3g %s, generated from set D = %.3g %s\n', ...
    Dest, [SPACE_UNITS '²/' TIME_UNITS], D, [SPACE_UNITS '²/' TIME_UNITS]);


end




