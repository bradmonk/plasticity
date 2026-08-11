function [varargout] = ExitNeck(varargin)
format compact;format short;
scsz = get(0,'ScreenSize');

%{

BIESS, KORKOTIAN, AND HOLCMAN 
PHYSICAL REVIEW E 76, 021922 2007

FIG. 2. Color online Molecular simulations in a three-dimensional 3D
model of a dendritic spine. The spine geometry is approximated by a
spherical head of radius R=0.5 m and a cylindrical neck with radius  and
length L. All particles were initially released from the center of the
spine head. a, b Mean sojourn time  standard of a diffusing molecule
in the spine neck and in the spine as a function of the spine neck length
a and the spine neck radius b assuming that no reentries into the spine
neck occur after the molecule entered the neck. The theoretical results as
predicted by formula 1 are shown as solid lines. c and d Same
settings as in a and b, except that multiple reentries into the spine
head were allowed. The theoretical prediction according to formula 2 are
shown as solid lines for =0.84 e Snapshots of one realization of the
Brownian trajectory as traced by one diffusing molecule in a 3D model of a
dendritic spine with possible reentries into the spine head R=0.5 m,
L=1.0 m, and =0.1 m. The snapshots are taken at different times t=0, 1,
5, 10, and 57.4 ms, where the particle reached the dendritic shaft.
Simulation parameters in a?e. Diffusion constant D=400 m2/s free
calcium, time step=1e-7 s, total simulation time: 250 ms.

--------------
the mean sojourn time of a diffusing molecule inside an empty spine of neck
length L, longer than the neck radius ? (L > ?) is (assuming no return)
given by:

? ? (V / 4D?) + (L² / 2D')		[1]	
--------------
where
?	mean sojourn time
L	spine neck length
?	neck radius
V	volume of spine head
D	diffusion coefficient

Formula [1] states the average spine escape time is the SUM of the 
mean first passage time to the head-to-neck portal (V / 4D?), plus,
mean time to travel through the neck (L² / 2D').
--------------
if...
?	???
L	.5
?	.1
V	.5
D	.5

tau =  (.5 / (4 *.5 * .1)) + ((.5^2) / (2 * .5))
tau = 2.75
(tau = 6.5 when D=.1)

? = 2.75 s
--------------

--------------
If a receptor can return to the head...
Po		probability of no return
?.h		mean sojourn time in spine head
t.n		mean time dot reaches dendrite before returning to head
t.r		mean time dot spends in neck before returning to head
x=??	where ?=1 ?=neck.radius ?? is exit portal location	

Po = ??L / 3D
t.n = L² / 6D
?.spine = (VL / 4D??²) + (L² / 2D)

thus ?.spine is the mean sojourn time of a diffusing molecule 
inside a spine (assuming it can return to the head).
--------------
if...
?	???
L = .5
e = .1
V = .5
D = .1

((V*L) / (4*L*e^2)) + ((L^2) / (2*D))

?.spine =  ((.5*.5) / (4*.5*.1^2)) + ((.5^2) / (2*.5))
?.spine = 12.75 s
(?.spine = 13.75 when D=.1)
--------------




--------------
if...
?	???
L = .5
e = .1
V = .5
D = .1

((V*L) / (4*L*e^2)) + ((L^2) / (2*D))

?.spine =  ((.5*.5) / (4*.5*.1^2)) + ((.5^2) / (2*.5))
?.spine = 12.75 s
(?.spine = 13.75 when D=.1)
--------------


L = .5
e = .1
V = .5
D = .1

neckL = [.2 .4 .6 .8 1.0]
neckR = [.05 .1 .15 .2 .25 .3]

tau1 = ((V) / (4*D*e)) + ((L^2) / (2*D))
tau2 = ((V*L) / (4*L*e^2)) + ((L^2) / (2*D))


%}



clc; close all; clear all;
%===========================================%
tau = [];	% mean sojourn time
L = .5;		% spine neck length
e = .1;		% neck radius
V = .5^3;	% volume of spine head
D = .5;		% diffusion coefficient

tauA = ((V) / (4*D*e)) + ((L^2) / (2*D));
tauB = ((V*L) / (4*D*(e^2))) + ((L^2) / (2*D));
%===========================================%
Ls = [.2 .4 .6 .8 1.0];
es = [.05 .1 .15 .2 .25 .3];

% Neck Length
tau1 = ((V) ./ (4*D*e)) + ((Ls.^2) ./ (2*D));
tau2 = ((V.*Ls) ./ (4*D .* (e^2))) + ((Ls.^2) ./ (2*D));

% Neck Radius
tau3 = ((V) ./ (4 * D .* es)) + ((L^2) ./ (2*D));
tau4 = ((V*L) ./ (4*D .* (es.^2))) + ((L^2) ./ (2*D));
%===========================================%

figure(1),
subplot(2,2,1),plot(Ls,tau1);
title('NECK LENGTH - NO RETURNS');
xlabel('Neck Length'); ylabel('Sojourn Time');
set(gca,'FontName','Helvetica','FontSize',12);
hold on;

subplot(2,2,2),plot(es,tau3)
title('NECK RADIUS - NO RETURNS');
xlabel('Neck Radius'); ylabel('Sojourn Time');

subplot(2,2,3),plot(Ls,tau2)
title('NECK LENGTH - WITH RETURNS');
xlabel('Neck Length'); ylabel('Sojourn Time');

subplot(2,2,4),plot(es,tau4)
title('NECK RADIUS - WITH RETURNS');
xlabel('Neck Radius'); ylabel('Sojourn Time');
















%===========================================%

figure(1), hold on;
subplot(2,2,1),plot(Ls,tau1);
tl = title('NECK LENGTH EFFECTS ON SPINE EXIT - NO RETURNS');
xl = xlabel('Neck Length'); yl = ylabel('Sojourn Time');
set([xl, yl],'FontSize',14,'FontName','Helvetica');
set([tl],'FontSize',12,'FontWeight','bold');

subplot(2,2,2),plot(es,tau3)
tl = title('NECK RADIUS EFFECTS ON SPINE EXIT - NO RETURNS');
xl = xlabel('Neck Radius'); yl = ylabel('Sojourn Time');
set([xl, yl],'FontSize',14,'FontName','Helvetica');
set([tl],'FontSize',12,'FontWeight','bold');

subplot(2,2,3),plot(Ls,tau2)
tl = title('NECK LENGTH EFFECTS ON SPINE EXIT - WITH RETURNS');
xl = xlabel('Neck Length'); yl = ylabel('Sojourn Time');
set([xl, yl],'FontSize',14,'FontName','Helvetica');
set([tl],'FontSize',12,'FontWeight','bold');

subplot(2,2,4),plot(es,tau4)
tl = title('NECK RADIUS EFFECTS ON SPINE EXIT - WITH RETURNS');
xl = xlabel('Neck Radius'); yl = ylabel('Sojourn Time');
set([xl, yl],'FontSize',14,'FontName','Helvetica');
set([tl],'FontSize',12,'FontWeight','bold');

set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, ['outputfigs/SpineNeckExit.png']);

%===========================================%
% USER ENTERED VALUES
NSteps = 200;			% number of steps (loop time)
Ndots = 10;				% number of particles
Scale = 1/10;			% scale of model (life:model)
TimeStep = 1000;		% time step (ms)
Delay = 0;				% slow animation rate by X%
POLYSz = [10 20];		% size of polygon enclosure (XY in µm)
DiffRateA = 0.2;		% diffusion rate coefficient A
DiffRateB = 0.1;		% diffusion rate coefficient B


% BASE DIFFUSION RATES EQUATIONS
Sc = Scale;					% scale of model (life:model)
t = TimeStep/1000;			% time step (ms)
dm = 2;                     % dimensions
Da = DiffRateA*t/Sc;		% Diffusion Rate A (D = L² / 2d*t)
Db = DiffRateB*t/Sc;		% Diffusion Rate B
Dr = Da/Db;					% Ratio of Da:Ds (1/Ls)^2;
Dn = Da/Dr;					% new D after scaling L
k = sqrt(dm*Da);			% stdev of D's step size distribution
L = sqrt(2*dm*Da);			% average diagonal (2D) step size
Lx = L/sqrt(2);             % average linear (1D) step size
Ls = 1/sqrt(Dr);			% scales Lx values for Dn
MSD = 2*dm*Da;				% mean squared displacement


XYL = zeros(2,Ndots);		% XY particle locations
XYS = zeros(2,Ndots);		% XY step sizes
XYLp = zeros(2,NSteps);		% preallocate matrix for trace dot

POLYSz = POLYSz./Scale;		% scale enclosures 
XWIDE = POLYSz(1)/2;		% half X enclosure size (will double below)
YHIGH = POLYSz(2)/2;		% half X enclosure size (will double below)
BOX = [-1 2 2 2] ./Scale;	% special area location [X Y W H]

fig1=figure(1);

%===========================================%
for Nt = 1:NSteps 

	XYS = (k * randn(2,Ndots));	% generates step sizes
	
	XYL = XYL+XYS;				% adds step to location
	
	XYLp(:,Nt) = XYL(:,1);		% save step of first dot (for trace)

	% Keep everything inside enclosure  %
	[XYL] = ENCLOSE(Nt,XYL,XWIDE,YHIGH,Ndots);

	% Plot live diffusion
	MAINPLOT(Nt,XYL,XWIDE,YHIGH,BOX,fig1);
	TRACEDOT(Nt,XYLp,XWIDE,YHIGH,BOX,fig1)

%-------------------------------%
if Delay>0; pause(Delay); end;
if mod(Nt,100)==0;disp(Nt); end;
end % for Nt = 1:Nsteps 
%===========================================%


%===============================%
varargout = {XYLp}; % export traced particle location
end
%===============================%



%-------------------------------------------%
% ENCLOSE: keep particles inside polygon
%-------------------------------------------%
function [XYL] = ENCLOSE(Nt,XYL,XWIDE,YHIGH,Ndots)

	for j = 1:Ndots 
		if XYL(1,j)>(XWIDE) || XYL(1,j)<(-XWIDE)
			XYL(1,j) = sign(XYL(1,j))*(XWIDE);
		elseif XYL(2,j)>(YHIGH) || XYL(2,j)<(-YHIGH)
			XYL(2,j) = sign(XYL(2,j))*(YHIGH);
		end
	end

end


%-------------------------------------------%
% MAINPLOT
%-------------------------------------------%
function [] = MAINPLOT(Nt,XYL,XWIDE,YHIGH,BOX,fig1)
%-------------------------------------------%
xlim = [-XWIDE XWIDE]; ylim = [-YHIGH YHIGH];
%---

figure(fig1)
subplot(1,2,1), 
scatter(XYL(1,:),XYL(2,:),5,[0 0 1]);
axis([xlim, ylim]);
rectangle('Position',BOX)
hold off;


%=================================%
%           3D PLOT
%---------------------------------%
%{
%=================================%
%           3D PLOT
%---------------------------------%
%     MAIN DOTS
%----------------------%
if do3DPLOT
figure(1);
subplot(5,5,[1 7]), 
%gscatter(GluR2xyl(1,:),GluR2xyl(2,:)); view(20, 30);
scatter3(GluR2xyl(1,:),GluR2xyl(2,:),G2Z(:),'.'), view(20, 30)
axis normal;
grid off
axis([xlim, ylim, zlim]);
set(gca, 'Box', 'on');
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
%-------%
%}
% if Nt == 201; keyboard; end
%=================================%

end


%-------------------------------------------%
% TRACEDOT
%-------------------------------------------%
function [] = TRACEDOT(Nt,XYLp,XWIDE,YHIGH,BOX,fig1)
%-------------------------------------------%
xlim = [-XWIDE XWIDE]; ylim = [-YHIGH YHIGH];
%---

LLg = 3;

if Nt > LLg
xp = (XYLp(1,(Nt-LLg):Nt));
yp = (XYLp(2,(Nt-LLg):Nt));
else
xp = (XYLp(1,Nt));
yp = (XYLp(2,Nt));
end


figure(fig1)
subplot(1,2,2),plot(xp,yp);
axis([xlim, ylim]);
rectangle('Position',BOX)
hold on;


 
end





