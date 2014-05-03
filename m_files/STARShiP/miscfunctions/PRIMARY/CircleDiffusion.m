function [varargout] = CircleDiffusion(varargin)
format compact;format short;close all;
scsz = get(0,'ScreenSize');
 
 
 
% USER ENTERED VALUES
NSteps = 500;                   % number of steps (loop time)
Ndots = 10;                     % number of particles
Scale = 1/10;                   % scale of model (life:model)
TimeStep = 10;	                % time step (ms)
Delay = .1;						% slow animation rate by X%
Rsphere = 5;
XYZSize = [20 20 (Rsphere*4)];  % size of enclosure (XYZ in µm)
DiffRateA = 0.2;                % diffusion rate coefficient A
DiffRateB = 0.1;                % diffusion rate coefficient B


 
% BASE DIFFUSION RATES EQUATIONS
Sc = Scale;             % scale of model (life:model)
t = TimeStep/1000;      % time step (ms)
dm = 2;                 % dimensions
Da = DiffRateA*t/Sc;    % Diffusion Rate A (D = L² / 2d*t)
Db = DiffRateB*t/Sc;    % Diffusion Rate B
Dr = Da/Db;             % Ratio of Da:Ds (1/Ls)^2;
Dn = Da/Dr;             % new D after scaling L
k = sqrt(dm*Da);        % stdev of D's step size distribution
L = sqrt(2*dm*Da);      % average diagonal (2D) step size
Lx = L/sqrt(2);         % average linear (1D) step size
Ls = 1/sqrt(Dr);        % scales Lx values for Dn
MSD = 2*dm*Da;          % mean squared displacement

 
TPL = zeros(2,Ndots);
TPS = zeros(2,Ndots);

% SCALE THINGS
Rsphere = Rsphere./Scale;
XYZSize = XYZSize./Scale;       % scale enclosures
XWIDE = XYZSize(1)/2;           % half X enclosure size (will double below)
YHIGH = XYZSize(2)/2;           % half X enclosure size (will double below)
ZHIGH = XYZSize(3)/2;
BOX = [-1 2 2 2] ./Scale;       % special area location [X Y W H]
r = Rsphere;
 
%===============================%
for Nt = 1:NSteps
   
        TPS = (k * randn(2,Ndots));     % generates step sizes
        TPL = TPL+TPS;					% adds step to location
        TPLp(:,Nt) = TPL(:,1);			% save step of first dot (for trace)
		

        %MAINPLOT(Nt,TPL,XWIDE,YHIGH,BOX);
		TRACEDOT2(Nt,TPLp,TPL,XWIDE,YHIGH,ZHIGH,BOX,Rsphere);
        %TRACEDOT(Nt,TPLp,TPL,XWIDE,YHIGH,ZHIGH,BOX,Rsphere);

% pause(Delay)
end % for Nt = 1:Nsteps


 
FINALTRACEDOT(Nt,TPLp,XWIDE,YHIGH,ZHIGH,BOX,Rsphere);
%===============================%
 
varargout = {TPLp}; % export traced particle location
end % end main function
%=====================================================%
 
 
 
%-------------------------------------------%
% TRACEDOT
%-------------------------------------------%
function [] = TRACEDOT(Nt,TPLp,TPL,XWIDE,YHIGH,ZHIGH,BOX,Rsphere)
%-------------------------------------------%
xlim = [-XWIDE XWIDE]; ylim = [-YHIGH YHIGH]; zlim = [-ZHIGH ZHIGH];
%---
 
r = Rsphere;
 
xp = r*sin(TPLp(1,:)).*cos(TPLp(2,:));
yp = r*sin(TPLp(1,:)).*sin(TPLp(2,:));
zp = r*cos(TPLp(1,:));

 
figure(1)
plot3(xp,yp,zp);
axis([xlim, ylim, zlim]);
grid on;
%---

% if mod(Nt,10)==0
% axis vis3d
% camorbit(10,0,'camera')
% drawnow
% end
 
end


%-------------------------------------------%
% TRACEDOT2
%-------------------------------------------%
function [] = TRACEDOT2(Nt,TPLp,TPL,XWIDE,YHIGH,ZHIGH,BOX,Rsphere)
%-------------------------------------------%
xlim = [-XWIDE XWIDE]; ylim = [-YHIGH YHIGH]; zlim = [-ZHIGH ZHIGH];
%---
 
r = Rsphere;
LLg = 3;

if Nt > LLg
xp = r*sin(TPLp(1,(Nt-LLg):Nt)).*cos(TPLp(2,(Nt-LLg):Nt));
yp = r*sin(TPLp(1,(Nt-LLg):Nt)).*sin(TPLp(2,(Nt-LLg):Nt));
zp = r*cos(TPLp(1,(Nt-LLg):Nt));
else
xp = r*sin(TPLp(1,Nt)).*cos(TPLp(2,Nt));
yp = r*sin(TPLp(1,Nt)).*sin(TPLp(2,Nt));
zp = r*cos(TPLp(1,Nt));
end

 
figure(1)
plot3(xp,yp,zp);
axis([xlim, ylim, zlim]);
grid on;
hold on;
%---

if mod(Nt,10)==0
axis vis3d
camorbit(10,0,'camera')
drawnow
end
 
end




%-------------------------------------------%
% FINALTRACEDOT
%-------------------------------------------%
function [] = FINALTRACEDOT(Nt,TPLp,XWIDE,YHIGH,ZHIGH,BOX,Rsphere)
%-------------------------------------------%
xlim = [-XWIDE XWIDE]; ylim = [-YHIGH YHIGH]; zlim = [-ZHIGH ZHIGH];
%---
 
r = Rsphere;
 
xp = r*sin(TPLp(1,:)).*cos(TPLp(2,:));
yp = r*sin(TPLp(1,:)).*sin(TPLp(2,:));
zp = r*cos(TPLp(1,:));
 
 
 
figure(1)
plot3(xp,yp,zp);
axis([xlim, ylim, zlim]);
rectangle('Position',BOX)
hold on;
rotate3d on
axis vis3d
for i=1:36
   camorbit(10,0,'data',[0 1 0])
   drawnow
   pause(.1)
end
axis vis3d
camorbit(10,0,'camera')
drawnow

 
end