function [varargout] = SimpleDiffusion(varargin)
format compact;format short;close all;clc
scsz = get(0,'ScreenSize');


% USER ENTERED VALUES
NSteps = 100;			% number of steps (loop time)
Ndots = 10;				% number of particles
Scale = 1/10;			% scale of model (life:model)
TimeStep = 1000;		% time step (ms)
Delay = 0;				% slow animation rate by X%
POLYSz = [10 15];		% size of polygon enclosure (XY in µm)
DiffRateA = 0.5;		% diffusion rate coefficient A
DiffRateB = 0.1;		% diffusion rate coefficient B


% BASE DIFFUSION RATES EQUATIONS
Sc = Scale;					% scale of model (life:model)
t = TimeStep/1000;			% time step (ms)
dm = 2;                     % dimensions
Da = DiffRateA*t/Sc;		% Diffusion Rate A (D = L² / 2d*t)
Db = DiffRateB*t/Sc;		% Diffusion Rate B
Dr = Da/Db;					% Ratio of Da:Ds (1/Ls)^2;
Dn = Da/Dr;					% new D after scaling L
k = sqrt(dm*Da);			% stdev of Ds step size distribution
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
BOARDER = [-5 -7.5 10 14.99] ./Scale; % [POLYSz(1)/2 POLYSz(2)/2 POLYSz(1) POLYSz(2)]


%================================================%
%               FIGURE SETUP
%------------------------------------------------%
Flh = figure(1);
set(Flh,'Units','pixels');  scnsize = get(0,'ScreenSize');
pos1 = [scnsize(3)/2  scnsize(4)/3  scnsize(3)/3.0  scnsize(4)/1.8];
set(Flh,'OuterPosition',pos1)
set(gcf,'Color',[.9,.9,.9])
%--------
xlim = [-XWIDE XWIDE]; ylim = [-YHIGH YHIGH];
%---
Ph1 = scatter(XYL(1,:),XYL(2,:),5,[0 0 1]);
axis([xlim, ylim]);
hold on
SPh2 = subplot(1,1,1);
Ph2 = plot(XYLp(1,:),XYLp(2,:));
axis([xlim, ylim]);
axis off
hold(SPh2,'on');
% rectangle('Position',BOX)
rectangle('Position',BOARDER)
%-------------------------------------------------%




%===========================================%
for Nt = 1:NSteps 

	XYS = (k * randn(2,Ndots));	% generates step sizes
	XYL = XYL+XYS;				% adds step to location
	XYLp(:,Nt) = XYL(:,1);		% save step of first dot (for trace)

	% Keep everything inside enclosure  %
	[XYL] = ENCLOSE(Nt,XYL,XWIDE,YHIGH,Ndots);

	% Plot live diffusion
	DUALPLOT(Nt,XYL,Ph1,XYLp,SPh2)
	
	
	
	% XYLS((Nt*2-1):(Nt*2),:) = XYL;
	XL(Nt,:) = XYL(1,:);
	YL(Nt,:) = XYL(2,:);
	
%-------------------------------%
% pause(Delay)
if mod(Nt,100)==0;disp(Nt); end;
end % for Nt = 1:Nsteps 
%===========================================%


%===============================%
% varargout = {XYLp}; % export traced particle location
varargout = {XL,YL};
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
% DUALPLOT: plots diffusion
%-------------------------------------------%
function [] = DUALPLOT(Nt,XYL,Ph1,XYLp,SPh2)
%-------------------------------------------%

LLg = 3;
if Nt > LLg
xp = (XYLp(1,(Nt-LLg):Nt));
yp = (XYLp(2,(Nt-LLg):Nt));
else
xp = (XYLp(1,Nt));
yp = (XYLp(2,Nt));
end

set(Ph1,'XData',XYL(1,:),'YData',XYL(2,:));
drawnow

plot(xp,yp,'Parent',SPh2);
drawnow

 
end


