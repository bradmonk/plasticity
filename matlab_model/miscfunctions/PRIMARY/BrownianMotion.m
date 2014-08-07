function [] = BrownianMotion(varargin)
close all;clear all;clc;

% USER ENTERED VALUES
NSteps = 400;			% number of steps (loop time)
Ndots = 10;				% number of particles
DiffusionRate = 0.5;	% diffusion rate coefficient (microns^2 / sec)
dT = 1;					% time step (seconds)
Sc = 1/10;				% scale of model
POLYSz = [10 15];		% size of enclosure (XY in microns)

% DIFFUSION RATE EQUATIONS
dm = 2;                     % dimensions
Dr = DiffusionRate*dT/Sc;	% Diffusion Rate A (D = L^2 / 2d*t)
k = sqrt(dm*Dr);			% stdev of Dr step size distribution

% Particle Location Matrix
XYL = zeros(2,Ndots);		% XY particle locations
XYLp = zeros(2,NSteps);		% preallocate matrix for trace dot

% Figure Setup (for animation)
[Ph1 SPh2 XWIDE YHIGH] = FIGSET(XYL,XYLp,POLYSz,Sc); 

%===============================%
for Nt = 1:NSteps 

	XYS = (k * randn(2,Ndots));	% generates step sizes
	XYL = XYL+XYS;				% adds step to location
	XYLp(:,Nt) = XYL(:,1);		% save step of first dot (for trace)

	% Keep everything inside enclosure  %
	[XYL] = ENCLOSE(XYL,XWIDE,YHIGH,Ndots);

	% Plot live diffusion
	DOTPLOT(Nt,XYL,Ph1,XYLp,SPh2)
		
%-------------------------------%
if mod(Nt,100)==0;disp(Nt); end;
end
%===============================%
end
%===========================================%



%-------------------------------------------%
% ENCLOSE: keep particles inside polygon
%-------------------------------------------%
function [XYL] = ENCLOSE(XYL,XWIDE,YHIGH,Ndots)

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
function [] = DOTPLOT(Nt,XYL,Ph1,XYLp,SPh2)
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


%-------------------------------------------%
% FIGSET: figure setup
%-------------------------------------------%
function [Ph1 SPh2 XWIDE YHIGH] = FIGSET(XYL,XYLp,POLYSz,Sc)
%-------------------------------------------%
Flh = figure(1);
set(Flh,'Units','pixels');  scsz = get(0,'ScreenSize');
set(Flh,'OuterPosition',[scsz(3)/2  scsz(4)/3  scsz(3)/3.0  scsz(4)/1.8])
set(gcf,'Color',[.9,.9,.9])
%--------
POLYSz = POLYSz./Sc;		% scale enclosures 
XWIDE = POLYSz(1)/2;		% half X enclosure size (will double below)
YHIGH = POLYSz(2)/2;		% half X enclosure size (will double below)
BOARDER = [-5 -7.5 10 14.99] ./Sc; % [POLYSz(1)/2 POLYSz(2)/2 POLYSz(1) POLYSz(2)]
%--------
xlim = [-XWIDE XWIDE]; ylim = [-YHIGH YHIGH];
Ph1 = scatter(XYL(1,:),XYL(2,:),5,[0 0 1]);
axis([xlim, ylim]); hold on;
SPh2 = subplot(1,1,1); plot(XYLp(1,:),XYLp(2,:)); 
axis([xlim, ylim]); axis off; hold(SPh2,'on'); 
rectangle('Position',BOARDER);
%-------------------------------------------------%

end
